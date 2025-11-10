function simulation_from_post(state, session_idx, parameters)
% Generate simulated data for the GLM model
% Use inferred post data
% Validate possible model structures
% UseGPU = gpuDeviceCount > 0;

if nargin < 3 || isempty(parameters)
    parameters = struct();
    parameters.parameter_tag = 'Default';
elseif ~isstruct(parameters)
    error('simulation_from_post:InvalidParameters', 'parameters must be a struct');
elseif ~isfield(parameters, 'parameter_tag')
    error('simulation_from_post:MissingField', 'parameters struct must include parameter_tag');
end
%% Get root folder
code_depth = 4;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
% include code folder and utils
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Main
do_plot = false;
UseGPU = false;
partial = true;

if isfield(parameters, 'do_plot')
    do_plot = logical(parameters.do_plot);
    parameters = rmfield(parameters, 'do_plot');
end

if isfield(parameters, 'UseGPU')
    UseGPU = logical(parameters.UseGPU);
    parameters = rmfield(parameters, 'UseGPU');
end

if isfield(parameters, 'partial')
    partial = logical(parameters.partial);
    parameters = rmfield(parameters, 'partial');
end

kernel_name = 'DeltaPure';
reg = 'L2=0_2';

n_trials = 30;
B = 30000;

% groups = {'Thal', 'ACC-E', 'ACC-I', 'VLPFC-E', 'VLPFC-I'};
groups = {'ACC', 'VLPFC', 'Thal-ACC', 'Thal-VLPFC'}; % Two groups of thalamic neurons. Project to one area only
n_groups = length(groups);

fprintf('State: %s\n', state);
session_name_full = sprintf('MuscimolPost%sCortexAlignLast', state);

%% default parameters
p = struct();
p.N_thalamus_per_area = 10;
p.Thalamus_h = -4.5;
p.TC_density = 0.1;
p.within_E = 2.7;
p.across_E = 1;
p.TC_E = 3;
p.TC_I = 3;
p.within_EI = 7.5;
p.within_IE = -5;
p.across_EI = 0.5;
p.across_IE = 0;
p.adjacent_distance = 4;
p.sync_wave_T = 100; % sync wave period
p.parameter_tag = parameters.parameter_tag;

override_fields = fieldnames(parameters);
for override_idx = 1:numel(override_fields)
    field = override_fields{override_idx};
    p.(field) = parameters.(field);
end

%% Load post-Muscimol experimental data
raster_folder = fullfile(root, 'Data', 'Working', 'raster');
raster_file = fullfile(raster_folder, sprintf('raster_%s_%d.mat', session_name_full, session_idx));
load(raster_file, 'cell_area');

model_folder = fullfile(root, 'Data', 'Working', 'GLM_models');
model_filename = sprintf('GLM_%s_s%d_shuffle0_%s_%s_epoch%d_fold0.mat', session_name_full, session_idx, kernel_name, reg, 2500);
model_path = fullfile(model_folder, model_filename);
load(model_path, 'N', 'model_par', "kernel");
n_PS_kernel = kernel.n_PS_kernel;
n_cortex = N;

n_cells = zeros(1, n_groups); 
for i = 1:n_groups
    n_cells(i) = sum(ismember(cell_area, groups{i}));
end

% set thalamus cell number
thalamus_per_area = round(p.N_thalamus_per_area);
if thalamus_per_area < 0
    error('simulation_from_post:InvalidThalamusCount', 'N_thalamus_per_area must be non-negative');
end
n_cells(3) = thalamus_per_area;
n_cells(4) = thalamus_per_area;

n_cells_cumulative = [0, cumsum(n_cells)];
N = sum(n_cells);
cell_areas = zeros(1, N);
for i = 1:n_groups
    cell_areas((n_cells_cumulative(i)+1):n_cells_cumulative(i+1)) = i;
end

% load model parameters into J and h
J = zeros(N, N, 3);
h = zeros(N, 1);
cortex_mask = zeros(1, N); % mask for cortex neurons
thalamus_mask = zeros(1, N); % mask for thalamus neurons
cortex_mask(1:n_cortex) = 1; % cortex neurons
thalamus_mask(n_cortex+1:end) = 1; % thalamus neurons

if partial
    partial_mask = cortex_mask; % only keep cortex neurons
    keep_flags = [1, 1, 0, 0]; % keep ACC and VLPFC only
else
    partial_mask = ones(1, N); % keep all neurons
    keep_flags = [1, 1, 1, 1]; % keep all groups
end

for kernel_idx = 1:3
    J_mat = model_par(:, (n_cortex*(kernel_idx-1) + n_PS_kernel + 2):(n_cortex*kernel_idx + n_PS_kernel + 1));
    J_mat(isnan(J_mat)) = 0;
    J(1:n_cortex, 1:n_cortex, kernel_idx) = J_mat;
end
h(1:n_cortex) = model_par(:, 1);
h(isnan(h)) = -20;

%% Generate basic connectivity matrix
seed = 137 + session_idx;
rng(seed);

% generate random TC connections
for i = 1:N
    if thalamus_mask(i) == 1
        continue; % skip thalamus neurons
    end
    for j = 1:N
        if thalamus_mask(j) == 0
            continue; % skip cortex neurons
        end
        if cell_areas(j) == 3 && cell_areas(i) ~= 1
            continue; % thalamus to ACC only
        end
        if cell_areas(j) == 4 && cell_areas(i) ~= 2
            continue; % thalamus to VLPFC only
        end

        for k = 1:3
            if rand() < p.TC_density
                J(i, j, k) = p.TC_E*max((1+0.2*randn()), 0.1);
            end
        end

    end
end

for k = 1:3
    J(:, :, k) = J(:, :, k) - diag(diag(J(:, :, k))); % remove self connection
end

% generate random h
for i = 1:N
    if thalamus_mask(i) == 1
        h(i) = p.Thalamus_h * max((1+0.1*randn()), 0.1);
    end
end

% load kernels
kernel_file = fullfile(root, 'Data', 'Working', 'kernel', sprintf('kernel_%s.mat', kernel_name));
load(kernel_file, 'conn_kernels', 'kernel_len');
conn_kernels_reversed = cell(1, 3);
for k = 1:3
    conn_kernels_reversed{k} = conn_kernels{k}(:, end:-1:1);
end

N_full = N;
J_base = J;
h_base = h;

%% generate random raster for three input types: Async, Sync, NoInput
session_names = {'FittedAsync', 'FittedSync', 'FittedNoInput'};

if do_plot
    J_ctr = J; %#ok<NASGU>
end
all_firing_rates = cell(1, 3);

Sync_T = p.sync_wave_T; % sync wave period
all_rasters = cell(1, 3);

for session_type_idx = 1:3
    N = N_full;
    J = J_base;
    h = h_base;
    session_base = session_names{session_type_idx};
    session_name_full = sprintf('%s%sCortexAlignLast%s', session_base, state, p.parameter_tag);

    if strcmp(session_base, 'FittedNoInput')
        % remove thalamus connections
        J(:, thalamus_mask==1, :) = 0;
        if do_plot
            J_Mus = J; %#ok<NASGU>
        end
    elseif strcmp(session_base, 'FittedSync')
        % % uniform thalamocortical connections
        % sync_J = mean(J(21:100, 1:20, :), 2);
        % J(21:100, 1:20, :) = repmat(sync_J, 1, 20, 1);
        if do_plot
            J_Sync = J; %#ok<NASGU>
        end
    end

    rasters = cell(1, n_trials);
    trial_len = ones(1, n_trials)*B;
    firing_rates = zeros(N, n_trials);

    % GPU if available
    if UseGPU
        fprintf("Using GPU\n");
        J = gpuArray(J);
        h = gpuArray(h);
        conn_kernels_reversed = cellfun(@gpuArray, conn_kernels_reversed, 'UniformOutput', false);
    end

    fprintf("%s  start...\n", session_name_full);

    % start parallel pool if not already started
    if isempty(gcp)
        parpool;
    end

    % progress bar
    fprintf("Progress: \n[");
    for i = 1:n_trials
        fprintf(".");
    end
    fprintf("]\n[");

    parfor trial_idx = 1:n_trials
        raster = zeros(N, B);

        if UseGPU
            raster = gpuArray(raster);
        end

        for t = 1:B
            htot = h;

            if strcmp(session_base, 'FittedSync')
                % sync wave
                if mod(t, Sync_T) < Sync_T/8
                    htot = htot + 1.5*thalamus_mask';
                else
                    htot = htot - 0.6*thalamus_mask';
                end
            end

            if t>=kernel_len
                for k = 1:3
                    kernel_vec = conn_kernels_reversed{k};
                    % raster_conv = raster(:, t-kernel_len+1:t);
                    % convolved = sum(raster_conv .* kernel, 2);
                    convolved = sum(raster(:, t-kernel_len+1:t) .* kernel_vec, 2);
                    h_k = J(:, :, k)*convolved;
                    htot = htot + h_k;
                end
            end
            lambda = exp(htot);
            p_i = lambda./(1+lambda);
            raster(:, t) = rand(N, 1) < p_i;
        end

        if UseGPU
            raster = gather(raster);
        end

        rasters{trial_idx} = raster;
        firing_rates(:, trial_idx) = sum(raster, 2)/B;
        fprintf("#");
    end

    fprintf("]\n%s end\n", session_name_full);

    all_rasters{session_type_idx} = rasters;
    all_firing_rates{session_type_idx} = firing_rates;

    if partial
        % remove thalamus and inhibitory neurons
        for i=1:n_trials
            rasters{i}=rasters{i}(partial_mask==1, :);
        end
        firing_rates = firing_rates(partial_mask==1, :);
        N = sum(partial_mask);
    end

    trial_num = n_trials;
    save_folder = fullfile(root, 'Data', 'Working', 'raster');
    check_path(save_folder);
    save_filename = sprintf('raster_%s_%d.mat', session_name_full, session_idx);
    save_path = fullfile(save_folder, save_filename);
    save(save_path, 'N', 'rasters', 'J', 'h', "trial_num", "session_name_full", "trial_len", "firing_rates");
    if partial
        kept_counts = n_cells(logical(keep_flags));
        kept_cumulative = cumsum(kept_counts);
        % starting index of each area
        borders = [1, kept_cumulative(1:end-1) + 1];
    else
        borders = n_cells_cumulative(2:end-1)+0.5;
    end
    border_folder = fullfile(root, 'Data', 'Working', 'border');
    check_path(border_folder);
    border_filename = sprintf('borders_%s_%d.mat', session_name_full, session_idx);
    border_path = fullfile(border_folder, border_filename);
    save(border_path, "borders");
end

%% Plot
% figure('Position', [100, 100, 800, 800]);
% hold on;
% for i = 1:3
%     plot(mean(all_firing_rates{i}, 2), 'o-');
% end
% hold off;
% title('Firing rates');
% xlabel('Neuron ID');
% ylabel('Firing rate (Hz)');
% legend(session_names, 'Location', 'Best');

% %% Plotting
% if do_plot
%     N = sum(n_cells);
%     model_names = {'Ctr', 'Sync', 'Mus'};
%     % plot J matrix
%     cmap = brewermap(256,'*RdBu');
%     all_J = {J_ctr, J_Sync, J_Mus};
%     figure('Position', [100, 100, 800, 800]);
%     tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
%     for i = 1:3
%         for k = 1:3
%             model_name = model_names{i};
%             nexttile;
%             hold on;
%             imagesc(squeeze(all_J{i}(:, :, k)));
%             % title(['J', num2str(k)]);
%             colorbar;
%             colormap(cmap);
%             clim([-7.5, 7.5])
%             axis square;
%             set(gca, 'XTick', [], 'YTick', []);
%             title([model_name, ', Kernel ', num2str(k)]);

%             % area borders
%             for j = 1:(n_groups+1)
%                 plot([n_cells_cumulative(j), n_cells_cumulative(j)]+0.5, [0, N]+0.5, 'k', 'LineWidth', 1);
%                 plot([0, N]+0.5, [n_cells_cumulative(j), n_cells_cumulative(j)]+0.5, 'k', 'LineWidth', 1);
%             end

%             % reverse y axis
%             set(gca, 'YDir', 'reverse');
%             xlim([0.5, N+0.5]);
%             ylim([0.5, N+0.5]);

%             hold off;
%         end
%     end

%     % plot connection counts
%     areas = zeros(n_groups, 2);
    
%     all_J = {J_ctr, J_Sync, J_Mus};
%     figure('Position', [100, 100, 800, 800]);
%     tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
%     for session_i = 1:3
%         for k = 1:3
%             nexttile;
%             hold on;
%             J = squeeze(all_J{session_i}(:, :, k));

%             % count connections
%             J_within1 = J((areas(1, 1)+1):areas(1, 2), (areas(1, 1)+1):areas(1, 2));
%             J_within2 = J((areas(2, 1)+1):areas(2, 2), (areas(2, 1)+1):areas(2, 2));
%             J_across1 = J((areas(1, 1)+1):areas(1, 2), (areas(2, 1)+1):areas(2, 2));
%             J_across2 = J((areas(2, 1)+1):areas(2, 2), (areas(1, 1)+1):areas(1, 2));
%             err = 1e-5;

%             within_pos = nnz(J_within1(:) > err)+nnz(J_within2(:) > err);
%             within_neg = nnz(J_within1(:) < -err)+nnz(J_within2(:) < -err);
%             across_pos = nnz(J_across1(:) > err)+nnz(J_across2(:) > err);
%             across_neg = nnz(J_across1(:) < -err)+nnz(J_across2(:) < -err);

%             J_counts = [within_pos, within_neg; across_pos, across_neg; within_pos+across_pos, within_neg+across_neg];

%             b = bar([1, 2], J_counts(1:2, :), 'grouped', 'BarWidth', 0.8,'FaceColor', 'flat');
%             % set colors for bars: red for positive, blue for negative
%             b(1).CData(1, :) = [0.9, 0, 0]; % red for positive
%             b(2).CData(1, :) = [0, 0, 0.9]; % blue for negative
%             b(1).CData(2, :) = [0.9, 0, 0]; % red for positive
%             b(2).CData(2, :) = [0, 0, 0.9]; % blue for negative

%             title([model_name, ', Kernel ', num2str(k)]);
%             xticks([1, 2]);
%             xticklabels({'Within', 'Across', 'All'});
%             ylabel('Count');
%             legend({'Pos', 'Neg'});
%             ylim([0, 1200]);
%             % text for bars
%             for i = 1:2
%                 for j = 1:2
%                     text(i+(j-1.5)*0.3, J_counts(i, j)+20, num2str(J_counts(i, j)),...
%                         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8);
%                 end
%             end

%             hold off;
%         end
%     end

%     % Plot firing rates
%     figure('Position', [100, 100, 800, 800]);
%     hold on;
%     for i = 1:3
%         plot(mean(all_firing_rates{i}, 2), 'o-');
%     end
%     hold off;
%     title('Firing rates');
%     xlabel('Neuron ID');
%     ylabel('Firing rate (Hz)');
%     legend(session_names, 'Location', 'Best');

%     %% Plot raster
%     names = {'Control (Task)', 'Sync (EyeClose)', 'Muscimol'};
%     figure('Position', [100, 100, 800, 800]);
%     tiledlayout(3, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
%     for i = 1:3
%         rasters = all_rasters{i};
%         nexttile;
%         % plot thalamus raster in first 10 trials, first 10 neurons
%         concatenated_raster = zeros(10, 1000);
%         for j = 1:1
%             concatenated_raster(:, :) = rasters{j}(1:10, 1:1000);
%         end
%         if i==3
%             concatenated_raster = concatenated_raster.*0;
%         end
%         imagesc(concatenated_raster, [0, 1]);
%         title(names{i});
%         xlim();
%     end
% end
%% Finished
% fprintf('Parameter: %s, Strength level %.1f all sessions done.\n', parameter_name, strength_levels(strength_idx));

end