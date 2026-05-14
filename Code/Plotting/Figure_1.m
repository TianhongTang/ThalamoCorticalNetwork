clear;
%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
% include code folder and utils
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));
% addpath(fullfile(root, 'Code', 'Utils', 'HELPER_GENERAL'));

%% Main

% Example session: Slayer Mus 6.
meta = struct();
meta.animal_name = 'Slayer';
meta.injection = 'Muscimol';
meta.align = 'Last';
meta.session_idx = 6;
meta.resting_dur_threshold = 15;

% areas    = {'Full',     'Full',      'Cortex',   'Cortex'};
% areas    = {'Cortex',   'Cortex',    'Cortex',   'Cortex'};
% preposts = {'Pre',      'Pre',       'Post',     'Post'};
% states   = {'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};
areas    = {'Cortex',   'Cortex',  };
preposts = {'Pre',      'Pre',     };
states   = {'RestOpen', 'RestClose'};
n_state = numel(areas);

% parameters
shuffle_N = 20;
err_multi = 2; % threshold for significant J, in multiples of the error estimate from GLM.

% figure
f = figure();
tiles = tiledlayout(4, n_state, "TileSpacing", "Compact", "Padding", "Compact");


%% Row 1: Raster of selected neurons

% Example session: Slayer Mus 6. Selected t_range. Justification: Best sleep period for Pre-eyeclose.
% Remove this range in population analysis, or add proper filters to select sleep periods.

% selected_neurons = 1:95;
% selected_neurons = [16:21, 39, 44, 61:67, 83:93];
% selected_neurons = [14:19, 33:40, 42, 44, 46, 49, 58:73, 83:95];
% selected_neurons = [33, 61];
selected_neurons = [2, 42];
% selected_neurons = [34, 92];
display_t_range = 20001:30000;

for i = 1:length(areas)
    % Load data
    meta.area = areas{i};
    meta.prepost = preposts{i};
    meta.state = states{i};
    meta.filename = generate_filename('raster', meta);

    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    % Plot
    % N = raster_data.meta.N;

    cell_area = raster_data.data.cell_area;
    % tile = nexttile(i*2-1);
    tile = nexttile(i);
    raster = raster_data.data.rasters{1}(:, display_t_range);

    cell_area = cell_area(selected_neurons);
    raster = raster(selected_neurons, :);

    N = numel(selected_neurons);
    colors = zeros(N, 3);
    for j = 1:N
        switch cell_area{j}
            case 'Thalamus'
                colors(j, :) = [1, 0, 1];
            case 'ACC'
                colors(j, :) = [0, 0, 1]; % blue
            case 'VLPFC'
                colors(j, :) = [1, 0, 0]; % red
            otherwise
                colors(j, :) = [0, 0, 0]; % black
        end
    end
    cla(tile);
    raster_visualization_plot(tile, raster, colors)
    title(sprintf('Example pair, %s', meta.state));
    xlabel('Time (ms)');
    % ylabel('Neuron No.');
    ylabel('');
    yticks([1, 2]);
    yticklabels({'ACC #2', 'VLPFC #1'});
    ytickangle(0);
    ylim([-1.5, 4.5]);
    % ylim([-1.5, 4.5]);
end

%% Row 2: Pairwise correlogram between selected neurons
% f = figure();
% tiles = tiledlayout(2, 2, "TileSpacing", "Compact", "Padding", "Compact");

% selected_neurons = [85, 86];
% selected_neurons = [33, 61];
% selected_neurons = [2, 42];
% selected_neurons = [35, 60];
t_range = 1:60000;

smooth_window = 30; % ms

for i = 1:length(areas)
    % Load data
    meta.area = areas{i};
    meta.prepost = preposts{i};
    meta.state = states{i};
    meta.filename = generate_filename('raster', meta);

    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    raster = raster_data.data.rasters{1};

    r1 = raster(selected_neurons(1), t_range);
    r2 = raster(selected_neurons(2), t_range);

    % Compute correlogram
    [correlogram, lags] = norm_xcorr(r1, r2, 1000);
    smooth_correlogram = movmean(correlogram, smooth_window);

    % Compute shuffled correlogram for control
    shuffle_correlograms = zeros(shuffle_N, length(correlogram));
    for j = 1:shuffle_N
        fprintf('%d/%d shuffles started\n', j, shuffle_N);
        shuffled_r2 = r2(randperm(length(r2)));
        [shuffle_corr, ~] = norm_xcorr(r1, shuffled_r2, 1000);
        shuffle_correlograms(j, :) = movmean(shuffle_corr, smooth_window);
        fprintf('%d/%d shuffles finished\n', j, shuffle_N);
    end
    shuffle_mean = mean(shuffle_correlograms, 1);
    shuffle_std = std(shuffle_correlograms, [], 1);
    % shuffle_mean = movmean(shuffle_mean, smooth_window);
    % shuffle_std = movmean(shuffle_std, smooth_window);

    % Plot
    nexttile(i+n_state);
    std_multiplier = 2;
    shuffle_upper = shuffle_mean + std_multiplier * shuffle_std;
    shuffle_lower = shuffle_mean - std_multiplier * shuffle_std;
    % plot shuffled correlogram as a shaded area
    fill([lags, fliplr(lags)], [shuffle_upper, fliplr(shuffle_lower)], [0.8, 0.8, 0.8], 'EdgeColor', 'none');

    % plot correlogram
    hold on;
    xline(0, 'k--');
    yline(1, 'k--');
    plot(lags, smooth_correlogram);
    hold off;

    title(sprintf('Correlogram: %s', meta.state));
    xlabel('Lag (ms)');
    ylabel('Normalized correlation');
    legend({'Shuffled 2SD', '', '', 'correlation'});
    ylim([-1, 5]);
end

%% Row 3: network plot
for i = 1:length(areas)
    % Load data
    meta.area = areas{i};
    meta.prepost = preposts{i};
    meta.state = states{i};
    meta.shuffle_idx = 0;
    meta.kernel_name = 'DeltaPure';
    meta.reg_name = 'L2=0_2';
    meta.epoch = 3000;
    meta.fold_idx = 0;
    
    % Load raster data for neuron area info
    meta.filename = generate_filename('raster', meta);
    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    % split neurons by area
    cell_area = raster_data.data.cell_area;
    filter1 = ismember(cell_area, {'ACC'}); % filter for ACC neurons
    filter2 = ismember(cell_area, {'VLPFC'}); % filter for VLPFC neurons

    % Load GLM data for connectivity info
    meta.filename = generate_filename('GLM', meta);
    GLM_data = load(fullfile(root, 'Data', 'Working', 'GLM', meta.filename));
    fprintf('Loaded GLM data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    N = GLM_data.meta.N;
    J = GLM_data.data.model_par(:, (2:N+1)); % kernel 1 weights
    err = GLM_data.data.model_err.total(:, (2:N+1));

    % Plot connectivity
    tile = nexttile(i+n_state*2);
    plot_network(tile, J(filter1, filter2), J(filter2, filter1), ...
        err(filter1, filter2), err(filter2, filter1), err_multi, ...
        selected_neurons(1), selected_neurons(2)-41);

    title(tile, sprintf('Network plot: %s', meta.state));
end

% %% Row 4: average CCG across all neuron pairs
% for i = 1:length(areas)
%     % Load data
%     meta.area = areas{i};
%     meta.prepost = preposts{i};
%     meta.state = states{i};
%     meta.filename = generate_filename('raster', meta);
% 
%     raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
%     fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
%     fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
%     fprintf('Trial_num: %d\n', raster_data.meta.trial_num);
% 
%     % split neurons by area
%     cell_area = raster_data.data.cell_area;
%     filter1 = ismember(cell_area, {'ACC'}); % filter for ACC neurons
%     filter2 = ismember(cell_area, {'VLPFC'}); % filter for VLPFC neurons
% 
%     raster = raster_data.data.rasters{1};
% 
%     r1 = raster(filter1, 1:60000);
%     r2 = raster(filter2, 1:60000);
%     N1 = sum(filter1);
%     N2 = sum(filter2);
% 
%     % Compute correlogram for each pair and average
%     all_correlograms = zeros(N1, N2, 2001);
%     for n1 = 1:N1
%         for n2 = 1:N2
%             fprintf('Computing correlogram for pair (%d, %d)\n', n1, n2);
%             [correlogram, lags] = norm_xcorr(r1(n1, :), r2(n2, :), 1000);
%             smooth_correlogram = movmean(correlogram, smooth_window);
%             all_correlograms(n1, n2, :) = smooth_correlogram;
%         end
%     end
% 
%     % Compute mean and sem across all pairs
%     mean_correlogram = squeeze(mean(all_correlograms, [1, 2]));
%     sem_correlogram = squeeze(std(all_correlograms, [], [1, 2])) / sqrt(N1 * N2);
% 
%     % Plot
%     nexttile(i+12);
%     se_multiplier = 2;
%     error_upper = mean_correlogram + se_multiplier * sem_correlogram;
%     error_lower = mean_correlogram - se_multiplier * sem_correlogram;
%     fill([lags, fliplr(lags)], [error_upper', fliplr(error_lower')], [0.8, 0.8, 0.8], 'EdgeColor', 'none');
% 
%     % plot correlogram
%     hold on;
%     xline(0, 'k--');
%     yline(1, 'k--');
%     plot(lags, mean_correlogram);
%     hold off;
% 
%     title(sprintf('Mean correlogram: %s-%s, %s', meta.prepost, meta.injection, meta.state));
%     xlabel('Lag (ms)');
%     ylabel('Normalized correlation');
%     legend({'Mean ± 2SEM', '', '', 'Mean correlation'});
%     ylim([0, 2]);
% end

%% Row 4: Population firing rate CCG
for i = 1:length(areas)
    % Load data
    meta.area = areas{i};
    meta.prepost = preposts{i};
    meta.state = states{i};
    meta.filename = generate_filename('raster', meta);

    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    % split neurons by area
    cell_area = raster_data.data.cell_area;
    filter1 = ismember(cell_area, {'ACC'}); % filter for ACC neurons
    filter2 = ismember(cell_area, {'VLPFC'}); % filter for VLPFC neurons

    raster = raster_data.data.rasters{1};

    r1 = raster(filter1, 1:60000);
    r2 = raster(filter2, 1:60000);
    N1 = sum(filter1);
    N2 = sum(filter2);

    % Compute mean firing rate for each area
    pop_r1 = mean(r1, 1);
    pop_r2 = mean(r2, 1);

    % Compute correlogram
    [correlogram, lags] = norm_xcorr(pop_r1, pop_r2, 1000);
    smooth_correlogram = movmean(correlogram, smooth_window);

    % Compute shuffled correlogram for control
    shuffle_correlograms = zeros(shuffle_N, length(correlogram));
    for j = 1:shuffle_N
        fprintf('%d/%d shuffles started\n', j, shuffle_N);
        % shuffled_r2 = pop_r2(randperm(length(pop_r2)));
        shuffled_r2 = circshift(pop_r2, randi(length(pop_r2))); 
        [shuffle_corr, ~] = norm_xcorr(pop_r1, shuffled_r2, 1000);
        shuffle_correlograms(j, :) = movmean(shuffle_corr, smooth_window);
        fprintf('%d/%d shuffles finished\n', j, shuffle_N);
    end
    shuffle_mean = mean(shuffle_correlograms, 1);
    shuffle_std = std(shuffle_correlograms, [], 1);
    % shuffle_mean = movmean(shuffle_mean, smooth_window);
    % shuffle_std = movmean(shuffle_std, smooth_window);

    % Plot
    nexttile(i+n_state*3);
    std_multiplier = 2;
    shuffle_upper = shuffle_mean + std_multiplier * shuffle_std;
    shuffle_lower = shuffle_mean - std_multiplier * shuffle_std;
    % plot shuffled correlogram as a shaded area
    fill([lags, fliplr(lags)], [shuffle_upper, fliplr(shuffle_lower)], [0.8, 0.8, 0.8], 'EdgeColor', 'none');

    % plot correlogram
    hold on;
    xline(0, 'k--');
    yline(1, 'k--');
    plot(lags, smooth_correlogram);
    hold off;

    title(sprintf('Population correlogram: %s', meta.state));
    xlabel('Lag (ms)');
    ylabel('Normalized correlation');
    legend({'Shuffled 2SD', '', '', 'correlation'});
    ylim([0, 2]);
end


%% Export to pdf
% Save current figure as vector PDF
fig = gcf;

% ----- Figure settings -----
save_folder = fullfile(root, 'Figures', 'Paper');
check_path(save_folder);
filename = fullfile(save_folder, 'Figure1.pdf');

figWidth  = 8.0;   % inches
figHeight = 12.0;   % inches

resolution = 300;  % dpi; mainly affects rasterized components
% -------------------------

% Set figure size on screen
set(fig, 'Units', 'inches');
fig.Position(3:4) = [figWidth, figHeight];

% Set paper size for PDF export
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [figWidth, figHeight]);
set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);

% Make background white
set(fig, 'Color', 'w');

% Export as vector PDF
exportgraphics(fig, filename, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);

%%
close(fig);


% %% Fig 2: J counting bar plot
% % counting
% pos_counts = zeros(1, length(areas));
% neg_counts = zeros(1, length(areas));
% max_counts = zeros(1, length(areas));
% for i = 1:length(areas)
%     % Load data
%     meta.area = areas{i};
%     meta.prepost = preposts{i};
%     meta.state = states{i};
%     meta.shuffle_idx = 0;
%     meta.kernel_name = 'DeltaPure';
%     meta.reg_name = 'L2=0_2';
%     meta.epoch = 3000;
%     meta.fold_idx = 0;

%     % Load raster data for neuron area info
%     meta.filename = generate_filename('raster', meta);
%     raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
%     fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
%     fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
%     fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

%     % split neurons by area
%     cell_area = raster_data.data.cell_area;
%     filter1 = ismember(cell_area, {'ACC'}); % filter for ACC neurons
%     filter2 = ismember(cell_area, {'VLPFC'}); % filter for VLPFC neurons

%     % Load GLM data for connectivity info
%     meta.filename = generate_filename('GLM', meta);
%     GLM_data = load(fullfile(root, 'Data', 'Working', 'GLM', meta.filename));
%     fprintf('Loaded GLM data for %s %s %s\n', meta.prepost, meta.state, meta.area);
%     N = GLM_data.meta.N;
%     J = GLM_data.data.model_par(:, (2:N+1)); % kernel 1 weights
%     err = GLM_data.data.model_err.total(:, (2:N+1));

%     % Count significant J
%     pos_count = sum(J(filter1, filter2) > err_multi*err(filter1, filter2), 'all') + ...
%     sum(J(filter2, filter1) > err_multi*err(filter2, filter1), 'all');
%     neg_count = sum(J(filter1, filter2) < -err_multi*err(filter1, filter2), 'all') + ...
%     sum(J(filter2, filter1) < -err_multi*err(filter2, filter1), 'all');
%     max_count = numel(J(filter1, filter2)) + numel(J(filter2, filter1));
%     pos_counts(i) = pos_count;
%     neg_counts(i) = neg_count;
%     max_counts(i) = max_count;
% end

% % Plot
% f = figure();
% % x axis: states. Groups: positive vs negative connections.
% pos_ratios = pos_counts ./ max_counts;
% neg_ratios = neg_counts ./ max_counts;




%% functions
function [correlogram, lags] = norm_xcorr(r1, r2, max_lag)
    % Compute normalized cross-correlogram, not only by the number of overlapping points, but also by the overall firing rates of the two neurons.
    % Compute raw cross-correlogram
    [raw_correlogram, lags] = xcorr(r1, r2, max_lag, 'unbiased');
    ave_rate = zeros(size(lags));
    for i = 1:length(lags)
        lag = lags(i);
        if lag < 0
            ave_rate(i) = mean(r1(1:end+lag)) * mean(r2(1-lag:end));
        else
            ave_rate(i) = mean(r1(1+lag:end)) * mean(r2(1:end-lag));
        end
    end
    valid_filter = ave_rate > 0;
    correlogram = zeros(size(raw_correlogram));
    correlogram(valid_filter) = raw_correlogram(valid_filter) ./ ave_rate(valid_filter);

    if any(~valid_filter)
        correlogram(~valid_filter) = 0;
        warning('Some lags have zero average firing rate product, resulting in 0 in the normalized correlogram.');
    end
end

function [p_low, p_high] = wilsonCI(M, N, alpha)
    % Wilson score interval
    % M = successes, N = trials
    % alpha = significance level (e.g., 0.05 for 95% CI)

    if nargin < 3
        alpha = 0.05;
    end

    p = M / N;
    z = norminv(1 - alpha/2); % Z for CI

    denominator = 1 + (z^2)/N;
    center = p + (z^2)/(2*N);
    radius = z * sqrt( (p*(1-p)/N) + (z^2)/(4*N^2) );

    p_low = (center - radius) / denominator;
    p_high = (center + radius) / denominator;
end

function p_val = twoProportionPValue(success1, trials1, success2, trials2)
    % Two-proportion z-test (two-tailed) between pre and post proportions

    if trials1 == 0 || trials2 == 0
        p_val = NaN;
        return;
    end

    p1 = success1 / trials1;
    p2 = success2 / trials2;
    pooled = (success1 + success2) / (trials1 + trials2);
    denom = sqrt(pooled * (1 - pooled) * (1 / trials1 + 1 / trials2));
    if denom == 0
        p_val = 1;
        return;
    end

    z = (p1 - p2) / denom;
    abs_z = abs(z);
    Phi = 0.5 * erfc(-abs_z / sqrt(2));
    p_val = 2 * (1 - Phi);
end
