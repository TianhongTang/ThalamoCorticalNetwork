%% Plot functional connectivity from trained model
% Matrix view and circular graph view

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

%% Main

% register tasks
controls = {'Muscimol', 'Saline'};
% controls = {'Muscimol'};
session_idxs_all = {1:10, 1:5}; % session indices for each control
% area_types = {'Full', 'Cortex'};
area_types = {'Cortex'};
prepost_types = {'Pre', 'Post'};
states = {'RestOpen', 'RestClose', 'Task'};
alignments = {'Last'};

kernel = 'DeltaPure';
reg = struct();
reg.l1=0;
reg.l2=0.2;
reg.name='L2=0_2';

tasks = {};
for control_idx = 1:length(controls)
    control = controls{control_idx};
    sessions = session_idxs_all{control_idx};
    for session_idx = sessions
        for area_idx = 1:length(area_types)
            area_type = area_types{area_idx};
            prepost_states = {};
            for prepost_idx = 1:length(prepost_types)
                prepost = prepost_types{prepost_idx};
                if strcmp(area_type, 'Full') && strcmp(prepost, 'Post')
                    % All thalamus data in post sessions are not available
                    continue;
                end
                for align_idx = 1:length(alignments)
                    alignment = alignments{align_idx};
                    for state_idx = 1:length(states)
                        state = states{state_idx};

                        task = struct();
                        task.control = control;
                        task.area_type = area_type;
                        task.prepost = prepost;
                        task.state = state;
                        task.alignment = alignment;
                        task.name = [control, prepost, state, area_type, 'Align', alignment];
                        task.session_idx = session_idx;
                        task.kernel = kernel;
                        task.reg = reg;
                        task.shuffle_size = 0; % number of shuffles
                        task.epoch = 2500;
                        tasks{end+1} = task;
                    end
                end
            end
        end
    end
end

% load kernel info
folder_name = fullfile(root, 'Data', 'Working', 'kernel');
file_name = sprintf('kernel_%s.mat', kernel);
file_path = fullfile(folder_name, file_name);
load(file_path, 'kernel_len', 'n_PS_kernel', 'n_conn_kernel');

cmap=brewermap(256,'*RdBu'); % colormap: blue-white-red

% run tasks: for each task = 1 model, plot matrix and circular graph
for task_idx = 1:length(tasks)
    task = tasks{task_idx};
    % load model data
    folder_name = fullfile(root, 'Data', 'Working', 'GLM_models');
    file_name = sprintf('GLM_%s_s%d_shuffle%d_%s_%s_epoch%d_fold0.mat', task.name, task.session_idx, 0, task.kernel, ...
        task.reg.name, task.epoch);
    file_path = fullfile(folder_name, file_name);
    load(file_path, 'model_par', 'model_err', 'N');

    % load border info, get area num. borders: starting index of each area
    folder_name = fullfile(root, 'Data', 'Working', 'border');
    file_name = sprintf('borders_%sPreCortex_%d.mat', task.control, task.session_idx);
    file_path = fullfile(folder_name, file_name);
    load(file_path, 'borders');
    area_num = length(borders);
    borders_raw = borders;

    % plot: 2 row, 3 col. Row 1: weight matrix, Row 2: circular graph. 3 kernels.
    f = figure('Position', [100, 100, 1200, 1200], 'Visible', 'off');
    t = tiledlayout(3, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for kernel_idx = 1:n_conn_kernel
        % extract full J matrix
        J_mat = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
        J_err = model_err.total(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));

        % weight matrix plot
        nexttile(kernel_idx);
        hold on;
        imagesc(J_mat);
        % imagesc(J_mat, 'AlphaData', ~isnan(J_mat));
        % set(gca, 'Color', [0, 0, 0]);
        colormap(cmap);
        clim([-5, 5]);
        colorbar;
        % border lines for areas
        for b = 2:length(borders_raw)
            border_pos = borders_raw(b)-0.5;
            plot([border_pos, border_pos], [0.5, N+0.5], 'k-', 'LineWidth', 1);
            plot([0.5, N+0.5], [border_pos, border_pos], 'k-', 'LineWidth', 1);
        end
        set(gca, 'YDir', 'reverse');
        xlim([0.5, N+0.5]);
        ylim([0.5, N+0.5]);
        axis square;
        box on;
        hold off;

        title(sprintf('Kernel %d', kernel_idx));
        xlabel('From Neuron');
        ylabel('To Neuron');

        
        % weight matrix plot, set non-significant to 0
        nexttile(kernel_idx + 3);
        hold on;
        J_mat_sig = J_mat;;
        nonsig_idx = abs(J_mat) <= J_err;
        J_mat_sig(nonsig_idx) = 0;
        imagesc(J_mat_sig);
        % imagesc(J_mat, 'AlphaData', ~isnan(J_mat));
        % set(gca, 'Color', [0, 0, 0]);
        colormap(cmap);
        clim([-5, 5]);
        colorbar;
        % border lines for areas
        for b = 2:length(borders_raw)
            border_pos = borders_raw(b)-0.5;
            plot([border_pos, border_pos], [0.5, N+0.5], 'k-', 'LineWidth', 1);
            plot([0.5, N+0.5], [border_pos, border_pos], 'k-', 'LineWidth', 1);
        end
        set(gca, 'YDir', 'reverse');
        xlim([0.5, N+0.5]);
        ylim([0.5, N+0.5]);
        axis square;
        box on;
        hold off;

        title(sprintf('Significant connections, Kernel %d', kernel_idx));
        xlabel('From Neuron');
        ylabel('To Neuron');

        % circular graph plot
        nexttile(kernel_idx + 6);
        % prepare data for circular graph
        % nodes: neurons
        % edges: significant connections
        ax = gca;
        node_colors = zeros(N, 3);
        area_names = {'ACC', 'VLPFC'};
        mode = 'across'; % 'within', 'across' or 'all'
        area_offset = 0.5;
        fill_ratio = 0.7;

    network_plot_hemi(ax, J_mat, J_err, [borders, N+1], node_colors, area_names, mode, area_offset, fill_ratio);
        title(sprintf('Network Graph, Kernel %d', kernel_idx));
    end

    % global title
    sgtitle(sprintf('Functional Connectivity: %s, Session %d', task.name, task.session_idx), 'FontSize', 16);

    % save figure
    folder_name = fullfile(root, 'Figures', 'Functional_Connectivity');
    check_path(folder_name);
    figure_name = sprintf('FC_%s_%d.png', task.name, task.session_idx);
    figure_path = fullfile(folder_name, figure_name);
    saveas(f, figure_path);
    set(f, 'PaperUnits', 'inches');
    set(f, 'PaperPosition', [0 0 15 10]);
    set(f, 'PaperSize', [15 10]);
    pdf_name = sprintf('FC_%s_%d.pdf', task.name, task.session_idx);
    pdf_path = fullfile(folder_name, pdf_name);
    print(f, '-painters', '-dpdf', pdf_path);

end