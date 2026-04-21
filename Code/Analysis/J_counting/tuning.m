%% tuning.m - J difference based on tuning

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
% model config
config.align                 = 'Last';
config.shuffle_idx           = 0;
config.kernel_name           = 'DeltaPure';
config.reg_name              = 'L2=0_2';
config.epoch                 = 3000;
config.fold_idx              = 1;
config.resting_dur_threshold = 15;

kernel_num = 3;

% Load metadata
metadata = load_meta(root);
tuning_metas = metadata.tuning;
tuning_num = numel(tuning_metas);

for i = 1:tuning_num
    % load tuning
    meta = tuning_metas{i};
    file_folder = fullfile(root, 'Data', 'Working', 'tuning');
    file_name = meta.file_name;
    file_path = fullfile(file_folder, file_name);
    load(file_path, 'meta', 'data');

    N = meta.N;

    J_all = zeros(N, N, kernel_num, 2, 2); % (N, N, n_kernels, pre/post, state);
    J_err_all = zeros(N, N, kernel_num, 2, 2);
    prepost_str = {'pre', 'post'};
    state_str = {'RestOpen', 'RestClosed'};
    for prepost_idx = 1:numel(prepost_str)
        for state_idx = 1:numel(state_str)
            % construct model meta
            model_meta = meta;
            model_meta.prepost = prepost_str{prepost_idx};
            model_meta.state = state_str{state_idx};
            for field = fieldnames(config)'
                model_meta.(field{1}) = config.(field{1});
            end

            % get model filename
            model_file_name = generate_filename(model_meta, 'GLM');
            model_file_path = fullfile(root, 'Data', 'Working', 'GLM', model_file_name);

            % load model
            model_file = load(model_file_path);
            model_par = model_file.data.model_par;
            model_err = model_file.data.model_err.total;
            model_N = model_file.meta.N;
            n_PS_kernel = model_file.data.kernel.meta.n_PS_kernel;
            assert(model_N == N, 'Model N does not match meta N');

            % extract J matrix
            for kernel_idx = 1:kernel_num
                start_idx = n_PS_kernel + 2 + (kernel_idx-1)*N;
                end_idx = start_idx + N - 1;
                J_mat = model_par(:, start_idx:end_idx);
                J_err = model_err(:, start_idx:end_idx);
                J_all(:, :, kernel_idx, prepost_idx, state_idx) = J_mat;
                J_err_all(:, :, kernel_idx, prepost_idx, state_idx) = J_err;
            end
        end
    end

    % define comparason groups
    % tuning format: (N, pre/post, offer1/offer2/BeforeCho/Cho2Info/Inof2Rew/AfterRew, ER/IxU/Unc/Sub)
    % 1. choice to info +/-, IxU
    filter_pls = data.tuning(:, 1, 4, 2);
    filter_mns = data.tuning(:, 1, 4, 3);

    % filter groups
    filters = {filter_pls, filter_mns};
    filter_names = {'Info +', 'Info -'};
    filter_num = numel(filter_names);

    % Plot bar graph for each state and each kernel
    for state_idx = 1:numel(state_str)
        state = state_str{state_idx};
        for kernel_idx = 1:kernel_num

            f = figure('Position', [100, 100, 400*filter_num, 400*filter_num], 'Visible', 'off');
            t = tiledlayout(filter_num, filter_num, 'TileSpacing', 'Compact', 'Padding', 'Compact');
            for i = 1:filter_num
                for j = 1:filter_num
                    nexttile;
                    filter_i = filters{i};
                    filter_j = filters{j};
                    J_ij = J_all(filter_i, filter_j, kernel_idx, :, state_idx);
                    J_err_ij = J_err_all(filter_i, filter_j, kernel_idx, :, state_idx);
                    pos_pre = sum(J_ij(:, :, 1, 1) > J_err_ij(:, :, 1, 1), 'all');
                    neg_pre = sum(J_ij(:, :, 1, 1) < -J_err_ij(:, :, 1, 1), 'all');
                    pos_post = sum(J_ij(:, :, 1, 2) > J_err_ij(:, :, 1, 2), 'all');
                    neg_post = sum(J_ij(:, :, 1, 2) < -J_err_ij(:, :, 1, 2), 'all');

                    % plot bar graph: x: pos/neg, color: pre/post
                    bar([1, 2]-0.15, [pos_pre, neg_pre], 0.3, 'r');
                    hold on;
                    bar([1, 2]+0.15, [pos_post, neg_post], 0.3, 'b');
                    xticks([1, 2]);
                    xticklabels({'Positive', 'Negative'});
                    title(sprintf('%s to %s', filter_names{j}, filter_names{i}));
                    legend('Pre', 'Post');
                end
            end

            % Subtitle and save figure
            suptitle(sprintf('%s, %s, Session %d, %s, Kernel %d', meta.animal_name, meta.injection, meta.session_idx, state, kernel_idx));
            save_folder = fullfile(root, 'Figures', 'Tuning');
            check_path(save_folder);
            file_name = sprintf('Tuning_J_%s_%s_S%d_%s_K%d.png', meta.animal_name, meta.injection, meta.session_idx, state, kernel_idx);
            save_path = fullfile(save_folder, file_name);
            saveas(f, save_path);
            close(f);
        end
    end
end

            


