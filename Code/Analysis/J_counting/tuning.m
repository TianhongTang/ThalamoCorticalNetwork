%% tuning.m - J difference based on tuning

clear;
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
% model config
config.align                 = 'Last';
config.shuffle_idx           = 0;
config.kernel_name           = 'DeltaPure';
config.reg_name              = 'L2=0_2';
config.epoch                 = 3000;
config.fold_idx              = 0;
config.resting_dur_threshold = 15;

kernel_num = 3;

% define comparason groups
% tuning format: (N, pre/post, offer1/offer2/BeforeCho/Cho2Info/Inof2Rew/AfterRew, ER/IxU/Unc/Sub)
% define tuning filters
filter_pls = @(tuning) logical(sum(tuning(:, 1, 2:5, 2), 3)); % Info +, Pre, choice to info
filter_mns = @(tuning) ~logical(sum(tuning(:, 1, 2:5, 2), 3)); % Info -, Pre, choice to info

% filter groups
filters = {filter_pls, filter_mns};
filter_names = {'Info +', 'Info -'};
filter_num = numel(filter_names);

% Load metadata
metadata = load_meta(root);
tuning_metas = metadata.tuning;
tuning_num = numel(tuning_metas);

all_J_mats = cell(tuning_num, 1);
all_J_err_mats = cell(tuning_num, 1);

%% Separate stats for each session
for i = 1:tuning_num
    % load tuning
    meta = tuning_metas(i);
    file_folder = fullfile(root, 'Data', 'Working', 'tuning');
    file_name = meta.file_name;
    file_path = fullfile(file_folder, file_name);
    load(file_path, 'meta', 'data');

    N = meta.N;
    tuning_session = data.tuning;

    J_all = zeros(N, N, kernel_num, 2, 2); % (N, N, n_kernels, pre/post, state);
    J_err_all = zeros(N, N, kernel_num, 2, 2);
    prepost_str = {'Pre', 'Post'};
    state_str = {'RestOpen', 'RestClose'};
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
            model_file_name = generate_filename('GLM', model_meta);
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
    all_J_mats{i} = J_all;
    all_J_err_mats{i} = J_err_all;

    % % Plot bar graph for each state and each kernel
    % for state_idx = 1:numel(state_str)
    %     state = state_str{state_idx};
    %     for kernel_idx = 1:kernel_num

    %         f = figure('Position', [100, 100, 400*filter_num, 400*filter_num], 'Visible', 'off');
    %         t = tiledlayout(filter_num, filter_num, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    %         for i = 1:filter_num
    %             for j = 1:filter_num
    %                 nexttile;
    %                 filter_i = filters{i}(tuning_session);
    %                 filter_j = filters{j}(tuning_session);
    %                 J_ij = J_all(filter_i, filter_j, kernel_idx, :, state_idx);
    %                 J_err_ij = J_err_all(filter_i, filter_j, kernel_idx, :, state_idx);
    %                 pos_pre = sum(J_ij(:, :, 1, 1) > J_err_ij(:, :, 1, 1), 'all');
    %                 neg_pre = sum(J_ij(:, :, 1, 1) < -J_err_ij(:, :, 1, 1), 'all');
    %                 pos_post = sum(J_ij(:, :, 1, 2) > J_err_ij(:, :, 1, 2), 'all');
    %                 neg_post = sum(J_ij(:, :, 1, 2) < -J_err_ij(:, :, 1, 2), 'all');

    %                 % plot bar graph: x: pos/neg, color: pre/post
    %                 bar([1, 2]-0.15, [pos_pre, neg_pre], 0.3, 'r');
    %                 hold on;
    %                 bar([1, 2]+0.15, [pos_post, neg_post], 0.3, 'b');
    %                 hold off;

    %                 max_y = max([pos_pre, neg_pre, pos_post, neg_post]) * 1.2;
    %                 if max_y == 0
    %                     max_y = 1;
    %                 end
    %                 ylim([0, max_y]);
    %                 xticks([1, 2]);
    %                 xticklabels({'Positive', 'Negative'});
    %                 title(sprintf('%s to %s', filter_names{j}, filter_names{i}));
    %                 legend('Pre', 'Post');
    %             end
    %         end

    %         % Subtitle and save figure
    %         sgtitle(sprintf('%s, %s, Session %d, %s, Kernel %d', meta.animal_name, meta.injection, meta.session_idx, state, kernel_idx));
    %         save_folder = fullfile(root, 'Figures', 'Tuning');
    %         check_path(save_folder);
    %         file_name = sprintf('Tuning_J_%s_%s_S%d_%s_K%d.png', meta.animal_name, meta.injection, meta.session_idx, state, kernel_idx);
    %         save_path = fullfile(save_folder, file_name);
    %         saveas(f, save_path);
    %         close(f);
    %     end
    % end
end

%% Summarize all sessions
injections = {'Saline', 'Muscimol'};
inj_kws = {'Sal', 'Mus'};
tuning_table = struct2table(tuning_metas);

for injection_idx = 1:numel(injections)
    injection = injections{injection_idx};
    inj_kw = inj_kws{injection_idx};
    selected_sessions = tuning_table(strcmp(tuning_table.injection, injection), :);
    selected_J = all_J_mats(strcmp(tuning_table.injection, injection));
    selected_J_err = all_J_err_mats(strcmp(tuning_table.injection, injection));
    selected_num = height(selected_sessions);
    fprintf('Processing injection: %s, with %d sessions.\n', injection, selected_num);

    for state_idx = 1:numel(state_str)
        state = state_str{state_idx};
        for kernel_idx = 1:kernel_num
            f = figure('Position', [100, 100, 400*filter_num, 400*filter_num], 'Visible', 'off');
            t = tiledlayout(filter_num, filter_num, 'TileSpacing', 'Compact', 'Padding', 'Compact');
            max_y_global = 0;
            max_conn_all = zeros(filter_num, filter_num); % save for text labelling
            p_val_all = zeros(filter_num, filter_num, 2); % save for text labelling
            for i = 1:filter_num
                for j = 1:filter_num    
                    pos_pre_all = 0;
                    neg_pre_all = 0;
                    pos_post_all = 0;
                    neg_post_all = 0;
                    max_conn = 0;

                    for session_idx = 1:selected_num
                        fprintf(' - Session %d/%d\n', session_idx, selected_num);
                        J_all = selected_J{session_idx};
                        J_err_all = selected_J_err{session_idx};

                        % load tuning
                        session_meta = selected_sessions(session_idx, :);
                        session_meta = table2struct(session_meta);
                        file_folder = fullfile(root, 'Data', 'Working', 'tuning');
                        file_name = session_meta.file_name;
                        file_path = fullfile(file_folder, file_name);
                        fprintf('   - Loading tuning from %s\n', file_path);
                        load(file_path, 'meta', 'data');
                        tuning_session = data.tuning;
                        filter_i = filters{i}(tuning_session);
                        filter_j = filters{j}(tuning_session);

                        % load model
                        prepost_str = {'Pre', 'Post'};
                        for prepost_idx = 1:numel(prepost_str)
                            J_filtered = J_all(filter_i, filter_j, kernel_idx, prepost_idx, state_idx);
                            J_err_filtered = J_err_all(filter_i, filter_j, kernel_idx, prepost_idx, state_idx);
                            pos_count = sum(J_filtered(:) > J_err_filtered(:));
                            neg_count = sum(J_filtered(:) < -J_err_filtered(:));
                            if prepost_idx == 1
                                pos_pre_all = pos_pre_all + pos_count;
                                neg_pre_all = neg_pre_all + neg_count;
                            else
                                pos_post_all = pos_post_all + pos_count;
                                neg_post_all = neg_post_all + neg_count;
                            end
                            if i==j
                                max_conn = max_conn + sum(filter_i) * (sum(filter_j)-1); % exclude self-connection
                            else
                                max_conn = max_conn + sum(filter_i) * sum(filter_j);
                            end
                        end
                    end
                    max_conn_all(i, j) = max_conn;

                    % plot bar graph of conn ratio: x: pos/neg, color: pre/post
                    if max_conn > 0
                        nexttile;
                        hold on;
                        bar([1, 2]-0.15, [pos_pre_all, neg_pre_all]/max_conn, 0.3, 'r');
                        bar([1, 2]+0.15, [pos_post_all, neg_post_all]/max_conn, 0.3, 'b');

                        % Error bar: Wilson confidence interval
                        [pre_pos_low, pre_pos_high]   = wilsonCI(pos_pre_all, max_conn);
                        [pre_neg_low, pre_neg_high]   = wilsonCI(neg_pre_all, max_conn);
                        [post_pos_low, post_pos_high] = wilsonCI(pos_post_all, max_conn);
                        [post_neg_low, post_neg_high] = wilsonCI(neg_post_all, max_conn);
                        x_all    = [[1, 2]-0.15, [1, 2]+0.15];
                        y_all    = [pos_pre_all, neg_pre_all, pos_post_all, neg_post_all]/max_conn;
                        err_low  = [pre_pos_low, pre_neg_low, post_pos_low, post_neg_low];
                        err_high = [pre_pos_high, pre_neg_high, post_pos_high, post_neg_high];
                        errorbar(x_all, y_all, y_all - err_low, err_high - y_all, 'k', 'LineStyle', 'none');

                        % Statistical test: two-proportion z-test between pre and post
                        p_val_pos = twoProportionPValue(pos_pre_all, max_conn, pos_post_all, max_conn);
                        p_val_neg = twoProportionPValue(neg_pre_all, max_conn, neg_post_all, max_conn);
                        p_val_all(i, j, 1) = p_val_pos;
                        p_val_all(i, j, 2) = p_val_neg;
                    end

                    max_y = max([pos_pre_all, neg_pre_all, pos_post_all, neg_post_all]) * 1.2;
                    max_y_global = max(max_y_global, max_y);

                    xticks([1, 2]);
                    xticklabels({'Positive', 'Negative'});
                    ylabel('Significant J Ratio');
                    title(sprintf('%s to %s', filter_names{j}, filter_names{i}));
                    legend('Pre', 'Post');
                end
            end

            % set same y limit for all subplots, and put text labels
            if max_y_global == 0
                max_y_global = 1;
            end

            max_y_global = 0.25; % for ratio plot
            
            for i = 1:filter_num
                for j = 1:filter_num
                    nexttile((i-1)*filter_num + j);
                    ylim([0, max_y_global]);
                    text(1.5, max_y_global * 0.9, sprintf('Total connection: %d', max_conn_all(i, j)), 'HorizontalAlignment', 'center');
                    for posneg_idx = 1:2
                        p = p_val_all(i, j, posneg_idx);
                        significant_text = 'N.S.';
                        if ~isnan(p) && p < 0.05
                            significant_text = '*';
                        end
                        if ~isnan(p) && p < 0.01
                            significant_text = '**';
                        end
                        if ~isnan(p) && p < 0.001
                            significant_text = '***';
                        end
                        text(posneg_idx, max_y_global * 0.75, sprintf('p = %.3f\n%s', p, significant_text), 'HorizontalAlignment', 'center');
                    end
                end
            end

            % save figure
            sgtitle(sprintf('%s, %s, Kernel %d', injection, state, kernel_idx));
            save_folder = fullfile(root, 'Figures', 'Tuning');
            check_path(save_folder);
            file_name = sprintf('Tuning_J_%s_%s_K%d.png', inj_kw, state, kernel_idx);
            save_path = fullfile(save_folder, file_name);
            saveas(f, save_path);
            close(f);

            % Relative post/pre ratio
        end
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