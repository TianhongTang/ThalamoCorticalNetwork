%% plot_hist_new.m - Plot histogram of J values for different kernels and distances.

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
resting_dur_threshold = 15;
mt = load_meta(root, 'table'); % metadata table

animals = {'Slayer', 'Emperor', 'Both'};
injection_types = {'Muscimol', 'Saline'};
states = {'RestOpen', 'RestClose'};
state_labels = {'Eyes Open', 'Eyes Closed'};
posneg_labels = {'Positive J', 'Negative J'};
prepost_labels = {'Pre', 'Post'};
alignment = 'Last';

%% Fig1: pre vs post. Fig2: open vs closed
for kernel_idx = 1:3
    for injection_idx = 1:2
        % fig 1: pre vs post
        f1 = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
        t1 = tiledlayout(f1, 2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
        injection = injection_types{injection_idx};
        for animal_idx = 1:3
            animal_name = animals{animal_idx};
            for state_idx = 1:2
                state = states{state_idx};
                state_label = state_labels{state_idx};

                J_counts = zeros(2, 2); % pos/neg x pre/post
                total_counts = zeros(2, 2); % total count for pos/neg in pre/post (should be the same)

                for prepost_idx = 1:2
                    prepost_str = prepost_labels{prepost_idx};

                    % filter metadata for current condition
                    % TODO: better condition filtering by a function
                    if strcmp(animal_name, 'Both')
                        animal_filter = strcmp(mt.GLM.animal_name, 'Slayer') | strcmp(mt.GLM.animal_name, 'Emperor');
                    else
                        animal_filter = strcmp(mt.GLM.animal_name, animal_name);
                    end
                    condition_filter = (animal_filter) & ...
                                        strcmp(mt.GLM.state, state) & ...
                                        strcmp(mt.GLM.kernel_name, "DeltaPure") & ...
                                        strcmp(mt.GLM.align, alignment) & ...
                                        strcmp(mt.GLM.area, "Cortex") & ...
                                        strcmp(mt.GLM.injection, injection) & ...
                                        (mt.GLM.epoch == 3000) & ...
                                        (mt.GLM.fold_idx == 0) & ...
                                        (mt.GLM.shuffle_idx == 0) & ...
                                        cellfun(@(x) ~isempty(x) && x == resting_dur_threshold, mt.GLM.resting_dur_threshold) & ...
                                        strcmp(mt.GLM.prepost, prepost_str);
                                        
                    condition_meta = mt.GLM(condition_filter, :);
                    if isempty(condition_meta)
                        warning('No data found for condition: %s, %s, %s', animal_name, state, prepost_str);
                        continue;
                    else
                        fprintf('Plotting condition: %s, %s, %s %s with %d sessions.\n', animal_name, state, prepost_str, injection, height(condition_meta));
                    end

                    % count significant J 
                    pos_count   = 0;
                    neg_count   = 0;
                    total_count = 0;
                    for session_idx = 1:height(condition_meta)
                        meta = condition_meta(session_idx, :);
                        meta = table2struct(meta);

                        % load model parameters and errors
                        data_folder = fullfile(root, 'Data', 'Working', 'GLM');
                        file_name = meta.file_name;
                        file_path = fullfile(data_folder, file_name);
                        if ~isfile(file_path)
                            warning('File not found: %s. Skipping.', file_path);
                            continue;
                        end
                        session_data = load(file_path, 'meta', 'data');
                        model_par = session_data.data.model_par;
                        model_err = session_data.data.model_err.total;

                        % load area borders
                        border_folder = fullfile(root, 'Data', 'Working', 'border');
                        border_file_name = generate_filename('border', meta);
                        border_file_path = fullfile(border_folder, border_file_name);
                        if ~isfile(border_file_path)
                            warning('Border file not found: %s. Skipping.', border_file_path);
                            continue;
                        end
                        border_data = load(border_file_path, 'data', 'meta');
                        borders = border_data.data.borders;
                        N = border_data.meta.N;
                        assert(numel(borders) == 2);
                        borders = [borders, N+1]; % add end border

                        selected_areas = {[1, 2], [2, 1]}; % between ACC and VLPFC

                        % count significant J values
                        [pos, neg, total] = J_count(model_par, model_err, kernel_idx, borders, selected_areas);
                        pos_count = pos_count + pos;
                        neg_count = neg_count + neg;
                        total_count = total_count + total;
                    end
                    J_counts(:, prepost_idx) = [pos_count; neg_count];
                    total_counts(:, prepost_idx) = [total_count; total_count];
                end

                % error bar and statistical test
                ratios = J_counts ./ total_counts;
                CIs = zeros(2, 2, 2); % pos/neg x pre/post x low/high
                p_vals = zeros(1, 2); % pos/neg
                sig_labels = cell(1, 2);
                for posneg_idx = 1:2
                    count = J_counts(posneg_idx, :);
                    total = total_counts(posneg_idx, :);
                    for prepost_idx = 1:2
                        [CI_low, CI_high] = wilsonCI(count(prepost_idx), total(prepost_idx), 0.05);
                        CIs(posneg_idx, prepost_idx, :) = [CI_low, CI_high];
                    end
                    p_val = twoProportionPValue(count(1), total(1), count(2), total(2));
                    p_vals(posneg_idx) = p_val;
                    if p_val < 0.001
                        sig_labels{posneg_idx} = '***';
                    elseif p_val < 0.01
                        sig_labels{posneg_idx} = '**';
                    elseif p_val < 0.05
                        sig_labels{posneg_idx} = '*';
                    else
                        sig_labels{posneg_idx} = 'N.S.';    
                    end
                end

                % Bar plot with error bars
                tile_idx = (state_idx-1)*3 + animal_idx;
                tile = nexttile(t1, tile_idx);
                ratios = ratios * 100; % convert to percentage, posneg x prepost
                low_errs = ratios - squeeze(CIs(:, :, 1))*100; % posneg x prepost
                high_errs = squeeze(CIs(:, :, 2))*100 - ratios; % posneg x prepost

                hold(tile, "on");
                % Pre
                bar(tile, (1:2)-0.15, ratios(:, 1), 'FaceColor', 'b', 'BarWidth', 0.3);
                % Post
                bar(tile, (1:2)+0.15, ratios(:, 2), 'FaceColor', 'r', 'BarWidth', 0.3);
                % Error bars
                errorbar(tile, (1:2)-0.15, ratios(:, 1), low_errs(:, 1), high_errs(:, 1), 'k', 'LineStyle', 'none');
                errorbar(tile, (1:2)+0.15, ratios(:, 2), low_errs(:, 2), high_errs(:, 2), 'k', 'LineStyle', 'none');
                % Pre vs Post significance
                text(tile, 1, max(ratios(:, 1))+5, sig_labels{1}, 'HorizontalAlignment', 'center', 'FontSize', 14);
                text(tile, 2, max(ratios(:, 2))+5, sig_labels{2}, 'HorizontalAlignment', 'center', 'FontSize', 14);
                hold(tile, "off");

                xticks(tile, 1:2);
                xticklabels(tile, posneg_labels);
                ylabel(tile, 'Significant J %');
                legend(tile, prepost_labels, 'Location', 'Best');
                title(tile, sprintf('%s, %s', animal_name, state_label));
                ylim([0, 30]);

                disp(J_counts);
                disp(total_counts);
            end
        end
        sgtitle(f1, sprintf("%s, Kernel %d, Pre vs Post", injection, kernel_idx), 'FontSize', 14);

        save_folder = fullfile(root, 'Figures', 'J_hist');
        check_path(save_folder);
        save_path = fullfile(save_folder, sprintf('J_hist_prepost_k%d_%s.png', kernel_idx, injection));
        saveas(f1, save_path);

        % fig 2: open vs closed, can combine with previous one
        f2 = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
        t2 = tiledlayout(f2, 2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

        injection = injection_types{injection_idx};
        for animal_idx = 1:3
            animal_name = animals{animal_idx};
            for prepost_idx = 1:2
                prepost_str = prepost_labels{prepost_idx};

                J_counts = zeros(2, 2); % pos/neg x open/closed
                total_counts = zeros(2, 2); % total count for pos/neg in open/closed

                for state_idx = 1:2
                    state = states{state_idx};
                    state_label = state_labels{state_idx};
                    % filter metadata for current condition 
                    % TODO: better condition filtering by a function
                    if strcmp(animal_name, 'Both')
                        animal_filter = strcmp(mt.GLM.animal_name, 'Slayer') | strcmp(mt.GLM.animal_name, 'Emperor');
                    else
                        animal_filter = strcmp(mt.GLM.animal_name, animal_name);
                    end
                    condition_filter = (animal_filter) & ...
                                        strcmp(mt.GLM.state, state) & ...
                                        strcmp(mt.GLM.kernel_name, "DeltaPure") & ...
                                        strcmp(mt.GLM.align, alignment) & ...
                                        strcmp(mt.GLM.area, "Cortex") & ...
                                        strcmp(mt.GLM.injection, injection) & ...
                                        (mt.GLM.epoch == 3000) & ...
                                        (mt.GLM.fold_idx == 0) & ...
                                        (mt.GLM.shuffle_idx == 0) & ...
                                        cellfun(@(x) ~isempty(x) && x == resting_dur_threshold, mt.GLM.resting_dur_threshold) & ...
                                        strcmp(mt.GLM.prepost, prepost_str);
                                        
                    condition_meta = mt.GLM(condition_filter, :);
                    if isempty(condition_meta)
                        warning('No data found for condition: %s, %s, %s', animal_name, state, prepost_str);
                        continue;
                    else
                        fprintf('Plotting condition: %s, %s, %s %s with %d sessions.\n', animal_name, state, prepost_str, injection, height(condition_meta));
                    end

                    % count significant J 
                    pos_count   = 0;
                    neg_count   = 0;
                    total_count = 0;
                    for session_idx = 1:height(condition_meta)
                        meta = condition_meta(session_idx, :);
                        meta = table2struct(meta);

                        % load model parameters and errors
                        data_folder = fullfile(root, 'Data', 'Working', 'GLM');
                        file_name = meta.file_name;
                        file_path = fullfile(data_folder, file_name);
                        if ~isfile(file_path)
                            warning('File not found: %s. Skipping.', file_path);
                            continue;
                        end
                        session_data = load(file_path, 'meta', 'data');
                        model_par = session_data.data.model_par;
                        model_err = session_data.data.model_err.total;

                        % load area borders
                        border_folder = fullfile(root, 'Data', 'Working', 'border');
                        border_file_name = generate_filename('border', meta);
                        border_file_path = fullfile(border_folder, border_file_name);
                        if ~isfile(border_file_path)
                            warning('Border file not found: %s. Skipping.', border_file_path);
                            continue;
                        end
                        border_data = load(border_file_path, 'data', 'meta');
                        borders = border_data.data.borders;
                        N = border_data.meta.N;
                        assert(numel(borders) == 2);
                        borders = [borders, N+1]; % add end border

                        selected_areas = {[1, 2], [2, 1]}; % between ACC and VLPFC

                        % count significant J values
                        [pos, neg, total] = J_count(model_par, model_err, kernel_idx, borders, selected_areas);
                        pos_count = pos_count + pos;
                        neg_count = neg_count + neg;
                        total_count = total_count + total;
                    end
                    J_counts(:, state_idx) = [pos_count; neg_count];
                    total_counts(:, state_idx) = [total_count; total_count];
                end

                % error bar and statistical test
                ratios = J_counts ./ total_counts;
                CIs = zeros(2, 2, 2); % pos/neg x open/closed x low/high
                p_vals = zeros(1, 2); % pos/neg
                sig_labels = cell(1, 2);
                for posneg_idx = 1:2
                    count = J_counts(posneg_idx, :);
                    total = total_counts(posneg_idx, :);
                    for state_idx = 1:2
                        [CI_low, CI_high] = wilsonCI(count(state_idx), total(state_idx), 0.05);
                        CIs(posneg_idx, state_idx, :) = [CI_low, CI_high];
                    end
                    p_val = twoProportionPValue(count(1), total(1), count(2), total(2));
                    p_vals(posneg_idx) = p_val;
                    if p_val < 0.001
                        sig_labels{posneg_idx} = '***';
                    elseif p_val < 0.01
                        sig_labels{posneg_idx} = '**';
                    elseif p_val < 0.05
                        sig_labels{posneg_idx} = '*';
                    else
                        sig_labels{posneg_idx} = 'N.S.';    
                    end
                end

                % Bar plot with error bars
                % fig 2: open vs closed
                tile_idx = (prepost_idx-1)*3 + animal_idx;
                tile = nexttile(t2, tile_idx);
                ratios = ratios * 100; % convert to percentage, posneg x prepost
                low_errs = ratios - squeeze(CIs(:, :, 1))*100; % posneg x prepost
                high_errs = squeeze(CIs(:, :, 2))*100 - ratios; % posneg x prepost

                hold(tile, "on");
                % Open
                bar(tile, (1:2)-0.15, ratios(:, 1), 'FaceColor', 'b', 'BarWidth', 0.3);
                % Closed
                bar(tile, (1:2)+0.15, ratios(:, 2), 'FaceColor', 'r', 'BarWidth', 0.3);
                % Error bars
                errorbar(tile, (1:2)-0.15, ratios(:, 1), low_errs(:, 1), high_errs(:, 1), 'k', 'LineStyle', 'none');
                errorbar(tile, (1:2)+0.15, ratios(:, 2), low_errs(:, 2), high_errs(:, 2), 'k', 'LineStyle', 'none');
                % Open vs Closed significance
                text(tile, 1, max(ratios(:, 1))+5, sig_labels{1}, 'HorizontalAlignment', 'center', 'FontSize', 14);
                text(tile, 2, max(ratios(:, 2))+5, sig_labels{2}, 'HorizontalAlignment', 'center', 'FontSize', 14);
                hold(tile, "off");

                xticks(tile, 1:2);
                xticklabels(tile, posneg_labels);
                ylabel(tile, 'Significant J %');
                legend(tile, state_labels, 'Location', 'Best');
                title(tile, sprintf('%s, %s', animal_name, prepost_str));
                ylim([0, 30]);

                disp(J_counts);
                disp(total_counts);
            end
        end
        sgtitle(f2, sprintf("%s, Kernel %d, Open vs Closed", injection, kernel_idx), 'FontSize', 14);

        save_folder = fullfile(root, 'Figures', 'J_hist');
        check_path(save_folder);
        save_path = fullfile(save_folder, sprintf('J_hist_openclosed_k%d_%s.png', kernel_idx, injection));
        saveas(f2, save_path);
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
