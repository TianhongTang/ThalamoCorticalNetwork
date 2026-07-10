%% J_bars_kernels_grouped_by_sign.m - Plot significant J percentage across kernels.
% Layout:
%   Columns = Positive J / Negative J.
%   X-axis within each tile = Kernel 1/2/3.
%   Bar groups within each kernel = compared conditions.
% Only All animals is used; individual Slayer/Emperor panels are removed.

clear;

%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end

addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Main parameters
resting_dur_threshold = 15;
mt = load_meta(root, 'table');

kernel_indices = 1:3;
kernel_labels = arrayfun(@(k) sprintf('Kernel %d', k), kernel_indices, 'UniformOutput', false);

injection_types = {'Muscimol', 'Saline'};
states = {'RestOpen', 'RestClose'};
state_labels = {'Eyes Open', 'Eyes Closed'};
posneg_labels = {'Positive J', 'Negative J'};
prepost_labels = {'Pre', 'Post'};
alignment = 'Last';

bar_width = 0.30;
bar_offset = 0.18;
y_limit = [0, 20];

% Pooled connection-level effect figures.
% Effects are in percentage points: 100 * (Post proportion - Pre proportion).
effect_injection_types = {'Saline', 'Muscimol'};
effect_injection_labels = {'Saline', 'Muscimol'};
effect_y_limit = [];      % [] = automatic. Example: [-10, 10].
effect_diff_y_limit = []; % [] = automatic. Example: [-10, 10].
effect_alpha = 0.05;

% Extended pooled effect plotting.
% 'difference' = Post - Pre in percentage points.
% 'ratio'      = Post / Pre significant-J ratio.
effect_metric_mode = 'difference';
make_combined_effect_figure = true;
make_single_injection_effect_figures = true;
make_prepost_percent_figures = true;
make_difference_only_figure = false;

% animal_filter = strcmp(mt.GLM.animal_name, 'Slayer') | strcmp(mt.GLM.animal_name, 'Emperor');
animal_filter = strcmp(mt.GLM.animal_name, 'Slayer');

%% Fig1: Pre vs Post. Fig2: Open vs Closed.
for injection_idx = 1:numel(injection_types)
    injection = injection_types{injection_idx};

    %% Figure 1: Pre vs Post.
    % Rows = states. Columns = Positive/Negative J.
    % X-axis = kernels. Bars = Pre/Post.
    f1 = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
    t1 = tiledlayout(f1, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for state_idx = 1:numel(states)
        state = states{state_idx};
        state_label = state_labels{state_idx};

        % kernel x sign x condition
        ratios_all = nan(numel(kernel_indices), 2, numel(prepost_labels));
        low_errs_all = nan(numel(kernel_indices), 2, numel(prepost_labels));
        high_errs_all = nan(numel(kernel_indices), 2, numel(prepost_labels));
        sig_labels_all = cell(numel(kernel_indices), 2);

        for kernel_i = 1:numel(kernel_indices)
            kernel_idx = kernel_indices(kernel_i);

            J_counts = zeros(2, numel(prepost_labels)); % pos/neg x pre/post
            total_counts = zeros(2, numel(prepost_labels));

            for prepost_idx = 1:numel(prepost_labels)
                prepost_str = prepost_labels{prepost_idx};

                condition_filter = animal_filter & ...
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
                    warning('No data found for condition: All animals, %s, %s %s, Kernel %d', ...
                        state, prepost_str, injection, kernel_idx);
                    continue;
                else
                    fprintf('Plotting condition: All animals, %s, %s %s, Kernel %d with %d sessions.\n', ...
                        state, prepost_str, injection, kernel_idx, height(condition_meta));
                end

                [pos_count, neg_count, total_count] = count_condition_J(root, condition_meta, kernel_idx);
                J_counts(:, prepost_idx) = [pos_count; neg_count];
                total_counts(:, prepost_idx) = [total_count; total_count];
            end

            [ratios, low_errs, high_errs, sig_labels] = compute_bar_stats(J_counts, total_counts);
            ratios_all(kernel_i, :, :) = ratios;
            low_errs_all(kernel_i, :, :) = low_errs;
            high_errs_all(kernel_i, :, :) = high_errs;
            sig_labels_all(kernel_i, :) = sig_labels;

            disp(J_counts);
            disp(total_counts);
        end

        for posneg_idx = 1:2
            tile_idx = (state_idx - 1) * 2 + posneg_idx;
            tile = nexttile(t1, tile_idx);

            plot_kernel_grouped_bars(tile, squeeze(ratios_all(:, posneg_idx, :)), ...
                squeeze(low_errs_all(:, posneg_idx, :)), squeeze(high_errs_all(:, posneg_idx, :)), ...
                sig_labels_all(:, posneg_idx), prepost_labels, kernel_labels, bar_offset, bar_width);

            ylabel(tile, 'Significant J %');
            title(tile, sprintf('%s, %s', state_label, posneg_labels{posneg_idx}));
            ylim(tile, y_limit);
        end
    end

    sgtitle(f1, sprintf("%s, Pre vs Post", injection), 'FontSize', 14);
    export_figure(root, f1, sprintf('J_bars_prepost_kernel_groups_by_sign_%s', injection));

    %% Figure 2: Open vs Closed.
    % Rows = Pre/Post. Columns = Positive/Negative J.
    % X-axis = kernels. Bars = Eyes Open/Eyes Closed.
    f2 = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
    t2 = tiledlayout(f2, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for prepost_idx = 1:numel(prepost_labels)
        prepost_str = prepost_labels{prepost_idx};

        % kernel x sign x condition
        ratios_all = nan(numel(kernel_indices), 2, numel(states));
        low_errs_all = nan(numel(kernel_indices), 2, numel(states));
        high_errs_all = nan(numel(kernel_indices), 2, numel(states));
        sig_labels_all = cell(numel(kernel_indices), 2);

        for kernel_i = 1:numel(kernel_indices)
            kernel_idx = kernel_indices(kernel_i);

            J_counts = zeros(2, numel(states)); % pos/neg x open/closed
            total_counts = zeros(2, numel(states));

            for state_idx = 1:numel(states)
                state = states{state_idx};

                condition_filter = animal_filter & ...
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
                    warning('No data found for condition: All animals, %s, %s %s, Kernel %d', ...
                        state, prepost_str, injection, kernel_idx);
                    continue;
                else
                    fprintf('Plotting condition: All animals, %s, %s %s, Kernel %d with %d sessions.\n', ...
                        state, prepost_str, injection, kernel_idx, height(condition_meta));
                end

                [pos_count, neg_count, total_count] = count_condition_J(root, condition_meta, kernel_idx);
                J_counts(:, state_idx) = [pos_count; neg_count];
                total_counts(:, state_idx) = [total_count; total_count];
            end

            [ratios, low_errs, high_errs, sig_labels] = compute_bar_stats(J_counts, total_counts);
            ratios_all(kernel_i, :, :) = ratios;
            low_errs_all(kernel_i, :, :) = low_errs;
            high_errs_all(kernel_i, :, :) = high_errs;
            sig_labels_all(kernel_i, :) = sig_labels;

            disp(J_counts);
            disp(total_counts);
        end

        for posneg_idx = 1:2
            tile_idx = (prepost_idx - 1) * 2 + posneg_idx;
            tile = nexttile(t2, tile_idx);

            plot_kernel_grouped_bars(tile, squeeze(ratios_all(:, posneg_idx, :)), ...
                squeeze(low_errs_all(:, posneg_idx, :)), squeeze(high_errs_all(:, posneg_idx, :)), ...
                sig_labels_all(:, posneg_idx), state_labels, kernel_labels, bar_offset, bar_width);

            ylabel(tile, 'Significant J %');
            title(tile, sprintf('%s, %s', prepost_str, posneg_labels{posneg_idx}));
            ylim(tile, y_limit);
        end
    end

    sgtitle(f2, sprintf("%s, Open vs Closed", injection), 'FontSize', 14);
    export_figure(root, f2, sprintf('J_bars_openclosed_kernel_groups_by_sign_%s', injection));
end

%% Extended pooled effect figures.
% Collect pooled Pre/Post counts once, then reuse them for all effect figures.
effect_stats_grid = collect_effect_stats_grid_pooled_ext(root, mt.GLM, animal_filter, states, ...
    effect_injection_types, kernel_indices, alignment, resting_dur_threshold);

if make_combined_effect_figure
    plot_effect_metric_combined_pooled_ext(root, effect_stats_grid, states, state_labels, posneg_labels, ...
        kernel_indices, kernel_labels, effect_injection_types, effect_injection_labels, ...
        bar_offset, bar_width, effect_metric_mode, effect_y_limit, effect_alpha);
end

if make_single_injection_effect_figures
    for effect_injection_idx = 1:numel(effect_injection_types)
        plot_effect_metric_single_injection_pooled_ext(root, effect_stats_grid, states, state_labels, posneg_labels, ...
            kernel_indices, kernel_labels, effect_injection_types, effect_injection_labels, effect_injection_idx, ...
            effect_metric_mode, effect_y_limit, effect_alpha);
    end
end

if make_prepost_percent_figures
    for effect_injection_idx = 1:numel(effect_injection_types)
        plot_prepost_percent_injection_pooled_ext(root, effect_stats_grid, states, state_labels, posneg_labels, ...
            kernel_indices, kernel_labels, effect_injection_types, effect_injection_labels, effect_injection_idx, ...
            bar_offset, bar_width, effect_alpha);
    end
end

if make_difference_only_figure
    plot_muscimol_minus_saline_effect_pooled_ext(root, effect_stats_grid, states, state_labels, posneg_labels, ...
        kernel_indices, kernel_labels, effect_injection_types, effect_metric_mode, effect_diff_y_limit, effect_alpha);
end

function plot_kernel_grouped_bars(tile, ratios, low_errs, high_errs, sig_labels, condition_labels, kernel_labels, bar_offset, bar_width)
    % ratios: kernel x condition
    n_kernel = size(ratios, 1);
    x = 1:n_kernel;

    hold(tile, "on");

    bar(tile, x - bar_offset, ratios(:, 1), 'FaceColor', 'b', 'BarWidth', bar_width);
    bar(tile, x + bar_offset, ratios(:, 2), 'FaceColor', 'r', 'BarWidth', bar_width);

    errorbar(tile, x - bar_offset, ratios(:, 1), low_errs(:, 1), high_errs(:, 1), ...
        'k', 'LineStyle', 'none');
    errorbar(tile, x + bar_offset, ratios(:, 2), low_errs(:, 2), high_errs(:, 2), ...
        'k', 'LineStyle', 'none');

    for kernel_i = 1:n_kernel
        y_text = max(ratios(kernel_i, :) + high_errs(kernel_i, :), [], 'omitnan') + 2;
        if ~isfinite(y_text)
            y_text = 2;
        end
        text(tile, x(kernel_i), y_text, sig_labels{kernel_i}, ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
    end

    hold(tile, "off");

    xticks(tile, x);
    xticklabels(tile, kernel_labels);
    legend(tile, condition_labels, 'Location', 'Best');
end

function [pos_count, neg_count, total_count] = count_condition_J(root, condition_meta, kernel_idx)
    pos_count = 0;
    neg_count = 0;
    total_count = 0;

    for session_idx = 1:height(condition_meta)
        meta = condition_meta(session_idx, :);
        meta = table2struct(meta);

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
        borders = [borders, N+1];

        selected_areas = {[1, 2], [2, 1]}; % between ACC and VLPFC

        [pos, neg, total] = J_count(model_par, model_err, kernel_idx, borders, selected_areas);
        pos_count = pos_count + pos;
        neg_count = neg_count + neg;
        total_count = total_count + total;
    end
end

function [ratios, low_errs, high_errs, sig_labels] = compute_bar_stats(J_counts, total_counts)
    ratios_raw = J_counts ./ total_counts;
    CIs = zeros(2, 2, 2); % pos/neg x condition x low/high
    sig_labels = cell(1, 2);

    for posneg_idx = 1:2
        count = J_counts(posneg_idx, :);
        total = total_counts(posneg_idx, :);

        for condition_idx = 1:2
            [CI_low, CI_high] = wilsonCI(count(condition_idx), total(condition_idx), 0.05);
            CIs(posneg_idx, condition_idx, :) = [CI_low, CI_high];
        end

        p_val = twoProportionPValue(count(1), total(1), count(2), total(2));
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

    ratios = ratios_raw * 100;
    low_errs = ratios - squeeze(CIs(:, :, 1)) * 100;
    high_errs = squeeze(CIs(:, :, 2)) * 100 - ratios;
end


function plot_prepost_effect_comparison_pooled(root, mt_glm, animal_filter, states, state_labels, ...
    posneg_labels, kernel_indices, kernel_labels, injection_types, injection_labels, ...
    alignment, resting_dur_threshold, bar_offset, bar_width, y_limit, alpha)

    if numel(injection_types) ~= 2 || ~strcmp(injection_types{1}, 'Saline') || ~strcmp(injection_types{2}, 'Muscimol')
        error('plot_prepost_effect_comparison_pooled expects injection_types = {''Saline'', ''Muscimol''}.');
    end

    f = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for state_idx = 1:numel(states)
        state = states{state_idx};
        state_label = state_labels{state_idx};

        for posneg_idx = 1:2
            effects = nan(numel(kernel_indices), numel(injection_types));
            low_errs = nan(numel(kernel_indices), numel(injection_types));
            high_errs = nan(numel(kernel_indices), numel(injection_types));
            sig_labels = cell(numel(kernel_indices), 1);

            for kernel_i = 1:numel(kernel_indices)
                kernel_idx = kernel_indices(kernel_i);
                stats_by_injection = cell(1, numel(injection_types));

                for injection_idx = 1:numel(injection_types)
                    injection = injection_types{injection_idx};
                    stats = collect_pooled_prepost_counts(root, mt_glm, animal_filter, state, injection, ...
                        kernel_idx, alignment, resting_dur_threshold);
                    stats_by_injection{injection_idx} = stats;

                    [effect_pct, ci_low_pct, ci_high_pct] = prepost_effect_ci(stats, posneg_idx, alpha);
                    effects(kernel_i, injection_idx) = effect_pct;
                    low_errs(kernel_i, injection_idx) = effect_pct - ci_low_pct;
                    high_errs(kernel_i, injection_idx) = ci_high_pct - effect_pct;
                end

                % Difference-in-differences test: (Muscimol Post-Pre) - (Saline Post-Pre) = 0.
                [~, ~, ~, p_val] = prepost_effect_difference_ci( ...
                    stats_by_injection{2}, stats_by_injection{1}, posneg_idx, alpha);
                sig_labels{kernel_i} = p_to_sig_label(p_val);
            end

            tile_idx = (state_idx - 1) * 2 + posneg_idx;
            tile = nexttile(t, tile_idx);
            plot_effect_grouped_bars(tile, effects, low_errs, high_errs, sig_labels, ...
                injection_labels, kernel_labels, bar_offset, bar_width, y_limit);
            ylabel(tile, 'Post - Pre significant J (percentage points)');
            title(tile, sprintf('%s, %s', state_label, posneg_labels{posneg_idx}));
            draw_zero_line(tile);
        end
    end

    sgtitle(f, 'Pooled Post - Pre effect: Saline vs Muscimol', 'FontSize', 14);
    export_figure(root, f, 'J_bars_post_minus_pre_effect_pooled_saline_muscimol');
end

function plot_muscimol_minus_saline_effect_pooled(root, mt_glm, animal_filter, states, state_labels, ...
    posneg_labels, kernel_indices, kernel_labels, injection_types, alignment, resting_dur_threshold, ...
    y_limit, alpha)

    if numel(injection_types) ~= 2 || ~strcmp(injection_types{1}, 'Saline') || ~strcmp(injection_types{2}, 'Muscimol')
        error('plot_muscimol_minus_saline_effect_pooled expects injection_types = {''Saline'', ''Muscimol''}.');
    end

    f = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for state_idx = 1:numel(states)
        state = states{state_idx};
        state_label = state_labels{state_idx};

        for posneg_idx = 1:2
            diff_effects = nan(numel(kernel_indices), 1);
            low_errs = nan(numel(kernel_indices), 1);
            high_errs = nan(numel(kernel_indices), 1);
            sig_labels = cell(numel(kernel_indices), 1);

            for kernel_i = 1:numel(kernel_indices)
                kernel_idx = kernel_indices(kernel_i);
                saline_stats = collect_pooled_prepost_counts(root, mt_glm, animal_filter, state, 'Saline', ...
                    kernel_idx, alignment, resting_dur_threshold);
                muscimol_stats = collect_pooled_prepost_counts(root, mt_glm, animal_filter, state, 'Muscimol', ...
                    kernel_idx, alignment, resting_dur_threshold);

                [diff_pct, ci_low_pct, ci_high_pct, p_val] = prepost_effect_difference_ci( ...
                    muscimol_stats, saline_stats, posneg_idx, alpha);

                diff_effects(kernel_i) = diff_pct;
                low_errs(kernel_i) = diff_pct - ci_low_pct;
                high_errs(kernel_i) = ci_high_pct - diff_pct;
                sig_labels{kernel_i} = p_to_sig_label(p_val);
            end

            tile_idx = (state_idx - 1) * 2 + posneg_idx;
            tile = nexttile(t, tile_idx);
            plot_effect_difference_bars(tile, diff_effects, low_errs, high_errs, sig_labels, kernel_labels, y_limit);
            ylabel(tile, '(Muscimol Post-Pre) - (Saline Post-Pre) (percentage points)');
            title(tile, sprintf('%s, %s', state_label, posneg_labels{posneg_idx}));
            draw_zero_line(tile);
        end
    end

    sgtitle(f, 'Pooled difference of effects: Muscimol - Saline', 'FontSize', 14);
    export_figure(root, f, 'J_bars_muscimol_minus_saline_prepost_effect_pooled');
end

function stats = collect_pooled_prepost_counts(root, mt_glm, animal_filter, state, injection, ...
    kernel_idx, alignment, resting_dur_threshold)

    pre_meta = filter_condition_meta(mt_glm, animal_filter, state, injection, 'Pre', alignment, resting_dur_threshold);
    post_meta = filter_condition_meta(mt_glm, animal_filter, state, injection, 'Post', alignment, resting_dur_threshold);

    [pre_pos, pre_neg, pre_total] = count_condition_J(root, pre_meta, kernel_idx);
    [post_pos, post_neg, post_total] = count_condition_J(root, post_meta, kernel_idx);

    stats = struct();
    stats.pre_counts = [pre_pos; pre_neg];
    stats.post_counts = [post_pos; post_neg];
    stats.pre_total = pre_total;
    stats.post_total = post_total;
end

function condition_meta = filter_condition_meta(mt_glm, animal_filter, state, injection, prepost, alignment, resting_dur_threshold)
    condition_filter = animal_filter & ...
                       strcmp(mt_glm.state, state) & ...
                       strcmp(mt_glm.kernel_name, "DeltaPure") & ...
                       strcmp(mt_glm.align, alignment) & ...
                       strcmp(mt_glm.area, "Cortex") & ...
                       strcmp(mt_glm.injection, injection) & ...
                       (mt_glm.epoch == 3000) & ...
                       (mt_glm.fold_idx == 0) & ...
                       (mt_glm.shuffle_idx == 0) & ...
                       cellfun(@(x) ~isempty(x) && x == resting_dur_threshold, mt_glm.resting_dur_threshold) & ...
                       strcmp(mt_glm.prepost, prepost);
    condition_meta = mt_glm(condition_filter, :);
end

function [effect_pct, ci_low_pct, ci_high_pct] = prepost_effect_ci(stats, posneg_idx, alpha)
    if stats.pre_total == 0 || stats.post_total == 0
        effect_pct = NaN;
        ci_low_pct = NaN;
        ci_high_pct = NaN;
        return;
    end

    p_pre = stats.pre_counts(posneg_idx) / stats.pre_total;
    p_post = stats.post_counts(posneg_idx) / stats.post_total;
    effect = p_post - p_pre;
    se = sqrt(p_post * (1 - p_post) / stats.post_total + ...
              p_pre * (1 - p_pre) / stats.pre_total);
    zcrit = norminv(1 - alpha / 2);

    effect_pct = 100 * effect;
    ci_low_pct = 100 * (effect - zcrit * se);
    ci_high_pct = 100 * (effect + zcrit * se);
end

function [diff_pct, ci_low_pct, ci_high_pct, p_val] = prepost_effect_difference_ci(stats_a, stats_b, posneg_idx, alpha)
    % stats_a - stats_b, where each effect is Post - Pre.
    if stats_a.pre_total == 0 || stats_a.post_total == 0 || stats_b.pre_total == 0 || stats_b.post_total == 0
        diff_pct = NaN;
        ci_low_pct = NaN;
        ci_high_pct = NaN;
        p_val = NaN;
        return;
    end

    a_pre = stats_a.pre_counts(posneg_idx) / stats_a.pre_total;
    a_post = stats_a.post_counts(posneg_idx) / stats_a.post_total;
    b_pre = stats_b.pre_counts(posneg_idx) / stats_b.pre_total;
    b_post = stats_b.post_counts(posneg_idx) / stats_b.post_total;

    diff = (a_post - a_pre) - (b_post - b_pre);
    se = sqrt(a_post * (1 - a_post) / stats_a.post_total + ...
              a_pre * (1 - a_pre) / stats_a.pre_total + ...
              b_post * (1 - b_post) / stats_b.post_total + ...
              b_pre * (1 - b_pre) / stats_b.pre_total);

    if se == 0
        if diff == 0
            p_val = 1;
        else
            p_val = 0;
        end
        ci_low = diff;
        ci_high = diff;
    else
        z = diff / se;
        p_val = twoSidedNormalPValue(z);
        zcrit = norminv(1 - alpha / 2);
        ci_low = diff - zcrit * se;
        ci_high = diff + zcrit * se;
    end

    diff_pct = 100 * diff;
    ci_low_pct = 100 * ci_low;
    ci_high_pct = 100 * ci_high;
end

function p_val = twoSidedNormalPValue(z)
    if ~isfinite(z)
        p_val = NaN;
        return;
    end
    abs_z = abs(z);
    Phi = 0.5 * erfc(-abs_z / sqrt(2));
    p_val = 2 * (1 - Phi);
end

function label = p_to_sig_label(p_val)
    if ~isfinite(p_val)
        label = 'n/a';
    elseif p_val < 0.001
        label = '***';
    elseif p_val < 0.01
        label = '**';
    elseif p_val < 0.05
        label = '*';
    else
        label = 'N.S.';
    end
end

function plot_effect_grouped_bars(tile, effects, low_errs, high_errs, sig_labels, condition_labels, kernel_labels, bar_offset, bar_width, y_limit)
    n_kernel = size(effects, 1);
    x = 1:n_kernel;

    hold(tile, 'on');
    bar(tile, x - bar_offset, effects(:, 1), 'FaceColor', 'b', 'BarWidth', bar_width);
    bar(tile, x + bar_offset, effects(:, 2), 'FaceColor', 'r', 'BarWidth', bar_width);

    errorbar(tile, x - bar_offset, effects(:, 1), low_errs(:, 1), high_errs(:, 1), ...
        'k', 'LineStyle', 'none');
    errorbar(tile, x + bar_offset, effects(:, 2), low_errs(:, 2), high_errs(:, 2), ...
        'k', 'LineStyle', 'none');

    set_effect_ylim(tile, effects, low_errs, high_errs, y_limit);
    yl = ylim(tile);
    y_offset = 0.04 * diff(yl);
    for kernel_i = 1:n_kernel
        y_text = max([effects(kernel_i, :) + high_errs(kernel_i, :), 0], [], 'omitnan') + y_offset;
        if ~isfinite(y_text)
            y_text = yl(2) - y_offset;
        end
        text(tile, x(kernel_i), y_text, sig_labels{kernel_i}, ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
    end

    hold(tile, 'off');
    xticks(tile, x);
    xticklabels(tile, kernel_labels);
    legend(tile, condition_labels, 'Location', 'Best');
end

function plot_effect_difference_bars(tile, diff_effects, low_errs, high_errs, sig_labels, kernel_labels, y_limit)
    n_kernel = numel(diff_effects);
    x = 1:n_kernel;

    hold(tile, 'on');
    bar(tile, x, diff_effects, 'FaceColor', [0.4, 0.4, 0.4], 'BarWidth', 0.45);
    errorbar(tile, x, diff_effects, low_errs, high_errs, 'k', 'LineStyle', 'none');

    set_effect_ylim(tile, diff_effects(:), low_errs(:), high_errs(:), y_limit);
    yl = ylim(tile);
    y_offset = 0.04 * diff(yl);
    for kernel_i = 1:n_kernel
        y_text = max([diff_effects(kernel_i) + high_errs(kernel_i), 0], [], 'omitnan') + y_offset;
        if ~isfinite(y_text)
            y_text = yl(2) - y_offset;
        end
        text(tile, x(kernel_i), y_text, sig_labels{kernel_i}, ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
    end

    hold(tile, 'off');
    xticks(tile, x);
    xticklabels(tile, kernel_labels);
end

function set_effect_ylim(tile, values, low_errs, high_errs, y_limit)
    if ~isempty(y_limit)
        ylim(tile, y_limit);
        return;
    end

    lower = values(:) - low_errs(:);
    upper = values(:) + high_errs(:);
    vals = [lower; upper; 0];
    vals = vals(isfinite(vals));
    if isempty(vals)
        ylim(tile, [-1, 1]);
        return;
    end
    yr = max(vals) - min(vals);
    if yr == 0
        yr = max(abs(vals));
    end
    if yr == 0
        yr = 1;
    end
    ylim(tile, [min(vals) - 0.15 * yr, max(vals) + 0.25 * yr]);
end

function draw_zero_line(tile)
    xl = xlim(tile);
    hold(tile, 'on');
    plot(tile, xl, [0, 0], 'k-', 'LineWidth', 0.75);
    xlim(tile, xl);
    hold(tile, 'off');
end

function export_figure(root, fig, output_stub)
    save_folder = fullfile(root, 'Figures', 'Paper', 'J_bars');
    check_path(save_folder);

    figWidth  = 12.0;
    figHeight = 8.0;
    resolution = 300;

    set(fig, 'Units', 'inches');
    fig.Position(3:4) = [figWidth, figHeight];

    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperSize', [figWidth, figHeight]);
    set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);
    set(fig, 'Color', 'w');

    preview_filename = fullfile(save_folder, [output_stub, '_preview.jpg']);
    exportgraphics(fig, preview_filename, ...
        'ContentType', 'image', ...
        'BackgroundColor', 'white', ...
        'Resolution', resolution);

    pdf_filename = fullfile(save_folder, [output_stub, '.pdf']);
    exportgraphics(fig, pdf_filename, ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white', ...
        'Resolution', resolution);

    close(fig);
end

function [p_low, p_high] = wilsonCI(M, N, alpha)
    if nargin < 3
        alpha = 0.05;
    end

    if N == 0
        p_low = NaN;
        p_high = NaN;
        return;
    end

    p = M / N;
    z = norminv(1 - alpha/2);

    denominator = 1 + (z^2)/N;
    center = p + (z^2)/(2*N);
    radius = z * sqrt((p*(1-p)/N) + (z^2)/(4*N^2));

    p_low = (center - radius) / denominator;
    p_high = (center + radius) / denominator;
end

function p_val = twoProportionPValue(success1, trials1, success2, trials2)
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


%% Extended pooled effect plotting functions.
function stats_grid = collect_effect_stats_grid_pooled_ext(root, mt_glm, animal_filter, states, injection_types, kernel_indices, alignment, resting_dur_threshold)
    stats_grid = cell(numel(states), numel(kernel_indices), numel(injection_types));
    for state_idx = 1:numel(states)
        state = states{state_idx};
        for kernel_i = 1:numel(kernel_indices)
            kernel_idx = kernel_indices(kernel_i);
            for injection_idx = 1:numel(injection_types)
                injection = injection_types{injection_idx};
                stats_grid{state_idx, kernel_i, injection_idx} = collect_pooled_prepost_counts(root, mt_glm, animal_filter, state, injection, ...
                    kernel_idx, alignment, resting_dur_threshold);
            end
        end
    end
end

function plot_effect_metric_combined_pooled_ext(root, stats_grid, states, state_labels, posneg_labels, ...
    kernel_indices, kernel_labels, injection_types, injection_labels, bar_offset, bar_width, metric_mode, y_limit, alpha)

    if numel(injection_types) ~= 2 || ~strcmp(injection_types{1}, 'Saline') || ~strcmp(injection_types{2}, 'Muscimol')
        error('Combined effect figure expects injection_types = {''Saline'', ''Muscimol''}.');
    end

    [values, low_errs, high_errs, prepost_sig_labels, did_values, did_low_errs, did_high_errs, did_sig_labels] = ...
        compute_metric_arrays_from_stats_grid_ext(stats_grid, states, posneg_labels, kernel_indices, injection_types, metric_mode, alpha);

    if isempty(y_limit)
        y_limit = common_ylim_for_metric_ext(values, low_errs, high_errs, metric_mode);
    end

    f = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for state_idx = 1:numel(states)
        for posneg_idx = 1:2
            tile = nexttile(t, (state_idx - 1) * 2 + posneg_idx);
            plot_metric_grouped_bars_ext(tile, squeeze(values(state_idx, posneg_idx, :, :)), ...
                squeeze(low_errs(state_idx, posneg_idx, :, :)), squeeze(high_errs(state_idx, posneg_idx, :, :)), ...
                squeeze(prepost_sig_labels(state_idx, posneg_idx, :, :)), squeeze(did_sig_labels(state_idx, posneg_idx, :)), ...
                injection_labels, kernel_labels, bar_offset, bar_width, y_limit, metric_mode);
            ylabel(tile, metric_axis_label_ext(metric_mode));
            title(tile, sprintf('%s, %s', state_labels{state_idx}, posneg_labels{posneg_idx}));
        end
    end

    sgtitle(f, sprintf('Pooled %s: Saline vs Muscimol', metric_axis_label_ext(metric_mode)), 'FontSize', 14, 'Interpreter', 'none');
    export_figure(root, f, sprintf('J_bars_%s_effect_pooled_saline_muscimol_combined', metric_output_stub_ext(metric_mode)));
end

function plot_effect_metric_single_injection_pooled_ext(root, stats_grid, states, state_labels, posneg_labels, ...
    kernel_indices, kernel_labels, injection_types, injection_labels, injection_idx, metric_mode, y_limit, alpha)

    [values, low_errs, high_errs, prepost_sig_labels] = ...
        compute_metric_arrays_from_stats_grid_ext(stats_grid, states, posneg_labels, kernel_indices, injection_types, metric_mode, alpha);

    values_one = values(:, :, :, injection_idx);
    low_one = low_errs(:, :, :, injection_idx);
    high_one = high_errs(:, :, :, injection_idx);
    labels_one = prepost_sig_labels(:, :, :, injection_idx);

    if isempty(y_limit)
        y_limit = common_ylim_for_metric_ext(values_one, low_one, high_one, metric_mode);
    end

    f = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    injection_label = injection_labels{injection_idx};

    for state_idx = 1:numel(states)
        for posneg_idx = 1:2
            tile = nexttile(t, (state_idx - 1) * 2 + posneg_idx);
            plot_metric_single_bars_ext(tile, squeeze(values_one(state_idx, posneg_idx, :)), ...
                squeeze(low_one(state_idx, posneg_idx, :)), squeeze(high_one(state_idx, posneg_idx, :)), ...
                squeeze(labels_one(state_idx, posneg_idx, :)), kernel_labels, y_limit, metric_mode, injection_label);
            ylabel(tile, metric_axis_label_ext(metric_mode));
            title(tile, sprintf('%s, %s', state_labels{state_idx}, posneg_labels{posneg_idx}));
        end
    end

    sgtitle(f, sprintf('Pooled %s: %s only', metric_axis_label_ext(metric_mode), injection_label), 'FontSize', 14, 'Interpreter', 'none');
    export_figure(root, f, sprintf('J_bars_%s_effect_pooled_%s_only', metric_output_stub_ext(metric_mode), lower(injection_types{injection_idx})));
end

function plot_prepost_percent_injection_pooled_ext(root, stats_grid, states, state_labels, posneg_labels, ...
    kernel_indices, kernel_labels, injection_types, injection_labels, injection_idx, bar_offset, bar_width, alpha)

    n_state = numel(states);
    n_sign = 2;
    n_kernel = numel(kernel_indices);
    ratios = nan(n_state, n_sign, n_kernel, 2);
    low_errs = nan(n_state, n_sign, n_kernel, 2);
    high_errs = nan(n_state, n_sign, n_kernel, 2);
    sig_labels = cell(n_state, n_sign, n_kernel);

    for state_idx = 1:n_state
        for kernel_i = 1:n_kernel
            stats = stats_grid{state_idx, kernel_i, injection_idx};
            for posneg_idx = 1:n_sign
                counts = [stats.pre_counts(posneg_idx), stats.post_counts(posneg_idx)];
                totals = [stats.pre_total, stats.post_total];
                for pp_idx = 1:2
                    [ci_low, ci_high] = wilsonCI(counts(pp_idx), totals(pp_idx), alpha);
                    p = counts(pp_idx) / totals(pp_idx);
                    ratios(state_idx, posneg_idx, kernel_i, pp_idx) = 100 * p;
                    low_errs(state_idx, posneg_idx, kernel_i, pp_idx) = 100 * (p - ci_low);
                    high_errs(state_idx, posneg_idx, kernel_i, pp_idx) = 100 * (ci_high - p);
                end
                p_val = twoProportionPValue(stats.pre_counts(posneg_idx), stats.pre_total, stats.post_counts(posneg_idx), stats.post_total);
                sig_labels{state_idx, posneg_idx, kernel_i} = p_to_sig_label(p_val);
            end
        end
    end

    y_limit = common_ylim_percent_ext(ratios, low_errs, high_errs);
    f = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    injection_label = injection_labels{injection_idx};

    for state_idx = 1:n_state
        for posneg_idx = 1:n_sign
            tile = nexttile(t, (state_idx - 1) * 2 + posneg_idx);
            plot_prepost_grouped_bars_ext(tile, squeeze(ratios(state_idx, posneg_idx, :, :)), ...
                squeeze(low_errs(state_idx, posneg_idx, :, :)), squeeze(high_errs(state_idx, posneg_idx, :, :)), ...
                squeeze(sig_labels(state_idx, posneg_idx, :)), {'Pre', 'Post'}, kernel_labels, bar_offset, bar_width, y_limit);
            ylabel(tile, 'Significant J %');
            title(tile, sprintf('%s, %s', state_labels{state_idx}, posneg_labels{posneg_idx}));
        end
    end

    sgtitle(f, sprintf('Pooled Pre vs Post: %s', injection_label), 'FontSize', 14, 'Interpreter', 'none');
    export_figure(root, f, sprintf('J_bars_prepost_percent_pooled_%s', lower(injection_types{injection_idx})));
end

function plot_muscimol_minus_saline_effect_pooled_ext(root, stats_grid, states, state_labels, posneg_labels, ...
    kernel_indices, kernel_labels, injection_types, metric_mode, y_limit, alpha)

    if numel(injection_types) ~= 2 || ~strcmp(injection_types{1}, 'Saline') || ~strcmp(injection_types{2}, 'Muscimol')
        error('DID effect figure expects injection_types = {''Saline'', ''Muscimol''}.');
    end

    [~, ~, ~, ~, did_values, did_low_errs, did_high_errs, did_sig_labels] = ...
        compute_metric_arrays_from_stats_grid_ext(stats_grid, states, posneg_labels, kernel_indices, injection_types, metric_mode, alpha);

    if isempty(y_limit)
        y_limit = common_ylim_for_metric_ext(did_values, did_low_errs, did_high_errs, metric_mode);
    end

    f = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    for state_idx = 1:numel(states)
        for posneg_idx = 1:2
            tile = nexttile(t, (state_idx - 1) * 2 + posneg_idx);
            plot_metric_single_bars_ext(tile, squeeze(did_values(state_idx, posneg_idx, :)), ...
                squeeze(did_low_errs(state_idx, posneg_idx, :)), squeeze(did_high_errs(state_idx, posneg_idx, :)), ...
                squeeze(did_sig_labels(state_idx, posneg_idx, :)), kernel_labels, y_limit, metric_mode, 'DID');
            ylabel(tile, did_axis_label_ext(metric_mode));
            title(tile, sprintf('%s, %s', state_labels{state_idx}, posneg_labels{posneg_idx}));
        end
    end
    sgtitle(f, sprintf('Pooled DID: Muscimol vs Saline (%s)', metric_axis_label_ext(metric_mode)), 'FontSize', 14, 'Interpreter', 'none');
    export_figure(root, f, sprintf('J_bars_did_%s_pooled_muscimol_vs_saline', metric_output_stub_ext(metric_mode)));
end

function [values, low_errs, high_errs, prepost_sig_labels, did_values, did_low_errs, did_high_errs, did_sig_labels] = ...
    compute_metric_arrays_from_stats_grid_ext(stats_grid, states, posneg_labels, kernel_indices, injection_types, metric_mode, alpha)

    n_state = numel(states);
    n_sign = numel(posneg_labels);
    n_kernel = numel(kernel_indices);
    n_injection = numel(injection_types);

    values = nan(n_state, n_sign, n_kernel, n_injection);
    low_errs = nan(n_state, n_sign, n_kernel, n_injection);
    high_errs = nan(n_state, n_sign, n_kernel, n_injection);
    prepost_sig_labels = cell(n_state, n_sign, n_kernel, n_injection);
    did_values = nan(n_state, n_sign, n_kernel);
    did_low_errs = nan(n_state, n_sign, n_kernel);
    did_high_errs = nan(n_state, n_sign, n_kernel);
    did_sig_labels = cell(n_state, n_sign, n_kernel);

    for state_idx = 1:n_state
        for posneg_idx = 1:n_sign
            for kernel_i = 1:n_kernel
                for injection_idx = 1:n_injection
                    stats = stats_grid{state_idx, kernel_i, injection_idx};
                    [v, ci_low, ci_high, p_val] = prepost_metric_ci_ext(stats, posneg_idx, metric_mode, alpha);
                    values(state_idx, posneg_idx, kernel_i, injection_idx) = v;
                    low_errs(state_idx, posneg_idx, kernel_i, injection_idx) = v - ci_low;
                    high_errs(state_idx, posneg_idx, kernel_i, injection_idx) = ci_high - v;
                    prepost_sig_labels{state_idx, posneg_idx, kernel_i, injection_idx} = p_to_sig_label(p_val);
                end
                [v, ci_low, ci_high, p_val] = did_metric_ci_ext(stats_grid{state_idx, kernel_i, 2}, stats_grid{state_idx, kernel_i, 1}, posneg_idx, metric_mode, alpha);
                did_values(state_idx, posneg_idx, kernel_i) = v;
                did_low_errs(state_idx, posneg_idx, kernel_i) = v - ci_low;
                did_high_errs(state_idx, posneg_idx, kernel_i) = ci_high - v;
                did_sig_labels{state_idx, posneg_idx, kernel_i} = p_to_sig_label(p_val);
            end
        end
    end
end

function [value, ci_low, ci_high, p_val] = prepost_metric_ci_ext(stats, posneg_idx, metric_mode, alpha)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            [value, ci_low, ci_high] = prepost_effect_ci(stats, posneg_idx, alpha);
            p_val = twoProportionPValue(stats.pre_counts(posneg_idx), stats.pre_total, stats.post_counts(posneg_idx), stats.post_total);
        case {'ratio', 'post/pre', 'post_over_pre'}
            [value, ci_low, ci_high, p_val] = prepost_ratio_ci_ext(stats, posneg_idx, alpha);
        otherwise
            error('Unknown metric_mode: %s', metric_mode);
    end
end

function [rr, ci_low, ci_high, p_val] = prepost_ratio_ci_ext(stats, posneg_idx, alpha)
    [post_count, post_total, pre_count, pre_total] = safe_counts_for_log_ratio_ext( ...
        stats.post_counts(posneg_idx), stats.post_total, stats.pre_counts(posneg_idx), stats.pre_total);
    if pre_total == 0 || post_total == 0 || pre_count <= 0 || post_count <= 0
        rr = NaN; ci_low = NaN; ci_high = NaN; p_val = NaN; return;
    end
    p_post = post_count / post_total;
    p_pre = pre_count / pre_total;
    rr = p_post / p_pre;
    log_rr = log(rr);
    se = sqrt(1 / post_count - 1 / post_total + 1 / pre_count - 1 / pre_total);
    if se == 0
        p_val = double(log_rr == 0);
        ci_low = rr; ci_high = rr; return;
    end
    zcrit = norminv(1 - alpha / 2);
    ci_low = exp(log_rr - zcrit * se);
    ci_high = exp(log_rr + zcrit * se);
    p_val = twoSidedNormalPValue(log_rr / se);
end

function [value, ci_low, ci_high, p_val] = did_metric_ci_ext(stats_a, stats_b, posneg_idx, metric_mode, alpha)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            [value, ci_low, ci_high, p_val] = prepost_effect_difference_ci(stats_a, stats_b, posneg_idx, alpha);
        case {'ratio', 'post/pre', 'post_over_pre'}
            [value, ci_low, ci_high, p_val] = did_ratio_ci_ext(stats_a, stats_b, posneg_idx, alpha);
        otherwise
            error('Unknown metric_mode: %s', metric_mode);
    end
end

function [ror, ci_low, ci_high, p_val] = did_ratio_ci_ext(stats_a, stats_b, posneg_idx, alpha)
    [a_post_count, a_post_total, a_pre_count, a_pre_total] = safe_counts_for_log_ratio_ext(stats_a.post_counts(posneg_idx), stats_a.post_total, stats_a.pre_counts(posneg_idx), stats_a.pre_total);
    [b_post_count, b_post_total, b_pre_count, b_pre_total] = safe_counts_for_log_ratio_ext(stats_b.post_counts(posneg_idx), stats_b.post_total, stats_b.pre_counts(posneg_idx), stats_b.pre_total);
    if any([a_post_count, a_pre_count, b_post_count, b_pre_count] <= 0) || any([a_post_total, a_pre_total, b_post_total, b_pre_total] == 0)
        ror = NaN; ci_low = NaN; ci_high = NaN; p_val = NaN; return;
    end
    log_rr_a = log((a_post_count / a_post_total) / (a_pre_count / a_pre_total));
    log_rr_b = log((b_post_count / b_post_total) / (b_pre_count / b_pre_total));
    log_ror = log_rr_a - log_rr_b;
    se = sqrt(1 / a_post_count - 1 / a_post_total + 1 / a_pre_count - 1 / a_pre_total + ...
              1 / b_post_count - 1 / b_post_total + 1 / b_pre_count - 1 / b_pre_total);
    ror = exp(log_ror);
    if se == 0
        p_val = double(log_ror == 0);
        ci_low = ror; ci_high = ror; return;
    end
    zcrit = norminv(1 - alpha / 2);
    ci_low = exp(log_ror - zcrit * se);
    ci_high = exp(log_ror + zcrit * se);
    p_val = twoSidedNormalPValue(log_ror / se);
end

function [num1, den1, num0, den0] = safe_counts_for_log_ratio_ext(num1, den1, num0, den0)
    if num1 == 0 || num0 == 0 || num1 == den1 || num0 == den0
        num1 = num1 + 0.5;
        num0 = num0 + 0.5;
        den1 = den1 + 1;
        den0 = den0 + 1;
    end
end

function plot_metric_grouped_bars_ext(tile, values, low_errs, high_errs, bar_sig_labels, did_sig_labels, condition_labels, kernel_labels, bar_offset, bar_width, y_limit, metric_mode)
    n_kernel = size(values, 1);
    x = 1:n_kernel;
    hold(tile, 'on');
    bar(tile, x - bar_offset, values(:, 1), 'FaceColor', 'b', 'BarWidth', bar_width, 'DisplayName', condition_labels{1});
    bar(tile, x + bar_offset, values(:, 2), 'FaceColor', 'r', 'BarWidth', bar_width, 'DisplayName', condition_labels{2});
    errorbar(tile, x - bar_offset, values(:, 1), low_errs(:, 1), high_errs(:, 1), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    errorbar(tile, x + bar_offset, values(:, 2), low_errs(:, 2), high_errs(:, 2), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    ylim(tile, y_limit);
    draw_metric_reference_line_ext(tile, metric_mode);
    yr = diff(y_limit);
    label_offset = 0.025 * yr;
    bracket_offset = 0.11 * yr;
    bracket_height = 0.025 * yr;
    for kernel_i = 1:n_kernel
        xpos = [x(kernel_i) - bar_offset, x(kernel_i) + bar_offset];
        for injection_idx = 1:2
            y_bar = values(kernel_i, injection_idx) + high_errs(kernel_i, injection_idx) + label_offset;
            text(tile, xpos(injection_idx), y_bar, bar_sig_labels{kernel_i, injection_idx}, ...
                'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'HandleVisibility', 'off');
        end
        y_bracket = max(values(kernel_i, :) + high_errs(kernel_i, :), [], 'omitnan') + bracket_offset;
        draw_sig_bracket_ext(tile, xpos(1), xpos(2), y_bracket, bracket_height, did_sig_labels{kernel_i});
    end
    hold(tile, 'off');
    xticks(tile, x);
    xticklabels(tile, kernel_labels);
    legend(tile, condition_labels, 'Location', 'Best');
end

function plot_metric_single_bars_ext(tile, values, low_errs, high_errs, sig_labels, kernel_labels, y_limit, metric_mode, bar_label)
    n_kernel = numel(values);
    x = 1:n_kernel;
    hold(tile, 'on');
    bar(tile, x, values, 'BarWidth', 0.45, 'DisplayName', char(bar_label));
    errorbar(tile, x, values, low_errs, high_errs, 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    ylim(tile, y_limit);
    draw_metric_reference_line_ext(tile, metric_mode);
    yr = diff(y_limit);
    label_offset = 0.035 * yr;
    for kernel_i = 1:n_kernel
        y_text = values(kernel_i) + high_errs(kernel_i) + label_offset;
        text(tile, x(kernel_i), y_text, sig_labels{kernel_i}, 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'HandleVisibility', 'off');
    end
    hold(tile, 'off');
    xticks(tile, x);
    xticklabels(tile, kernel_labels);
end

function plot_prepost_grouped_bars_ext(tile, ratios, low_errs, high_errs, sig_labels, condition_labels, kernel_labels, bar_offset, bar_width, y_limit)
    n_kernel = size(ratios, 1);
    x = 1:n_kernel;
    hold(tile, 'on');
    bar(tile, x - bar_offset, ratios(:, 1), 'FaceColor', 'b', 'BarWidth', bar_width, 'DisplayName', condition_labels{1});
    bar(tile, x + bar_offset, ratios(:, 2), 'FaceColor', 'r', 'BarWidth', bar_width, 'DisplayName', condition_labels{2});
    errorbar(tile, x - bar_offset, ratios(:, 1), low_errs(:, 1), high_errs(:, 1), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    errorbar(tile, x + bar_offset, ratios(:, 2), low_errs(:, 2), high_errs(:, 2), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    ylim(tile, y_limit);
    yr = diff(y_limit);
    bracket_height = 0.025 * yr;
    bracket_offset = 0.06 * yr;
    for kernel_i = 1:n_kernel
        y_bracket = max(ratios(kernel_i, :) + high_errs(kernel_i, :), [], 'omitnan') + bracket_offset;
        draw_sig_bracket_ext(tile, x(kernel_i) - bar_offset, x(kernel_i) + bar_offset, y_bracket, bracket_height, sig_labels{kernel_i});
    end
    hold(tile, 'off');
    xticks(tile, x);
    xticklabels(tile, kernel_labels);
    legend(tile, condition_labels, 'Location', 'Best');
end

function draw_sig_bracket_ext(tile, x1, x2, y, h, label)
    plot(tile, [x1, x1, x2, x2], [y, y + h, y + h, y], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(tile, mean([x1, x2]), y + h, label, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold', 'HandleVisibility', 'off');
end

function y_limit = common_ylim_for_metric_ext(values, low_errs, high_errs, metric_mode)
    lower = values(:) - low_errs(:);
    upper = values(:) + high_errs(:);
    switch lower(metric_mode)
        case {'ratio', 'post/pre', 'post_over_pre'}
            ref = 1;
        otherwise
            ref = 0;
    end
    vals = [lower; upper; ref];
    vals = vals(isfinite(vals));
    if isempty(vals)
        y_limit = [ref - 1, ref + 1];
        return;
    end
    y_min = min(vals);
    y_max = max(vals);
    yr = y_max - y_min;
    if yr == 0
        yr = max(abs([y_min, y_max, 1]));
    end
    if yr == 0
        yr = 1;
    end
    y_limit = [y_min - 0.18 * yr, y_max + 0.45 * yr];
    if strcmpi(metric_mode, 'ratio') || strcmpi(metric_mode, 'post/pre')
        y_limit(1) = max(0, y_limit(1));
    end
end

function y_limit = common_ylim_percent_ext(values, low_errs, high_errs)
    vals = [values(:) - low_errs(:); values(:) + high_errs(:); 0];
    vals = vals(isfinite(vals));
    if isempty(vals)
        y_limit = [0, 1];
        return;
    end
    y_min = min(vals);
    y_max = max(vals);
    yr = y_max - y_min;
    if yr == 0, yr = max(abs(vals)); end
    if yr == 0, yr = 1; end
    y_limit = [max(0, y_min - 0.10 * yr), y_max + 0.35 * yr];
end

function draw_metric_reference_line_ext(tile, metric_mode)
    xl = xlim(tile);
    hold(tile, 'on');
    switch lower(metric_mode)
        case {'ratio', 'post/pre', 'post_over_pre'}
            plot(tile, xl, [1, 1], 'k--', 'LineWidth', 0.9, 'HandleVisibility', 'off');
        otherwise
            plot(tile, xl, [0, 0], 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    end
    xlim(tile, xl);
end

function label = metric_axis_label_ext(metric_mode)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            label = 'Post - Pre significant J (percentage points)';
        case {'ratio', 'post/pre', 'post_over_pre'}
            label = 'Post / Pre significant J ratio';
        otherwise
            label = metric_mode;
    end
end

function label = did_axis_label_ext(metric_mode)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            label = '(Muscimol Post-Pre) - (Saline Post-Pre) (percentage points)';
        case {'ratio', 'post/pre', 'post_over_pre'}
            label = '(Muscimol Post/Pre) / (Saline Post/Pre)';
        otherwise
            label = metric_mode;
    end
end

function stub = metric_output_stub_ext(metric_mode)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            stub = 'post_minus_pre';
        case {'ratio', 'post/pre', 'post_over_pre'}
            stub = 'post_over_pre';
        otherwise
            stub = regexprep(metric_mode, '[^A-Za-z0-9_-]', '_');
    end
end
