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

%% Figure 3: pooled Post - Pre effect. Bars = Saline/Muscimol.
plot_prepost_effect_comparison_pooled(root, mt.GLM, animal_filter, states, state_labels, ...
    posneg_labels, kernel_indices, kernel_labels, effect_injection_types, effect_injection_labels, ...
    alignment, resting_dur_threshold, bar_offset, bar_width, effect_y_limit, effect_alpha);

%% Figure 4: pooled difference of effects: Muscimol effect - Saline effect.
plot_muscimol_minus_saline_effect_pooled(root, mt.GLM, animal_filter, states, state_labels, ...
    posneg_labels, kernel_indices, kernel_labels, effect_injection_types, ...
    alignment, resting_dur_threshold, effect_diff_y_limit, effect_alpha);

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
