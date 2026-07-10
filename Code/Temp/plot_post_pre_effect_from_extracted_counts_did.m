%% plot_post_pre_effect_from_extracted_counts_did.m
% Hard-coded counts extracted from the provided bar plots.
% Plot pooled Post - Pre effect for Saline vs Muscimol.
%
% Statistics match J_bars_kernels_grouped_by_sign_effects_pooled.m:
%   Bar height = 100 * (Post proportion - Pre proportion)
%   Bar error  = normal-approx CI for a difference of two pooled proportions
%   Star label = difference-in-differences z-test:
%                (Muscimol Post-Pre) - (Saline Post-Pre) = 0
%
% Dimensions for counts/totals:
%   state   x sign   x prepost x kernel
%   1 Open    1 Pos     1 Pre     1 K1
%   2 Close   2 Neg     2 Post    2 K2, 3 K3

clear;

%% Data
kernels = 1:3;
kernel_labels = arrayfun(@(k) sprintf('Kernel %d', k), kernels, 'UniformOutput', false);
injection_labels = {'Saline', 'Muscimol'};
state_labels = {'Eyes Open', 'Eyes Closed'};
sign_labels = {'Positive J', 'Negative J'};
alpha = 0.05;

make_difference_figure = true;
export_figures = true;
save_folder = fullfile(pwd, 'Figures');

counts = struct();
totals = struct();

% Saline counts extracted from the Saline figure.
totals.Saline = 7442 * ones(2, 2, 2, 3);
counts.Saline = nan(2, 2, 2, 3);

% Saline, Positive J
counts.Saline(1, 1, 1, :) = [526, 443, 339]; % Open Pre
counts.Saline(1, 1, 2, :) = [632, 519, 398]; % Open Post
counts.Saline(2, 1, 1, :) = [876, 783, 489]; % Close Pre
counts.Saline(2, 1, 2, :) = [857, 753, 526]; % Close Post

% Saline, Negative J
counts.Saline(1, 2, 1, :) = [324, 307, 219]; % Open Pre
counts.Saline(1, 2, 2, :) = [396, 337, 226]; % Open Post
counts.Saline(2, 2, 1, :) = [308, 272, 291]; % Close Pre
counts.Saline(2, 2, 2, :) = [418, 407, 337]; % Close Post

% Muscimol counts extracted from the Muscimol figure.
totals.Muscimol = 10834 * ones(2, 2, 2, 3);
counts.Muscimol = nan(2, 2, 2, 3);

% Muscimol, Positive J
counts.Muscimol(1, 1, 1, :) = [619, 488, 286];  % Open Pre
counts.Muscimol(1, 1, 2, :) = [900, 745, 500];  % Open Post
counts.Muscimol(2, 1, 1, :) = [1055, 848, 556]; % Close Pre
counts.Muscimol(2, 1, 2, :) = [1154, 959, 685]; % Close Post

% Muscimol, Negative J
counts.Muscimol(1, 2, 1, :) = [331, 278, 130]; % Open Pre
counts.Muscimol(1, 2, 2, :) = [496, 398, 252]; % Open Post
counts.Muscimol(2, 2, 1, :) = [437, 341, 238]; % Close Pre
counts.Muscimol(2, 2, 2, :) = [685, 606, 454]; % Close Post

%% Compute effects and difference-in-differences statistics
n_state = numel(state_labels);
n_sign = numel(sign_labels);
n_kernel = numel(kernels);
n_injection = numel(injection_labels);

effect_pct = nan(n_state, n_sign, n_kernel, n_injection);
low_err_pct = nan(n_state, n_sign, n_kernel, n_injection);
high_err_pct = nan(n_state, n_sign, n_kernel, n_injection);

diff_effect_pct = nan(n_state, n_sign, n_kernel);
diff_low_err_pct = nan(n_state, n_sign, n_kernel);
diff_high_err_pct = nan(n_state, n_sign, n_kernel);
diff_p = nan(n_state, n_sign, n_kernel);
diff_stars = cell(n_state, n_sign, n_kernel);

for state_idx = 1:n_state
    for sign_idx = 1:n_sign
        for kernel_idx = 1:n_kernel
            stats_by_injection = cell(1, n_injection);

            for inj_idx = 1:n_injection
                inj = injection_labels{inj_idx};
                stats = make_prepost_stats_from_counts(counts.(inj), totals.(inj), state_idx, sign_idx, kernel_idx);
                stats_by_injection{inj_idx} = stats;

                [effect, ci_low, ci_high] = prepost_effect_ci(stats, alpha);
                effect_pct(state_idx, sign_idx, kernel_idx, inj_idx) = effect;
                low_err_pct(state_idx, sign_idx, kernel_idx, inj_idx) = effect - ci_low;
                high_err_pct(state_idx, sign_idx, kernel_idx, inj_idx) = ci_high - effect;
            end

            % Difference-in-differences: Muscimol effect - Saline effect.
            [diff_effect, diff_ci_low, diff_ci_high, p_val] = prepost_effect_difference_ci( ...
                stats_by_injection{2}, stats_by_injection{1}, alpha);

            diff_effect_pct(state_idx, sign_idx, kernel_idx) = diff_effect;
            diff_low_err_pct(state_idx, sign_idx, kernel_idx) = diff_effect - diff_ci_low;
            diff_high_err_pct(state_idx, sign_idx, kernel_idx) = diff_ci_high - diff_effect;
            diff_p(state_idx, sign_idx, kernel_idx) = p_val;
            diff_stars{state_idx, sign_idx, kernel_idx} = p_to_sig_label(p_val);
        end
    end
end

%% Figure 1: Post - Pre effect, Saline vs Muscimol
bar_width = 0.30;
bar_offset = 0.18;

effect_fig = figure('Color', 'w', 'Position', [100, 100, 1300, 900]);
t = tiledlayout(effect_fig, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for state_idx = 1:n_state
    for sign_idx = 1:n_sign
        tile_idx = (state_idx - 1) * 2 + sign_idx;
        ax = nexttile(t, tile_idx);

        effects = squeeze(effect_pct(state_idx, sign_idx, :, :));
        low_errs = squeeze(low_err_pct(state_idx, sign_idx, :, :));
        high_errs = squeeze(high_err_pct(state_idx, sign_idx, :, :));
        sig_labels = squeeze(diff_stars(state_idx, sign_idx, :));

        plot_effect_grouped_bars(ax, effects, low_errs, high_errs, sig_labels, ...
            injection_labels, kernel_labels, bar_offset, bar_width, []);
        ylabel(ax, 'Post - Pre significant J (percentage points)');
        title(ax, sprintf('%s, %s', state_labels{state_idx}, sign_labels{sign_idx}), 'FontWeight', 'bold');
        draw_zero_line(ax);
    end
end

sgtitle(effect_fig, 'Pooled Post - Pre effect: Saline vs Muscimol');

%% Figure 2: Muscimol effect - Saline effect
if make_difference_figure
    diff_fig = figure('Color', 'w', 'Position', [100, 100, 1300, 900]);
    t2 = tiledlayout(diff_fig, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for state_idx = 1:n_state
        for sign_idx = 1:n_sign
            tile_idx = (state_idx - 1) * 2 + sign_idx;
            ax = nexttile(t2, tile_idx);

            values = squeeze(diff_effect_pct(state_idx, sign_idx, :));
            low_errs = squeeze(diff_low_err_pct(state_idx, sign_idx, :));
            high_errs = squeeze(diff_high_err_pct(state_idx, sign_idx, :));
            sig_labels = squeeze(diff_stars(state_idx, sign_idx, :));

            plot_effect_difference_bars(ax, values, low_errs, high_errs, sig_labels, kernel_labels, []);
            ylabel(ax, '(Muscimol Post-Pre) - (Saline Post-Pre) (percentage points)');
            title(ax, sprintf('%s, %s', state_labels{state_idx}, sign_labels{sign_idx}), 'FontWeight', 'bold');
            draw_zero_line(ax);
        end
    end

    sgtitle(diff_fig, 'Pooled difference of effects: Muscimol - Saline');
else
    diff_fig = [];
end

%% Summary table
summary_rows = {};
for state_idx = 1:n_state
    for sign_idx = 1:n_sign
        for kernel_idx = 1:n_kernel
            summary_rows(end+1, :) = { ... %#ok<SAGROW>
                state_labels{state_idx}, ...
                sign_labels{sign_idx}, ...
                kernel_idx, ...
                effect_pct(state_idx, sign_idx, kernel_idx, 1), ...
                low_err_pct(state_idx, sign_idx, kernel_idx, 1), ...
                high_err_pct(state_idx, sign_idx, kernel_idx, 1), ...
                effect_pct(state_idx, sign_idx, kernel_idx, 2), ...
                low_err_pct(state_idx, sign_idx, kernel_idx, 2), ...
                high_err_pct(state_idx, sign_idx, kernel_idx, 2), ...
                diff_effect_pct(state_idx, sign_idx, kernel_idx), ...
                diff_low_err_pct(state_idx, sign_idx, kernel_idx), ...
                diff_high_err_pct(state_idx, sign_idx, kernel_idx), ...
                diff_p(state_idx, sign_idx, kernel_idx), ...
                diff_stars{state_idx, sign_idx, kernel_idx}};
        end
    end
end

summary_table = cell2table(summary_rows, 'VariableNames', ...
    {'State', 'Sign', 'Kernel', ...
     'Saline_PostMinusPre_pctpt', 'Saline_LowErr_pctpt', 'Saline_HighErr_pctpt', ...
     'Muscimol_PostMinusPre_pctpt', 'Muscimol_LowErr_pctpt', 'Muscimol_HighErr_pctpt', ...
     'MuscimolMinusSaline_pctpt', 'MuscimolMinusSaline_LowErr_pctpt', ...
     'MuscimolMinusSaline_HighErr_pctpt', 'DifferenceInDifference_P', 'Stars'});
disp(summary_table);

%% Optional export
if export_figures
    if ~isfolder(save_folder)
        mkdir(save_folder);
    end
    exportgraphics(effect_fig, fullfile(save_folder, 'post_pre_effect_extracted_counts_did_preview.jpg'), ...
        'ContentType', 'image', 'Resolution', 300, 'BackgroundColor', 'white');
    exportgraphics(effect_fig, fullfile(save_folder, 'post_pre_effect_extracted_counts_did.pdf'), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    if make_difference_figure && ~isempty(diff_fig)
        exportgraphics(diff_fig, fullfile(save_folder, 'muscimol_minus_saline_effect_extracted_counts_preview.jpg'), ...
            'ContentType', 'image', 'Resolution', 300, 'BackgroundColor', 'white');
        exportgraphics(diff_fig, fullfile(save_folder, 'muscimol_minus_saline_effect_extracted_counts.pdf'), ...
            'ContentType', 'vector', 'BackgroundColor', 'white');
    end
end

function stats = make_prepost_stats_from_counts(count_arr, total_arr, state_idx, sign_idx, kernel_idx)
    stats = struct();
    stats.pre_counts = count_arr(state_idx, sign_idx, 1, kernel_idx);
    stats.post_counts = count_arr(state_idx, sign_idx, 2, kernel_idx);
    stats.pre_total = total_arr(state_idx, sign_idx, 1, kernel_idx);
    stats.post_total = total_arr(state_idx, sign_idx, 2, kernel_idx);
end

function [effect_pct, ci_low_pct, ci_high_pct] = prepost_effect_ci(stats, alpha)
    if stats.pre_total == 0 || stats.post_total == 0
        effect_pct = NaN;
        ci_low_pct = NaN;
        ci_high_pct = NaN;
        return;
    end

    p_pre = stats.pre_counts / stats.pre_total;
    p_post = stats.post_counts / stats.post_total;
    effect = p_post - p_pre;
    se = sqrt(p_post * (1 - p_post) / stats.post_total + ...
              p_pre * (1 - p_pre) / stats.pre_total);
    zcrit = norminv(1 - alpha / 2);

    effect_pct = 100 * effect;
    ci_low_pct = 100 * (effect - zcrit * se);
    ci_high_pct = 100 * (effect + zcrit * se);
end

function [diff_pct, ci_low_pct, ci_high_pct, p_val] = prepost_effect_difference_ci(stats_a, stats_b, alpha)
    % stats_a - stats_b, where each effect is Post - Pre.
    % In this script: stats_a = Muscimol, stats_b = Saline.
    if stats_a.pre_total == 0 || stats_a.post_total == 0 || stats_b.pre_total == 0 || stats_b.post_total == 0
        diff_pct = NaN;
        ci_low_pct = NaN;
        ci_high_pct = NaN;
        p_val = NaN;
        return;
    end

    a_pre = stats_a.pre_counts / stats_a.pre_total;
    a_post = stats_a.post_counts / stats_a.post_total;
    b_pre = stats_b.pre_counts / stats_b.pre_total;
    b_post = stats_b.post_counts / stats_b.post_total;

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

function plot_effect_grouped_bars(ax, effects, low_errs, high_errs, sig_labels, condition_labels, kernel_labels, bar_offset, bar_width, y_limit)
    n_kernel = size(effects, 1);
    x = 1:n_kernel;

    hold(ax, 'on');
    bar(ax, x - bar_offset, effects(:, 1), 'FaceColor', 'b', 'BarWidth', bar_width, 'DisplayName', condition_labels{1});
    bar(ax, x + bar_offset, effects(:, 2), 'FaceColor', 'r', 'BarWidth', bar_width, 'DisplayName', condition_labels{2});

    errorbar(ax, x - bar_offset, effects(:, 1), low_errs(:, 1), high_errs(:, 1), ...
        'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    errorbar(ax, x + bar_offset, effects(:, 2), low_errs(:, 2), high_errs(:, 2), ...
        'k', 'LineStyle', 'none', 'HandleVisibility', 'off');

    set_effect_ylim(ax, effects, low_errs, high_errs, y_limit);
    yl = ylim(ax);
    y_offset = 0.04 * diff(yl);
    for kernel_i = 1:n_kernel
        y_text = max([effects(kernel_i, :) + high_errs(kernel_i, :), 0], [], 'omitnan') + y_offset;
        if ~isfinite(y_text)
            y_text = yl(2) - y_offset;
        end
        text(ax, x(kernel_i), y_text, sig_labels{kernel_i}, ...
            'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    end

    hold(ax, 'off');
    xticks(ax, x);
    xticklabels(ax, kernel_labels);
    legend(ax, condition_labels, 'Location', 'Best');
    box(ax, 'off');
end

function plot_effect_difference_bars(ax, values, low_errs, high_errs, sig_labels, kernel_labels, y_limit)
    n_kernel = numel(values);
    x = 1:n_kernel;

    hold(ax, 'on');
    bar(ax, x, values, 'FaceColor', [0.4, 0.4, 0.4], 'BarWidth', 0.45);
    errorbar(ax, x, values, low_errs, high_errs, 'k', 'LineStyle', 'none');

    set_effect_ylim(ax, values(:), low_errs(:), high_errs(:), y_limit);
    yl = ylim(ax);
    y_offset = 0.04 * diff(yl);
    for kernel_i = 1:n_kernel
        y_text = max([values(kernel_i) + high_errs(kernel_i), 0], [], 'omitnan') + y_offset;
        if ~isfinite(y_text)
            y_text = yl(2) - y_offset;
        end
        text(ax, x(kernel_i), y_text, sig_labels{kernel_i}, ...
            'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    end

    hold(ax, 'off');
    xticks(ax, x);
    xticklabels(ax, kernel_labels);
    box(ax, 'off');
end

function set_effect_ylim(ax, values, low_errs, high_errs, y_limit)
    if ~isempty(y_limit)
        ylim(ax, y_limit);
        return;
    end

    lower = values(:) - low_errs(:);
    upper = values(:) + high_errs(:);
    vals = [lower; upper; 0];
    vals = vals(isfinite(vals));
    if isempty(vals)
        ylim(ax, [-1, 1]);
        return;
    end
    yr = max(vals) - min(vals);
    if yr == 0
        yr = max(abs(vals));
    end
    if yr == 0
        yr = 1;
    end
    ylim(ax, [min(vals) - 0.15 * yr, max(vals) + 0.25 * yr]);
end

function draw_zero_line(ax)
    xl = xlim(ax);
    hold(ax, 'on');
    plot(ax, xl, [0, 0], 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    xlim(ax, xl);
    hold(ax, 'off');
end
