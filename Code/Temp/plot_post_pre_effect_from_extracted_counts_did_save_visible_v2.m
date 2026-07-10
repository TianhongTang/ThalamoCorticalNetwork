%% plot_post_pre_effect_from_extracted_counts_extended.m
% Hard-coded counts extracted from bar plots.
% Generates:
%   1) Saline + Muscimol Post-Pre or Post/Pre metric with per-bar Pre-vs-Post tests
%      and higher Muscimol-vs-Saline DID bracket.
%   2) Saline-only metric figure.
%   3) Muscimol-only metric figure.
%   4) Saline Pre vs Post significant J % figure.
%   5) Muscimol Pre vs Post significant J % figure.
%
% metric_mode:
%   'difference' : Post - Pre in percentage points.
%                  Bar CI uses normal approximation for difference of two proportions.
%                  Per-bar p uses two-proportion z-test.
%                  DID compares (Muscimol Post-Pre) - (Saline Post-Pre).
%   'ratio'      : Post / Pre proportion ratio.
%                  Bar CI and p use log risk-ratio approximation.
%                  DID compares ratio-of-ratios: (Muscimol Post/Pre) / (Saline Post/Pre).

clear;

%% Get root folder
script_path = mfilename('fullpath');
if isempty(script_path)
    root = pwd;
else
    root = fileparts(script_path);
end
for root_search_i = 1:8
    if isfolder(fullfile(root, 'Data'))
        break;
    end
    parent_root = fileparts(root);
    if strcmp(parent_root, root)
        break;
    end
    root = parent_root;
end
if ~isfolder(fullfile(root, 'Data')) && isfolder(fullfile(pwd, 'Data'))
    root = pwd;
end

%% Parameters
metric_mode = 'ratio'; % 'difference' or 'ratio'.
alpha = 0.05;

make_combined_figure = true;
make_single_injection_figures = true;
make_prepost_percent_figures = true;
make_difference_only_figure = false; % optional diagnostic figure for DID values only.

figure_output_mode = 'save'; % 'visible' or 'save'.
[figure_visible, export_figures] = parse_figure_output_mode(figure_output_mode);
set(0, 'DefaultFigureVisible', figure_visible);
save_folder = fullfile(root, 'Figures', 'Paper', 'J_bars_extracted');

bar_width = 0.30;
bar_offset = 0.18;

kernels = 1:3;
kernel_labels = arrayfun(@(k) sprintf('Kernel %d', k), kernels, 'UniformOutput', false);
injection_labels = {'Saline', 'Muscimol'};
state_labels = {'Eyes Open', 'Eyes Closed'};
sign_labels = {'Positive J', 'Negative J'};
prepost_labels = {'Pre', 'Post'};

%% Data
% Dimensions: state x sign x prepost x kernel.
% state: 1 Open, 2 Close.
% sign: 1 Positive, 2 Negative.
% prepost: 1 Pre, 2 Post.
% kernel: 1/2/3.
counts = struct();
totals = struct();

% Saline counts extracted from the Saline figure.
totals.Saline = 7442 * ones(2, 2, 2, 3);
counts.Saline = nan(2, 2, 2, 3);

counts.Saline(1, 1, 1, :) = [526, 443, 339]; % Open Positive Pre
counts.Saline(1, 1, 2, :) = [632, 519, 398]; % Open Positive Post
counts.Saline(2, 1, 1, :) = [876, 783, 489]; % Close Positive Pre
counts.Saline(2, 1, 2, :) = [857, 753, 526]; % Close Positive Post

counts.Saline(1, 2, 1, :) = [324, 307, 219]; % Open Negative Pre
counts.Saline(1, 2, 2, :) = [396, 337, 226]; % Open Negative Post
counts.Saline(2, 2, 1, :) = [308, 272, 291]; % Close Negative Pre
counts.Saline(2, 2, 2, :) = [418, 407, 337]; % Close Negative Post

% Muscimol counts extracted from the Muscimol figure.
totals.Muscimol = 10834 * ones(2, 2, 2, 3);
counts.Muscimol = nan(2, 2, 2, 3);

counts.Muscimol(1, 1, 1, :) = [619, 488, 286];  % Open Positive Pre
counts.Muscimol(1, 1, 2, :) = [900, 745, 500];  % Open Positive Post
counts.Muscimol(2, 1, 1, :) = [1055, 848, 556]; % Close Positive Pre
counts.Muscimol(2, 1, 2, :) = [1154, 959, 685]; % Close Positive Post

counts.Muscimol(1, 2, 1, :) = [331, 278, 130]; % Open Negative Pre
counts.Muscimol(1, 2, 2, :) = [496, 398, 252]; % Open Negative Post
counts.Muscimol(2, 2, 1, :) = [437, 341, 238]; % Close Negative Pre
counts.Muscimol(2, 2, 2, :) = [685, 606, 454]; % Close Negative Post

%% Compute metric and statistics
n_state = numel(state_labels);
n_sign = numel(sign_labels);
n_kernel = numel(kernels);
n_injection = numel(injection_labels);

metric_value = nan(n_state, n_sign, n_kernel, n_injection);
low_err = nan(n_state, n_sign, n_kernel, n_injection);
high_err = nan(n_state, n_sign, n_kernel, n_injection);
prepost_p = nan(n_state, n_sign, n_kernel, n_injection);
prepost_stars = cell(n_state, n_sign, n_kernel, n_injection);

did_value = nan(n_state, n_sign, n_kernel);
did_low_err = nan(n_state, n_sign, n_kernel);
did_high_err = nan(n_state, n_sign, n_kernel);
did_p = nan(n_state, n_sign, n_kernel);
did_stars = cell(n_state, n_sign, n_kernel);

for state_idx = 1:n_state
    for sign_idx = 1:n_sign
        for kernel_idx = 1:n_kernel
            stats_by_injection = cell(1, n_injection);
            for inj_idx = 1:n_injection
                inj = injection_labels{inj_idx};
                stats = make_prepost_stats_from_counts(counts.(inj), totals.(inj), state_idx, sign_idx, kernel_idx);
                stats_by_injection{inj_idx} = stats;

                [v, ci_low, ci_high, p_val] = prepost_metric_ci(stats, metric_mode, alpha);
                metric_value(state_idx, sign_idx, kernel_idx, inj_idx) = v;
                low_err(state_idx, sign_idx, kernel_idx, inj_idx) = v - ci_low;
                high_err(state_idx, sign_idx, kernel_idx, inj_idx) = ci_high - v;
                prepost_p(state_idx, sign_idx, kernel_idx, inj_idx) = p_val;
                prepost_stars{state_idx, sign_idx, kernel_idx, inj_idx} = p_to_sig_label(p_val);
            end

            % DID: Muscimol metric vs Saline metric.
            [v, ci_low, ci_high, p_val] = did_metric_ci(stats_by_injection{2}, stats_by_injection{1}, metric_mode, alpha);
            did_value(state_idx, sign_idx, kernel_idx) = v;
            did_low_err(state_idx, sign_idx, kernel_idx) = v - ci_low;
            did_high_err(state_idx, sign_idx, kernel_idx) = ci_high - v;
            did_p(state_idx, sign_idx, kernel_idx) = p_val;
            did_stars{state_idx, sign_idx, kernel_idx} = p_to_sig_label(p_val);
        end
    end
end

%% Figures
metric_label = metric_axis_label(metric_mode);
metric_stub = metric_output_stub(metric_mode);

if make_combined_figure
    y_lim_combined = common_ylim_for_metric(metric_value, low_err, high_err, metric_mode);
    combined_fig = figure('Color', 'w', 'Position', [100, 100, 1300, 900]);
    t = tiledlayout(combined_fig, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for state_idx = 1:n_state
        for sign_idx = 1:n_sign
            ax = nexttile(t, (state_idx - 1) * 2 + sign_idx);
            values = squeeze(metric_value(state_idx, sign_idx, :, :));
            lows = squeeze(low_err(state_idx, sign_idx, :, :));
            highs = squeeze(high_err(state_idx, sign_idx, :, :));
            bar_labels = squeeze(prepost_stars(state_idx, sign_idx, :, :));
            bracket_labels = squeeze(did_stars(state_idx, sign_idx, :));

            plot_metric_grouped_bars(ax, values, lows, highs, bar_labels, bracket_labels, ...
                injection_labels, kernel_labels, bar_offset, bar_width, y_lim_combined, metric_mode);
            ylabel(ax, metric_label);
            title(ax, sprintf('%s, %s', state_labels{state_idx}, sign_labels{sign_idx}), 'FontWeight', 'bold');
        end
    end
    sgtitle(combined_fig, sprintf('Pooled %s: Saline vs Muscimol', metric_label), 'Interpreter', 'none');
else
    combined_fig = [];
end

single_figs = cell(1, n_injection);
if make_single_injection_figures
    for inj_idx = 1:n_injection
        inj = injection_labels{inj_idx};
        y_lim_single = common_ylim_for_metric(metric_value(:, :, :, inj_idx), low_err(:, :, :, inj_idx), high_err(:, :, :, inj_idx), metric_mode);
        f = figure('Color', 'w', 'Position', [100, 100, 1300, 900]);
        t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
        for state_idx = 1:n_state
            for sign_idx = 1:n_sign
                ax = nexttile(t, (state_idx - 1) * 2 + sign_idx);
                values = squeeze(metric_value(state_idx, sign_idx, :, inj_idx));
                lows = squeeze(low_err(state_idx, sign_idx, :, inj_idx));
                highs = squeeze(high_err(state_idx, sign_idx, :, inj_idx));
                labels = squeeze(prepost_stars(state_idx, sign_idx, :, inj_idx));
                plot_metric_single_bars(ax, values, lows, highs, labels, kernel_labels, y_lim_single, metric_mode, inj);
                ylabel(ax, metric_label);
                title(ax, sprintf('%s, %s', state_labels{state_idx}, sign_labels{sign_idx}), 'FontWeight', 'bold');
            end
        end
        sgtitle(f, sprintf('Pooled %s: %s only', metric_label, inj), 'Interpreter', 'none');
        single_figs{inj_idx} = f;
    end
end

prepost_figs = cell(1, n_injection);
if make_prepost_percent_figures
    for inj_idx = 1:n_injection
        inj = injection_labels{inj_idx};
        prepost_figs{inj_idx} = plot_prepost_percent_figure(counts.(inj), totals.(inj), inj, ...
            state_labels, sign_labels, kernel_labels, prepost_labels, alpha, bar_offset, bar_width);
    end
end

if make_difference_only_figure
    diff_fig = figure('Color', 'w', 'Position', [100, 100, 1300, 900]);
    t2 = tiledlayout(diff_fig, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    y_lim_did = common_ylim_for_metric(did_value, did_low_err, did_high_err, metric_mode);
    for state_idx = 1:n_state
        for sign_idx = 1:n_sign
            ax = nexttile(t2, (state_idx - 1) * 2 + sign_idx);
            values = squeeze(did_value(state_idx, sign_idx, :));
            lows = squeeze(did_low_err(state_idx, sign_idx, :));
            highs = squeeze(did_high_err(state_idx, sign_idx, :));
            labels = squeeze(did_stars(state_idx, sign_idx, :));
            plot_metric_single_bars(ax, values, lows, highs, labels, kernel_labels, y_lim_did, metric_mode, 'DID');
            ylabel(ax, did_axis_label(metric_mode));
            title(ax, sprintf('%s, %s', state_labels{state_idx}, sign_labels{sign_idx}), 'FontWeight', 'bold');
        end
    end
    sgtitle(diff_fig, sprintf('Pooled DID: Muscimol vs Saline (%s)', metric_label), 'Interpreter', 'none');
else
    diff_fig = [];
end

%% Summary table
summary_rows = {};
for state_idx = 1:n_state
    for sign_idx = 1:n_sign
        for kernel_idx = 1:n_kernel
            summary_rows(end+1, :) = { ... %#ok<SAGROW>
                state_labels{state_idx}, sign_labels{sign_idx}, kernel_idx, metric_mode, ...
                metric_value(state_idx, sign_idx, kernel_idx, 1), ...
                low_err(state_idx, sign_idx, kernel_idx, 1), ...
                high_err(state_idx, sign_idx, kernel_idx, 1), ...
                prepost_p(state_idx, sign_idx, kernel_idx, 1), ...
                prepost_stars{state_idx, sign_idx, kernel_idx, 1}, ...
                metric_value(state_idx, sign_idx, kernel_idx, 2), ...
                low_err(state_idx, sign_idx, kernel_idx, 2), ...
                high_err(state_idx, sign_idx, kernel_idx, 2), ...
                prepost_p(state_idx, sign_idx, kernel_idx, 2), ...
                prepost_stars{state_idx, sign_idx, kernel_idx, 2}, ...
                did_value(state_idx, sign_idx, kernel_idx), ...
                did_low_err(state_idx, sign_idx, kernel_idx), ...
                did_high_err(state_idx, sign_idx, kernel_idx), ...
                did_p(state_idx, sign_idx, kernel_idx), ...
                did_stars{state_idx, sign_idx, kernel_idx}};
        end
    end
end
summary_table = cell2table(summary_rows, 'VariableNames', ...
    {'State', 'Sign', 'Kernel', 'MetricMode', ...
     'Saline_Metric', 'Saline_LowErr', 'Saline_HighErr', 'Saline_PrePost_P', 'Saline_PrePost_Stars', ...
     'Muscimol_Metric', 'Muscimol_LowErr', 'Muscimol_HighErr', 'Muscimol_PrePost_P', 'Muscimol_PrePost_Stars', ...
     'DID_Metric', 'DID_LowErr', 'DID_HighErr', 'DID_P', 'DID_Stars'});
disp(summary_table);

%% Optional export
if export_figures
    if ~isfolder(save_folder)
        mkdir(save_folder);
    end
    if ~isempty(combined_fig)
        export_figure_local(combined_fig, save_folder, sprintf('post_pre_metric_combined_%s_extracted_counts', metric_stub));
    end
    for inj_idx = 1:n_injection
        if ~isempty(single_figs{inj_idx})
            export_figure_local(single_figs{inj_idx}, save_folder, sprintf('post_pre_metric_%s_%s_extracted_counts', lower(injection_labels{inj_idx}), metric_stub));
        end
        if ~isempty(prepost_figs{inj_idx})
            export_figure_local(prepost_figs{inj_idx}, save_folder, sprintf('pre_vs_post_percent_%s_extracted_counts', lower(injection_labels{inj_idx})));
        end
    end
    if make_difference_only_figure && ~isempty(diff_fig)
        export_figure_local(diff_fig, save_folder, sprintf('did_metric_%s_extracted_counts', metric_stub));
    end
end

%% Local functions
function stats = make_prepost_stats_from_counts(count_arr, total_arr, state_idx, sign_idx, kernel_idx)
    stats = struct();
    stats.pre_counts = count_arr(state_idx, sign_idx, 1, kernel_idx);
    stats.post_counts = count_arr(state_idx, sign_idx, 2, kernel_idx);
    stats.pre_total = total_arr(state_idx, sign_idx, 1, kernel_idx);
    stats.post_total = total_arr(state_idx, sign_idx, 2, kernel_idx);
end

function [value, ci_low, ci_high, p_val] = prepost_metric_ci(stats, metric_mode, alpha)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            [value, ci_low, ci_high] = prepost_difference_ci(stats, alpha);
            p_val = twoProportionPValue(stats.pre_counts, stats.pre_total, stats.post_counts, stats.post_total);
        case {'ratio', 'post/pre', 'post_over_pre'}
            [value, ci_low, ci_high, p_val] = prepost_ratio_ci(stats, alpha);
        otherwise
            error('Unknown metric_mode: %s', metric_mode);
    end
end

function [value_pct, ci_low_pct, ci_high_pct] = prepost_difference_ci(stats, alpha)
    if stats.pre_total == 0 || stats.post_total == 0
        value_pct = NaN; ci_low_pct = NaN; ci_high_pct = NaN; return;
    end
    p_pre = stats.pre_counts / stats.pre_total;
    p_post = stats.post_counts / stats.post_total;
    effect = p_post - p_pre;
    se = sqrt(p_post * (1 - p_post) / stats.post_total + p_pre * (1 - p_pre) / stats.pre_total);
    zcrit = norminv(1 - alpha / 2);
    value_pct = 100 * effect;
    ci_low_pct = 100 * (effect - zcrit * se);
    ci_high_pct = 100 * (effect + zcrit * se);
end

function [rr, ci_low, ci_high, p_val] = prepost_ratio_ci(stats, alpha)
    [post_count, post_total, pre_count, pre_total] = safe_counts_for_log_ratio( ...
        stats.post_counts, stats.post_total, stats.pre_counts, stats.pre_total);
    if pre_total == 0 || post_total == 0 || pre_count <= 0 || post_count <= 0
        rr = NaN; ci_low = NaN; ci_high = NaN; p_val = NaN; return;
    end
    p_post = post_count / post_total;
    p_pre = pre_count / pre_total;
    rr = p_post / p_pre;
    log_rr = log(rr);
    se = sqrt(1 / post_count - 1 / post_total + 1 / pre_count - 1 / pre_total);
    if se == 0
        p_val = double(log_rr ~= 0) * 0 + double(log_rr == 0) * 1;
        ci_low = rr; ci_high = rr; return;
    end
    zcrit = norminv(1 - alpha / 2);
    ci_low = exp(log_rr - zcrit * se);
    ci_high = exp(log_rr + zcrit * se);
    p_val = twoSidedNormalPValue(log_rr / se);
end

function [value, ci_low, ci_high, p_val] = did_metric_ci(stats_a, stats_b, metric_mode, alpha)
    % stats_a = Muscimol, stats_b = Saline.
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            [value, ci_low, ci_high, p_val] = did_difference_ci(stats_a, stats_b, alpha);
        case {'ratio', 'post/pre', 'post_over_pre'}
            [value, ci_low, ci_high, p_val] = did_ratio_ci(stats_a, stats_b, alpha);
        otherwise
            error('Unknown metric_mode: %s', metric_mode);
    end
end

function [diff_pct, ci_low_pct, ci_high_pct, p_val] = did_difference_ci(stats_a, stats_b, alpha)
    if stats_a.pre_total == 0 || stats_a.post_total == 0 || stats_b.pre_total == 0 || stats_b.post_total == 0
        diff_pct = NaN; ci_low_pct = NaN; ci_high_pct = NaN; p_val = NaN; return;
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
        p_val = double(diff ~= 0) * 0 + double(diff == 0) * 1;
        ci_low = diff; ci_high = diff;
    else
        zcrit = norminv(1 - alpha / 2);
        ci_low = diff - zcrit * se;
        ci_high = diff + zcrit * se;
        p_val = twoSidedNormalPValue(diff / se);
    end
    diff_pct = 100 * diff;
    ci_low_pct = 100 * ci_low;
    ci_high_pct = 100 * ci_high;
end

function [ror, ci_low, ci_high, p_val] = did_ratio_ci(stats_a, stats_b, alpha)
    [a_post_count, a_post_total, a_pre_count, a_pre_total] = safe_counts_for_log_ratio(stats_a.post_counts, stats_a.post_total, stats_a.pre_counts, stats_a.pre_total);
    [b_post_count, b_post_total, b_pre_count, b_pre_total] = safe_counts_for_log_ratio(stats_b.post_counts, stats_b.post_total, stats_b.pre_counts, stats_b.pre_total);
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
        p_val = double(log_ror ~= 0) * 0 + double(log_ror == 0) * 1;
        ci_low = ror; ci_high = ror; return;
    end
    zcrit = norminv(1 - alpha / 2);
    ci_low = exp(log_ror - zcrit * se);
    ci_high = exp(log_ror + zcrit * se);
    p_val = twoSidedNormalPValue(log_ror / se);
end

function [num1, den1, num0, den0] = safe_counts_for_log_ratio(num1, den1, num0, den0)
    % Continuity correction only when log-ratio would be undefined.
    if num1 == 0 || num0 == 0 || num1 == den1 || num0 == den0
        num1 = num1 + 0.5;
        num0 = num0 + 0.5;
        den1 = den1 + 1;
        den0 = den0 + 1;
    end
end

function f = plot_prepost_percent_figure(count_arr, total_arr, injection_label, state_labels, sign_labels, kernel_labels, prepost_labels, alpha, bar_offset, bar_width)
    n_state = numel(state_labels);
    n_sign = numel(sign_labels);
    n_kernel = numel(kernel_labels);
    ratios = nan(n_state, n_sign, n_kernel, 2);
    low_errs = nan(n_state, n_sign, n_kernel, 2);
    high_errs = nan(n_state, n_sign, n_kernel, 2);
    p_labels = cell(n_state, n_sign, n_kernel);

    for state_idx = 1:n_state
        for sign_idx = 1:n_sign
            for kernel_idx = 1:n_kernel
                for pp_idx = 1:2
                    M = count_arr(state_idx, sign_idx, pp_idx, kernel_idx);
                    N = total_arr(state_idx, sign_idx, pp_idx, kernel_idx);
                    [lo, hi] = wilsonCI(M, N, alpha);
                    p = M / N;
                    ratios(state_idx, sign_idx, kernel_idx, pp_idx) = 100 * p;
                    low_errs(state_idx, sign_idx, kernel_idx, pp_idx) = 100 * (p - lo);
                    high_errs(state_idx, sign_idx, kernel_idx, pp_idx) = 100 * (hi - p);
                end
                p_val = twoProportionPValue(count_arr(state_idx, sign_idx, 1, kernel_idx), total_arr(state_idx, sign_idx, 1, kernel_idx), ...
                                            count_arr(state_idx, sign_idx, 2, kernel_idx), total_arr(state_idx, sign_idx, 2, kernel_idx));
                p_labels{state_idx, sign_idx, kernel_idx} = p_to_sig_label(p_val);
            end
        end
    end

    y_lim = common_ylim_percent(ratios, low_errs, high_errs);
    f = figure('Color', 'w', 'Position', [100, 100, 1300, 900]);
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    for state_idx = 1:n_state
        for sign_idx = 1:n_sign
            ax = nexttile(t, (state_idx - 1) * 2 + sign_idx);
            plot_prepost_grouped_bars(ax, squeeze(ratios(state_idx, sign_idx, :, :)), ...
                squeeze(low_errs(state_idx, sign_idx, :, :)), squeeze(high_errs(state_idx, sign_idx, :, :)), ...
                squeeze(p_labels(state_idx, sign_idx, :)), prepost_labels, kernel_labels, bar_offset, bar_width, y_lim);
            ylabel(ax, 'Significant J %');
            title(ax, sprintf('%s, %s', state_labels{state_idx}, sign_labels{sign_idx}), 'FontWeight', 'bold');
        end
    end
    sgtitle(f, sprintf('Pooled Pre vs Post: %s', injection_label), 'Interpreter', 'none');
end

function plot_metric_grouped_bars(ax, values, low_errs, high_errs, bar_sig_labels, did_sig_labels, condition_labels, kernel_labels, bar_offset, bar_width, y_lim, metric_mode)
    n_kernel = size(values, 1);
    x = 1:n_kernel;
    hold(ax, 'on');
    bar(ax, x - bar_offset, values(:, 1), 'FaceColor', 'b', 'BarWidth', bar_width, 'DisplayName', condition_labels{1});
    bar(ax, x + bar_offset, values(:, 2), 'FaceColor', 'r', 'BarWidth', bar_width, 'DisplayName', condition_labels{2});
    errorbar(ax, x - bar_offset, values(:, 1), low_errs(:, 1), high_errs(:, 1), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    errorbar(ax, x + bar_offset, values(:, 2), low_errs(:, 2), high_errs(:, 2), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    ylim(ax, y_lim);
    draw_metric_reference_line(ax, metric_mode);
    yr = diff(y_lim);
    label_offset = 0.025 * yr;
    bracket_offset = 0.11 * yr;
    bracket_height = 0.025 * yr;
    for kernel_i = 1:n_kernel
        xpos = [x(kernel_i) - bar_offset, x(kernel_i) + bar_offset];
        for inj_idx = 1:2
            y_bar = values(kernel_i, inj_idx) + high_errs(kernel_i, inj_idx) + label_offset;
            text(ax, xpos(inj_idx), y_bar, bar_sig_labels{kernel_i, inj_idx}, ...
                'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'HandleVisibility', 'off');
        end
        y_bracket = max(values(kernel_i, :) + high_errs(kernel_i, :), [], 'omitnan') + bracket_offset;
        draw_sig_bracket(ax, xpos(1), xpos(2), y_bracket, bracket_height, did_sig_labels{kernel_i});
    end
    hold(ax, 'off');
    xticks(ax, x);
    xticklabels(ax, kernel_labels);
    legend(ax, condition_labels, 'Location', 'Best');
    box(ax, 'off');
end

function plot_metric_single_bars(ax, values, low_errs, high_errs, sig_labels, kernel_labels, y_lim, metric_mode, bar_label)
    n_kernel = numel(values);
    x = 1:n_kernel;
    hold(ax, 'on');
    bar(ax, x, values, 'BarWidth', 0.45, 'DisplayName', char(bar_label));
    errorbar(ax, x, values, low_errs, high_errs, 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    ylim(ax, y_lim);
    draw_metric_reference_line(ax, metric_mode);
    yr = diff(y_lim);
    label_offset = 0.035 * yr;
    for kernel_i = 1:n_kernel
        y_text = values(kernel_i) + high_errs(kernel_i) + label_offset;
        text(ax, x(kernel_i), y_text, sig_labels{kernel_i}, 'HorizontalAlignment', 'center', ...
            'FontSize', 12, 'FontWeight', 'bold', 'HandleVisibility', 'off');
    end
    hold(ax, 'off');
    xticks(ax, x);
    xticklabels(ax, kernel_labels);
    box(ax, 'off');
end

function plot_prepost_grouped_bars(ax, ratios, low_errs, high_errs, sig_labels, condition_labels, kernel_labels, bar_offset, bar_width, y_lim)
    n_kernel = size(ratios, 1);
    x = 1:n_kernel;
    hold(ax, 'on');
    bar(ax, x - bar_offset, ratios(:, 1), 'FaceColor', 'b', 'BarWidth', bar_width, 'DisplayName', condition_labels{1});
    bar(ax, x + bar_offset, ratios(:, 2), 'FaceColor', 'r', 'BarWidth', bar_width, 'DisplayName', condition_labels{2});
    errorbar(ax, x - bar_offset, ratios(:, 1), low_errs(:, 1), high_errs(:, 1), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    errorbar(ax, x + bar_offset, ratios(:, 2), low_errs(:, 2), high_errs(:, 2), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    ylim(ax, y_lim);
    yr = diff(y_lim);
    bracket_height = 0.025 * yr;
    bracket_offset = 0.06 * yr;
    for kernel_i = 1:n_kernel
        y_bracket = max(ratios(kernel_i, :) + high_errs(kernel_i, :), [], 'omitnan') + bracket_offset;
        draw_sig_bracket(ax, x(kernel_i) - bar_offset, x(kernel_i) + bar_offset, y_bracket, bracket_height, sig_labels{kernel_i});
    end
    hold(ax, 'off');
    xticks(ax, x);
    xticklabels(ax, kernel_labels);
    legend(ax, condition_labels, 'Location', 'Best');
    box(ax, 'off');
end

function draw_sig_bracket(ax, x1, x2, y, h, label)
    plot(ax, [x1, x1, x2, x2], [y, y + h, y + h, y], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(ax, mean([x1, x2]), y + h, label, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold', 'HandleVisibility', 'off');
end

function y_lim = common_ylim_for_metric(values, low_errs, high_errs, metric_mode)
    y_lower = values(:) - low_errs(:);
    upper = values(:) + high_errs(:);
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            ref = 0;
        case {'ratio', 'post/pre', 'post_over_pre'}
            ref = 1;
        otherwise
            ref = 0;
    end
    vals = [y_lower; upper; ref];
    vals = vals(isfinite(vals));
    if isempty(vals)
        y_lim = [ref - 1, ref + 1];
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
    y_lim = [y_min - 0.18 * yr, y_max + 0.45 * yr];
    if strcmpi(metric_mode, 'ratio') || strcmpi(metric_mode, 'post/pre')
        y_lim(1) = max(0, y_lim(1));
    end
end

function y_lim = common_ylim_percent(values, low_errs, high_errs)
    vals = [values(:) - low_errs(:); values(:) + high_errs(:); 0];
    vals = vals(isfinite(vals));
    if isempty(vals)
        y_lim = [0, 1];
        return;
    end
    y_min = min(vals);
    y_max = max(vals);
    yr = y_max - y_min;
    if yr == 0, yr = max(abs(vals)); end
    if yr == 0, yr = 1; end
    y_lim = [max(0, y_min - 0.10 * yr), y_max + 0.35 * yr];
end

function draw_metric_reference_line(ax, metric_mode)
    xl = xlim(ax);
    hold(ax, 'on');
    switch lower(metric_mode)
        case {'ratio', 'post/pre', 'post_over_pre'}
            plot(ax, xl, [1, 1], 'k--', 'LineWidth', 0.9, 'HandleVisibility', 'off');
        otherwise
            plot(ax, xl, [0, 0], 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    end
    xlim(ax, xl);
end

function label = metric_axis_label(metric_mode)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            label = 'Post - Pre significant J (percentage points)';
        case {'ratio', 'post/pre', 'post_over_pre'}
            label = 'Post / Pre significant J ratio';
        otherwise
            label = metric_mode;
    end
end

function label = did_axis_label(metric_mode)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            label = '(Muscimol Post-Pre) - (Saline Post-Pre) (percentage points)';
        case {'ratio', 'post/pre', 'post_over_pre'}
            label = '(Muscimol Post/Pre) / (Saline Post/Pre)';
        otherwise
            label = metric_mode;
    end
end

function stub = metric_output_stub(metric_mode)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            stub = 'post_minus_pre';
        case {'ratio', 'post/pre', 'post_over_pre'}
            stub = 'post_over_pre';
        otherwise
            stub = regexprep(metric_mode, '[^A-Za-z0-9_-]', '_');
    end
end

function [p_low, p_high] = wilsonCI(M, N, alpha)
    if nargin < 3
        alpha = 0.05;
    end
    if N == 0
        p_low = NaN; p_high = NaN; return;
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
    p_val = twoSidedNormalPValue(z);
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

function [figure_visible, export_figures] = parse_figure_output_mode(output_mode)
    switch lower(string(output_mode))
        case {"visible", "show", "on"}
            figure_visible = 'on';
            export_figures = false;
        case {"save", "export", "off"}
            figure_visible = 'off';
            export_figures = true;
        otherwise
            error('Unknown figure_output_mode: %s. Use ''visible'' or ''save''.', string(output_mode));
    end
end

function export_figure_local(fig, save_folder, output_stub)
    if ~isfolder(save_folder)
        mkdir(save_folder);
    end
    exportgraphics(fig, fullfile(save_folder, [output_stub, '_preview.jpg']), ...
        'ContentType', 'image', 'Resolution', 300, 'BackgroundColor', 'white');
    exportgraphics(fig, fullfile(save_folder, [output_stub, '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    close(fig);
end
