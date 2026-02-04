%% Plot J count from trained model
% 

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

% session_types = {'Muscimol', 'Saline', 'Simulated'};
% session_types = {'Muscimol', 'Saline'};
% session_types = {'EmperorMus', 'EmperorSal'};
session_types = {'SlayerSal', 'SlayerMus', 'EmperorSal', 'EmperorMus'};
session_type_num = length(session_types);
kernel_num = 3;

% state_labels = {'Task', 'RestOpen', 'RestClose'};
% state_names = {'Eyes Open', 'Eyes Open', 'Eye Closed'};
state_labels = {'RestOpen', 'RestClose'};
state_names = {'Eyes Open', 'Eye Closed'};
selected_states = {'RestOpen', 'RestClose'};
% selected_states = {'Task', 'RestClose'};
[~, selected_state_idx] = intersect(state_labels, selected_states, 'stable');

if isempty(selected_state_idx)
    error('plot_J_count:NoStatesSelected', 'selected_states must include at least one valid state');
end

for session_type_idx = 1:session_type_num
    session_type = session_types{session_type_idx};
    
    % load J count data
    folder_name = fullfile(root, 'Data', 'Working', 'J_count');
    file_name = sprintf('Jcount_%s.mat', session_type);
    file_path = fullfile(folder_name, file_name);
    % J_count: (area i, area j, session, posneg, kernel, state, prepost)
    % J_count_by_area: (within/across, session, posneg, kernel, state, prepost)
    load(file_path, 'J_count', 'J_count_by_area', 'J_ratio', 'J_ratio_by_area', 'max_count', 'max_count_by_area', 'session_num'); 

    % bar plot
    f = figure('Position', [100, 100, 800, 800], 'Visible', 'off');
    t = tiledlayout(4, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    area_types = {'Within Area', 'Across Area'};
    posneg_types = {'Positive', 'Negative'};
    
    for area_type_idx = 1:2
        area_type = area_types{area_type_idx};
        for posneg_idx = 1:2
            posneg = posneg_types{posneg_idx};
            for kernel_idx = 1:kernel_num
                ax = nexttile;

                % merge all sessions
                counts = squeeze(J_count_by_area(area_type_idx, :, posneg_idx, kernel_idx, :, :)); % (session, state, prepost)
                max_counts = squeeze(max_count_by_area(area_type_idx, :, posneg_idx, kernel_idx, :, :)); % (session, state, prepost)
                counts = squeeze(sum(counts, 1)); % (state, prepost)
                max_counts = squeeze(sum(max_counts, 1)); % (state, prepost)
                counts = counts(selected_state_idx, :);
                max_counts = max_counts(selected_state_idx, :);

                num_states = size(counts, 1);
                ratios = zeros(num_states, 1);
                ratio_ci = zeros(num_states, 2);
                p_vals = zeros(num_states, 1);
                pre_successes = zeros(num_states, 1);
                pre_trials = zeros(num_states, 1);
                post_successes = zeros(num_states, 1);
                post_trials = zeros(num_states, 1);

                for state_idx = 1:num_states
                    successes_pre = counts(state_idx, 1);
                    trials_pre = max_counts(state_idx, 1);
                    successes_post = counts(state_idx, 2);
                    trials_post = max_counts(state_idx, 2);

                    [ratio_val, ci_high, ci_low, p_val] = ratio_of_proportions(successes_post, trials_post, successes_pre, trials_pre);
                    ratios(state_idx) = ratio_val;
                    ratio_ci(state_idx, :) = [ci_low, ci_high];
                    p_vals(state_idx) = p_val;

                    pre_successes(state_idx) = successes_pre;
                    pre_trials(state_idx) = trials_pre;
                    post_successes(state_idx) = successes_post;
                    post_trials(state_idx) = trials_post;
                end

                b = bar(ratios, 'FaceColor', [1, 1, 1], 'BarWidth', 0.6);
                hold on;
                lower_err = ratios - ratio_ci(:, 1);
                upper_err = ratio_ci(:, 2) - ratios;
                errorbar(1:num_states, ratios, lower_err, upper_err, 'k', 'LineStyle', 'none', 'LineWidth', 1);

                for state_idx = 1:num_states
                    stars = significanceStars(p_vals(state_idx));
                    if ~isempty(stars)
                        star_y = ratios(state_idx) + upper_err(state_idx) + 0.05;
                        text(state_idx, star_y, stars, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
                    end
                end

                yline(1, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
                y_upper = ratios + upper_err;
                y_upper = y_upper(isfinite(y_upper));
                if isempty(y_upper)
                    y_upper = ratios(isfinite(ratios));
                end
                if isempty(y_upper)
                    y_upper = 1;
                end
                max_y = max(y_upper);

                % if num_states == 2
                %     [p_pair, pair_valid] = ratio_ttest(post_successes, post_trials, pre_successes, pre_trials);
                %     if pair_valid
                %         y_line = max_y + 0.1;
                %         plot([1, 2], [y_line, y_line], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
                %         star_pair = significanceStars(p_pair);
                %         if isempty(star_pair)
                %             label_text = 'N.S.';
                %         else
                %             label_text = star_pair;
                %         end
                %         text(1.5, y_line + 0.05, label_text, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
                %         max_y = max(max_y, y_line + 0.05);
                %     end
                % end

                hold off;
                xticks(1:num_states);
                xticklabels(state_names(selected_state_idx));
                ylabel('Post / Pre ratio');
                title(sprintf('%s - %s - Kernel %d', area_type, posneg, kernel_idx));
                ylim([0, max(3, max_y + 0.2)]);
            end
        end
    end

    % save figure
    folder_name = fullfile(root, 'Figures', 'J_count');
    check_path(folder_name);
    figure_name = sprintf('J_count_%s.png', session_type);
    figure_path = fullfile(folder_name, figure_name);
    saveas(f, figure_path);
    set(f, 'PaperUnits', 'inches');
    set(f, 'PaperPosition', [0 0 12 12]);
    set(f, 'PaperSize', [12 12]);
    pdf_name = sprintf('J_count_%s.pdf', session_type);
    pdf_path = fullfile(folder_name, pdf_name);
    print(f, '-painters', '-dpdf', pdf_path);

    % figure 2: post/pre ratios with confidence intervals and tests
    f = figure('Position', [100, 100, 800, 800], 'Visible', 'off');
    t = tiledlayout(4, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    set(f, 'PaperUnits', 'inches');
    set(f, 'PaperPosition', [0 0 12 12]);
    set(f, 'PaperSize', [12 12]);

    area_types = {'Within Area', 'Across Area'};
    posneg_types = {'Positive', 'Negative'};
    
    for area_type_idx = 1:2
        area_type = area_types{area_type_idx};
        for posneg_idx = 1:2
            posneg = posneg_types{posneg_idx};
            for kernel_idx = 1:kernel_num
                nexttile;

                % merge all sessions
                counts = squeeze(J_count_by_area(area_type_idx, :, posneg_idx, kernel_idx, :, :)); % (session, state, prepost)
                max_counts = squeeze(max_count_by_area(area_type_idx, :, posneg_idx, kernel_idx, :, :)); % (session, state, prepost)
                counts = squeeze(sum(counts, 1)); % (state, prepost)
                max_counts = squeeze(sum(max_counts, 1)); % (state, prepost)
                counts = counts(selected_state_idx, :);
                max_counts = max_counts(selected_state_idx, :);
                ratio = counts ./ max_counts; % (state, prepost)
                CI = zeros(2, size(ratio, 1), size(ratio, 2)); % (low/high, state, prepost)
                for state_idx = 1:size(ratio, 1)
                    for prepost_idx = 1:size(ratio, 2)
                        M = counts(state_idx, prepost_idx);
                        N = max_counts(state_idx, prepost_idx);
                        [p_low, p_high] = wilsonCI(M, N, 0.05);
                        CI(1, state_idx, prepost_idx) = p_low;
                        CI(2, state_idx, prepost_idx) = p_high;
                    end
                end
                CI_low = squeeze(CI(1, :, :)) * 100; % lower CI bound (% units)
                CI_high = squeeze(CI(2, :, :)) * 100; % upper CI bound (% units)

                % bar plot of ratios, with error bars
                bar_data = ratio * 100;
                b = bar(bar_data);
                set(b, 'BarWidth', 0.6);
                % set(b, 'BarWidth', 0.6, 'FaceColor', [0.7, 0.7, 0.7]);
                hold on;
                lower_err = bar_data - CI_low;
                upper_err = CI_high - bar_data;
                x_points = zeros(size(bar_data));
                for prepost_idx = 1:size(bar_data, 2)
                    x = b(prepost_idx).XEndPoints;
                    errorbar(x, bar_data(:, prepost_idx), lower_err(:, prepost_idx), upper_err(:, prepost_idx), 'k', 'LineStyle', 'none', 'LineWidth', 1);
                    x_points(:, prepost_idx) = x;
                end

                % significance annotations for pre vs post using two-proportion z-test
                for state_idx = 1:size(bar_data, 1)
                    success_pre = counts(state_idx, 1);
                    trials_pre = max_counts(state_idx, 1);
                    success_post = counts(state_idx, 2);
                    trials_post = max_counts(state_idx, 2);
                    p_val = twoProportionPValue(success_pre, trials_pre, success_post, trials_post);
                    stars = significanceStars(p_val);
                    if ~isempty(stars)
                        x_mid = mean(x_points(state_idx, :));
                        y_candidates = bar_data(state_idx, :) + upper_err(state_idx, :);
                        y_candidates = y_candidates(isfinite(y_candidates));
                        if isempty(y_candidates)
                            y_candidates = bar_data(state_idx, isfinite(bar_data(state_idx, :)));
                        end
                        if isempty(y_candidates)
                            y_candidates = 0;
                        end
                        star_y = max(y_candidates) + 1;
                        text(x_mid, star_y, stars, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
                    end
                end

                max_y = max(bar_data + upper_err, [], 'all');

                if size(bar_data, 1) == 2
                    % statistical comparison between the two states (e.g., Task vs RestClose)
                    success_state1_pre = counts(1, 1);
                    trials_state1_pre = max_counts(1, 1);
                    success_state1_post = counts(1, 2);
                    trials_state1_post = max_counts(1, 2);
                    success_state2_pre = counts(2, 1);
                    trials_state2_pre = max_counts(2, 1);
                    success_state2_post = counts(2, 2);
                    trials_state2_post = max_counts(2, 2);
                    % Pre comparison
                    p_val_pre = twoProportionPValue(success_state1_pre, trials_state1_pre, success_state2_pre, trials_state2_pre);
                    stars_pre = significanceStars(p_val_pre);
                    if ~isempty(stars_pre)
                        y_line = max_y + 3;
                        x1 = x_points(1, 1);
                        x2 = x_points(2, 1);
                        plot([x1, x2], [y_line, y_line], 'k-', 'LineWidth', 1);
                        text(mean([x1, x2]), y_line + 1, ['Pre ', stars_pre], 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
                        max_y = max(max_y, y_line + 1);
                    end
                    % Post comparison
                    p_val_post = twoProportionPValue(success_state1_post, trials_state1_post, success_state2_post, trials_state2_post);
                    stars_post = significanceStars(p_val_post);
                    if ~isempty(stars_post)
                        y_line = max_y + 3;
                        x1 = x_points(1, 2);
                        x2 = x_points(2, 2);
                        plot([x1, x2], [y_line, y_line], 'k-', 'LineWidth', 1);
                        text(mean([x1, x2]), y_line + 1, ['Post ', stars_post], 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
                        max_y = max(max_y, y_line + 1);
                    end
                end

                hold off;
                xticks(1:length(selected_state_idx));
                xticklabels(state_names(selected_state_idx));
                % xlabel('State');
                ylabel('Significant connection Ratio (%)');
                title(sprintf('%s - %s - Kernel %d', area_type, posneg, kernel_idx));
                legend({'Pre', 'Post'}, 'Location', 'Best');
                y_limits = bar_data + upper_err;
                y_limits = y_limits(isfinite(y_limits));
                if isempty(y_limits)
                    y_limits = bar_data(isfinite(bar_data));
                end
                if isempty(y_limits)
                    y_limits = 0;
                end
                ylim([0, max(30, max(y_limits) + 5)]);
            end
        end
    end

    % save figure
    folder_name = fullfile(root, 'Figures', 'J_count');
    check_path(folder_name);
    figure_name = sprintf('J_count_ratio_%s.png', session_type);
    figure_path = fullfile(folder_name, figure_name);
    saveas(f, figure_path);
    pdf_name = sprintf('J_count_ratio_%s.pdf', session_type);
    pdf_path = fullfile(folder_name, pdf_name);
    print(f, '-painters', '-dpdf', pdf_path);

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

function [ratio, ci_high, ci_low, p] = ratio_of_proportions(significant_post, trials_post, significant_pre, trials_pre)
% RATIO_OF_PROPORTIONS
%   [ratio, ci_high, ci_low, p] = ratio_of_proportions(a, n1, b, n0)
%
% Inputs:
%   significant_post : a  (successes in the POST group)
%   trials_post      : n1 (trials in the POST group)
%   significant_pre  : b  (successes in the PRE  group)
%   trials_pre       : n0 (trials in the PRE  group)
%
% Outputs:
%   ratio   : Relative Risk (RR) = (a/n1) / (b/n0)
%   ci_high : Upper 95% CI bound for RR (asymmetric, log-scale)
%   ci_low  : Lower 95% CI bound for RR (asymmetric, log-scale)
%   p       : Two-sided p-value (Fisher exact if small counts; else z-test)
%
% Notes:
%   - CI is computed on the log scale and exponentiated (Katz method).
%   - If any cell count is zero, applies a Haldane–Anscombe 0.5 correction.
%   - The p-value comes from Fisher's exact test when any cell < 5 and
%     fishertest is available; otherwise from a pooled two-proportion z-test.
%
% Example:
%   [RR, hi, lo, p] = ratio_of_proportions(30, 120, 18, 100)

    % ---- Input checks
    a  = double(significant_post);
    n1 = double(trials_post);
    b  = double(significant_pre);
    n0 = double(trials_pre);

    if any([n1, n0] <= 0) || any([a, b] < 0) || a > n1 || b > n0 || ...
       any(~isfinite([a, n1, b, n0]))
        error('Inputs must be finite, with 0 <= successes <= trials and trials > 0.');
    end

    % ---- Basic proportions & RR (point estimate)
    p1 = a / n1;
    p0 = b / n0;

    % Handle zero denominator for RR gracefully
    if p0 == 0 && p1 == 0
        ratio = NaN;  % 0/0 undefined
    elseif p0 == 0
        ratio = Inf;  % division by zero -> infinite RR
    else
        ratio = p1 / p0;
    end

    % ---- 2x2 table cells
    c = n1 - a;  % failures post
    d = n0 - b;  % failures pre

    % ---- Confidence interval for RR (log scale, Katz); guard zeros
    aa = a; bb = b; nn1 = n1; nn0 = n0;
    use_haldene = (a == 0) || (b == 0) || (c == 0) || (d == 0);

    if use_haldene
        % Haldane–Anscombe correction (add 0.5 to each cell)
        aa = a + 0.5; bb = b + 0.5;
        cc = c + 0.5; dd = d + 0.5;
        nn1 = aa + cc; nn0 = bb + dd;
        p1_adj = aa / nn1;
        p0_adj = bb / nn0;
        RR_for_CI = p1_adj / p0_adj;
        se_logRR = sqrt(1/aa - 1/nn1 + 1/bb - 1/nn0);
    else
        RR_for_CI = ratio;
        se_logRR = sqrt(1/a - 1/n1 + 1/b - 1/n0);
    end

    z = 1.96; % 95% CI
    if ~isfinite(RR_for_CI) || ~isfinite(se_logRR)
        ci_low = NaN;
        ci_high = NaN;
    else
        ci_low  = exp(log(RR_for_CI) - z*se_logRR);
        ci_high = exp(log(RR_for_CI) + z*se_logRR);
    end

    % ---- Significance test
    % If any cell < 5 AND fishertest exists, use Fisher's exact test
    use_fisher = any([a, b, c, d] < 5) && exist('fishertest','file') == 2;

    if use_fisher
        % Contingency table: rows=POST/ PRE, cols=success/failure
        tbl = [a, c; b, d];
        try
            [~, p] = fishertest(tbl, 'Tail', 'two');
        catch
            % Fallback to z-test if fishertest errors for any reason
            p = twoProportionPValue(a, n1, b, n0);
        end
    else
        p = twoProportionPValue(a, n1, b, n0);
    end
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

function stars = significanceStars(p_val)
    % Convert p-value to significance stars

    if isnan(p_val) || p_val >= 0.05
        stars = '';
    elseif p_val < 0.001
        stars = '***';
    elseif p_val < 0.01
        stars = '**';
    else
        stars = '*';
    end
end

function [p_value, is_valid] = ratio_ttest(success_post, trials_post, success_pre, trials_pre)
    valid = (trials_post > 0) & (trials_pre > 0);
    if nnz(valid) < 2
        p_value = NaN;
        is_valid = false;
        return;
    end

    p_post = success_post(valid) ./ trials_post(valid);
    p_pre = success_pre(valid) ./ trials_pre(valid);
    ratios = p_post ./ p_pre;

    [~, p_value] = ttest(ratios, ones(size(ratios)), 'Alpha', 0.05);
    is_valid = true;
end