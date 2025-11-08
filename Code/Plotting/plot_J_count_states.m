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
session_types = {'Muscimol', 'Saline'};
session_type_num = length(session_types);
kernel_num = 3;

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
                nexttile;

                % merge all sessions
                counts = squeeze(J_count_by_area(area_type_idx, :, posneg_idx, kernel_idx, :, :)); % (session, state, prepost)
                max_counts = squeeze(max_count_by_area(area_type_idx, :, posneg_idx, kernel_idx, :, :)); % (session, state, prepost)
                counts = squeeze(sum(counts, 1)); % (state, prepost)
                max_counts = squeeze(sum(max_counts, 1)); % (state, prepost)
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
                CI_low = squeeze(CI(1, :, :))'; % (prepost, state)
                CI_high = squeeze(CI(2, :, :))';
                CI_low = CI_low * 100; % lower CI bound (% units)
                CI_high = CI_high * 100; % upper CI bound (% units)

                % bar plot of ratios, with error bars (prepost on x-axis, state as series)
                bar_data = (ratio') * 100; % (prepost, state)
                b = bar(bar_data, 'grouped');
                set(b, 'BarWidth', 0.6);
                hold on;
                lower_err = bar_data - CI_low;
                upper_err = CI_high - bar_data;
                x_points = zeros(size(bar_data));
                for state_idx = 1:size(bar_data, 2)
                    x = b(state_idx).XEndPoints;
                    errorbar(x, bar_data(:, state_idx), lower_err(:, state_idx), upper_err(:, state_idx), 'k', 'LineStyle', 'none', 'LineWidth', 1);
                    x_points(:, state_idx) = x;
                end

                % significance annotations comparing Task (state 1) vs RestClose (state 3)
                task_state_idx = 1;
                restclose_state_idx = 3;
                for prepost_idx = 1:size(bar_data, 1)
                    success_task = counts(task_state_idx, prepost_idx);
                    trials_task = max_counts(task_state_idx, prepost_idx);
                    success_restclose = counts(restclose_state_idx, prepost_idx);
                    trials_restclose = max_counts(restclose_state_idx, prepost_idx);
                    p_val = twoProportionPValue(success_task, trials_task, success_restclose, trials_restclose);
                    stars = significanceStars(p_val);
                    if ~isempty(stars)
                        x_task = x_points(prepost_idx, task_state_idx);
                        x_restclose = x_points(prepost_idx, restclose_state_idx);
                        x_mid = mean([x_task, x_restclose]);
                        y_candidates = [bar_data(prepost_idx, task_state_idx) + upper_err(prepost_idx, task_state_idx), ...
                                        bar_data(prepost_idx, restclose_state_idx) + upper_err(prepost_idx, restclose_state_idx)];
                        y_candidates = y_candidates(isfinite(y_candidates));
                        if isempty(y_candidates)
                            y_candidates = [bar_data(prepost_idx, task_state_idx), bar_data(prepost_idx, restclose_state_idx)];
                            y_candidates = y_candidates(isfinite(y_candidates));
                        end
                        if isempty(y_candidates)
                            y_candidates = 0;
                        end
                        star_y = max(y_candidates) + 1;
                        plot([x_task, x_restclose], [star_y - 0.5, star_y - 0.5], 'k-', 'LineWidth', 1);
                        text(x_mid, star_y, stars, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'k');
                    end
                end
                hold off;
                xticks(1:2);
                xticklabels({'Pre', 'Post'});
                xlabel('Pre/Post');
                ylabel('J Ratio (%)');
                title(sprintf('%s - %s - Kernel %d', area_type, posneg, kernel_idx));
                legend({'Task', 'RestOpen', 'RestClose'}, 'Location', 'Best');
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
    figure_name = sprintf('J_count_%s.png', session_type);
    figure_path = fullfile(folder_name, figure_name);
    saveas(f, figure_path);
    set(f, 'PaperUnits', 'inches');
    set(f, 'PaperPosition', [0 0 12 12]);
    set(f, 'PaperSize', [12 12]);
    pdf_name = sprintf('J_count_%s.pdf', session_type);
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

function p_val = twoProportionPValue(success1, trials1, success2, trials2)
    % Two-proportion z-test (two-tailed) between two Bernoulli samples

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