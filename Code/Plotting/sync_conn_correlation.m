%% resting_state_scatter.m - Scatter plot of connection density, eyes open vs. closed
% Used dataset: SlayerMus-Pre, SlayerSal-Pre, ZeppelinNoinj-Pre, ZeppelinMus-Pre, ZeppelinSal-Pre,
% EmperorMus-Pre, EmperorSal-Pre. (Slayer Noinj does not have eyeID)

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
session_types = {'SlayerMus', 'SlayerSal', 'ZeppelinNoinj', 'ZeppelinMus', 'ZeppelinSal', 'EmperorMus', 'EmperorSal'};
% session_types = {'SlayerMus', 'SlayerSal', 'EmperorMus', 'EmperorSal'};
session_type_num = length(session_types);
for kernel_idx = 1:3
    f = figure('Position', [100, 100, 800, 800], 'Visible', 'on');
    t = tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    data_limits = zeros(4, 2); % (tile, [min, max])
    data_limits(:, 1) = inf;  data_limits(:, 2) = -inf; 

    total_count = zeros(4, 2); % (tile, x/y)
    total_max_count = zeros(4); % (tile)
    total_disagreement = zeros(4, 2); % (tile, n01/n10)

    for session_type_idx = 1:session_type_num
        session_type = session_types{session_type_idx};
        folder_name = fullfile(root, 'Data', 'Working', 'J_count');
        file_name = sprintf('Jcount_Full_%s.mat', session_type);
        file_path = fullfile(folder_name, file_name);
        load(file_path, 'J_count', 'J_count_by_area', 'J_ratio', 'J_ratio_by_area', 'max_count', 'max_count_by_area',...
        'disagreement_resting', 'disagreement_prepost', 'disagreement_resting_by_area', 'disagreement_prepost_by_area',...
        'session_num', 'kernel', 'reg', 'epoch');

        if session_num == 0
            fprintf('No sessions for %s, skipping...\n', session_type);
            continue; % skip if no sessions
        end
        
        area_types = {'Within Area', 'Across Area'};
        posneg_types = {'Positive', 'Negative'};
        for area_type_idx = 1:2
            area_type = area_types{area_type_idx};
            for posneg_idx = 1:2
                posneg = posneg_types{posneg_idx};
                tile_idx = posneg_idx + (area_type_idx - 1) * 2;
                nexttile(tile_idx); % row: area type, col: posneg
                hold on;

                % marker style by animal
                if contains(session_type, 'Slayer') % Slayer: red circle
                    marker = 'o';
                    color = 'r';
                elseif contains(session_type, 'Zeppelin') % Zeppelin: blue square
                    marker = 's';
                    color = 'b';
                elseif contains(session_type, 'Emperor') % Emperor: black triangle
                    marker = '^';
                    color = 'k';
                else 
                    marker = 'x';
                    color = 'm';
                end

                x = squeeze(J_ratio(area_type_idx, :, posneg_idx, kernel_idx, 1, 1)); % Eyes Open
                y = squeeze(J_ratio(area_type_idx, :, posneg_idx, kernel_idx, 2, 1)); % Eyes Closed
                xy_max_count = squeeze(max_count(area_type_idx, :, posneg_idx, kernel_idx, 1, 1)); % (session, state, prepost)
                filter = xy_max_count > 0;
                if sum(filter) == 0
                    continue; % skip if no connections
                end
                x = x(filter);
                y = y(filter);
                xy_max_count = xy_max_count(filter);
                marker_size = log(xy_max_count) * 10; % marker size by log of max count
            
                scatter(x, y, marker_size, color, marker, 'filled', 'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');

                % accumulate totals for this tile
                x_count = sum(squeeze(J_count(area_type_idx, :, posneg_idx, kernel_idx, 1, 1))); % total count across sessions
                y_count = sum(squeeze(J_count(area_type_idx, :, posneg_idx, kernel_idx, 2, 1)));
                xy_max_count = sum(squeeze(max_count(area_type_idx, :, posneg_idx, kernel_idx, 1, 1)));
                disagree_01 = sum(squeeze(disagreement_resting(1, area_type_idx, :, posneg_idx, kernel_idx, 1))); % count of connections present in Open but not Closed
                disagree_10 = sum(squeeze(disagreement_resting(2, area_type_idx, :, posneg_idx, kernel_idx, 1))); % count of connections present in Closed but not Open
                total_count(tile_idx, 1) = total_count(tile_idx, 1) + x_count;
                total_count(tile_idx, 2) = total_count(tile_idx, 2) + y_count;
                total_max_count(tile_idx) = total_max_count(tile_idx) + xy_max_count;
                total_disagreement(tile_idx, 1) = total_disagreement(tile_idx, 1) + disagree_01;
                total_disagreement(tile_idx, 2) = total_disagreement(tile_idx, 2) + disagree_10;
                
                % update data limits
                data_limits(tile_idx, 1) = min([data_limits(tile_idx, 1), min(x), min(y)]);
                data_limits(tile_idx, 2) = max([data_limits(tile_idx, 2), max(x), max(y)]);
            end
        end
    end
    sgtitle(sprintf('Connection Density Scatter Plot (Kernel %d)', kernel_idx));

    % Axes settings for each tile
    for area_type_idx = 1:2
        for posneg_idx = 1:2
            tile_idx = posneg_idx + (area_type_idx - 1) * 2;

            % calculate limits
            lim_range = data_limits(tile_idx, 2) - data_limits(tile_idx, 1);
            if lim_range == 0
                lim_min = 0;
                lim_max = 1;
                fprintf('Warning: No data in tile %d, setting default limits.\n', tile_idx);
            else
                lim_min = data_limits(tile_idx, 1) - 0.1 * lim_range;
                lim_max = data_limits(tile_idx, 2) + 0.1 * lim_range;
                % lim_min = max(lim_min, 0);
                % lim_max = min(lim_max, 1);
            end

            % text
            nexttile(tile_idx);
            title_str = sprintf('%s, %s', area_types{area_type_idx}, posneg_types{posneg_idx});
            title(title_str);
            xlabel('Eyes Open Connection Density');
            ylabel('Eyes Closed Connection Density');
            axis equal;

            % Statistics
            px = total_count(tile_idx, 1)/total_max_count(tile_idx);
            py = total_count(tile_idx, 2)/total_max_count(tile_idx);
            max_count = total_max_count(tile_idx);

            % % statistics: two-proportion z-test
            % p_pool = (total_count(tile_idx, 1) + total_count(tile_idx, 2)) / (2 * max_count);
            % se = sqrt(p_pool * (1 - p_pool) * (2 / max_count));
            % if se == 0
            %     fprintf('Warning: Standard error is zero in tile %d, skipping p-value calculation.\n', tile_idx);
            %     p_val = 1;
            % else
            %     z = (px - py) / se;
            %     p_val = 2 * (1 - normcdf(abs(z))); % two-tailed
            %     text_str = sprintf('p = %.3f', p_val);
            %     text_color = 'k';
            %     if p_val < 0.001
            %         text_str = [text_str, ' (***)']; %#ok<AGROW>
            %         text_color = 'r';
            %     elseif p_val < 0.01
            %         text_str = [text_str, ' (**)']; %#ok<AGROW>
            %         text_color = 'r';
            %     elseif p_val < 0.05
            %         text_str = [text_str, ' (*)']; %#ok<AGROW>
            %         text_color = 'r';
            %     end
            %     text_x = lim_min + 0.05 * lim_range;
            %     text_y = lim_max - 0.15 * lim_range;
            %     text(text_x, text_y, sprintf('Open: %.2f%% (%d / %d)', px * 100, total_count(tile_idx, 1), max_count), 'HorizontalAlignment', 'left', 'FontSize', 10);
            %     text(text_x, text_y - 0.05 * lim_range, sprintf('Closed: %.2f%% (%d / %d)', py * 100, total_count(tile_idx, 2), max_count), 'HorizontalAlignment', 'left', 'FontSize', 10);
            %     text(text_x, text_y - 0.1 * lim_range, text_str, 'HorizontalAlignment', 'left', 'FontSize', 10, 'Color', text_color);
            % end

            % Alternative statistics: McNemar's test
            n01 = total_disagreement(tile_idx, 1); % Open yes, Closed no
            n10 = total_disagreement(tile_idx, 2); % Open no, Closed yes

            p_val = mcnemar_test(n01, n10);
            text_str = sprintf('McNemar p = %.3f', p_val);
            text_color = 'k';
            if p_val < 0.001
                text_str = [text_str, ' (***)']; %#ok<AGROW>
                text_color = 'r';
            elseif p_val < 0.01
                text_str = [text_str, ' (**)']; %#ok<AGROW>
                text_color = 'r';
            elseif p_val < 0.05
                text_str = [text_str, ' (*)']; %#ok<AGROW>
                text_color = 'r';
            end
            text_x = lim_min + 0.05 * lim_range;
            text_y = lim_max - 0.15 * lim_range;
            text(text_x, text_y, sprintf('Open: %.2f%% (%d / %d)', px * 100, total_count(tile_idx, 1), max_count), 'HorizontalAlignment', 'left', 'FontSize', 10);
            text(text_x, text_y - 0.05 * lim_range, sprintf('Closed: %.2f%% (%d / %d)', py * 100, total_count(tile_idx, 2), max_count), 'HorizontalAlignment', 'left', 'FontSize', 10);
            text(text_x, text_y - 0.1 * lim_range, text_str, 'HorizontalAlignment', 'left', 'FontSize', 10, 'Color', text_color);

            % legend
            hold on;
            scatter(nan, nan, 50, 'r', 'o', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Slayer');
            scatter(nan, nan, 50, 'b', 's', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Zeppelin');
            scatter(nan, nan, 50, 'k', '^', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', 'Emperor');
            legend('Location', 'southeast');
            hold off;

            % limits and reference line
            xlim([lim_min, lim_max]);
            ylim([lim_min, lim_max]);
            hold on;
            plot([lim_min, lim_max], [lim_min, lim_max], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off'); % x=y line
            hold off;
        end
    end

    % Save figure
    save_path = fullfile(root, 'Figures', 'RestingStateScatters');
    check_path(save_path);
    save_file_name = sprintf('RestingStateScatter_Kernel%d.png', kernel_idx);
    saveas(f, fullfile(save_path, save_file_name));
    close(f);
end

function p_val = mcnemar_test(n01, n10)
    % McNemar's test for paired nominal data
    % n01: count of pairs where condition 1 is yes and condition 2 is no
    % n10: count of pairs where condition 1 is no and condition 2 is yes

    if n01 + n10 == 0
        p_val = 1; % no disagreement, p-value is 1
        return;
    end

    % % Method 1: Chi-square with continuity correction
    % chi2 = (abs(n01 - n10) - 1)^2 / (n01 + n10);
    % p_val = 1 - chi2cdf(chi2, 1); % chi-squared distribution with 1 degree of freedom

    % Method 2: Exact binomial test
    n = n01 + n10;
    k = min(n01, n10);
    p_val = 2 * binocdf(k, n, 0.5); % two-tailed
    p_val = min(p_val, 1); % cap at 1

end