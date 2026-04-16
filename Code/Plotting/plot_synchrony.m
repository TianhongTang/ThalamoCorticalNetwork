% plot_synchrony.m - Synchrony plot for eyes-open vs eyes-closed states

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

result_folder = fullfile(root, 'Data', 'Working', 'Analysis', 'Synchrony');
result_name = sprintf('synchrony_results.mat');
result_path = fullfile(result_folder, result_name);
load(result_path, 'tasks');

% fig 1: thalamus within area synchrony
f = figure("Position", [100, 100, 1200, 600], "Visible", "off");
t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
% task_filter = false(1, length(tasks));

for plot_idx = 1:2 % 1: chi, 2: Pearson's r
    nexttile;
    hold on;
    max_max = -Inf;
    min_min = Inf;
    
    % dummy plot to create legend
    x = [NaN, NaN];
    y = [NaN, NaN];
    scatter(x, y, 50, 'o', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', 'none', 'DisplayName', 'Slayer');
    scatter(x, y, 50, 's', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', 'none', 'DisplayName', 'Zeppelin');
    scatter(x, y, 50, '^', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', 'none', 'DisplayName', 'Emperor');
    % scatter(x, y, 50, 's', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', 'none', 'DisplayName', 'Emperor Saline');
    % scatter(x, y, 50, 'o', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', 'none', 'DisplayName', 'Emperor Muscimol');


    for i = 1:length(tasks)
        task = tasks{i};
        if strcmp(task.filter_type, 'Area') && ...
        length(task.area_filter) == 1 && ...
        strcmp(task.area_filter{1}, 'Thalamus') && ...
        ~task.across_area 
        % && strcmp(task.dataset(1:6), 'Empero')
            task_filter(i) = true;
            % plot on figure
            if plot_idx == 1
                x = task.sync_chi_x;
                y = task.sync_chi_y;
                title_str = 'Chi';
            else
                x = task.sync_corr_x;
                y = task.sync_corr_y;
                x_err = task.sync_corr_x_err;
                y_err = task.sync_corr_y_err;
                title_str = "Pearson's r";
            end
            nan_filter = ~isnan(x) & ~isnan(y);
            x = x(nan_filter);
            y = y(nan_filter);

            % Error bars for Pearson's r
            % if plot_idx == 2
            %     x_err = x_err(nan_filter);
            %     y_err = y_err(nan_filter);
            %     x_neg = x - x_err;
            %     x_pos = x + x_err;
            %     y_neg = y - y_err;
            %     y_pos = y + y_err;
            %     errorbar(x, y, y_err, y_err, x_err, x_err, 'Color',...
            %      [0, 0, 0, 0.3], 'LineStyle', 'none', 'CapSize', 5, 'HandleVisibility','off');
            % end

            % markers for different datasets
            switch task.dataset(1:6)
                case 'Slayer' % Slayer
                    color = [1, 0, 0]; % red
                    marker = 'o'; % circle
                case 'Zeppel' % Zeppelin
                    color = [0, 0, 1]; % blue
                    marker = 's'; % square
                case 'Empero' % Emperor
                    color = [0, 0, 0]; % black
                    marker = '^'; % triangle
                otherwise 
                    color = [0.5, 0.5, 0.5]; % gray
            end
            % switch task.dataset(end-2:end)
            %     case 'Mus' % Muscimol
            %         marker = 'o'; % circle
            %     case 'Sal' % Saline
            %         marker = 's'; % square
            %     case 'inj' % No injection
            %         marker = 'd'; % diamond
            %     otherwise
            %         marker = '^'; % triangle
            % end
            % marker = 'o'; % for now, use circle for all

            scatter(x, y, 50, marker, 'filled', 'MarkerFaceAlpha', 0.6,...
             'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'HandleVisibility','off');
            if ~isempty(x)
                max_max = max(max_max, max(x));
                min_min = min(min_min, min(x));
                if plot_idx == 2
                    max_max = max(max_max, max(x + x_err));
                    min_min = min(min_min, min(x - x_err));
                end
            end
            if ~isempty(y)
                max_max = max(max_max, max(y));
                min_min = min(min_min, min(y));
                if plot_idx == 2
                    max_max = max(max_max, max(y + y_err));
                    min_min = min(min_min, min(y - y_err));
                end
            end
        end
    end

    lim_range = max_max - min_min;
    max_lim = max_max + 0.1 * lim_range;
    min_lim = min_min - 0.1 * lim_range;

    axis equal;
    xlim([min_lim, max_lim]);
    ylim([min_lim, max_lim]);
    plot([min_lim, max_lim], [min_lim, max_lim], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
    title(title_str);
    xlabel('Eyes Open');
    ylabel('Eyes Closed');
    legend();
end
sgtitle('Thalamus Within Area Synchrony: Eyes Open vs Eyes Closed');
fig_save_folder = fullfile(root, 'Figures', 'Synchrony');
check_path(fig_save_folder);
fig_save_path = fullfile(fig_save_folder, 'synchrony_thalamus_within_area.png');
saveas(f, fig_save_path);
close(f);

% fig 2: cortex single area synchrony
controls = {'Mus', 'Sal'};
for control_idx = 1:length(controls)
    control = controls{control_idx};

    f = figure("Position", [100, 100, 1200, 600], "Visible", "off");
    t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    area_list = {'ACC', 'VLPFC'};
    prepost_list = {'Pre', 'Post'};

    colors = [0, 0, 1; 1, 0, 0];
    markers = {'+', 'x'};
    for plot_idx = 1:2 % 1: chi, 2: Pearson's
        nexttile;
        hold on;
        max_max = -Inf;
        min_min = Inf;

        % dummy plot to create legend
        x = [NaN, NaN];
        y = [NaN, NaN];
        scatter(x, y, 50, '+', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', colors(1, :), 'MarkerEdgeColor', 'none', 'DisplayName', 'ACC Pre');
        scatter(x, y, 50, '+', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', colors(2, :), 'MarkerEdgeColor', 'none', 'DisplayName', 'ACC Post');
        scatter(x, y, 50, 'x', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', colors(1, :), 'MarkerEdgeColor', 'none', 'DisplayName', 'VLPFC Pre');
        scatter(x, y, 50, 'x', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', colors(2, :), 'MarkerEdgeColor', 'none', 'DisplayName', 'VLPFC Post');

        for area_idx = 1:length(area_list)
            area = area_list{area_idx};
            marker = markers{area_idx};
            for prepost_idx = 1:length(prepost_list)
                prepost = prepost_list{prepost_idx};
                color = colors(prepost_idx, :);
                
                for i = 1:length(tasks)
                    task = tasks{i};
                    if strcmp(task.filter_type, 'Area') && ...
                    length(task.area_filter) == 1 && ...
                    strcmp(task.area_filter{1}, area) && ...
                    ~task.across_area && ...
                    strcmp(task.prepost, prepost) && ...
                    contains(task.dataset, control)
                        % plot on figure
                        if plot_idx == 1
                            x = task.sync_chi_x;
                            y = task.sync_chi_y;
                            title_str = 'Chi';
                        else
                            x = task.sync_corr_x;
                            y = task.sync_corr_y;
                            x_err = task.sync_corr_x_err;
                            y_err = task.sync_corr_y_err;
                            title_str = "Pearson's r";
                        end
                        nan_filter = ~isnan(x) & ~isnan(y);
                        x = x(nan_filter);
                        y = y(nan_filter);
                        if plot_idx == 2
                            x_err = x_err(nan_filter);
                            y_err = y_err(nan_filter);
                            % errorbar(x, y, y_err, y_err, x_err, x_err, 'o', 'Color', [0, 0, 1, 0.3], 'LineStyle', 'none', 'CapSize', 0);
                        end
                        scatter(x, y, 50, marker, 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'HandleVisibility','off');
                        
                        if ~isempty(x)
                            max_max = max(max_max, max(x));
                            min_min = min(min_min, min(x));
                        end
                        if ~isempty(y)
                            max_max = max(max_max, max(y));
                            min_min = min(min_min, min(y));
                        end
                    end
                end
            end
        end
        axis equal;
        lim_range = max_max - min_min;
        max_lim = max_max + 0.1 * lim_range;
        min_lim = min_min - 0.1 * lim_range;
        xlim([min_lim, max_lim]);
        ylim([min_lim, max_lim]);
        plot([min_lim, max_lim], [min_lim, max_lim], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
        title(title_str);
        xlabel('Eyes Open');
        ylabel('Eyes Closed');
    end
    sgtitle(sprintf('Cortex Single Area Synchrony (%s): Eyes Open vs Eyes Closed', control));
    fig_save_folder = fullfile(root, 'Figures', 'Synchrony');
    check_path(fig_save_folder);
    fig_save_path = fullfile(fig_save_folder, sprintf('synchrony_cortex_single_area_%s.png', control));
    saveas(f, fig_save_path);
    close(f);
end

% fig 3: cortex across area synchrony
fprintf('Fig 3: Plotting cortex across area synchrony...\n');
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums');
prepost_types = {'Pre', 'Post'};

controls = {'Mus', 'Sal'};
for control_idx = 1:length(controls)
    control = controls{control_idx};

    f = figure("Position", [100, 100, 600, 600], "Visible", "off");
    t = tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    prepost_list = {'Pre', 'Post'};

    colors = [0, 0, 1; 1, 0, 0];
    for plot_idx = 2:2 % 1: chi, 2: Pearson's
        nexttile;
        hold on;
        max_max = -Inf;
        min_min = Inf;

        % dummy plot to create legend
        x = [NaN, NaN];
        y = [NaN, NaN];
        scatter(x, y, 50, 'o', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', 'none', 'DisplayName', 'Pre');
        scatter(x, y, 50, '^', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', 'none', 'DisplayName', 'Post');
        plot(x, y, '-', 'Color', [1,0,0], 'LineWidth', 0.5, 'DisplayName', 'Slayer');
        plot(x, y, '-', 'Color', [0,0,1], 'LineWidth', 0.5, 'DisplayName', 'Zeppelin');
        plot(x, y, '-', 'Color', [0,0,0], 'LineWidth', 0.5, 'DisplayName', 'Emperor');

        for dataset_idx = 1:dataset_num
            dataset_name = dataset_names{dataset_idx};
            if ~contains(dataset_name, control)
                continue;
            end
            session_num = session_nums(dataset_idx);
            pre_task = NaN;
            post_task = NaN;

            for i = 1:length(tasks)
                task = tasks{i};
                if strcmp(task.filter_type, 'Area') && ...
                length(task.area_filter) == 2 && ...
                task.across_area && ...
                strcmp(task.dataset, dataset_name)
                    fprintf('Plotting task: Dataset: %s, Prepost: %s, session num: %d\n', task.dataset, task.prepost, numel(task.sessions));
                    
                    if strcmp(task.prepost, 'Pre')
                        assert(isa(pre_task, 'double') && isnan(pre_task), 'Pre task already assigned.');
                        pre_task = task;
                    elseif strcmp(task.prepost, 'Post')
                        assert(isa(post_task, 'double') && isnan(post_task), 'Post task already assigned.');
                        post_task = task;
                    end
                end
            end
            if ~(isa(pre_task, 'struct') && isa(post_task, 'struct'))
                fprintf('Skipping dataset %s due to missing pre or post task.\n', dataset_name);
                continue;
            end

            % plot on figure
            if plot_idx == 1
                pre_x = pre_task.sync_chi_x;
                pre_y = pre_task.sync_chi_y;
                post_x = post_task.sync_chi_x;
                post_y = post_task.sync_chi_y;
                title_str = 'Chi';
            else
                pre_x = pre_task.sync_corr_x;
                pre_y = pre_task.sync_corr_y;
                post_x = post_task.sync_corr_x;
                post_y = post_task.sync_corr_y;
                % x_err = post_task.sync_corr_x_err;
                % y_err = post_task.sync_corr_y_err;
                title_str = "Pearson's r";
            end
            nan_filter = ~isnan(pre_x) & ~isnan(pre_y) & ~isnan(post_x) & ~isnan(post_y);
            point_num = sum(nan_filter);
            fprintf(' - Dataset %s: Plotting %d pairs.\n', dataset_name, point_num);
            pre_x = pre_x(nan_filter);
            pre_y = pre_y(nan_filter);
            post_x = post_x(nan_filter);
            post_y = post_y(nan_filter);

            switch dataset_name(1:6)
                case 'Slayer' % Slayer
                    color = [1, 0, 0]; % red
                case 'Zeppel' % Zeppelin
                    color = [0, 0, 1]; % blue
                case 'Empero' % Emperor
                    color = [0, 0, 0]; % black
                otherwise
                    color = [0.5, 0.5, 0.5]; % gray
            end

            scatter(pre_x, pre_y, 50, 'o', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none', 'HandleVisibility','off');
            scatter(post_x, post_y, 50, '^', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'none','HandleVisibility','off');
            for j = 1:length(pre_x)
                plot([pre_x(j), post_x(j)], [pre_y(j), post_y(j)], 'Marker', 'none', 'Color', color, 'LineWidth', 0.5, 'HandleVisibility','off');
            end
            
            if point_num > 0
                max_max = max([max_max, max(pre_x), max(pre_y), max(post_x), max(post_y)]);
                min_min = min([min_min, min(pre_x), min(pre_y), min(post_x), min(post_y)]);
            end
        end
        axis equal;
        lim_range = max_max - min_min;
        max_lim = max_max + 0.1 * lim_range;
        min_lim = min_min - 0.1 * lim_range;
        xlim([min_lim, max_lim]);
        ylim([min_lim, max_lim]);
        plot([min_lim, max_lim], [min_lim, max_lim], 'k--', 'LineWidth', 1, 'HandleVisibility','off');
        title(title_str);
        xlabel('Eyes Open');
        ylabel('Eyes Closed');
        legend('Location', 'best');
    end
    sgtitle(sprintf('Cortex ACC-VLPFC Synchrony (%s): Eyes Open vs Eyes Closed', control));
    fig_save_folder = fullfile(root, 'Figures', 'Synchrony');
    check_path(fig_save_folder);
    fig_save_path = fullfile(fig_save_folder, sprintf('synchrony_cortex_across_area_%s.png', control));
    saveas(f, fig_save_path);
    close(f);
end