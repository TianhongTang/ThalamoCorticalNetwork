% synchrony.m - Synchrony plot for eyes-open vs eyes-closed states

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
USE_EXISTING_RESULTS = true; % if true, load existing results if available
result_folder = fullfile(root, 'Data', 'Working', 'Analysis', 'Synchrony');
check_path(result_folder);
result_name = sprintf('synchrony_results.mat');
result_path = fullfile(result_folder, result_name);
skip = isfile(result_path) && USE_EXISTING_RESULTS;

if skip
    % load existing results
    load(result_path, 'tasks');
else
    %% register tasks
    tasks = {};

    % register MY dataset metadata
    metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
    metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
    load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums');
    prepost_types = {'Pre', 'Post'};
    
    % Remove 'SlayerNoinj' from dataset_names and dataset_num
    noinj_idx = find(contains(dataset_names, 'SlayerNoinj'));
    fprintf('Removing dataset: %s\n', dataset_names{noinj_idx});
    dataset_names(noinj_idx) = [];
    dataset_num = length(dataset_names);
    session_nums(noinj_idx) = [];

    fprintf('Total datasets after removal: %d\n', dataset_num);
    for i = 1:dataset_num
        fprintf(' - Dataset %d: %s, Sessions: %d\n', i, dataset_names{i}, session_nums(i));
    end
    fprintf('-------------------------\n');

    % thalamus within area synchrony comparison tasks
    for dataset_idx = 1:dataset_num
        dataset_name = dataset_names{dataset_idx};
        session_num = session_nums(dataset_idx);

        prepost = 'Pre';
        task = struct();
        task.dataset = dataset_name;
        task.prepost = prepost;
        task.xlabel = 'Eyes Open';
        task.ylabel = 'Eyes Closed';
        task.xdata = [dataset_name, prepost, 'RestOpen', 'Full'];
        task.ydata = [dataset_name, prepost, 'RestClose', 'Full'];
        task.across_area = false;
        task.filter_type = 'Area';
        task.area_filter = {'Thalamus'};
        task.sessions = 1:session_num;

        tasks{end+1} = task; %#ok<SAGROW>
    end

    % cortical single area synchrony comparison tasks
    single_areas = {'ACC', 'VLPFC'};
    for area_idx = 1:length(single_areas)
        for dataset_idx = 1:dataset_num
            for prepost_idx = 1:length(prepost_types)
                prepost = prepost_types{prepost_idx};
                dataset_name = dataset_names{dataset_idx};
                session_num = session_nums(dataset_idx);

                if contains(dataset_name, 'Noinj') && strcmp(prepost, 'Post')
                    % No injection session does not have post sessions 
                    continue;
                end

                task = struct();
                task.dataset = dataset_name;
                task.prepost = prepost;
                task.xlabel = 'Eyes Open';
                task.ylabel = 'Eyes Closed';
                task.xdata = [dataset_name, prepost, 'RestOpen', 'Cortex'];
                task.ydata = [dataset_name, prepost, 'RestClose', 'Cortex'];
                task.filter_type = 'Area';
                task.area_filter = single_areas(area_idx);
                task.across_area = false;
                task.sessions = 1:session_num;

                tasks{end+1} = task; %#ok<SAGROW>
            end
        end
    end

    % cortical across area synchrony comparison tasks
    for dataset_idx = 1:dataset_num
        for prepost_idx = 1:length(prepost_types)
            prepost = prepost_types{prepost_idx};
            dataset_name = dataset_names{dataset_idx};
            session_num = session_nums(dataset_idx);

            if contains(dataset_name, 'Noinj') && strcmp(prepost, 'Post')
                % No injection session does not have post sessions 
                continue;
            end

            task = struct();
            task.dataset = dataset_name;
            task.prepost = prepost;
            task.xlabel = 'Eyes Open';
            task.ylabel = 'Eyes Closed';
            task.xdata = [dataset_name, prepost, 'RestOpen', 'Cortex'];
            task.ydata = [dataset_name, prepost, 'RestClose', 'Cortex'];
            task.across_area = true;
            task.filter_type = 'Area';
            task.area_filter = {'ACC', 'VLPFC'};
            task.sessions = 1:session_num;

            tasks{end+1} = task; %#ok<SAGROW>
        end
    end

    % cortical multi-area synchrony comparison tasks
    for dataset_idx = 1:dataset_num
        for prepost_idx = 1:length(prepost_types)
            prepost = prepost_types{prepost_idx};
            dataset_name = dataset_names{dataset_idx};
            session_num = session_nums(dataset_idx);

            if contains(dataset_name, 'Noinj') && strcmp(prepost, 'Post')
                % No injection session does not have post sessions 
                continue;
            end

            task = struct();
            task.dataset = dataset_name;
            task.prepost = prepost;
            task.xlabel = 'Eyes Open';
            task.ylabel = 'Eyes Closed';
            task.xdata = [dataset_name, prepost, 'RestOpen', 'Cortex'];
            task.ydata = [dataset_name, prepost, 'RestClose', 'Cortex'];
            task.across_area = false;
            task.filter_type = 'Area';
            task.area_filter = {'ACC', 'VLPFC'};
            task.sessions = 1:session_num;

            tasks{end+1} = task; %#ok<SAGROW>
        end
    end

    % register KZ dataset metadata
    task = struct();
    task.dataset = 'KZ';
    task.prepost = 'NA';
    task.xlabel = 'Eyes Open';
    task.ylabel = 'Eyes Closed';
    task.xdata = 'KZRestOpen';
    task.ydata = 'KZRestClose';
    task.filter_type = 'Area';
    task.area_filter = {'ThalamusAnterior'};
    task.across_area = false;

    % load session filters
    filter_path = fullfile(root, 'Data', 'Working', 'Meta', 'session_filters_KZ.mat');
    load(filter_path, 'property_filters');
    % add session indices for filtered sessions
    filter = [property_filters.combined];
    task.sessions = find(filter);
    tasks{end+1} = task;


    %% run tasks: calculate synchrony and save to task struct
    myGaussian = @(x, mu, sigma) exp(-((x - mu).^2) / (2*sigma^2)) / (sigma*sqrt(2*pi));
    smooth_kernel = myGaussian(-200:200, 0, 50); 
    smooth_kernel = smooth_kernel / sum(smooth_kernel);

    task_num = length(tasks);
    for task_idx = 1:task_num
        task = tasks{task_idx};
        fprintf('=========================\n');
        fprintf('Task %d/%d: Dataset: %s, Prepost: %s\n', task_idx, task_num, task.dataset, task.prepost);
        fprintf(' X: %s, Y: %s\n', task.xdata, task.ydata);
        fprintf(' Filter Type: %s, Area Filter: %s, Across Area: %d\n', task.filter_type, strjoin(task.area_filter, ','), task.across_area);
        fprintf(' Sessions to process: %d\n', length(task.sessions));
        fprintf('-------------------------\n');

        % target stats
        sync_chi_x = NaN(length(task.sessions), 1);
        sync_chi_y = NaN(length(task.sessions), 1);
        sync_corr_x = NaN(length(task.sessions), 1);
        sync_corr_y = NaN(length(task.sessions), 1);
        sync_corr_x_err = NaN(length(task.sessions), 1);
        sync_corr_y_err = NaN(length(task.sessions), 1);

        for i = 1:length(task.sessions)
            session_idx = task.sessions(i);
            fprintf(' - Session %d/%d: %d\n', i, length(task.sessions), session_idx);
            % load x and y data
            ax_data_list = {task.xdata, task.ydata};
            for ax_idx = 1:2
                raster_folder = fullfile(root, 'Data', 'Working', 'raster');
                data_name = sprintf('raster_%s_%d.mat', ax_data_list{ax_idx}, session_idx);
                data_path = fullfile(raster_folder, data_name);
                load(data_path, 'rasters', 'N', 'trial_num', 'cell_area');

                if ax_idx == 1
                    fprintf('   - X Data: %s\n', data_name);
                else
                    fprintf('   - Y Data: %s\n', data_name);
                end

                if numel(rasters) == 0
                    fprintf('   - No trials found. Skipping...\n');
                    if ax_idx == 1
                        sync_chi_x(i) = NaN;
                        sync_corr_x(i) = NaN;
                        sync_corr_x_err(i) = NaN;
                    else
                        sync_chi_y(i) = NaN;
                        sync_corr_y(i) = NaN;
                        sync_corr_y_err(i) = NaN;
                    end

                    continue;
                end
                fprintf('   - %d trials found, %d rasters.\n', trial_num, numel(rasters));

                % smooth and concatenate trials
                smoothed_rasters = cellfun(@(x) conv2(x, smooth_kernel, 'valid'), rasters, 'UniformOutput', false);
                total_length = sum(cellfun(@(x) size(x, 2), rasters));
                total_length_smooth = sum(cellfun(@(x) size(x, 2), smoothed_rasters));
                raster_concat = zeros(N, total_length);
                raster_smooth_concat = zeros(N, total_length_smooth);
                pointer = 1;
                pointer_smooth = 1;
                for t = 1:trial_num
                    t_len = size(rasters{t}, 2);
                    raster_concat(:, pointer:(pointer + t_len - 1)) = rasters{t};
                    pointer = pointer + t_len;

                    t_len_smooth = size(smoothed_rasters{t}, 2);
                    raster_smooth_concat(:, pointer_smooth:(pointer_smooth + t_len_smooth - 1)) = smoothed_rasters{t};
                    pointer_smooth = pointer_smooth + t_len_smooth;
                end

                if ~task.across_area
                    % within area synchrony
                    filter = true(1, N);
                    if strcmp(task.filter_type, 'Area')
                        area_filter = task.area_filter;
                        % if cell_area is in area_filter
                        filter_area = cellfun(@(x) any(strcmp(x, area_filter)), cell_area); 
                        filter = filter & filter_area;
                    end

                    raster_concat = raster_concat(filter, :);
                    raster_smooth_concat = raster_smooth_concat(filter, :);
                    N_filtered = size(raster_concat, 1);

                    if N_filtered == 0
                        chi = NaN;
                        corr = NaN;
                        corr_err = NaN;
                        fprintf(' - Session %d, Ax %d: No neurons after filtering. Skipping...\n', session_idx, ax_idx);
                    else
                        chi = calc_synchrony(raster_smooth_concat);
                        [corr, corr_err] = calc_synchrony_corr_within(raster_smooth_concat);
                    end
                    
                else
                    % across area synchrony
                    assert(length(task.area_filter) == 2, 'Area filter must have two areas for across area synchrony.');
                    area1 = task.area_filter{1};
                    area2 = task.area_filter{2};
                    filter_area1 = cellfun(@(x) strcmp(x, area1), cell_area);
                    filter_area2 = cellfun(@(x) strcmp(x, area2), cell_area);
                    raster_area1 = raster_smooth_concat(filter_area1, :);
                    raster_area2 = raster_smooth_concat(filter_area2, :);
                    N_filtered1 = size(raster_area1, 1);
                    N_filtered2 = size(raster_area2, 1);
                    if N_filtered1 == 0 || N_filtered2 == 0
                        corr = NaN;
                        corr_err = NaN;
                        chi = NaN;
                        fprintf(' - Session %d, Ax %d: No neurons after filtering. Skipping...\n', session_idx, ax_idx);
                    else
                        [corr, corr_err] = calc_synchrony_corr_across(raster_area1, raster_area2);
                        chi = NaN; % not defined for across area synchrony
                    end
                end

                if ax_idx == 1
                    sync_chi_x(i) = chi;
                    sync_corr_x(i) = corr;
                    sync_corr_x_err(i) = corr_err;
                else
                    sync_chi_y(i) = chi;
                    sync_corr_y(i) = corr;
                    sync_corr_y_err(i) = corr_err;
                end
            end
        end
        % save to task struct
        task.sync_chi_x = sync_chi_x;
        task.sync_chi_y = sync_chi_y;
        task.sync_corr_x = sync_corr_x;
        task.sync_corr_y = sync_corr_y;
        task.sync_corr_x_err = sync_corr_x_err;
        task.sync_corr_y_err = sync_corr_y_err;
        tasks{task_idx} = task;
    end

    % save results
    result_folder = fullfile(root, 'Data', 'Working', 'Analysis', 'Synchrony');
    check_path(result_folder);
    result_name = sprintf('synchrony_results.mat');
    result_path = fullfile(result_folder, result_name);
    save(result_path, 'tasks', '-v7.3');
end

%% Plot figures

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
    % scatter(x, y, 50, 'o', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [1,0,0], 'MarkerEdgeColor', 'none', 'DisplayName', 'Slayer');
    % scatter(x, y, 50, 'o', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [0,0,1], 'MarkerEdgeColor', 'none', 'DisplayName', 'Zeppelin');
    scatter(x, y, 50, 's', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', 'none', 'DisplayName', 'Emperor Saline');
    scatter(x, y, 50, 'o', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', 'none', 'DisplayName', 'Emperor Muscimol');

    for i = 1:length(tasks)
        task = tasks{i};
        if strcmp(task.filter_type, 'Area') && ...
        length(task.area_filter) == 1 && ...
        strcmp(task.area_filter{1}, 'Thalamus') && ...
        ~task.across_area && strcmp(task.dataset(1:6), 'Empero')
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
            if plot_idx == 2
                x_err = x_err(nan_filter);
                y_err = y_err(nan_filter);
                x_neg = x - x_err;
                x_pos = x + x_err;
                y_neg = y - y_err;
                y_pos = y + y_err;
                errorbar(x, y, y_err, y_err, x_err, x_err, 'Color',...
                 [0, 0, 0, 0.3], 'LineStyle', 'none', 'CapSize', 5, 'HandleVisibility','off');
            end

            % markers for different datasets
            switch task.dataset(1:6)
                case 'Slayer' % Slayer
                    color = [1, 0, 0]; % red
                case 'Zeppel' % Zeppelin
                    color = [0, 0, 1]; % blue
                case 'Empero' % Emperor
                    color = [0, 0, 0]; % black
                otherwise 
                    color = [0.5, 0.5, 0.5]; % gray
            end
            switch task.dataset(end-2:end)
                case 'Mus' % Muscimol
                    marker = 'o'; % circle
                case 'Sal' % Saline
                    marker = 's'; % square
                case 'inj' % No injection
                    marker = 'd'; % diamond
                otherwise
                    marker = '^'; % triangle
            end
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