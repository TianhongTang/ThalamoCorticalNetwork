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