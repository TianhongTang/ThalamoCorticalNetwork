%% Merge session data from multiple areas

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
ALLOW_MISSING_AREA = true;
SKIP_EXISTING = false;

% load data
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums', 'cortex_files', 'thalamus_files', 'eyeID_files');

areas = {'ACC', 'VLPFC', 'Thalamus'};

%% register tasks
merge_types = {'Full', 'Cortex'};
prepost_types = {'Pre', 'Post'};
states = {'RestOpen', 'RestClose', 'Task'};

tasks = {};
for dataset_idx = 1:dataset_num
    dataset_name = dataset_names{dataset_idx};
    session_num = session_nums(dataset_idx);
    for session_idx = 1:session_num
        date = thalamus_files{dataset_idx}{session_idx}(1:8);
        for merge_idx = 1:length(merge_types)
            merge_type = merge_types{merge_idx};
            for prepost_idx = 1:length(prepost_types)
                prepost = prepost_types{prepost_idx};
                if contains(dataset_name, 'Noinj') && strcmp(prepost, 'Post')
                    % No injection session does not have post sessions 
                    continue;
                end
                if strcmp(merge_type, 'Full') && strcmp(prepost, 'Post')
                    % All thalamus data in post sessions are not available
                    continue;
                end
                for state_idx = 1:length(states)
                    state = states{state_idx};
                    if strcmp(dataset_name, 'SlayerNoinj') && (strcmp(state, 'RestOpen') || strcmp(state, 'RestClose'))
                        % SlayerNoinj does not have RestOpen and RestClose sessions
                        continue;
                    end

                    [animal_name, injection] = split_dataset_name(dataset_name);

                    % construct tasks
                    task              = struct();
                    task.dataset_name = dataset_name;
                    task.animal_name  = animal_name;
                    task.injection    = injection;
                    task.session_idx  = session_idx;
                    task.merge_type   = merge_type;
                    task.prepost      = prepost;
                    task.state        = state;
                    task.date         = date;
                    if strcmp(merge_type, 'Cortex')
                        task.areas = {'ACC', 'VLPFC'};
                    elseif strcmp(merge_type, 'Full')
                        task.areas = {'ACC', 'VLPFC', 'Thalamus'};
                    end
                    tasks{end+1} = task; %#ok<SAGROW>
                end
            end
        end
    end
end

%% run tasks
task_num = numel(tasks);
total_tic = tic;
for task_idx = 1:task_num
    task         = tasks{task_idx};
    dataset_name = task.dataset_name;
    session_idx  = task.session_idx;
    merge_type   = task.merge_type;
    prepost      = task.prepost;
    state        = task.state;
    areas        = task.areas;
    area_num     = numel(areas);

    meta             = struct();
    meta.animal_name = task.animal_name;
    meta.injection   = task.injection;
    meta.prepost     = prepost;
    meta.state       = state;
    meta.area        = merge_type;
    meta.align       = 'None';
    meta.session_idx = session_idx;
    meta.date        = task.date;

    % check if merged file already exists
    save_folder = fullfile(root, 'Data', 'Working', 'raster');
    check_path(save_folder);
    save_name = generate_filename('raster', meta);
    % session_name = sprintf('%s_%d', [dataset_name, prepost, state, merge_type], session_idx);
    save_path = fullfile(save_folder, save_name);

    if isfile(save_path) && SKIP_EXISTING
        fprintf('Merged file already exists for Task %d/%d: session%d, %s. Skipping...\n', ...
            task_idx, task_num, session_idx, [dataset_name, prepost, state, merge_type]);
        continue;
    end

    fprintf('===================\n');
    fprintf('Merging Task %d/%d: session%d, %s...\n\n', ...
        task_idx, task_num, session_idx, [dataset_name, prepost, state, merge_type]);

    tic;
    N           = 0;
    dt          = NaN;
    cell_id     = cell(1, 0);   % (1, N)
    cell_area   = cell(1, 0);   % (1, N)
    cuetype     = zeros(1, 0);  % (1, trial_num)
    trial_num   = NaN;
    trial_len   = NaN;
    borders     = [];           % area borders. Marking the beginning of each area.
    sort_ranges = zeros(0, 2);
    sort_idx    = zeros(1, 0);

    for area_idx = 1:area_num
        area = areas{area_idx};
        area_meta = meta;
        area_meta.area = area;

        % load area data
        folder_name = fullfile(root, 'Data', 'Working', 'raster');
        file_name = generate_filename('raster', area_meta);
        area_file = fullfile(folder_name, file_name);

        borders = [borders, N+1]; %#ok<AGROW>
        % check if data exists
        if ~isfile(area_file)
            if ALLOW_MISSING_AREA
                fprintf('Warning: Missing file %s. Skipping area %s.\n', area_file, area);
                continue;
            else
                error('File not found: %s', area_file); %#ok<UNRCH>
            end
        end
        area_data = load(area_file);

        if area_data.meta.N==0
            continue;
        end

        if isnan(trial_num)
            dt           = area_data.meta.dt;
            trial_num    = area_data.meta.trial_num;
            trial_len    = area_data.data.trial_len;
            cuetype      = area_data.data.cuetype;
            rasters      = cell(1, trial_num);
            firing_rates = cell(1, trial_num);
            spikes       = cell(0, trial_num);
            channel      = zeros(1, 0);
            for i=1:trial_num
                rasters{1, i} = zeros(0, trial_len(i));
                firing_rates{1, i} = zeros(0, 1);
            end
        else
            assert(trial_num == area_data.meta.trial_num, ...
                'trial num not match! Area: %s, Dataset: %s, state: %s', area, dataset_name, state);
            if strcmp(state, 'Task')
                assert(~((any(cuetype~=area_data.data.cuetype & ~isnan(cuetype))) || any(trial_len ~= area_data.data.trial_len)), ...
                    'trial info not match! Area: %s, Dataset: %s, state: %s', area, dataset_name, state);
            end
        end

        cell_area   = [cell_area, area_data.data.cell_area];             %#ok<AGROW>
        cell_id     = [cell_id, area_data.data.cell_id];                 %#ok<AGROW>
        spikes      = [spikes;area_data.data.spikes];                    %#ok<AGROW>
        channel     = [channel, area_data.data.channel];                 %#ok<AGROW>
        sort_ranges = [sort_ranges; N+1, N+area_data.meta.N];            %#ok<AGROW>
        [~, sort_idx(N+1:N+area_data.meta.N)] = sort(area_data.data.channel);
        sort_idx(N+1:N+area_data.meta.N)      = sort_idx(N+1:N+area_data.meta.N) + N;
        
        N = N+area_data.meta.N;
        for i=1:trial_num
            rasters{1, i} = [rasters{1, i}; area_data.data.rasters{1, i}];
            firing_rates{1, i} = [firing_rates{1, i}; area_data.data.firing_rates{1, i}];
        end
    end

    % construct meta and data struct for merged session
    meta.N           = N;
    meta.dt          = dt;
    meta.trial_num   = trial_num;
    meta.file_name   = generate_filename('raster', meta);
    meta.trial_len   = sort(trial_len, 'descend');
    meta.max_len     = max(trial_len);
    meta.min_len     = min(trial_len);
    meta.total_len   = sum(trial_len);

    data = struct();
    data.rasters      = rasters;
    data.spikes       = spikes;
    data.trial_len    = trial_len;
    data.cell_id      = cell_id;
    data.cell_area    = cell_area;
    data.channel      = channel;
    data.cuetype      = cuetype;
    data.firing_rates = firing_rates;

    % save data
    save_folder = fullfile(root, 'Data', 'Working', 'raster');
    check_path(save_folder);
    save_name = meta.file_name;
    save_path = fullfile(save_folder, save_name);
    save(save_path, "meta", "data", "-v7.3");

    % construct border and sort index data for merged session
    border_meta = struct();
    border_meta.animal_name = task.animal_name;
    border_meta.injection   = task.injection;
    border_meta.prepost     = prepost;
    border_meta.area        = merge_type;
    border_meta.align       = 'None';
    border_meta.session_idx = session_idx;
    border_meta.N           = N;
    border_meta.area_num    = area_num;
    border_meta.file_name   = generate_filename('border', border_meta);

    border_data = struct();
    border_data.borders = borders;

    sortidx_meta = struct();
    sortidx_meta.animal_name = task.animal_name;
    sortidx_meta.injection   = task.injection;
    sortidx_meta.prepost     = prepost;
    sortidx_meta.state       = state;
    sortidx_meta.area        = merge_type;
    sortidx_meta.align       = 'None';
    sortidx_meta.session_idx = session_idx;
    sortidx_meta.kernel      = 'None';
    sortidx_meta.criterion   = 'channel';
    sortidx_meta.N           = N;
    sortidx_meta.file_name   = generate_filename('sortidx', sortidx_meta);

    sortidx_data = struct();
    sortidx_data.sort_idx = sort_idx;

    % save border file
    meta = border_meta; data = border_data;
    save_folder = fullfile(root, 'Data', 'Working', 'border');
    check_path(save_folder);
    save_name = meta.file_name;
    save_path = fullfile(save_folder, save_name);
    save(save_path, "meta", "data", "-v7.3");

    % save sort index file
    meta = sortidx_meta; data = sortidx_data;
    save_folder = fullfile(root, 'Data', 'Working', 'sortidx');
    check_path(save_folder);
    save_name = meta.file_name;
    save_path = fullfile(save_folder, save_name);
    save(save_path, "meta", "data", "-v7.3");
    
    toc;
end
fprintf('Total merging time: %.2f seconds.\n', toc(total_tic));
