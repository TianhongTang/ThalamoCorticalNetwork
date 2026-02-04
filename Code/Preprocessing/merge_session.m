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

                    % construct tasks
                    task = struct();
                    task.dataset_name = dataset_name;
                    task.session_idx = session_idx;
                    task.merge_type = merge_type;
                    task.prepost = prepost;
                    task.state = state;
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
    task = tasks{task_idx};
    dataset_name = task.dataset_name;
    session_idx = task.session_idx;
    merge_type = task.merge_type;
    prepost = task.prepost;
    state = task.state;
    areas = task.areas;
    area_num = numel(areas);

    fprintf('===================\n');
    fprintf('Merging Task %d/%d: session%d, %s...\n\n', ...
        task_idx, task_num, session_idx, [dataset_name, prepost, state, merge_type]);

    tic;
    N=0;
    cell_id = cell(1, 0); % (1, N)
    cell_area = cell(1, 0); % (1, N)
    cuetype = zeros(1, 0); % (1, trial_num)
    trial_num = NaN;
    trial_len = NaN;
    borders = []; % area borders. Marking the beginning of each area.
    sort_ranges = zeros(0, 2);
    sort_idx = zeros(1, 0);

    for area_idx = 1:area_num
        area = areas{area_idx};
        % load area data
        folder_name = fullfile(root, 'Data', 'Working', 'raster');
        file_name = sprintf('raster_%s%s%s%s_%d.mat', ...
            dataset_name, prepost, state, area, session_idx);
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
        data = load(area_file);
        % "rasters", "spikes", "firing_rates", "trial_num", "trial_len", ...
        % "session_name_full", "N", "cuetype", "cell_id", "channel", 'dt'

        if data.N==0
            continue;
        end

        if isnan(trial_num)
            trial_num = data.trial_num;
            trial_len = data.trial_len;
            cuetype = data.cuetype;
            rasters = cell(1, trial_num);
            firing_rates = cell(1, trial_num);
            spikes = cell(0, trial_num);
            channel = zeros(1, 0);
            for i=1:trial_num
                rasters{1, i} = zeros(0, trial_len(i));
                firing_rates{1, i} = zeros(0, 1);
            end
        else
            assert(trial_num == data.trial_num, ...
                'trial num not match! Area: %s, Dataset: %s, state: %s', area, dataset_name, state);
            if strcmp(state, 'Task')
                assert(~((any(cuetype~=data.cuetype & ~isnan(cuetype))) || any(trial_len ~= data.trial_len)), ...
                    'trial info not match! Area: %s, Dataset: %s, state: %s', area, dataset_name, state);
            end
        end
        cell_area = [cell_area, repmat({area}, 1, data.N)]; %#ok<AGROW>
        cell_id = [cell_id, data.cell_id]; %#ok<AGROW>
        spikes = [spikes;data.spikes]; %#ok<AGROW>
        channel = [channel, data.channel]; %#ok<AGROW>
        sort_ranges = [sort_ranges; N+1, N+data.N]; %#ok<AGROW>
        [~, sort_idx(N+1:N+data.N)] = sort(data.channel);
        
        N = N+data.N;
        for i=1:trial_num
            rasters{1, i} = [rasters{1, i}; data.rasters{1, i}];
            firing_rates{1, i} = [firing_rates{1, i}; data.firing_rates{1, i}];
        end
    end

    % save data
    save_folder = fullfile(root, 'Data', 'Working', 'raster');
    check_path(save_folder);
    save_name = sprintf('raster_%s_%d.mat', [dataset_name, prepost, state, merge_type], session_idx);
    session_name = sprintf('%s_%d', [dataset_name, prepost, state, merge_type], session_idx);
    save_path = fullfile(save_folder, save_name);
    save(save_path, 'N', "cell_id", "cuetype", "firing_rates", "trial_num",...
        "rasters", "spikes", "session_name", "trial_len", "cell_area", "channel", '-v7.3');
    % save border file
    save_folder = fullfile(root, 'Data', 'Working', 'border');
    check_path(save_folder);
    save_name = sprintf('borders_%s%s_%d.mat', dataset_name, merge_type, session_idx);
    save_path = fullfile(save_folder, save_name);
    save(save_path, "borders");
    % save sort index file
    save_folder = fullfile(root, 'Data', 'Working', 'sort_idx');
    check_path(save_folder);
    save_name = sprintf('sortidx_%s%s_%d.mat', dataset_name, merge_type, session_idx);
    save_path = fullfile(save_folder, save_name);
    save(save_path, "area_num", "areas", "sort_ranges", "sort_idx", "-v7.3");
    
    toc;
end
fprintf('Total merging time: %.2f seconds.\n', toc(total_tic));
