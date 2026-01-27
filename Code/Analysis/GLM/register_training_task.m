% register_training_task.m - Register training tasks for GLM analysis


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
%% Slayer sessions
task_name = 'Slayer';
% load metadata
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums', 'cortex_files', 'thalamus_files', 'eyeID_files');

merge_types = {'Full', 'Cortex'};
prepost_types = {'Pre', 'Post'};
% states = {'RestOpen', 'RestClose', 'Task'};
states = {'RestOpen', 'RestClose'};


tasks = cell(0, 1);
task_idx = 0;
for dataset_idx = [1, 2, 3]
    dataset_name = dataset_names{dataset_idx};
    session_num = session_nums(dataset_idx);
    for session_idx = 1:session_num
        for merge_idx = 1:length(merge_types)
            merge_type = merge_types{merge_idx};
            for prepost_idx = 1:length(prepost_types)
                prepost = prepost_types{prepost_idx};
                for state_idx = 1:length(states)
                    state = states{state_idx};

                    if contains(dataset_name, 'Noinj') && strcmp(prepost, 'Post')
                        % No injection session does not have post sessions 
                        continue;
                    end
                    if strcmp(merge_type, 'Full') && strcmp(prepost, 'Post')
                        % All thalamus data in post sessions are not available
                        continue;
                    end
                    if strcmp(dataset_name, 'SlayerNoinj') && contains(state, 'Rest')
                        % SlayerNoinj does not have RestOpen and RestClose sessions
                        continue;
                    end

                    task = struct();
                    task_idx = task_idx + 1;

                    task.state = state;
                    task.dataset_name = sprintf('%s%s%s%sAlignLast', dataset_name, prepost, state, merge_type);
                    task.border_name = sprintf('%s%s', dataset_name, merge_type);
                    task.session_name = sprintf('%s_%d', task.dataset_name, session_idx);
                    task.session_idx = session_idx;

                    % confirm raster file exists
                    raster_file_folder = fullfile(root, 'Data', 'Working', 'raster');
                    raster_file_name = sprintf('raster_%s_%d.mat', task.dataset_name, session_idx);
                    raster_file_path = fullfile(raster_file_folder, raster_file_name);
                    if ~isfile(raster_file_path)
                        throw(MException('RegisterTrainingTask:RasterFileNotFound', ...
                            'Raster file does not exist: %s. Skip this task.', raster_file_path));
                        % fprintf('Raster file does not exist: %s. Skip this task.\n', raster_file_path);
                        % task_idx = task_idx - 1;
                        % continue;
                    end

                    % training configurations
                    config = struct();
                    config.max_epochs = 3000;
                    config.save_interval = 100;
                    config.kernel = 'DeltaPure';
                    config.crossval_fold_num = 3;
                    reg = struct();
                    reg.l1 = 0;
                    reg.l2 = 1;
                    reg.name = 'L2=0.2';
                    config.reg = reg;
                    config.shuffle_size = 0;
                    task.config = config;

                    tasks{task_idx} = task;
                end
            end
        end
    end
end

task_file_folder = fullfile(root, 'Data', 'Working', 'Meta', 'training_tasks');
check_path(task_file_folder);
task_file_name = sprintf('training_task_%s.mat', task_name);
task_file_path = fullfile(task_file_folder, task_file_name);
save(task_file_path, 'tasks');



%% Zeppelin sessions
task_name = 'Zeppelin';
% load metadata
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums', 'cortex_files', 'thalamus_files', 'eyeID_files');

merge_types = {'Full', 'Cortex'};
prepost_types = {'Pre', 'Post'};
% states = {'RestOpen', 'RestClose', 'Task'};
states = {'RestOpen', 'RestClose'};


tasks = cell(0, 1);
task_idx = 0;
for dataset_idx = [4, 5, 6]
    dataset_name = dataset_names{dataset_idx};
    session_num = session_nums(dataset_idx);
    for session_idx = 1:session_num
        for merge_idx = 1:length(merge_types)
            merge_type = merge_types{merge_idx};
            for prepost_idx = 1:length(prepost_types)
                prepost = prepost_types{prepost_idx};
                for state_idx = 1:length(states)
                    state = states{state_idx};

                    if contains(dataset_name, 'Noinj') && strcmp(prepost, 'Post')
                        % No injection session does not have post sessions 
                        continue;
                    end
                    if strcmp(merge_type, 'Full') && strcmp(prepost, 'Post')
                        % All thalamus data in post sessions are not available
                        continue;
                    end

                    task = struct();
                    task_idx = task_idx + 1;

                    task.state = state;
                    task.dataset_name = sprintf('%s%s%s%sAlignLast', dataset_name, prepost, state, merge_type);
                    task.border_name = sprintf('%s%s', dataset_name, merge_type);
                    task.session_name = sprintf('%s_%d', task.dataset_name, session_idx);
                    task.session_idx = session_idx;

                    % confirm raster file exists
                    raster_file_folder = fullfile(root, 'Data', 'Working', 'raster');
                    raster_file_name = sprintf('raster_%s_%d.mat', task.dataset_name, session_idx);
                    raster_file_path = fullfile(raster_file_folder, raster_file_name);
                    if ~isfile(raster_file_path)
                        throw(MException('RegisterTrainingTask:RasterFileNotFound', ...
                            'Raster file does not exist: %s. Skip this task.', raster_file_path));
                        % fprintf('Raster file does not exist: %s. Skip this task.\n', raster_file_path);
                        % task_idx = task_idx - 1;
                        % continue;
                    end

                    % training configurations
                    config = struct();
                    config.max_epochs = 3000;
                    config.save_interval = 100;
                    config.kernel = 'DeltaPure';
                    config.crossval_fold_num = 3;
                    reg = struct();
                    reg.l1 = 0;
                    reg.l2 = 1;
                    reg.name = 'L2=0.2';
                    config.reg = reg;
                    config.shuffle_size = 0;
                    task.config = config;

                    tasks{task_idx} = task;
                end
            end
        end
    end
end

task_file_folder = fullfile(root, 'Data', 'Working', 'Meta', 'training_tasks');
check_path(task_file_folder);
task_file_name = sprintf('training_task_%s.mat', task_name);
task_file_path = fullfile(task_file_folder, task_file_name);
save(task_file_path, 'tasks');



%% Emperor sessions
task_name = 'Emperor';
% load metadata
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums', 'cortex_files', 'thalamus_files', 'eyeID_files');

merge_types = {'Full', 'Cortex'};
prepost_types = {'Pre', 'Post'};
% states = {'RestOpen', 'RestClose', 'Task'};
states = {'RestOpen', 'RestClose'};


tasks = cell(0, 1);
task_idx = 0;
for dataset_idx = [7, 8, 9]
    dataset_name = dataset_names{dataset_idx};
    session_num = session_nums(dataset_idx);
    for session_idx = 1:session_num
        for merge_idx = 1:length(merge_types)
            merge_type = merge_types{merge_idx};
            for prepost_idx = 1:length(prepost_types)
                prepost = prepost_types{prepost_idx};
                for state_idx = 1:length(states)
                    state = states{state_idx};

                    if contains(dataset_name, 'Noinj') && strcmp(prepost, 'Post')
                        % No injection session does not have post sessions 
                        continue;
                    end
                    if strcmp(merge_type, 'Full') && strcmp(prepost, 'Post')
                        % All thalamus data in post sessions are not available
                        continue;
                    end

                    task = struct();
                    task_idx = task_idx + 1;

                    task.state = state;
                    task.dataset_name = sprintf('%s%s%s%sAlignLast', dataset_name, prepost, state, merge_type);
                    task.border_name = sprintf('%s%s', dataset_name, merge_type);
                    task.session_name = sprintf('%s_%d', task.dataset_name, session_idx);
                    task.session_idx = session_idx;

                    % confirm raster file exists
                    raster_file_folder = fullfile(root, 'Data', 'Working', 'raster');
                    raster_file_name = sprintf('raster_%s_%d.mat', task.dataset_name, session_idx);
                    raster_file_path = fullfile(raster_file_folder, raster_file_name);
                    if ~isfile(raster_file_path)
                        throw(MException('RegisterTrainingTask:RasterFileNotFound', ...
                            'Raster file does not exist: %s. Skip this task.', raster_file_path));
                        % fprintf('Raster file does not exist: %s. Skip this task.\n', raster_file_path);
                        % task_idx = task_idx - 1;
                        % continue;
                    end

                    % training configurations
                    config = struct();
                    config.max_epochs = 3000;
                    config.save_interval = 100;
                    config.kernel = 'DeltaPure';
                    config.crossval_fold_num = 3;
                    reg = struct();
                    reg.l1 = 0;
                    reg.l2 = 1;
                    reg.name = 'L2=0.2';
                    config.reg = reg;
                    config.shuffle_size = 0;
                    task.config = config;

                    tasks{task_idx} = task;
                end
            end
        end
    end
end

task_file_folder = fullfile(root, 'Data', 'Working', 'Meta', 'training_tasks');
check_path(task_file_folder);
task_file_name = sprintf('training_task_%s.mat', task_name);
task_file_path = fullfile(task_file_folder, task_file_name);
save(task_file_path, 'tasks');

% %% KZ dataset

% task_name = 'KZ';
% % load metadata
% meta_folder = fullfile(root, 'Data', 'Working', 'Meta');
% check_path(meta_folder);
% meta_file_name = 'all_session_info_KZ.mat';
% meta_file_path = fullfile(meta_folder, meta_file_name);
% load(meta_file_path, 'all_session_info', 'segmentNames');

% session_num = length(all_session_info);
% area_num = length(segmentNames);
% fprintf('Total number of sessions: %d\n', session_num);

% % load session filter
% session_filter_file = fullfile(root, 'Data', 'Working', 'Meta', 'session_filters_KZ.mat');
% load(session_filter_file, 'property_filters');
% applied_filter = property_filters.combined_filter;

% filtered_sessions = all_session_info(applied_filter);
% filtered_session_num = length(filtered_sessions);
% filtered_session_idxs = find(applied_filter);
% fprintf('Number of sessions after filtering: %d\n', filtered_session_num);

% tasks = cell(filtered_session_num, 1);
% task_idx = 0;

% states = {'RestOpen', 'RestClose'};

% for session_idx = 1:filtered_session_num
%     % data information
%     session_info = filtered_sessions(session_idx);

%     for state_idx = 1:length(states)
%         state = states{state_idx};
%         task = struct();
%         task_idx = task_idx + 1;

%         task.session_name = session_info.sessionname;
%         task.monkey_name = session_info.monkeyName;
%         task.session_length = session_info.sessionlength;
%         task.state = state;
%         task.dataset_name = sprintf('KZ%s', state);
%         task.session_idx = filtered_session_idxs(session_idx);

%         % training configurations
%         config = struct();
%         config.max_epochs = 5000;
%         config.save_interval = 100;
%         config.kernel = 'DeltaPure';
%         config.crossval_fold_num = 3;
%         reg = struct();
%         reg.l1 = 0;
%         reg.l2 = 0.2;
%         reg.name = 'L2=0_2';
%         config.reg = reg;
%         config.shuffle_size = 0;
%         task.config = config;

%         tasks{task_idx} = task;
%     end
% end

% task_file_folder = fullfile(root, 'Data', 'Working', 'Meta', 'training_tasks');
% check_path(task_file_folder);
% task_file_name = sprintf('training_task_%s.mat', task_name);
% task_file_path = fullfile(task_file_folder, task_file_name);
% save(task_file_path, 'tasks');
