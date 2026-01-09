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


%% Zeppelin sessions


%% Emperor sessions


%% KZ dataset

task_name = 'KZ';
% load metadata
meta_folder = fullfile(root, 'Data', 'Working', 'Meta');
check_path(meta_folder);
meta_file_name = 'all_session_info_KZ.mat';
meta_file_path = fullfile(meta_folder, meta_file_name);
load(meta_file_path, 'all_session_info', 'segmentNames');

session_num = length(all_session_info);
area_num = length(segmentNames);
fprintf('Total number of sessions: %d\n', session_num);

% load session filter
session_filter_file = fullfile(root, 'Data', 'Working', 'Meta', 'session_filters_KZ.mat');
load(session_filter_file, 'session_filters');
applied_filter = session_filters.combined_filter;

filtered_sessions = all_session_info(applied_filter);
filtered_session_num = length(filtered_sessions);
filtered_session_idxs = find(applied_filter);
fprintf('Number of sessions after filtering: %d\n', filtered_session_num);

tasks = cell(filtered_session_num, 1);
task_idx = 0;

states = {'RestOpen', 'RestClose'};

for session_idx = 1:filtered_session_num
    % data information
    session_info = filtered_sessions(session_idx);

    for state_idx = 1:length(states)
        state = states{state_idx};
        task = struct();
        task_idx = task_idx + 1;

        task.session_name = session_info.sessionname;
        task.monkey_name = session_info.monkeyName;
        task.session_length = session_info.sessionlength;
        task.state = state;
        task.dataset_name = sprintf('KZ%s', state);
        task.session_idx = filtered_session_idxs(session_idx);

        % training configurations
        config = struct();
        config.max_epochs = 2500;
        config.save_interval = 100;
        config.kernel = 'DeltaPure';
        reg = struct();
        reg.l1 = 0;
        reg.l2 = 0.2;
        reg.name = 'L2=0_2';
        config.reg = reg;
        config.shuffle_size = 0;
        task.config = config;

        tasks{task_idx} = task;
    end
end

task_file_folder = fullfile(root, 'Data', 'Working', 'Meta', 'training_tasks');
check_path(task_file_folder);
task_file_name = sprintf('training_task_%s.mat', task_name);
task_file_path = fullfile(task_file_folder, task_file_name);
save(task_file_path, 'tasks');
