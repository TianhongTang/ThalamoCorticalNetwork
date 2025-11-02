%% Align all states to the same duration

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

% Load data
controls = {'Muscimol', 'Saline'};
session_idxs_all = {1:10, 1:5}; % session indices for each control
area_types = {'Full', 'Cortex'};
prepost_types = {'Pre', 'Post'};
states = {'RestOpen', 'RestClose', 'Task'};
kernel = 'DeltaPure';

tasks = {};
% register tasks
for control_idx = 1:length(controls)
    control = controls{control_idx};
    sessions = session_idxs_all{control_idx};
    for session_idx = sessions
        for area_idx = 1:length(area_types)
            area_type = area_types{area_idx};
            task = struct();
            task.control = control;
            task.area_type = area_type;
            task.session_idx = session_idx;

            prepost_states = {};
            for prepost_idx = 1:length(prepost_types)
                prepost = prepost_types{prepost_idx};
                if strcmp(area_type, 'Full') && strcmp(prepost, 'Post')
                    % All thalamus data in post sessions are not available
                    continue;
                end
                for state_idx = 1:length(states)
                    state = states{state_idx};
                    prepost_states{end+1} = struct('prepost', prepost, 'state', state);
                end
            end
            task.prepost_states = prepost_states;
            tasks{end+1} = task;
        end
    end
end

% run tasks
task_num = length(tasks);
for task_idx = 1:task_num
    task = tasks{task_idx};
    control = task.control;
    session_idx = task.session_idx;
    area_type = task.area_type;
    prepost_states = task.prepost_states;
    fprintf('Task %d/%d: %s Session %d Area %s with %d prepost_states\n', ...
        task_idx, task_num, control, session_idx, area_type, length(prepost_states));
    
    % find min duration among prepost_states
    min_duration = 10000000;
    for ps_idx = 1:length(prepost_states)
        prepost = prepost_states{ps_idx}.prepost;
        state = prepost_states{ps_idx}.state;
        folder_name = fullfile(root, 'Data', 'Working', 'raster');
    end
end

% check data size
for session = sessions
    min_duration = 10000000;
    for task_idx = 1:length(tasks)
        task = tasks{task_idx};

        file_path = ['../GLM_data/', control, task, '/GLMdata_', control, task, '_', int2str(session), '_', kernel, '_0.mat'];
        load(file_path, 'N', 'B');
        file_path = ['../GLM_data/', control, task, '/raster_', control, task, '_', int2str(session), '_0.mat'];
        load(file_path, 'n_trial');

        min_duration = min(min_duration, B);
        
        fprintf('Session %d, task %s, N = %d, B = %d, n_trial = %d\n', session, task, N, B, n_trial);
    end
    seconds = floor(min_duration/1000);
    fprintf('Session %d, min_duration = %d, %d:%d\n', session, min_duration, floor(seconds/60), mod(seconds, 60));

    % align data to min_duration
    % modes = {'First', 'Last', 'Random'};
    modes = {'Last'};
    for task_idx = 1:length(tasks)
        task = tasks{task_idx};
        fprintf('Aligning %s %s\n', control, task);
        fprintf('Loading...');
        file_path = ['../GLM_data/', control, task, '/GLMdata_', control, task, '_', int2str(session), '_', kernel, '_0.mat'];
        load(file_path, 'N', 'B', 'PS_kernels', 'conn_kernels', 'kernel_len', 'n_PS_kernel',...
            'n_conn_kernel', 'predjs_PS', 'predjs_conn', 'raster');
        raster_full = raster;
        predjs_PS_full = predjs_PS;
        predjs_conn_full = predjs_conn;
        B_full = B;

        % raster and border file 
        file_path = ['../GLM_data/', control, task, '/raster_', control, task, '_', int2str(session), '_0.mat'];
        load(file_path, 'n_trial', 'cell_id', 'cell_area', 'channel', 'session_name_full');
        file_path = ['../GLM_data/', control, task, '/borders_', control, task, '_', int2str(session), '.mat'];
        load(file_path, 'borders');
        fprintf('Done\n');

        for mode_idx = 1:length(modes)
            mode = modes{mode_idx};
            fprintf('Mode: %s\n', mode);
            switch mode
                case 'First'
                    selected = 1:min_duration;
                    % fprintf('First\n');
                case 'Last'
                    selected = (B_full-min_duration+1):B_full;
                    % fprintf('Last\n');
                case 'Random'
                    selected = randperm(B_full, min_duration);
                    % fprintf('Random\n');
            end
            % if strcmp(mode, 'First')
            %     selected = 1:min_duration;
            %     fprintf('First\n');
            % elseif strcmp(mode, 'Last')
            %     selected = (B-min_duration+1):B;
            %     fprintf('Last\n');
            % elseif strcmp(mode, 'Random')
            %     selected = randperm(B, min_duration);
            %     fprintf('Random\n');
            % end

            raster = raster_full(:, selected);
            predjs_PS = predjs_PS_full(:, selected, :);
            predjs_conn = predjs_conn_full(:, selected, :);
            % if n_PS_kernel > 0
            %     predjs_PS = predjs_PS_full(:, selected, :);
            % else
            %     predjs_PS = [];
            % end
            % if n_conn_kernel > 0
            %     predjs_conn = predjs_conn_full(:, selected, :);
            % else
            %     predjs_conn = [];
            % end
            B = min_duration;

            full_name = [control, task, '_Align', mode];
            fprintf('Saving...');
            file_path = ['../GLM_data/', full_name, '/GLMdata_', full_name, '_', int2str(session), '_', kernel, '_0.mat'];
            check_path(['../GLM_data/', full_name]);
            save(file_path, 'N', 'B', 'PS_kernels', 'conn_kernels', 'kernel_len', 'n_PS_kernel',...
                'n_conn_kernel', 'predjs_PS', 'predjs_conn', 'raster');

            % raster file and border file
            file_path = ['../GLM_data/', full_name, '/raster_', full_name, '_', int2str(session), '_0.mat'];
            save(file_path, 'n_trial', 'cell_id', 'cell_area', 'channel', 'session_name_full');
            file_path = ['../GLM_data/', full_name, '/borders_', full_name, '_', int2str(session), '.mat'];
            save(file_path, 'borders');

            fprintf('Done.\n');
            fprintf('CheckSum: %f\n', sum(raster, 'all'));
        end
    end
end
