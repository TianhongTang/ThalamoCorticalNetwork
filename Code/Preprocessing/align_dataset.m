%% Align raster files of all states to the same duration

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

% load data
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums', 'cortex_files', 'thalamus_files', 'eyeID_files');

% Load data
area_types = {'Full', 'Cortex'};
prepost_types = {'Pre', 'Post'};
states = {'RestOpen', 'RestClose'};
kernel = 'DeltaPure';

% load kernel length
folder_name = fullfile(root, 'Data', 'Working', 'kernel');
file_name = sprintf('kernel_%s.mat', kernel);
file_path = fullfile(folder_name, file_name);
load(file_path, 'kernel_len');
align_kernel_name = kernel;
align_kernel_len = kernel_len;

tasks = {};
% register tasks
for dataset_idx = 1:6
    dataset_name = dataset_names{dataset_idx};
    session_num = session_nums(dataset_idx);
    for session_idx = 1:session_num
        for area_idx = 1:length(area_types)
            area_type = area_types{area_idx};
            task = struct();
            task.dataset_name = dataset_name;
            task.area_type = area_type;
            task.session_idx = session_idx;

            prepost_states = {};
            for prepost_idx = 1:length(prepost_types)
                prepost = prepost_types{prepost_idx};
                if strcmp(area_type, 'Full') && strcmp(prepost, 'Post')
                    % All thalamus data in post sessions are not available
                    continue;
                end
                if contains(dataset_name, 'Noinj') && strcmp(prepost, 'Post')
                    % No injection session does not have post sessions 
                    continue;
                end

                for state_idx = 1:length(states)
                    state = states{state_idx};
                    if strcmp(dataset_name, 'SlayerNoinj') && contains(state, 'Rest')
                        % SlayerNoinj does not have RestOpen and RestClose sessions
                        continue;
                    end
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
    dataset_name = task.dataset_name;
    session_idx = task.session_idx;
    area_type = task.area_type;
    prepost_states = task.prepost_states;
    fprintf('----------------------------------------\n');
    fprintf('Task %d/%d: %s Session %d Area %s with %d prepost_states\n', ...
        task_idx, task_num, dataset_name, session_idx, area_type, length(prepost_states));
    
    % find min duration among prepost_states
    min_duration = 10000000;
    for ps_idx = 1:length(prepost_states)
        prepost = prepost_states{ps_idx}.prepost;
        state = prepost_states{ps_idx}.state;

        % load length info
        folder_name = fullfile(root, 'Data', 'Working', 'raster');
        file_name = sprintf('raster_%s_%d.mat', [dataset_name, prepost, state, area_type], session_idx);
        raster_path = fullfile(folder_name, file_name);
        load(raster_path, 'trial_num', 'trial_len', 'N');

        B = sum(trial_len - align_kernel_len + 1); % effective length

        min_duration = min(min_duration, B);
        fprintf('%s%s%s%s Session %d, N = %d, B = %d, trial_num = %d\n', dataset_name, prepost, state, area_type, session_idx, N, B, trial_num);
    end
    seconds = floor(min_duration/1000);
    fprintf('%s%s Session %d, min_duration = %d (%d:%d)\n', dataset_name, area_type, session_idx, min_duration, floor(seconds/60), mod(seconds, 60));

    % align data to min_duration
    modes = {'First', 'Last'};
    for ps_idx = 1:length(prepost_states)
        prepost = prepost_states{ps_idx}.prepost;
        state = prepost_states{ps_idx}.state;
        fprintf('Aligning %s ...\n', [prepost, state]);

        for mode_idx = 1:length(modes)
            % raster file and border file
            fprintf('Loading...');
            folder_name = fullfile(root, 'Data', 'Working', 'raster');
            file_name = sprintf('raster_%s_%d.mat', [dataset_name, prepost, state, area_type], session_idx);
            raster_path = fullfile(folder_name, file_name);
            load(raster_path, 'trial_num', 'cell_id', 'cell_area', 'channel', 'session_name', 'trial_len', 'rasters');
            % folder_name = fullfile(root, 'Data', 'Working', 'border');
            % file_name = sprintf('borders_%s_%d.mat', [dataset_name, prepost, area_type], session_idx);
            % raster_path = fullfile(folder_name, file_name);
            % load(raster_path, 'borders');
            fprintf('Done\n');

            mode = modes{mode_idx};
            fprintf('Mode: %s\n', mode);

            raster_aligned = cell(1, 0);
            trial_len_aligned = zeros(1, 0);
            current_eff_dur = 0;

            switch mode
                case 'First'
                    for trial = 1:trial_num
                        trial_eff_dur = trial_len(trial) - kernel_len + 1;
                        if current_eff_dur + trial_eff_dur >= min_duration
                            cut_len = min_duration - current_eff_dur + kernel_len - 1;
                            raster_aligned{end+1} = rasters{trial}(:, 1:cut_len);
                            trial_len_aligned(end+1) = cut_len;
                            break;
                        else
                            raster_aligned{end+1} = rasters{trial};
                            current_eff_dur = current_eff_dur + trial_eff_dur;
                            trial_len_aligned(end+1) = trial_eff_dur + kernel_len - 1;
                        end
                    end
                case 'Last'
                    for trial = trial_num:-1:1
                        trial_eff_dur = trial_len(trial) - kernel_len + 1;
                        if current_eff_dur + trial_eff_dur >= min_duration
                            cut_len = min_duration - current_eff_dur + kernel_len - 1;
                            raster_aligned{end+1} = rasters{trial}(:, end-cut_len+1:end);
                            trial_len_aligned(end+1) = cut_len;
                            break;
                        else
                            raster_aligned{end+1} = rasters{trial};
                            current_eff_dur = current_eff_dur + trial_eff_dur;
                            trial_len_aligned(end+1) = trial_eff_dur + kernel_len - 1;
                        end
                    end
            end

            if strcmp(mode, 'Last')
                raster_aligned = fliplr(raster_aligned);
                trial_len_aligned = fliplr(trial_len_aligned);
            end

            rasters = raster_aligned;
            trial_len = trial_len_aligned;
            trial_num = length(rasters);

            firing_rates = cell(1, trial_num);
            for i = 1:trial_num
                firing_rates{i} = mean(rasters{i}, 2);
            end

            full_name = [dataset_name, prepost, state, area_type, 'Align', mode];
            fprintf('Saving...');
            save_folder = fullfile(root, 'Data', 'Working', 'raster');
            check_path(save_folder);
            file_name = sprintf('raster_%s_%d.mat', full_name, session_idx);
            file_path = fullfile(save_folder, file_name);
            save(file_path, 'N', 'trial_num', 'cell_id', 'cell_area', 'channel', 'session_name', 'trial_len', 'rasters', 'firing_rates');

            fprintf('Done.\n');
        end
    end
end

