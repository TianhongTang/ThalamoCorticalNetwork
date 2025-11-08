%% GLM_main.m - GLM raster fitting
% Input: raw raster file "../GLM_data/[dataset_name]/raster_[dataset_name]_[session]_0.mat"

% Required variables in the input file:
% "n_trial": integer, trial number.
% "trial_len": int(1, n_trial), time bin number of each trial.
% "rasters": cell(1, n_trial), each element is a trial.
% Each raster is a (N, trial_len(i)) binary matrix, N is number of neurons.

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
force_rebuild = false;
force_retrain = false;
debug = true;

% register tasks
controls = {'Muscimol', 'Saline'};
% controls = {'Muscimol'};
session_idxs_all = {1:10, 1:5}; % session indices for each control
% area_types = {'Full', 'Cortex'};
area_types = {'Cortex'};
prepost_types = {'Pre', 'Post'};
states = {'RestOpen', 'RestClose', 'Task'};
alignments = {'Last'};

kernel = 'DeltaPure';
reg = struct();
reg.l1=0;
reg.l2=0.2;
reg.name='L2=0_2';

tasks = {};
for control_idx = 1:length(controls)
    control = controls{control_idx};
    sessions = session_idxs_all{control_idx};
    for session_idx = sessions
        for area_idx = 1:length(area_types)
            area_type = area_types{area_idx};
            prepost_states = {};
            for prepost_idx = 1:length(prepost_types)
                prepost = prepost_types{prepost_idx};
                if strcmp(area_type, 'Full') && strcmp(prepost, 'Post')
                    % All thalamus data in post sessions are not available
                    continue;
                end
                for align_idx = 1:length(alignments)
                    alignment = alignments{align_idx};
                    for state_idx = 1:length(states)
                        state = states{state_idx};

                        task = struct();
                        task.control = control;
                        task.area_type = area_type;
                        task.prepost = prepost;
                        task.state = state;
                        task.alignment = alignment;
                        task.name = [control, prepost, state, area_type, 'Align', alignment];
                        task.session_idx = session_idx;
                        task.kernel = kernel;
                        task.reg = reg;
                        task.shuffle_size = 0; % number of shuffles
                        task.epoch = 2500;
                        tasks{end+1} = task;
                    end
                end
            end
        end
    end
end


% run tasks
failed_list = {};
skipped = 0;
failed= 0;
success = 0;

total_time = 0;
task_num = length(tasks);
for task_idx=1:task_num
    task = tasks{task_idx};
    % try
        tick_session = tic;
        fprintf("Task %d/%d: %s, session%d\n", task_idx, task_num, task.name, task.session_idx);
        skip_flag = true;
        dataset_name = task.name;
        session_idx = task.session_idx;
        kernel_name = task.kernel;
        reg = task.reg;
        shuffle_size=task.shuffle_size;
        max_epoch=task.epoch;
        
        %% generate shuffled raster
        fprintf("Shuffle rasters\n");
        tic;
        for shuffle_id=0:shuffle_size
            % skip if already exists
            target_folder = fullfile(root, 'Data', 'Working', 'raster');
            target_file = sprintf('shuffled_%s_%d_%d.mat', dataset_name, session_idx, shuffle_id);
            target_path = fullfile(target_folder, target_file);
            if isfile(target_path) && ~force_rebuild
                fprintf("Skip %d. \n", shuffle_id);
                continue;
            end

            skip_flag = false;
            if shuffle_id==0
                % original data
                shuffle_type = "None";
            else
                shuffle_type = "Across trial";
            end
            shuffle(dataset_name, session_idx, shuffle_id, shuffle_id, shuffle_type);
        end
        toc;

        %% split cross validation folds
        fprintf("Cross-validation split\n");
        tic;
        for shuffle_id=0:shuffle_size % seed=0: original data (no shuffle)
            % skip if already exists
            target_folder = fullfile(root, 'Data', 'Working', 'crossval_split');
            target_file = sprintf('crossval_%s_%d_%d.mat', dataset_name, session_idx, shuffle_id);
            target_path = fullfile(target_folder, target_file);
            if isfile(target_path) && ~force_rebuild
                fprintf("Skip %d. \n", shuffle_id);
                continue;
            end

            skip_flag = false;
            if strcmp(task.state, 'Task')
                split_type = 'Trial';
            else
                split_type = 'Time';
            end
            crossval_split(dataset_name, session_idx, shuffle_id, 3, split_type);
        end

        %% convolve predj and combine trials
        fprintf("Convolution\n");
        tic;
        for shuffle_id=0:shuffle_size % seed=0: original data (no shuffle)
            % skip if already exists
            target_folder = fullfile(root, 'Data', 'Working', 'GLM_data');
            target_file = sprintf('GLMdata_%s_%d_%d_%s.mat', dataset_name, session_idx, shuffle_id, kernel_name);
            target_path = fullfile(target_folder, target_file);
            if isfile(target_path) && ~force_rebuild
                fprintf("Skip %d. \n", shuffle_id);
                continue;
            end
            skip_flag = false;
            convolution(dataset_name, session_idx, shuffle_id, kernel_name);
        end
        toc;

        %% GLM inference
        for shuffle_seed=0:shuffle_size
            fprintf("Training %d\n", shuffle_seed);
            skip_flag = false;
            tic;
            for fold_idx = 0:0
                % skip if already exists
                target_folder = fullfile(root, 'Data', 'Working', 'GLM_models');
                target_file = sprintf('GLM_%s_s%d_shuffle%d_%s_%s_epoch%d_fold%d.mat', dataset_name, session_idx, shuffle_id, kernel_name, ...
                    reg.name, max_epoch, fold_idx);
                target_path = fullfile(target_folder, target_file);
                if isfile(target_path) && ~force_retrain
                    fprintf("Skip. \n");
                    continue;
                end
                GLM_multi_kernel_crossval(dataset_name, session_idx, kernel_name, shuffle_id, max_epoch, reg, 1, 5e-3, fold_idx);
            end
            toc;
        end

        %% plot
        % fprintf("Plotting\n");
        % tic;
        % channel_file = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', int2str(session_idx), ...
        % '_0.mat'];
        % load(channel_file, "channel");
        % plot_GLM_sorted(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size, "idx", channel);
        % toc;

        if skip_flag
            skipped = skipped + 1;
        else
            success = success + 1;
        end

    % catch ME
    %     fprintf("Failed: %s\n", ME.message);
    %     failed = failed + 1;
    %     failed_list{end+1} = {dataset_name, int2str(session_idx), ME.message};
    %     if debug
    %         throw(ME);
    %     end
    % end
end

fprintf("Total: %d, Success: %d, Skipped: %d, Failed: %d\n", total_training, success, skipped, failed);
% save failed_list
save('../GLM_data/failed_list.mat', 'failed_list');
if failed>0
    for i=1:length(failed_list)
        fprintf("Failed: %s, %s, %s\n", failed_list{i}{1}, failed_list{i}{2}, failed_list{i}{3});
    end
end