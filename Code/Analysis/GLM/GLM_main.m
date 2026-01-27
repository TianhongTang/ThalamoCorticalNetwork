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
gpuDeviceTable

force_rebuild = false;
force_retrain = false;
debug = true;

training_tasks = {'Slayer'};

kernel = 'DeltaPure';
reg = struct();
reg.l1=0;
reg.l2=0.2;
reg.name='L2=0_2';

tasks = {};

for training_idx = 1:length(training_tasks)
    training_task = training_tasks{training_idx};

    % load training task
    training_task_folder = fullfile(root, 'Data', 'Working', 'Meta', 'training_tasks');
    training_task_name = sprintf('training_task_%s.mat', training_task);
    training_task_path = fullfile(training_task_folder, training_task_name);
    if ~isfile(training_task_path)
        fprintf('Training task file not found: %s\n', training_task_path);
        continue;
    end
    load(training_task_path, 'tasks');
    task_num = length(tasks);
    fprintf('Loaded training task: %s, with %d sessions.\n', training_task, task_num);

    % run tasks
    failed_list = {};
    skipped = 0;
    failed = 0;
    success = 0;

    total_time = 0;

    % run each session
    for task_idx = 1:task_num
        task = tasks{task_idx};
        try
            tick_session = tic;
            fprintf("Task %d/%d: %s, session%d\n", task_idx, task_num, task.dataset_name, task.session_idx);
            skip_flag = true;
            dataset_name = task.dataset_name;
            border_name = task.border_name;
            session_idx = task.session_idx;

            config = task.config;
            kernel_name = config.kernel;
            reg = config.reg;
            shuffle_size=config.shuffle_size;
            max_epoch=config.max_epochs;
            fold_num = config.crossval_fold_num;
            
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
                crossval_split(dataset_name, session_idx, shuffle_id, fold_num, split_type);
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
            fprintf("Plotting\n");
            tic;
            channel_folder = fullfile(root, 'Data', 'Working', 'raster');
            channel_file = sprintf('raster_%s_%d.mat', dataset_name, session_idx);
            channel_path = fullfile(channel_folder, channel_file);
            load(channel_path, "channel");
            plot_GLM_sorted(dataset_name, border_name, session_idx, kernel_name, max_epoch, reg, shuffle_size, "idx", channel);
            toc;

            if skip_flag
                skipped = skipped + 1;
            else
                success = success + 1;
            end

        catch ME
            fprintf("Failed: %s\n", ME.message);
            failed = failed + 1;
            % failed_list{end+1} = {dataset_name, int2str(session_idx), ME.message};
            if debug
                rethrow(ME);
            end
        end
    end
    
    % fprintf("Total: %d, Success: %d, Skipped: %d, Failed: %d\n", total_training, success, skipped, failed);
    % % save failed_list
    % save('../GLM_data/failed_list.mat', 'failed_list');
    % if failed>0
    %     for i=1:length(failed_list)
    %         fprintf("Failed: %s, %s, %s\n", failed_list{i}{1}, failed_list{i}{2}, failed_list{i}{3});
    %     end
    % end
end
