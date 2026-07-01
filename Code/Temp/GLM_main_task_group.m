%% GLM_main_task_group.m - GLM raster fitting for one or more registered task groups.
% Select task groups using environment variables:
%   GLM_TRAINING_TASK=LongExpParams
%   GLM_TRAINING_TASKS=LongExpParams,LongExpMix
%   GLM_TASK_GROUP=LongExpParams
% If none is set, SLURM_ARRAY_TASK_ID maps 1..8 to default kernel groups.
% If neither env var nor array id is set, all default kernel groups are run sequentially.

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
try
    gpuDeviceTable
    gpuDevice(1); % Use the first GPU visible to this Slurm task.
catch ME
    warning('GPU setup failed or no GPU is visible: %s', ME.message);
end

% Defaults are conservative for parallel array jobs.
% Override from the submit script with GLM_FORCE_REBUILD / GLM_FORCE_RETRAIN / GLM_FORCE_REPLOT / GLM_DEBUG.
force_rebuild = get_env_bool('GLM_FORCE_REBUILD', false);
force_retrain = get_env_bool('GLM_FORCE_RETRAIN', false);
force_replot  = get_env_bool('GLM_FORCE_REPLOT',  false);
debug         = get_env_bool('GLM_DEBUG',         false);

shuffle_seed_shift = 34;

available_training_tasks = { ...
    'LongExpParams', ...
    'LongExpMix', ...
    'LongGaussParams', ...
    'LongGaussMix', ...
    'LongGDerivParams', ...
    'LongGDerivMix', ...
    'LongStepParams', ...
    'LongStepAdaptive' ...
};
training_tasks = resolve_training_tasks(available_training_tasks);
fprintf('Running training task group(s): %s\n', strjoin(training_tasks, ', '));
fprintf('force_rebuild=%d, force_retrain=%d, force_replot=%d, debug=%d\n', ...
    force_rebuild, force_retrain, force_replot, debug);

failed_list = {};
skipped = 0;
failed = 0;
success = 0;

total_task_num = 0;

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
    fprintf('Loaded training task: %s, with %d tasks.\n', training_task, task_num);
    total_task_num = total_task_num + task_num;

    % run each registered task
    for task_idx = 1:task_num
        task = tasks(task_idx);
        dataset_name = 'UNKNOWN';
        session_idx = NaN;
        try
            fprintf('Task group %s, task %d/%d\n', training_task, task_idx, task_num);
            skip_flag = true;
            dataset_name = task.dataset_name;
            border_name  = task.border_name; %#ok<NASGU>
            session_idx  = task.session_idx;

            config       = task.config;
            kernel_name  = config.kernel;
            reg          = config.reg;
            shuffle_size = config.shuffle_size;
            max_epoch    = config.max_epochs;
            fold_num     = config.crossval_fold_num;

            task.kernel_name = kernel_name;
            task.reg_name = reg.name;

            fprintf('  Dataset: %s, session %d, kernel %s\n', dataset_name, session_idx, kernel_name);

            %% generate shuffled raster
            fprintf('  Shuffle rasters\n');
            tic;
            for shuffle_id = 0:shuffle_size
                shuffle_meta = task;
                shuffle_meta.shuffle_idx = shuffle_id;
                target_folder = fullfile(root, 'Data', 'Working', 'shuffled');
                target_file = generate_filename('shuffled', shuffle_meta);
                target_path = fullfile(target_folder, target_file);
                if isfile(target_path) && ~force_rebuild
                    fprintf('    Skip shuffle %d.\n', shuffle_id);
                    continue;
                end

                skip_flag = false;
                if shuffle_id == 0
                    shuffle_type = "None";
                else
                    shuffle_type = "Across trial";
                end
                shuffle(shuffle_meta, shuffle_id, shuffle_id + shuffle_seed_shift, shuffle_type);
            end
            toc;

            %% split cross validation folds
            fprintf('  Cross-validation split\n');
            tic;
            for shuffle_id = 0:shuffle_size
                crossval_meta = task;
                crossval_meta.shuffle_idx = shuffle_id;

                target_folder = fullfile(root, 'Data', 'Working', 'crossval');
                target_file = generate_filename('crossval', crossval_meta);
                target_path = fullfile(target_folder, target_file);
                if isfile(target_path) && ~force_rebuild
                    fprintf('    Skip crossval %d.\n', shuffle_id);
                    continue;
                end

                skip_flag = false;
                if strcmp(task.state, 'Task')
                    split_type = 'Trial';
                else
                    split_type = 'Time';
                end
                crossval_split(crossval_meta, fold_num, split_type);
            end
            toc;

            %% convolve predictors and combine trials
            fprintf('  Convolution\n');
            tic;
            for shuffle_id = 0:shuffle_size
                conv_meta = task;
                conv_meta.shuffle_idx = shuffle_id;

                target_folder = fullfile(root, 'Data', 'Working', 'GLMdata');
                target_file = generate_filename('GLMdata', conv_meta);
                target_path = fullfile(target_folder, target_file);
                if isfile(target_path) && ~force_rebuild
                    fprintf('    Skip convolution %d.\n', shuffle_id);
                    continue;
                end
                skip_flag = false;
                convolution(conv_meta, kernel_name);
            end
            toc;

            %% GLM inference
            for shuffle_id = 0:shuffle_size
                fprintf('  Training shuffle %d\n', shuffle_id);
                tic;
                for fold_idx = 0:0
                    model_meta = task;
                    model_meta.shuffle_idx = shuffle_id;
                    model_meta.fold_idx = fold_idx;
                    model_meta.fold_num = fold_num;
                    model_meta.epoch = max_epoch;

                    target_folder = fullfile(root, 'Data', 'Working', 'GLM');
                    target_file = generate_filename('GLM', model_meta);
                    target_path = fullfile(target_folder, target_file);
                    if isfile(target_path) && ~force_retrain
                        fprintf('    Skip training fold %d.\n', fold_idx);
                        continue;
                    end
                    skip_flag = false;
                    GLM_multi_kernel_crossval(model_meta, max_epoch, reg, 1, 5e-3, fold_idx);
                end
                toc;
            end

            % Optional plotting block can be re-enabled here if needed.
            if force_replot
                warning('force_replot is set, but plotting block is disabled in this script.');
            end

            if skip_flag
                skipped = skipped + 1;
            else
                success = success + 1;
            end

        catch ME
            fprintf('Failed: %s\n', ME.message);
            failed = failed + 1;
            failed_list{end + 1} = {training_task, dataset_name, int2str(session_idx), ME.message}; %#ok<SAGROW>
            if debug
                rethrow(ME);
            end
        end
    end
end

fprintf('Total: %d, Success: %d, Skipped: %d, Failed: %d\n', total_task_num, success, skipped, failed);

% save failed_list, one file per job/task selection to avoid parallel overwrite.
folder = fullfile(root, 'Data', 'Working', 'log');
check_path(folder);
run_label = sanitize_filename(strjoin(training_tasks, '_'));
job_id = getenv('SLURM_JOB_ID');
array_id = getenv('SLURM_ARRAY_TASK_ID');
if ~isempty(job_id) && ~isempty(array_id)
    log_name = sprintf('failed_list_%s_job%s_array%s.mat', run_label, job_id, array_id);
else
    log_name = sprintf('failed_list_%s.mat', run_label);
end
save(fullfile(folder, log_name), 'failed_list');
for i = 1:numel(failed_list)
    fprintf('Failed: %s, %s, %s, %s\n', failed_list{i}{1}, failed_list{i}{2}, failed_list{i}{3}, failed_list{i}{4});
end

%% Local functions
function value = get_env_bool(name, default_value)
    raw = strtrim(getenv(name));
    if isempty(raw)
        value = default_value;
        return;
    end
    raw = lower(raw);
    if any(strcmp(raw, {'1', 'true', 't', 'yes', 'y', 'on'}))
        value = true;
    elseif any(strcmp(raw, {'0', 'false', 'f', 'no', 'n', 'off'}))
        value = false;
    else
        error('Invalid boolean environment variable %s=%s', name, raw);
    end
end

function training_tasks = resolve_training_tasks(available_training_tasks)
    env_value = first_nonempty_env({'GLM_TASK_GROUP', 'GLM_TRAINING_TASK', 'GLM_TRAINING_TASKS'});
    if ~isempty(env_value)
        env_value = strrep(env_value, ';', ',');
        parts = strsplit(env_value, ',');
        parts = strtrim(parts);
        training_tasks = parts(~cellfun(@isempty, parts));
        validate_training_tasks(training_tasks, available_training_tasks);
        return;
    end

    array_id = strtrim(getenv('SLURM_ARRAY_TASK_ID'));
    if ~isempty(array_id)
        idx = str2double(array_id);
        if ~isfinite(idx) || idx < 1 || idx > numel(available_training_tasks) || idx ~= round(idx)
            error('SLURM_ARRAY_TASK_ID=%s is invalid. Expected an integer from 1 to %d.', ...
                  array_id, numel(available_training_tasks));
        end
        training_tasks = available_training_tasks(idx);
        return;
    end

    training_tasks = available_training_tasks;
end

function value = first_nonempty_env(names)
    value = '';
    for i = 1:numel(names)
        candidate = strtrim(getenv(names{i}));
        if ~isempty(candidate)
            value = candidate;
            return;
        end
    end
end

function validate_training_tasks(training_tasks, available_training_tasks)
    for i = 1:numel(training_tasks)
        if ~ismember(training_tasks{i}, available_training_tasks)
            error('Unknown training task group: %s. Available groups: %s', ...
                training_tasks{i}, strjoin(available_training_tasks, ', '));
        end
    end
end

function out = sanitize_filename(in)
    out = regexprep(char(in), '[^A-Za-z0-9_=-]', '_');
end
