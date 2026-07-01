% register_training_task_refactored.m - Register training tasks for GLM analysis.
%
% This version removes the duplicated Slayer/Zeppelin/Emperor blocks.
% Edit the parameter block below to change datasets, states, areas, alignments,
% kernels, resting threshold, training config, and missing-file behavior.

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

%% Parameters

% Metadata source.
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');

% Output.
task_file_folder = fullfile(root, 'Data', 'Working', 'Meta', 'training_tasks');

% Task groups. Indices refer to dataset_names/session_nums in PDS_dataset_info.mat.
% All animals are included in each group. Groups are split by kernel family:
% one parameter-sweep group and one composite-kernel group per family.
all_dataset_indices = 1:9;

long_time_suffixes = {'5', '10', '20', '40', '80', '160', '320', '640', '1000'};
long_step_suffixes = {'0_5', '5_10', '10_20', '20_40', '40_80', ...
                      '80_160', '160_320', '320_640', '640_1000', '1000_3000'};

long_exp_param_kernels = strcat('LongExp', long_time_suffixes);
long_gauss_param_kernels = strcat('LongGaussC', long_time_suffixes);
long_gderiv_param_kernels = strcat('LongGDeriv', long_time_suffixes);
long_step_param_kernels = strcat('LongStepB', long_step_suffixes);

task_groups = struct( ...
    'task_name', { ...
        'LongExpParams', 'LongExpMix', ...
        'LongGaussCParams', 'LongGaussCMix', ...
        'LongGDerivParams', 'LongGDerivMix', ...
        'LongStepBParams', 'LongStepBAdaptive' ...
    }, ...
    'dataset_indices', { ...
        all_dataset_indices, all_dataset_indices, ...
        all_dataset_indices, all_dataset_indices, ...
        all_dataset_indices, all_dataset_indices, ...
        all_dataset_indices, all_dataset_indices ...
    }, ...
    'kernel_types', { ...
        long_exp_param_kernels, {'LongExpMix'}, ...
        long_gauss_param_kernels, {'LongGaussCMix'}, ...
        long_gderiv_param_kernels, {'LongGDerivMix'}, ...
        long_step_param_kernels, {'LongStepBAdaptive'} ...
    } ...
);

% Task dimensions.
task_params = struct();
% task_params.merge_types = {'Cortex'};
task_params.merge_types = {'Full', 'Cortex'};

task_params.prepost_types = {'Pre', 'Post'};
task_params.states = {'RestOpen', 'RestClose'};
% task_params.states = {'RestOpen', 'RestClose', 'Task'};

% task_params.align_types = {'Longest'};
task_params.align_types = {'Longest', 'Last'};

task_params.resting_dur_threshold = 15;

% Skip rules.
skip_rules = struct();
skip_rules.skip_noinj_post = true;        % No-injection sessions do not have Post.
skip_rules.skip_full_post = true;         % Post sessions do not have thalamus/Full.
skip_rules.skip_slayer_noinj_rest = true; % SlayerNoinj does not have RestOpen/RestClose.
skip_rules.skip_missing_raster = false;   % false = error; true = skip missing raster.

% Training configuration.
training_config = struct();
training_config.max_epochs = 3000;
training_config.save_interval = 100;
% Kernel list is assigned per task group above.
% Keep this field for compatibility with build_training_tasks/make_task_struct.
training_config.kernel_types = {};

training_config.crossval_fold_num = 3;
training_config.shuffle_size = 0;

training_config.reg = struct();
training_config.reg.l1 = 0;
training_config.reg.l2 = 0.2;
training_config.reg.name = 'L2=0_2';

%% Load metadata once.
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums', ...
    'cortex_files', 'thalamus_files', 'eyeID_files'); %#ok<ASGLU>

%% Register each task group.
check_path(task_file_folder);

for group_idx = 1:numel(task_groups)
    group = task_groups(group_idx);

    group_training_config = training_config;
    group_training_config.kernel_types = group.kernel_types;

    tasks = build_training_tasks( ...
        root, ...
        dataset_names, ...
        session_nums, ...
        group.dataset_indices, ...
        task_params, ...
        group_training_config, ...
        skip_rules);

    task_file_name = sprintf('training_task_%s.mat', group.task_name);
    task_file_path = fullfile(task_file_folder, task_file_name);
    save(task_file_path, 'tasks');

    fprintf('Saved %d tasks for %s to:\n%s\n', numel(tasks), group.task_name, task_file_path);
end

%% Local functions

function tasks = build_training_tasks(root, dataset_names, session_nums, dataset_indices, task_params, training_config, skip_rules)
    tasks = struct();
    task_idx = 0;

    for dataset_idx = dataset_indices
        dataset_name = dataset_names{dataset_idx};
        session_num = session_nums(dataset_idx);

        for session_idx = 1:session_num
            for align_idx = 1:numel(task_params.align_types)
                align = task_params.align_types{align_idx};

                for kernel_idx = 1:numel(training_config.kernel_types)
                    kernel = training_config.kernel_types{kernel_idx};

                    for merge_idx = 1:numel(task_params.merge_types)
                        merge_type = task_params.merge_types{merge_idx};

                        for prepost_idx = 1:numel(task_params.prepost_types)
                            prepost = task_params.prepost_types{prepost_idx};

                            for state_idx = 1:numel(task_params.states)
                                state = task_params.states{state_idx};

                                if should_skip_task(dataset_name, merge_type, prepost, state, skip_rules)
                                    continue;
                                end

                                task = make_task_struct( ...
                                    dataset_name, ...
                                    session_idx, ...
                                    merge_type, ...
                                    prepost, ...
                                    state, ...
                                    align, ...
                                    kernel, ...
                                    task_params, ...
                                    training_config);

                                raster_file_path = get_raster_file_path(root, task);
                                if ~isfile(raster_file_path)
                                    if skip_rules.skip_missing_raster
                                        fprintf('Raster file does not exist: %s. Skipping this task.\n', raster_file_path);
                                        continue;
                                    else
                                        throw(MException('RegisterTrainingTask:RasterFileNotFound', ...
                                            'Raster file does not exist: %s. Stop registration.', raster_file_path));
                                    end
                                end

                                task_idx = task_idx + 1;
                                if task_idx == 1
                                    tasks = task; % initialize struct array
                                else
                                    tasks(task_idx) = task; %#ok<AGROW>
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    if task_idx == 0
        tasks = struct([]);
    end
end

function tf = should_skip_task(dataset_name, merge_type, prepost, state, skip_rules)
    tf = false;

    if skip_rules.skip_noinj_post && contains(dataset_name, 'Noinj') && strcmp(prepost, 'Post')
        tf = true;
        return;
    end

    if skip_rules.skip_full_post && strcmp(merge_type, 'Full') && strcmp(prepost, 'Post')
        tf = true;
        return;
    end

    if skip_rules.skip_slayer_noinj_rest && strcmp(dataset_name, 'SlayerNoinj') && contains(state, 'Rest')
        tf = true;
        return;
    end
end

function task = make_task_struct(dataset_name, session_idx, merge_type, prepost, state, align, kernel, task_params, training_config)
    [animal_name, injection] = split_dataset_name(dataset_name);

    task = struct();
    task.animal_name = animal_name;
    task.injection = injection;
    task.prepost = prepost;
    task.state = state;
    task.area = merge_type;
    task.align = align;
    task.session_idx = session_idx;
    task.dataset_name = sprintf('%s%s%s%s%s', dataset_name, prepost, state, merge_type, align);
    task.border_name = sprintf('%s%s', dataset_name, merge_type);
    task.session_name = sprintf('%s_%d', task.dataset_name, session_idx);
    task.resting_dur_threshold = task_params.resting_dur_threshold;
    config = training_config;
    config.kernel = kernel;
    if isfield(config, 'kernel_types')
        config = rmfield(config, 'kernel_types');
    end
    task.config = config;
end

function raster_file_path = get_raster_file_path(root, task)
    raster_file_folder = fullfile(root, 'Data', 'Working', 'raster');
    raster_file_name = generate_filename('raster', task);
    raster_file_path = fullfile(raster_file_folder, raster_file_name);
end
