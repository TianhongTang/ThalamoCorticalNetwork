% register_training_task_kernel_groups_refactored.m - Register GLM training tasks by kernel group.
%
% Task grouping:
%   - All animals/datasets are included in every task group.
%   - Each kernel family creates two task groups:
%       1) parameter kernels, e.g. LongGaussParamsC -> LongGaussC5 ... LongGaussC1000
%       2) composite kernel,  e.g. LongGaussMix    -> LongGaussMix
%
% The family prefix does not include B/C.  B/C are represented only as the
% parameter-group infix for LongGauss and LongStep.

clear;

%% Get root folder
code_depth = 4;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end

addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Parameters
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');

task_file_folder = fullfile(root, 'Data', 'Working', 'Meta', 'training_tasks');

% Dataset groups. Indices refer to dataset_names/session_nums in PDS_dataset_info.mat.
% 1:3 Slayer, 4:6 Zeppelin, 7:9 Emperor in the current metadata convention.
all_dataset_indices = 1:9;

% Kernel-group definitions.
% family_prefix: base kernel family name, without B/C.
% param_infix:   extra infix used only for the parameter-kernel group.
% composite_tag: suffix used by the composite kernel name.
long_time_suffixes = {'5', '10', '20', '40', '80', '160', '320', '640', '1000'};
long_step_suffixes = {'0_5', '5_10', '10_20', '20_40', '40_80', ...
                      '80_160', '160_320', '320_640', '640_1000', '1000_3000'};

kernel_families = [ ...
    make_kernel_family('LongExp',    '',  long_time_suffixes, 'Mix'), ...
    make_kernel_family('LongGauss',  'C', long_time_suffixes, 'Mix'), ...
    make_kernel_family('LongGDeriv', '',  long_time_suffixes, 'Mix'), ...
    make_kernel_family('LongStep',   'B', long_step_suffixes, 'Adaptive') ...
];

task_groups = build_kernel_task_groups(all_dataset_indices, kernel_families);

% Task dimensions.
task_params = struct();
task_params.merge_types = {'Cortex'};
task_params.prepost_types = {'Pre', 'Post'};
task_params.states = {'RestOpen', 'RestClose'};
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
training_config.kernel_types = {}; % Assigned per task group below.
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

function spec = make_kernel_family(family_prefix, param_infix, param_suffixes, composite_tag)
    spec = struct();
    spec.family_prefix = family_prefix;
    spec.param_infix = param_infix;
    spec.param_suffixes = param_suffixes;
    spec.composite_tag = composite_tag;
end

function task_groups = build_kernel_task_groups(dataset_indices, kernel_families)
    task_groups = struct('task_name', {}, 'dataset_indices', {}, 'kernel_types', {});

    for spec_idx = 1:numel(kernel_families)
        spec = kernel_families(spec_idx);

        param_task_name = make_param_task_name(spec.family_prefix, spec.param_infix);
        composite_task_name = make_composite_task_name(spec.family_prefix, spec.composite_tag);

        param_kernel_prefix = make_param_kernel_prefix(spec.family_prefix, spec.param_infix);
        param_kernel_types = strcat(param_kernel_prefix, spec.param_suffixes);

        base_idx = numel(task_groups);

        task_groups(base_idx + 1).task_name = param_task_name;
        task_groups(base_idx + 1).dataset_indices = dataset_indices;
        task_groups(base_idx + 1).kernel_types = param_kernel_types;

        task_groups(base_idx + 2).task_name = composite_task_name;
        task_groups(base_idx + 2).dataset_indices = dataset_indices;
        task_groups(base_idx + 2).kernel_types = {composite_task_name};
    end
end

function name = make_param_task_name(family_prefix, ~)
    % Examples:
    %   LongExp + ''  -> LongExpParams
    %   LongGauss + C -> LongGaussParamsC
    %   LongStep + B  -> LongStepParamsB
    name = sprintf('%sParams', family_prefix);
end

function name = make_composite_task_name(family_prefix, composite_tag)
    % Examples:
    %   LongGauss + Mix      -> LongGaussMix
    %   LongStep  + Adaptive -> LongStepAdaptive
    name = sprintf('%s%s', family_prefix, composite_tag);
end

function prefix = make_param_kernel_prefix(family_prefix, param_infix)
    % Examples:
    %   LongGauss + C -> LongGaussC
    %   LongStep  + B -> LongStepB
    prefix = sprintf('%s%s', family_prefix, param_infix);
end

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
                                    tasks = task;
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
