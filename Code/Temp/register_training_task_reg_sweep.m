% register_training_task_reg_sweep.m
% Register GLM training tasks for a regularization-strength sweep.
%
% Fixed kernels:
%   LongGaussC5, LongGaussC80, LongGaussC320, LongGaussC1000
%
% Regularization groups:
%   L1:  l1 = value, l2 = 0
%   L2:  l1 = 0,     l2 = value
%   L12: l1 = value, l2 = value
%
% Each regularization group is saved as one task file containing all
% selected regularization values and all four kernels.

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

%% Paths
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');
task_file_folder = fullfile(metadata_folder, 'training_tasks');

%% Dataset selection
% Indices refer to dataset_names/session_nums in PDS_dataset_info.mat.
% 1:3 Slayer, 4:6 Zeppelin, 7:9 Emperor in the current convention.
all_dataset_indices = 1:9;

%% Fixed kernels and regularization sweep
fixed_kernel_types = { ...
    'LongGaussC5', ...
    'LongGaussC80', ...
    'LongGaussC320', ...
    'LongGaussC1000' ...
};

reg_values = [0, 0.01, 0.05, 0.1, 0.2, 0.4, ...
              0.6, 0.8, 1, 2, 5, 10];

reg_types = {'L1', 'L2', 'L12'};

task_groups = build_regularization_task_groups( ...
    all_dataset_indices, fixed_kernel_types, reg_values, reg_types);

%% Task dimensions
task_params = struct();
task_params.merge_types = {'Cortex'};
task_params.prepost_types = {'Pre', 'Post'};
task_params.states = {'RestOpen', 'RestClose'};
task_params.align_types = {'Longest', 'Last'};
task_params.resting_dur_threshold = 15;

%% Skip rules
skip_rules = struct();
skip_rules.skip_noinj_post = true;        % No-injection sessions do not have Post.
skip_rules.skip_full_post = true;         % Post sessions do not have thalamus/Full.
skip_rules.skip_slayer_noinj_rest = true; % SlayerNoinj does not have RestOpen/RestClose.
skip_rules.skip_missing_raster = false;   % false = error; true = skip missing raster.

%% Shared training configuration
training_config = struct();
training_config.max_epochs = 3000;
training_config.save_interval = 100;
training_config.kernel_types = {}; % Assigned per task group.
training_config.reg_configs = struct([]); % Assigned per task group.
training_config.crossval_fold_num = 3;
training_config.shuffle_size = 0;

%% Load metadata once
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums', ...
    'cortex_files', 'thalamus_files', 'eyeID_files'); %#ok<ASGLU>

%% Register each regularization task group
check_path(task_file_folder);

fprintf('Fixed kernels (%d): %s\n', ...
    numel(fixed_kernel_types), strjoin(fixed_kernel_types, ', '));
fprintf('Regularization values (%d): %s\n', ...
    numel(reg_values), strjoin(arrayfun(@num2str, reg_values, 'UniformOutput', false), ', '));
fprintf('Task groups: %s\n\n', strjoin(reg_types, ', '));

for group_idx = 1:numel(task_groups)
    group = task_groups(group_idx);

    fprintf('===== Registering group %d/%d: %s =====\n', ...
        group_idx, numel(task_groups), group.task_name);

    group_training_config = training_config;
    group_training_config.kernel_types = group.kernel_types;
    group_training_config.reg_configs = group.reg_configs;

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

    fprintf('Saved %d tasks for %s to:\n%s\n\n', ...
        numel(tasks), group.task_name, task_file_path);
end

fprintf('Finished registering %d regularization task groups.\n', numel(task_groups));

%% Local functions

function task_groups = build_regularization_task_groups( ...
    dataset_indices, kernel_types, reg_values, reg_types)

    task_groups = struct( ...
        'task_name', {}, ...
        'dataset_indices', {}, ...
        'kernel_types', {}, ...
        'reg_configs', {});

    for type_idx = 1:numel(reg_types)
        reg_type = upper(char(reg_types{type_idx}));

        task_groups(type_idx).task_name = sprintf('LongGaussReg%s', reg_type);
        task_groups(type_idx).dataset_indices = dataset_indices;
        task_groups(type_idx).kernel_types = kernel_types;
        task_groups(type_idx).reg_configs = ...
            build_regularization_configs(reg_type, reg_values);
    end
end

function reg_configs = build_regularization_configs(reg_type, reg_values)
    reg_configs = repmat(struct('l1', 0, 'l2', 0, 'name', ''), ...
        1, numel(reg_values));

    for value_idx = 1:numel(reg_values)
        value = reg_values(value_idx);

        if ~isscalar(value) || ~isfinite(value) || value < 0
            error('Regularization values must be finite non-negative scalars.');
        end

        switch upper(reg_type)
            case 'L1'
                l1 = value;
                l2 = 0;
            case 'L2'
                l1 = 0;
                l2 = value;
            case 'L12'
                l1 = value;
                l2 = value;
            otherwise
                error('Unknown regularization type: %s', reg_type);
        end

        reg_configs(value_idx).l1 = l1;
        reg_configs(value_idx).l2 = l2;
        reg_configs(value_idx).name = sprintf( ...
            '%s=%s', upper(reg_type), format_reg_value(value));
    end
end

function text = format_reg_value(value)
    % Follow the existing metadata naming convention, e.g. 0.2 -> 0_2.
    text = sprintf('%.15g', value);
    text = strrep(text, '.', '_');
    text = strrep(text, '-', 'm');
    text = strrep(text, '+', '');
end

function tasks = build_training_tasks( ...
    root, dataset_names, session_nums, dataset_indices, ...
    task_params, training_config, skip_rules)

    tasks = struct([]);
    task_idx = 0;

    if isempty(training_config.kernel_types)
        error('training_config.kernel_types is empty.');
    end
    if isempty(training_config.reg_configs)
        error('training_config.reg_configs is empty.');
    end

    for dataset_idx = dataset_indices
        dataset_name = dataset_names{dataset_idx};
        session_num = session_nums(dataset_idx);

        for session_idx = 1:session_num
            for align_idx = 1:numel(task_params.align_types)
                align = task_params.align_types{align_idx};

                % Keep each regularization value as a contiguous task block.
                for reg_idx = 1:numel(training_config.reg_configs)
                    reg_config = training_config.reg_configs(reg_idx);

                    for kernel_idx = 1:numel(training_config.kernel_types)
                        kernel = training_config.kernel_types{kernel_idx};

                        for merge_idx = 1:numel(task_params.merge_types)
                            merge_type = task_params.merge_types{merge_idx};

                            for prepost_idx = 1:numel(task_params.prepost_types)
                                prepost = task_params.prepost_types{prepost_idx};

                                for state_idx = 1:numel(task_params.states)
                                    state = task_params.states{state_idx};

                                    if should_skip_task( ...
                                            dataset_name, merge_type, ...
                                            prepost, state, skip_rules)
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
                                        reg_config, ...
                                        task_params, ...
                                        training_config);

                                    raster_file_path = get_raster_file_path(root, task);
                                    if ~isfile(raster_file_path)
                                        if skip_rules.skip_missing_raster
                                            fprintf(['Raster file does not exist: %s. ' ...
                                                'Skipping this task.\n'], raster_file_path);
                                            continue;
                                        else
                                            throw(MException( ...
                                                'RegisterTrainingTask:RasterFileNotFound', ...
                                                ['Raster file does not exist: %s. ' ...
                                                 'Stop registration.'], ...
                                                raster_file_path));
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
    end
end

function tf = should_skip_task( ...
    dataset_name, merge_type, prepost, state, skip_rules)

    tf = false;

    if skip_rules.skip_noinj_post && ...
            contains(dataset_name, 'Noinj') && strcmp(prepost, 'Post')
        tf = true;
        return;
    end

    if skip_rules.skip_full_post && ...
            strcmp(merge_type, 'Full') && strcmp(prepost, 'Post')
        tf = true;
        return;
    end

    if skip_rules.skip_slayer_noinj_rest && ...
            strcmp(dataset_name, 'SlayerNoinj') && contains(state, 'Rest')
        tf = true;
    end
end

function task = make_task_struct( ...
    dataset_name, session_idx, merge_type, prepost, state, align, ...
    kernel, reg_config, task_params, training_config)

    [animal_name, injection] = split_dataset_name(dataset_name);

    task = struct();
    task.animal_name = animal_name;
    task.injection = injection;
    task.prepost = prepost;
    task.state = state;
    task.area = merge_type;
    task.align = align;
    task.session_idx = session_idx;
    task.dataset_name = sprintf('%s%s%s%s%s', ...
        dataset_name, prepost, state, merge_type, align);
    task.border_name = sprintf('%s%s', dataset_name, merge_type);
    task.session_name = sprintf('%s_%d', task.dataset_name, session_idx);
    task.resting_dur_threshold = task_params.resting_dur_threshold;

    config = training_config;
    config.kernel = kernel;
    config.reg = reg_config;

    removable_fields = {'kernel_types', 'reg_configs'};
    for field_idx = 1:numel(removable_fields)
        field_name = removable_fields{field_idx};
        if isfield(config, field_name)
            config = rmfield(config, field_name);
        end
    end

    task.config = config;
end

function raster_file_path = get_raster_file_path(root, task)
    raster_file_folder = fullfile(root, 'Data', 'Working', 'raster');
    raster_file_name = generate_filename('raster', task);
    raster_file_path = fullfile(raster_file_folder, raster_file_name);
end
