%% Data Summary
% Summarize current available data in Data/Working

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
data_folder = fullfile(root, 'Data', 'Working');
subfolders = dir(data_folder);

%% models
model_folder = fullfile(data_folder, 'GLM_models');
model_files = dir(fullfile(model_folder, 'GLM_*.mat'));
all_files = dir(model_folder);
unexpected_files = setdiff({all_files.name}, {model_files.name});
if ~isempty(unexpected_files)
    fprintf('****Unexpected files in GLM_models folder****:\n');
    for i = 1:length(unexpected_files)
        fprintf(' - %s\n', unexpected_files{i});
    end
    fprintf('****End of unexpected files****\n\n');
end

all_tasks = struct();
for i = 1:length(model_files)
    model_file = model_files(i);
    file_name = model_file.name;
    pattern = "^GLM_(\w+)_s(\d+)_shuffle(\d+)_(\w+)_(.+?)_epoch(\d+)_fold(\d+)\.mat$";

    t = regexp(file_name, pattern, "tokens", "once");

    task_name   = t{1};
    session_idx = str2double(t{2});
    shuffle_idx = str2double(t{3});
    Kernel      = t{4};
    reg         = t{5};
    epoch       = str2double(t{6});
    fold_idx    = str2double(t{7});

    if ~isfield(all_tasks, task_name)
        all_tasks.(task_name) = struct();
        all_tasks.(task_name).config_keys = {};
        all_tasks.(task_name).config_contents = {};
    end

    config_key = sprintf('kernel%s_reg%s', Kernel, reg);
    if ~ismember(config_key, all_tasks.(task_name).config_keys)
        all_tasks.(task_name).config_keys{end+1} = config_key;
        config = struct();
        config.kernel = Kernel;
        config.reg = reg;
        config.max_shuffle = -1;
        config.max_epoch = -1;
        config.max_fold = -1;
        config.sessions = [];
        config.model_count = 0;
        all_tasks.(task_name).config_contents{end+1} = config;
    end
    config_idx = find(strcmp(all_tasks.(task_name).config_keys, config_key));
    config = all_tasks.(task_name).config_contents{config_idx};
    config.max_shuffle = max(config.max_shuffle, shuffle_idx);
    config.max_epoch = max(config.max_epoch, epoch);
    config.max_fold = max(config.max_fold, fold_idx);
    config.model_count = config.model_count + 1;
    if ~ismember(session_idx, config.sessions)
        config.sessions = [config.sessions, session_idx];
    end
    all_tasks.(task_name).config_contents{config_idx} = config;
end

%% Print summary
task_names = fieldnames(all_tasks);
for i = 1:length(task_names)
    task_name = task_names{i};
    fprintf('Task: %s\n', task_name);
    configs = all_tasks.(task_name).config_keys;
    for j = 1:length(configs)
        config_key = configs{j};
        config = all_tasks.(task_name).config_contents{j};
        sessions = sort(config.sessions);
        session_num = length(sessions);

        fprintf(' - Kernel: %s, Reg: %s\n', config.kernel, config.reg);
        fprintf('   - Max Shuffle: %d, Max Epoch: %d, Max Fold: %d, Model Count: %d\n', ...
            config.max_shuffle, config.max_epoch, config.max_fold, config.model_count);
        fprintf('   - Sessions: ');
        % compress consecutive sessions into ranges
        start_idx = 1;
        for k = 2:length(sessions)
            if sessions(k) ~= sessions(k-1) + 1
                if start_idx == k-1
                    fprintf('%d, ', sessions(start_idx));
                else
                    fprintf('%d-%d, ', sessions(start_idx), sessions(k-1));
                end
                start_idx = k;
            end
        end
        if start_idx == length(sessions)
            fprintf('%d, ', sessions(start_idx));
        else
            fprintf('%d-%d, ', sessions(start_idx), sessions(end));
        end
        fprintf('\n');
        fprintf('   - Total sessions: %d\n\n', session_num);
    end
end