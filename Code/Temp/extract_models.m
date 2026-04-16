%% Copy all used models to a separate folder for easier access
% keep original folder structure

clear;
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

target_root = fullfile(root, 'Copy');
check_path(target_root);

%% Main

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

% load kernel info
folder_name = fullfile(root, 'Data', 'Working', 'kernel');
file_name = sprintf('kernel_%s.mat', kernel);
file_path = fullfile(folder_name, file_name);
% copy kernel file
target_folder = fullfile(target_root, 'Data', 'Working', 'kernel');
check_path(target_folder);
copyfile(file_path, target_folder);

for task_idx = 1:length(tasks)
    task = tasks{task_idx};
    % load model data
    folder_name = fullfile(root, 'Data', 'Working', 'GLM_models');
    file_name = sprintf('GLM_%s_s%d_shuffle%d_%s_%s_epoch%d_fold0.mat', task.name, task.session_idx, 0, task.kernel, ...
        task.reg.name, task.epoch);
    file_path = fullfile(folder_name, file_name);
    % copy model file
    target_folder = fullfile(target_root, 'Data', 'Working', 'GLM_models');
    check_path(target_folder);
    copyfile(file_path, target_folder);

    % load border info, get area num. borders: starting index of each area
    folder_name = fullfile(root, 'Data', 'Working', 'border');
    file_name = sprintf('borders_%sPreCortex_%d.mat', task.control, task.session_idx);
    file_path = fullfile(folder_name, file_name);
    % copy border file
    target_folder = fullfile(target_root, 'Data', 'Working', 'border');
    check_path(target_folder);
    copyfile(file_path, target_folder);

end