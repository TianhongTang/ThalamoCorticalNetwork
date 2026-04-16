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
reg_levels = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5];
% reg_names = {'Mix', 'L1', 'L2'};
reg_names = {'L2'};
regs = cell(length(reg_levels) * length(reg_names), 1);
for i = 1:length(reg_levels)
    reg_level = reg_levels(i);
    for j = 1:length(reg_names)
        reg_name = reg_names{j};
        reg = struct();
        switch reg_name
            case 'Mix'
                reg.l1 = reg_level;
                reg.l2 = reg_level;
            case 'L1'
                reg.l1 = reg_level;
                reg.l2 = 0;
            case 'L2'
                reg.l1 = 0;
                reg.l2 = reg_level;
        end
        reg.name = sprintf('%s=%d', reg_name, reg_level*100);
        regs{(i-1)*length(reg_names) + j} = reg;
    end
end
controls = {'Muscimol'};
session_idxs_all = {1:10}; % session indices for each control
area_types = {'Cortex'};
prepost_types = {'Pre', 'Post'};
states = {'RestOpen', 'RestClose', 'Task'};
alignments = {'Last'};
kernel = 'DeltaPure';

tasks = {};
found_count = 0;
total_count = 0;
% check valid files
for control_idx = 1:length(controls)
    control = controls{control_idx};
    sessions = session_idxs_all{control_idx};
    for session_idx = sessions
        for area_idx = 1:length(area_types)
            area_type = area_types{area_idx};
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

                        for reg_idx = 1:length(regs)
                            for fold_idx = 1:3
                                reg = regs{reg_idx};

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

                                folder_name = fullfile(root, 'Data', 'Working', 'GLM_models');
                                file_name = sprintf('GLM_%s_s%d_shuffle%d_%s_%s_epoch%d_fold%d.mat', ...
                                    task.name, task.session_idx, 0, task.kernel, task.reg.name, task.epoch, fold_idx);
                                file_path = fullfile(folder_name, file_name);
                                total_count = total_count + 1;
                                if isfile(file_path)
                                    found_count = found_count + 1;
                                else
                                    fprintf('Missing file: %s\n', file_path);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
fprintf('Found %d / %d files.\n', found_count, total_count);