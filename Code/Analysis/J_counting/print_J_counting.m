%% print_J_counting.m - Print J counting results in text files
% 

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
session_types = {'SlayerSal', 'SlayerMus', 'EmperorSal', 'EmperorMus'};
states = {'RestOpen', 'RestClose'};
area_types = {'Within', 'Across'};
prepost_strs = {'Pre ', 'Post'};

for st_idx = 1:length(session_types)
    session_type = session_types{st_idx};
    % load J counting results
    folder_name = fullfile(root, 'Data', 'Working', 'J_count');
    file_name = sprintf('Jcount_%s.mat', session_type);
    file_path = fullfile(folder_name, file_name);
    load(file_path, 'J_count', 'J_count_by_area', 'J_ratio', 'J_ratio_by_area', 'max_count', 'max_count_by_area', 'session_num');
    % J_count_by_area/max_count_by_area: (within/across, session, posneg, kernel, state, prepost)

    for session_idx = 1:session_num
        for area_type_idx = 1:2 % within/across areas
            fprintf('=========================\n');
            fprintf('%s, Session: %d/%d, %s areas.\n', session_type, session_idx, session_num, area_types{area_type_idx});
            for kernel_idx = 1:3
                for state_idx = 1:length(states)
                    fprintf(' - State: %s, Kernel: %d\n', states{state_idx}, kernel_idx);
                    for prepost_idx = 1:2
                        prepost_str = prepost_strs{prepost_idx};
                        pos_count = J_count_by_area(area_type_idx, session_idx, 1, kernel_idx, state_idx, prepost_idx);
                        neg_count = J_count_by_area(area_type_idx, session_idx, 2, kernel_idx, state_idx, prepost_idx);
                        max_count_curr = max_count_by_area(area_type_idx, session_idx, 1, kernel_idx, state_idx, prepost_idx);
                        pos_ratio = J_ratio_by_area(area_type_idx, session_idx, 1, kernel_idx, state_idx, prepost_idx);
                        neg_ratio = J_ratio_by_area(area_type_idx, session_idx, 2, kernel_idx, state_idx, prepost_idx);
                        fprintf('   - %s: Pos %.2f%%(%d/%d), Neg %.2f%%(%d/%d)\n', ......
                            prepost_str, pos_ratio * 100, pos_count, max_count_curr, neg_ratio * 100, neg_count, max_count_curr);
                    end
                end
            end
        end
    end
end
