%% J_counting.m - Count significant J numbers in each area.
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
% load metadata
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums', 'cortex_files', 'thalamus_files', 'eyeID_files');

% session_types = {'EmperorSal', 'EmperorMus'};
session_types = {'SlayerSal', 'SlayerMus', 'EmperorSal', 'EmperorMus'};
kernel = 'DeltaPure';
reg = 'L2=0_2';
epoch = '3000';
area_names = {'ACC', 'VLPFC', 'Thalamus'};
filter_threshold = 1;
% states = {'RestOpen', 'RestClose'};
% state_num = length(states);
prepost_types = {'Pre', 'Post'};

for session_type_idx = 1:length(session_types)
    session_type = session_types{session_type_idx};
    switch session_type
        case 'Muscimol'
            sessions = 1:10;
            states = {'Task', 'RestOpen', 'RestClose'};
        case 'Saline'
            sessions = 1:5;
            states = {'Task', 'RestOpen', 'RestClose'};
        case 'Simulated'
            sessions = 1:10;
            states = {'Task', 'RestClose'};
        case 'KZ'
            % load filtered sessions for KZ
            folder_name = fullfile(root, 'Data', 'Working', 'filtered_sessions');
            states = {'RestOpen', 'RestClose'};
        case 'SlayerSal'
            sessions = [1,3,4,5];
            states = {'RestOpen', 'RestClose'};
        case 'SlayerMus'
            sessions = 1:8;
            states = {'RestOpen', 'RestClose'};
        case 'EmperorSal'
            sessions = 1:3;
            states = {'RestOpen', 'RestClose'};
        case 'EmperorMus'
            sessions = 1:2;
            states = {'RestOpen', 'RestClose'};
    end

    session_num = length(sessions);
    state_num = length(states);

    % load kernel info
    folder_name = fullfile(root, 'Data', 'Working', 'kernel');
    file_name = sprintf('kernel_%s.mat', kernel);
    file_path = fullfile(folder_name, file_name);
    load(file_path, 'kernel_len', 'n_PS_kernel', 'n_conn_kernel');

    % (area i, area j, session, posneg, kernel, state, prepost)
    J_count = zeros(3, 3, session_num, 2, n_conn_kernel, state_num, 2);
    max_count = zeros(3, 3, session_num, 2, n_conn_kernel, state_num, 2); % maximum possible connections
    J_count_by_area = zeros(2, session_num, 2, n_conn_kernel, state_num, 2); % (within/across, session, posneg, kernel, state, prepost)
    max_count_by_area = zeros(2, session_num, 2, n_conn_kernel, state_num, 2); % (within/across, session, posneg, kernel, state, prepost)

    for prepost_idx = 1:2
        prepost = prepost_types{prepost_idx};
        for state_idx = 1:state_num
            state = states{state_idx};
            for session_idx = 1:session_num
                session_idx_file = sessions(session_idx);
                fprintf('Processing session %d/%d, state %s, prepost %s\n', session_idx, session_num, state, prepost);
                
                % load border info, get area num. borders: starting index of each area
                folder_name = fullfile(root, 'Data', 'Working', 'border');
                file_name = sprintf('borders_%sPreCortex_%d.mat', session_type, session_idx_file);
                file_path = fullfile(folder_name, file_name);
                load(file_path, 'borders');
                area_num = length(borders);
                borders_raw = borders;

                % load model data
                folder_name = fullfile(root, 'Data', 'Working', 'GLM_models');
                file_name = sprintf('GLM_%s%s%sCortexAlignLongest_s%d_shuffle0_%s_%s_epoch%s_fold0.mat', ...
                    session_type, prepost, state, session_idx_file, kernel, reg, epoch); 
                file_path = fullfile(folder_name, file_name);
                load(file_path, 'model_par', 'model_err', 'N');
                borders = [borders_raw, N+1]; % add end border

                for kernel_idx = 1:n_conn_kernel
                    % extract full J matrix
                    J_mat = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
                    J_err = model_err.total(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));

                    % count significant connections
                    for i = 1:area_num
                        for j = 1:area_num
                            i_filter = borders(i):borders(i+1)-1;
                            j_filter = borders(j):borders(j+1)-1;
                            i_num = borders(i+1) - borders(i);
                            j_num = borders(j+1) - borders(j);
                            data_mat = J_mat(i_filter, j_filter);
                            error_mat = J_err(i_filter, j_filter);
                            % count significant positive and negative connections
                            pos_count = sum(data_mat > filter_threshold * error_mat, "all");
                            neg_count = sum(data_mat < -filter_threshold * error_mat, "all");
                            J_count(i, j, session_idx, 1, kernel_idx, state_idx, prepost_idx) = pos_count;
                            J_count(i, j, session_idx, 2, kernel_idx, state_idx, prepost_idx) = neg_count;
                            if i == j
                                max_count(i, j, session_idx, 1, kernel_idx, state_idx, prepost_idx) = i_num * (i_num - 1); % exclude self-connection
                                max_count(i, j, session_idx, 2, kernel_idx, state_idx, prepost_idx) = i_num * (i_num - 1);
                                J_count_by_area(1, session_idx, 1, kernel_idx, state_idx, prepost_idx) = ...
                                    J_count_by_area(1, session_idx, 1, kernel_idx, state_idx, prepost_idx) + pos_count;
                                J_count_by_area(1, session_idx, 2, kernel_idx, state_idx, prepost_idx) = ...
                                    J_count_by_area(1, session_idx, 2, kernel_idx, state_idx, prepost_idx) + neg_count;
                                max_count_by_area(1, session_idx, 1, kernel_idx, state_idx, prepost_idx) = ...
                                    max_count_by_area(1, session_idx, 1, kernel_idx, state_idx, prepost_idx) + i_num * (i_num - 1);
                                max_count_by_area(1, session_idx, 2, kernel_idx, state_idx, prepost_idx) = ...
                                    max_count_by_area(1, session_idx, 2, kernel_idx, state_idx, prepost_idx) + i_num * (i_num - 1);
                            else
                                max_count(i, j, session_idx, 1, kernel_idx, state_idx, prepost_idx) = i_num * j_num;
                                max_count(i, j, session_idx, 2, kernel_idx, state_idx, prepost_idx) = i_num * j_num;
                                J_count_by_area(2, session_idx, 1, kernel_idx, state_idx, prepost_idx) = ...
                                    J_count_by_area(2, session_idx, 1, kernel_idx, state_idx, prepost_idx) + pos_count;
                                J_count_by_area(2, session_idx, 2, kernel_idx, state_idx, prepost_idx) = ...
                                    J_count_by_area(2, session_idx, 2, kernel_idx, state_idx, prepost_idx) + neg_count;
                                max_count_by_area(2, session_idx, 1, kernel_idx, state_idx, prepost_idx) = ...
                                    max_count_by_area(2, session_idx, 1, kernel_idx, state_idx, prepost_idx) + i_num * j_num;
                                max_count_by_area(2, session_idx, 2, kernel_idx, state_idx, prepost_idx) = ...
                                    max_count_by_area(2, session_idx, 2, kernel_idx, state_idx, prepost_idx) + i_num * j_num;
                            end
                        end
                    end
                end
            end
        end
    end

    J_ratio = J_count ./ max_count;
    J_ratio(isnan(J_ratio)) = 0;
    J_ratio_by_area = J_count_by_area ./ max_count_by_area;
    J_ratio_by_area(isnan(J_ratio_by_area)) = 0;

    % save results
    folder_name = fullfile(root, 'Data', 'Working', 'J_count');
    check_path(folder_name);
    file_name = sprintf('Jcount_%s.mat', session_type);
    file_path = fullfile(folder_name, file_name);
    save(file_path, 'J_count', 'J_count_by_area', 'J_ratio', 'J_ratio_by_area', 'max_count', 'max_count_by_area', 'session_num', 'kernel', 'reg', 'epoch');
    
    % % Export to Excel
    % excel_mat = zeros(session_num*state_num*2, n_conn_kernel*2*2); % (session*state*prepost, kernel*within/across*posneg)
    % row_idx = 1;
    % for session_idx = 1:session_num
    %     for state_idx = 1:state_num
    %         for prepost_idx = 1:2
    %             for kernel_idx = 1:n_conn_kernel
    %                 % within-area positive
    %                 excel_mat(row_idx, (kernel_idx-1)*4 + 1) = J_ratio_by_area(1, session_idx, 1, kernel_idx, state_idx, prepost_idx);
    %                 % within-area negative
    %                 excel_mat(row_idx, (kernel_idx-1)*4 + 2) = J_ratio_by_area(1, session_idx, 2, kernel_idx, state_idx, prepost_idx);
    %                 % across-area positive
    %                 excel_mat(row_idx, (kernel_idx-1)*4 + 3) = J_ratio_by_area(2, session_idx, 1, kernel_idx, state_idx, prepost_idx);
    %                 % across-area negative
    %                 excel_mat(row_idx, (kernel_idx-1)*4 + 4) = J_ratio_by_area(2, session_idx, 2, kernel_idx, state_idx, prepost_idx);
    %             end
    %             row_idx = row_idx + 1;
    %         end
    %     end
    % end

    % excel_folder = fullfile(root, 'Data', 'Working', 'J_count', 'Excel');
    % check_path(excel_folder);
    % excel_file = sprintf('Jcount_%s.xlsx', session_type);
    % excel_path = fullfile(excel_folder, excel_file);
    % writematrix(excel_mat, excel_path, 'Sheet', 'J_ratio_by_area');

end