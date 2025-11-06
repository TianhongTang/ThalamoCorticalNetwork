%% J_counting.m - Count significant J numbers in each area.
% 

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
clear;

session_types = {'Muscimol', 'Saline', 'Simulated'};

kernel = 'DeltaPure';
reg = 'L2=20';
epoch = '2500';
area_names = {'ACC', 'VLPFC', 'Thalamus'};
filter_threshold = 1;
states = {'Task', 'RestOpen', 'RestClose'};
state_num = length(states);

for session_type_idx = 1:length(session_types)
    session_type = session_types{session_type_idx};

    switch session_type
        case 'Muscimol'
            session_num = 10;
            states = {'Task', 'RestOpen', 'RestClose'};
        case 'Saline'
            session_num = 5;
            states = {'Task', 'RestOpen', 'RestClose'};
        case 'Simulated'
            session_num = 10;
            states = {'Task', 'RestClose'};
    end

    % load border info, get area_num. borders: starting index of each area
    folder_name = fullfile(root, 'Data', 'Working', 'border');
    file_name = sprintf('borders_%sPreCortex_1.mat', session_type);
    file_path = fullfile(folder_name, file_name);
    load(file_path, 'borders');
    area_num = length(borders);
    borders_raw = borders;

    % load kernel info
    folder_name = fullfile(root, 'Data', 'Working', 'kernel');
    file_name = sprintf('kernel_%s.mat', kernel);
    file_path = fullfile(folder_name, file_name);
    load(file_path, 'kernel_len', 'n_PS_kernel', 'n_conn_kernel');

    % (area i, area j, session, posneg, kernel, state, prepost)
    J_count = zeros(3, 3, session_num, 2, n_conn_kernel, state_num, 2);

    for prepost_idx = 1:2
        prepost = prepost_types{prepost_idx};
        for state_idx = 1:state_num
            state = states{state_idx};
            for session_idx = 1:session_num
                fprintf('Processing session %d/%d, state %s, prepost %s\n', session_idx, session_num, state, prepost);

                % load model data
                folder_name = fullfile(root, 'Data', 'Working', 'GLM_models');
                file_name = sprintf('GLM_%s%s%sCortexAlignLast_s%d_shuffle0_%s_%s_epoch%s_fold0.mat', ...
                    session_type, prepost, state, session_idx, kernel, reg, epoch); 
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
                            data_mat = J_mat(i_filter, j_filter);
                            error_mat = J_err(i_filter, j_filter);
                            % count significant positive and negative connections
                            pos_count = sum(data_mat > filter_threshold * error_mat);
                            neg_count = sum(data_mat < -filter_threshold * error_mat);
                            J_count(i, j, session_idx, 1, kernel_idx, state_idx, prepost_idx) = pos_count;
                            J_count(i, j, session_idx, 2, kernel_idx, state_idx, prepost_idx) = neg_count;
                        end
                    end
                end
            end
        end
    end
    % save results
    folder_name = fullfile(root, 'Data', 'Working', 'J_count');
    check_path(folder_name);
    file_name = sprintf('Jcount_%s.mat', session_type);
    file_path = fullfile(folder_name, file_name);
    save(file_path, 'J_count', 'session_num', 'kernel', 'reg', 'epoch');
end