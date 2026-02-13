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

mode = 'Full'; % Cortex: Pre and Post; Full: only Pre

kernel = 'DeltaPure';
reg = 'L2=0_2';
epoch = '3000';
area_names = {'ACC', 'VLPFC', 'Thalamus'};
filter_threshold = 1;
% states = {'RestOpen', 'RestClose'};
% state_num = length(states);
if strcmp(mode, 'Cortex') %#ok<*UNRCH>
    session_types = {'SlayerSal', 'SlayerMus', 'ZeppelinMus', 'ZeppelinSal', 'EmperorSal', 'EmperorMus'};
    prepost_types = {'Pre', 'Post'};
    area_type = 'Cortex';
elseif strcmp(mode, 'Full') %#ok<*UNRCH>
    session_types = {'SlayerSal', 'SlayerMus', 'ZeppelinNoinj', 'ZeppelinMus', 'ZeppelinSal', 'EmperorSal', 'EmperorMus'};
    prepost_types = {'Pre'};
    area_type = 'Full';
end
align = 'AlignLast';

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
            sessions = [1,2,3,4,5];
            states = {'RestOpen', 'RestClose'};
        case 'SlayerMus'
            sessions = 1:8;
            states = {'RestOpen', 'RestClose'};
        case 'ZeppelinNoinj'
            sessions = 1:8;
            states = {'RestOpen', 'RestClose'};
        case 'ZeppelinMus'
            sessions = 1:1;
            states = {'RestOpen', 'RestClose'};
        case 'ZeppelinSal'
            sessions = 1:1;
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
    J_count_by_area = zeros(3, 3, session_num, 2, n_conn_kernel, state_num, numel(prepost_types));
    max_count_by_area = zeros(3, 3, session_num, 2, n_conn_kernel, state_num, numel(prepost_types)); % maximum possible connections
    J_count = zeros(2, session_num, 2, n_conn_kernel, state_num, numel(prepost_types)); % (within/across, session, posneg, kernel, state, prepost)
    max_count = zeros(2, session_num, 2, n_conn_kernel, state_num, numel(prepost_types)); % (within/across, session, posneg, kernel, state, prepost)
    
    % disagreement between open/close and pre/post, for McNemar's test.
    disagreement_resting = zeros(2, 2, session_num, 2, n_conn_kernel, numel(prepost_types)); % (n10/n01, within/across, session, posneg, kernel, prepost)
    disagreement_prepost = zeros(2, 2, session_num, 2, n_conn_kernel, state_num); % (n10/n01, within/across, session, posneg, kernel, state)
    disagreement_mats = cell(2, session_num, 2, n_conn_kernel, 2, 2); % (within/across, session, posneg, kernel, state, prepost)
    disagreement_resting_by_area = zeros(2, 3, 3, session_num, 2, n_conn_kernel, numel(prepost_types)); % (n10/n01, i, j, session, posneg, kernel, prepost)
    disagreement_prepost_by_area = zeros(2, 3, 3, session_num, 2, n_conn_kernel, state_num); % (n10/n01, i, j, session, posneg, kernel, state)
    disagreement_mats_by_area = cell(3, 3, session_num, 2, n_conn_kernel, 2, 2); % (i, j, session, posneg, kernel, state, prepost)

    for prepost_idx = 1:numel(prepost_types)
        prepost = prepost_types{prepost_idx};
        for state_idx = 1:state_num
            state = states{state_idx};
            for session_idx = 1:session_num
                session_idx_file = sessions(session_idx);
                fprintf('Processing session %d/%d, state %s, prepost %s\n', session_idx, session_num, state, prepost);
                
                % load border info, get area num. borders: starting index of each area
                folder_name = fullfile(root, 'Data', 'Working', 'border');
                file_name = sprintf('borders_%sPre%s_%d.mat', session_type, area_type, session_idx_file);
                file_path = fullfile(folder_name, file_name);
                load(file_path, 'borders');
                area_num = length(borders);
                borders_raw = borders;

                % load model data
                folder_name = fullfile(root, 'Data', 'Working', 'GLM_models');
                file_name = sprintf('GLM_%s%s%s%s%s_s%d_shuffle0_%s_%s_epoch%s_fold0.mat', ...
                    session_type, prepost, state, area_type, align, session_idx_file, kernel, reg, epoch); 
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
                            pos_mat = data_mat > filter_threshold * error_mat;
                            neg_mat = data_mat < -filter_threshold * error_mat;
                            pos_count = sum(pos_mat, "all");
                            neg_count = sum(neg_mat, "all");
                            J_count_by_area(i, j, session_idx, 1, kernel_idx, state_idx, prepost_idx) = pos_count;
                            J_count_by_area(i, j, session_idx, 2, kernel_idx, state_idx, prepost_idx) = neg_count;
                            if i == j
                                max_count_by_area(i, j, session_idx, 1, kernel_idx, state_idx, prepost_idx) = i_num * (i_num - 1); % exclude self-connection
                                max_count_by_area(i, j, session_idx, 2, kernel_idx, state_idx, prepost_idx) = i_num * (i_num - 1);
                                J_count(1, session_idx, 1, kernel_idx, state_idx, prepost_idx) = ...
                                    J_count(1, session_idx, 1, kernel_idx, state_idx, prepost_idx) + pos_count;
                                J_count(1, session_idx, 2, kernel_idx, state_idx, prepost_idx) = ...
                                    J_count(1, session_idx, 2, kernel_idx, state_idx, prepost_idx) + neg_count;
                                max_count(1, session_idx, 1, kernel_idx, state_idx, prepost_idx) = ...
                                    max_count(1, session_idx, 1, kernel_idx, state_idx, prepost_idx) + i_num * (i_num - 1);
                                max_count(1, session_idx, 2, kernel_idx, state_idx, prepost_idx) = ...
                                    max_count(1, session_idx, 2, kernel_idx, state_idx, prepost_idx) + i_num * (i_num - 1);
                            else
                                max_count_by_area(i, j, session_idx, 1, kernel_idx, state_idx, prepost_idx) = i_num * j_num;
                                max_count_by_area(i, j, session_idx, 2, kernel_idx, state_idx, prepost_idx) = i_num * j_num;
                                J_count(2, session_idx, 1, kernel_idx, state_idx, prepost_idx) = ...
                                    J_count(2, session_idx, 1, kernel_idx, state_idx, prepost_idx) + pos_count;
                                J_count(2, session_idx, 2, kernel_idx, state_idx, prepost_idx) = ...
                                    J_count(2, session_idx, 2, kernel_idx, state_idx, prepost_idx) + neg_count;
                                max_count(2, session_idx, 1, kernel_idx, state_idx, prepost_idx) = ...
                                    max_count(2, session_idx, 1, kernel_idx, state_idx, prepost_idx) + i_num * j_num;
                                max_count(2, session_idx, 2, kernel_idx, state_idx, prepost_idx) = ...
                                    max_count(2, session_idx, 2, kernel_idx, state_idx, prepost_idx) + i_num * j_num;
                            end

                            % record disagreement for McNemar's test
                            if strcmp(state, 'RestOpen')
                                state_idx_disagree = 1;
                            elseif strcmp(state, 'RestClose')
                                state_idx_disagree = 2;
                            else
                                state_idx_disagree = 0; % not included in disagreement analysis
                            end
                            if strcmp(prepost, 'Pre')
                                prepost_idx_disagree = 1;
                            elseif strcmp(prepost, 'Post')
                                prepost_idx_disagree = 2;
                            end
                            if i==j
                                area_type_idx = 1; % within-area
                            else
                                area_type_idx = 2; % across-area
                            end
                            
                            % record mat contents for disagreement analysis (for McNemar's test)
                            if state_idx_disagree > 0 && prepost_idx_disagree > 0
                                pos_mat_all = disagreement_mats{area_type_idx, session_idx, 1, kernel_idx, state_idx_disagree, prepost_idx_disagree}; % positive
                                neg_mat_all = disagreement_mats{area_type_idx, session_idx, 2, kernel_idx, state_idx_disagree, prepost_idx_disagree}; % negative
                                pos_mat_by_area = disagreement_mats_by_area{i, j, session_idx, 1, kernel_idx, state_idx_disagree, prepost_idx_disagree};
                                neg_mat_by_area = disagreement_mats_by_area{i, j, session_idx, 2, kernel_idx, state_idx_disagree, prepost_idx_disagree};
                                pos_mat_all = [pos_mat_all; pos_mat(:)];
                                neg_mat_all = [neg_mat_all; neg_mat(:)];
                                pos_mat_by_area = [pos_mat_by_area; pos_mat(:)];
                                neg_mat_by_area = [neg_mat_by_area; neg_mat(:)];
                                disagreement_mats{area_type_idx, session_idx, 1, kernel_idx, state_idx_disagree, prepost_idx_disagree} = pos_mat_all;
                                disagreement_mats{area_type_idx, session_idx, 2, kernel_idx, state_idx_disagree, prepost_idx_disagree} = neg_mat_all;
                                disagreement_mats_by_area{i, j, session_idx, 1, kernel_idx, state_idx_disagree, prepost_idx_disagree} = pos_mat_by_area;
                                disagreement_mats_by_area{i, j, session_idx, 2, kernel_idx, state_idx_disagree, prepost_idx_disagree} = neg_mat_by_area;
                            end
                        end
                    end
                end
            end
        end
    end

    % disagreement counts for McNemar's test
    % There must be a better way to organize this...
    for session_idx = 1:session_num
        for posneg_idx = 1:2
            for kernel_idx = 1:n_conn_kernel
                % for reference
                % disagreement_resting (n10/n01, within/across, session, posneg, kernel, prepost)
                % disagreement_prepost (n10/n01, within/across, session, posneg, kernel, state)

                % within/across
                for area_type_idx = 1:2
                    open_pre = disagreement_mats{area_type_idx, session_idx, posneg_idx, kernel_idx, 1, 1}; % state 1, prepost 1
                    open_post = disagreement_mats{area_type_idx, session_idx, posneg_idx, kernel_idx, 1, 2}; % state 1, prepost 2
                    close_pre = disagreement_mats{area_type_idx, session_idx, posneg_idx, kernel_idx, 2, 1}; % state 2, prepost 1
                    close_post = disagreement_mats{area_type_idx, session_idx, posneg_idx, kernel_idx, 2, 2}; % state 2, prepost 2
                    
                    % resting state disagreement
                    if numel(open_pre) == numel(close_pre) && numel(open_pre) > 0
                        n10_pre = sum(open_pre & ~close_pre); % positive in open but not close
                        n01_pre = sum(~open_pre & close_pre); % positive in close but not open
                        disagreement_resting(1, area_type_idx, session_idx, posneg_idx, kernel_idx, 1) = n10_pre;
                        disagreement_resting(2, area_type_idx, session_idx, posneg_idx, kernel_idx, 1) = n01_pre;
                    else
                        fprintf('Warning: Disagreement analysis for session %d, posneg %d, kernel %d, area type %d in pre state: open %d, close %d\n',...
                         session_idx, posneg_idx, kernel_idx, area_type_idx, numel(open_pre), numel(close_pre));
                    end
                    if numel(open_post) == numel(close_post) && numel(open_post) > 0
                        n10_post = sum(open_post & ~close_post); % positive in open but not close
                        n01_post = sum(~open_post & close_post); % positive in close but not open
                        disagreement_resting(1, area_type_idx, session_idx, posneg_idx, kernel_idx, 2) = n10_post;
                        disagreement_resting(2, area_type_idx, session_idx, posneg_idx, kernel_idx, 2) = n01_post;
                    else
                        fprintf('Warning: Disagreement analysis for session %d, posneg %d, kernel %d, area type %d in post state: open %d, close %d\n',...
                         session_idx, posneg_idx, kernel_idx, area_type_idx, numel(open_post), numel(close_post));
                    end

                    % prepost disagreement
                    if numel(open_pre) == numel(open_post) && numel(open_pre)
                        n10_open = sum(open_pre & ~open_post); % positive in pre but not post
                        n01_open = sum(~open_pre & open_post); % positive in post but not pre
                        disagreement_prepost(1, area_type_idx, session_idx, posneg_idx, kernel_idx, 1) = n10_open;
                        disagreement_prepost(2, area_type_idx, session_idx, posneg_idx, kernel_idx, 1) = n01_open;
                    else
                        fprintf('Warning: Disagreement analysis for session %d, posneg %d, kernel %d, area type %d in open state: pre %d, post %d\n',...
                         session_idx, posneg_idx, kernel_idx, area_type_idx, numel(open_pre), numel(open_post));
                    end
                    if numel(close_pre) == numel(close_post) && numel(close_pre) > 0
                        n10_close = sum(close_pre & ~close_post); % positive in pre but not post
                        n01_close = sum(~close_pre & close_post); % positive in post but not pre
                        disagreement_prepost(1, area_type_idx, session_idx, posneg_idx, kernel_idx, 2) = n10_close;
                        disagreement_prepost(2, area_type_idx, session_idx, posneg_idx, kernel_idx, 2) = n01_close;
                    else
                        fprintf('Warning: Disagreement analysis for session %d, posneg %d, kernel %d, area type %d in close state: pre %d, post %d\n',...
                         session_idx, posneg_idx, kernel_idx, area_type_idx, numel(close_pre), numel(close_post));
                    end
                end

                % by area
                for i = 1:area_num
                    for j = 1:area_num
                        open_pre = disagreement_mats_by_area{i, j, session_idx, posneg_idx, kernel_idx, 1, 1}; % state 1, prepost 1
                        open_post = disagreement_mats_by_area{i, j, session_idx, posneg_idx, kernel_idx, 1, 2}; % state 1, prepost 2
                        close_pre = disagreement_mats_by_area{i, j, session_idx, posneg_idx, kernel_idx, 2, 1}; % state 2, prepost 1
                        close_post = disagreement_mats_by_area{i, j, session_idx, posneg_idx, kernel_idx, 2, 2}; % state 2, prepost 2
                        
                        % resting state disagreement
                        if numel(open_pre) == numel(close_pre) && numel(open_pre) > 0
                            n10_pre = sum(open_pre & ~close_pre); % positive in open but not close
                            n01_pre = sum(~open_pre & close_pre); % positive in close but not open
                            disagreement_resting_by_area(1, i, j, session_idx, posneg_idx, kernel_idx, 1) = n10_pre;
                            disagreement_resting_by_area(2, i, j, session_idx, posneg_idx, kernel_idx, 1) = n01_pre;
                        else
                            fprintf('Warning: Disagreement analysis for session %d, posneg %d, kernel %d, area %d-%d in pre state: open %d, close %d\n',...
                            session_idx, posneg_idx, kernel_idx, i, j, numel(open_pre), numel(close_pre));
                        end
                        if numel(open_post) == numel(close_post) && numel(open_post) > 0
                            n10_post = sum(open_post & ~close_post); % positive in open but not close
                            n01_post = sum(~open_post & close_post); % positive in close but not open
                            disagreement_resting_by_area(1, i, j, session_idx, posneg_idx, kernel_idx, 2) = n10_post;
                            disagreement_resting_by_area(2, i, j, session_idx, posneg_idx, kernel_idx, 2) = n01_post;
                        else
                            fprintf('Warning: Disagreement analysis for session %d, posneg %d, kernel %d, area %d-%d in post state: open %d, close %d\n',...
                            session_idx, posneg_idx, kernel_idx, i, j, numel(open_post), numel(close_post));
                        end

                        % prepost disagreement
                        if numel(open_pre) == numel(open_post) && numel(open_pre)
                            n10_open = sum(open_pre & ~open_post); % positive in pre but not post
                            n01_open = sum(~open_pre & open_post); % positive in post but not pre
                            disagreement_prepost_by_area(1, i, j, session_idx, posneg_idx, kernel_idx, 1) = n10_open;
                            disagreement_prepost_by_area(2, i, j, session_idx, posneg_idx, kernel_idx, 1) = n01_open;
                        else
                            fprintf('Warning: Disagreement analysis for session %d, posneg %d, kernel %d, area %d-%d in open state: pre %d, post %d\n',...
                            session_idx, posneg_idx, kernel_idx, i, j, numel(open_pre), numel(open_post));
                        end
                        if numel(close_pre) == numel(close_post) && numel(close_pre) > 0
                            n10_close = sum(close_pre & ~close_post); % positive in pre but not post
                            n01_close = sum(~close_pre & close_post); % positive in post but not pre
                            disagreement_prepost_by_area(1, i, j, session_idx, posneg_idx, kernel_idx, 2) = n10_close;
                            disagreement_prepost_by_area(2, i, j, session_idx, posneg_idx, kernel_idx, 2) = n01_close;
                        else
                            fprintf('Warning: Disagreement analysis for session %d, posneg %d, kernel %d, area %d-%d in close state: pre %d, post %d\n',...
                            session_idx, posneg_idx, kernel_idx, i, j, numel(close_pre), numel(close_post));
                        end
                    end
                end
            end
        end
    end

    
    % ratios
    J_ratio = J_count ./ max_count;
    J_ratio(isnan(J_ratio)) = 0;
    J_ratio_by_area = J_count_by_area ./ max_count_by_area;
    J_ratio_by_area(isnan(J_ratio_by_area)) = 0;

    % save results
    folder_name = fullfile(root, 'Data', 'Working', 'J_count');
    check_path(folder_name);
    file_name = sprintf('Jcount_%s_%s.mat', mode, session_type);
    file_path = fullfile(folder_name, file_name);
    save(file_path, 'J_count', 'J_count_by_area', 'J_ratio', 'J_ratio_by_area', 'max_count', 'max_count_by_area',...
        'disagreement_resting', 'disagreement_prepost', 'disagreement_resting_by_area', 'disagreement_prepost_by_area',...
        'session_num', 'kernel', 'reg', 'epoch');
    
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