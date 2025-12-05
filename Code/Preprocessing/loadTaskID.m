% load cell encodings 

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
% load task ID

% session_types = {'Muscimol', 'Saline', 'SimRec'};
% session_types = {'Muscimol_pre', 'Saline_pre', 'SimRec'};
session_types = {'Muscimol', 'Saline'};
area_names = {'ACC', 'Thalamus', 'VLPFC'};
% area_names = {'ACC', 'VLPFC'};

unique_sessions_all = ...
    {{'10272023', '11012023', '11102023', '11172023', '12012023',...
    '12082023', '12152023', '12292023', '01052024', '01122024'},...
    {'01302024', '02022024', '02092024', '02162024', '02292024'},...
    {'08112023', '08142023', '08152023', '08162023', '08172023'}}; % Muscimol, Saline, SimRec

n_session_types = length(session_types);
n_area_names = length(area_names);

cell_ID_all = cell(n_session_types, n_area_names);
task_ID_all = cell(n_session_types, n_area_names);
beta_all = cell(n_session_types, n_area_names);

% construct taskID labels
prepost = {'Pre', 'Post'};
stages = {'Offer1', 'Offer2', 'BeforeCho', 'Cho2Info', 'Info2Rew', 'AfterRew'};
regressors = {'ER', 'IxU', 'UNC', 'Sub'};
taskID_labels = cell(length(prepost), length(stages), length(regressors));
for prepost_idx = 1:length(prepost)
    for stage_idx = 1:length(stages)
        for reg_idx = 1:length(regressors)
            taskID_labels{prepost_idx, stage_idx, reg_idx} = ...
                [prepost{prepost_idx}, stages{stage_idx}, '_', regressors{reg_idx}];
        end
    end
end

for session_type_idx = 1:n_session_types
    for area_name_idx = 1:n_area_names
        fprintf('-------------------\n');

        session_type = session_types{session_type_idx};
        area_name = area_names{area_name_idx};
        unique_sessions = unique_sessions_all{session_type_idx};
        n_sessions = length(unique_sessions);

        if strcmp(area_name, 'ACC') || strcmp(area_name, 'VLPFC')
            folder_name = fullfile(root, 'Data', 'Experimental', 'taskID');
            file_name = ['Sl_', area_name, '_', session_type, '_prevspost_popChoice_uSelect_04112025_trM.mat'];
            file_path = fullfile(folder_name, file_name);
            load(file_path, 'uSelect');
            unitID = uSelect.unitID;
            N = size(unitID, 1);

            % load task ID
            taskID = cell(length(prepost), length(stages));
            taskID{1, 1} = uSelect.p_glmOff1_pre; taskID{1, 2} = uSelect.p_glmOff2_pre;
            taskID{2, 1} = uSelect.p_glmOff1_post; taskID{2, 2} = uSelect.p_glmOff2_post;
            taskID(1, 3:6) = uSelect.p_glmCho_pre;
            taskID(2, 3:6) = uSelect.p_glmCho_post;

            beta = cell(length(prepost), length(stages));
            beta{1, 1} = uSelect.beta_glmOff1_pre{1}; beta{1, 2} = uSelect.beta_glmOff2_pre{1};
            beta{2, 1} = uSelect.beta_glmOff1_post{1}; beta{2, 2} = uSelect.beta_glmOff2_post{1};
            beta(1, 3:6) = uSelect.beta_glmCho_pre;
            beta(2, 3:6) = uSelect.beta_glmCho_post;

            % unpack taskID
            taskID_unpacked = zeros(N, length(prepost), length(stages), length(regressors));
            beta_unpacked = zeros(N, length(prepost), length(stages), length(regressors));
            for prepost_idx = 1:length(prepost)
                for stage_idx = 1:length(stages)
                    taskID_unpacked(:, prepost_idx, stage_idx, :) = taskID{prepost_idx, stage_idx};
                    beta_unpacked(:, prepost_idx, stage_idx, :) = beta{prepost_idx, stage_idx};
                end
            end

            taskID = taskID_unpacked;
            beta = beta_unpacked;

        elseif strcmp(area_name, 'Thalamus')
            folder_name = fullfile(root, 'Data', 'Experimental', 'taskID');
            file_name = ['Sl_', area_name, '_', session_type, '_pre_popPav.mat'];
            file_path = fullfile(folder_name, file_name);
            load(file_path, 'pop');
            unitID = pop.unitID;
        else
            error('Invalid area name: %s', area_name);
        end

        % store cell ID and task ID
        cell_ID_all{session_type_idx, area_name_idx} = unitID;
        task_ID_all{session_type_idx, area_name_idx} = taskID;
        beta_all{session_type_idx, area_name_idx} = beta;

        % print session info
        fprintf('Session type: %s, Area: %s\n', session_type, area_name);
        fprintf('Cell number: %d\n', N);

        available_sessions = {};

        for i = 1:N
            ID = unitID{i};
            % split by '_' and '-'
            ID_parts = strsplit(ID, {'_', '-'});
            if length(ID_parts) ~= 3
                error('Invalid ID format: %s', ID);
            end

            unit_session = ID_parts{1};
            unit_ID = ID_parts{3};
            % add session to available_sessions if not already present
            if ~ismember(unit_session, available_sessions)
                available_sessions{end+1} = unit_session;
            end
        end

        n_available = length(available_sessions);
        fprintf('Number of available sessions: %d\n', n_available);
        fprintf('Available sessions: %s\n', strjoin(available_sessions, ', '));
    end
end

%% Load current data, match cell ID and task ID

session_types = {'Muscimol', 'Saline'};
% states = {'Task', 'RestClose', 'RestOpen'};
states = {'PreTask', 'PostTask', 'PreRestClose', 'PostRestClose', 'PreRestOpen', 'PostRestOpen'};
n_states = length(states);
% area_names = {'ACC', 'VLPFC'};
area_names = {'ACC', 'Thalamus', 'VLPFC'};


for session_type_idx = 1:length(session_types)
    session_type = session_types{session_type_idx};
    
    if strcmp(session_type, 'Muscimol')
        n_sessions = 10;
    else
        n_sessions = 5;
    end

    for state_idx = 1:n_states
        state = states{state_idx};
        model_name_full = [session_type, state, 'CortexAlignLast'];

        for session_idx = 1:n_sessions

            if strcmp(session_type, 'Muscimol') && (session_idx == 2 || session_idx == 3)
                continue; % skip sessions 2 and 3 for Muscimol
            end
            folder_name = fullfile(root, 'Data', 'Working', 'raster');
            file_name = ['raster_', model_name_full, '_', num2str(session_idx), '.mat'];
            file_path = fullfile(folder_name, file_name);
            load(file_path, 'cell_area', 'cell_id', 'session_name');
            N = size(cell_area, 2);

            tuning_session = zeros(N, length(prepost), length(stages), length(regressors));
            beta_session = zeros(N, length(prepost), length(stages), length(regressors));
            found = 0;
            
            for i = 1:N
                ID = [session_name, '-008_', cell_id{i}];
                area = cell_area{i};
                area_idx = find(strcmp(area_names, area));

                cell_ID = cell_ID_all{session_type_idx, area_idx};
                task_ID = task_ID_all{session_type_idx, area_idx};
                beta = beta_all{session_type_idx, area_idx};

                cell_idx = find(strcmp(cell_ID, ID));
                if isempty(cell_idx)
                    fprintf('Cell ID %s not found in session %s, area %s\n', ID, session_name, area);
                    tuning_session(i, :, :, :) = NaN;
                    beta_session(i, :, :, :) = NaN;
                else
                    tuning_session(i, :, :, :) = task_ID(cell_idx, :, :, :);
                    beta_session(i, :, :, :) = beta(cell_idx, :, :, :);
                    found = found + 1;
                end
            end
            fprintf('Session %s, Area %s, Found %d/%d cells\n', session_name, area_names{area_idx}, found, N);

            % save task ID
            save_folder = fullfile(root, 'Data', 'Working', 'tuning');
            check_path(save_folder);
            save_file = ['tuning_', model_name_full, '_', num2str(session_idx), '.mat'];
            save_path = fullfile(save_folder, save_file);
            save(save_path, 'tuning_session', 'beta_session', 'taskID_labels', 'cell_area', 'cell_id', 'session_name_full');
        end
    end
end