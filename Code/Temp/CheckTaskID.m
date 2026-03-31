% session_types = {'Muscimol', 'Saline', 'SimRec'};
% session_types = {'Muscimol_pre', 'Saline_pre', 'SimRec'};
session_types = {'Muscimol', 'Saline'};
% area_names = {'ACC', 'Thalamus', 'VLPFC'};
area_names = {'ACC', 'VLPFC'};

unique_sessions_all = ...
    {{'10272023', '11012023', '11102023', '11172023', '12012023',...
    '12082023', '12152023', '12292023', '01052024', '01122024'},...
    {'01302024', '02022024', '02092024', '02162024', '02292024'},...
    {'08112023', '08142023', '08152023', '08162023', '08172023'}};% Muscimol, Saline, SimRec

n_session_types = length(session_types);
n_area_names = length(area_names);

cell_ID_all = cell(n_session_types, n_area_names);
task_ID_all = cell(n_session_types, n_area_names);

% construct taskID labels
prepost = {'Pre', 'Post'};
stages = {'Offer1', 'Offer2', 'BeforeCho', 'Cho2Info', 'Info2Rew', 'AfterRew'};
regressors = {'ER', 'IxU', 'UNC', 'Sub'};
taskID_types = cell(length(prepost), length(stages), length(regressors));
for prepost_idx = 1:length(prepost)
    for stage_idx = 1:length(stages)
        for reg_idx = 1:length(regressors)
            taskID_types{prepost_idx, stage_idx, reg_idx} = ...
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
            file_name = ['../fromMengxi/Sl_', area_name, '_', session_type, '_prevspost_popChoice_uSelect_04112025_trM.mat'];
            load(file_name, 'uSelect');
            unitID = uSelect.unitID;
            N = size(unitID, 1);

            % load task ID
            taskID = cell(length(prepost), length(stages));
            taskID{1, 1} = uSelect.p_glmOff1_pre; taskID{1, 2} = uSelect.p_glmOff2_pre;
            taskID{2, 1} = uSelect.p_glmOff1_post; taskID{2, 2} = uSelect.p_glmOff2_post;
            taskID(1, 3:6) = uSelect.p_glmCho_pre;
            taskID(2, 3:6) = uSelect.p_glmCho_post;

            % unpack taskID
            taskID_unpacked = zeros(N, length(prepost), length(stages), length(regressors));
            for prepost_idx = 1:length(prepost)
                for stage_idx = 1:length(stages)
                    taskID_unpacked(:, prepost_idx, stage_idx, :) = taskID{prepost_idx, stage_idx};
                end
            end

            taskID = taskID_unpacked;

        elseif strcmp(area_name, 'Thalamus')
            file_name = ['../fromMengxi/Sl_', area_name, '_', session_type, '_pre_popPav.mat'];
            load(file_name, 'pop');
            unitID = pop.unitID;
        else
            error('Invalid area name: %s', area_name);
        end

        % store cell ID and task ID
        cell_ID_all{session_type_idx, area_name_idx} = unitID;
        task_ID_all{session_type_idx, area_name_idx} = taskID;

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