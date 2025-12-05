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
unique_sessions_all = ...
    {{'10272023', '11012023', '11102023', '11172023', '12012023',...
    '12082023', '12152023', '12292023', '01052024', '01122024'},...
    {'01302024', '02022024', '02092024', '02162024', '02292024'},...
    {'08112023', '08142023', '08152023', '08162023', '08172023'}};

% load data
controls = {'Muscimol', 'Saline', 'SimRec'};
areas = {'ACC', 'Thalamus', 'VLPFC'};
area_num = length(areas);

for control_idx = 1:length(controls)
    control = controls{control_idx};
    unique_sessions = unique_sessions_all{control_idx};
    session_num = length(unique_sessions);

    PDS_data_all = cell(1, area_num);

    % load index files
    for area_idx = 1:length(areas)
        area = areas{area_idx};

        idx_file_folder = fullfile(root, 'Data', 'Experimental', 'taskID');
        idx_file = fullfile(idx_file_folder, 'Sl_Thalamus_Muscimol.mat');
        
        load(idx_file, "idxUncI_posi", "idxUncN_posi", "idxUncI_nega", "idxUncN_nega", ...
            "unitID_all");
        idx_all = [idxUncI_posi, idxUncN_posi, idxUncI_nega, idxUncN_nega];
    end


    for session = 1:10
        session_suffix = '-001_';

        raster_file = ['../GLM_data/MuscimolPre_full/raster_MuscimolPre_full_', int2str(session), ...
            '_0.mat'];
        load(raster_file, "cell_area", "cell_id", "session_name_full");
        N = length(cell_area);
        cell_type = zeros(1, N);
        for i=1:N
            if strcmp(cell_area{i},"Thalamus")
                cell_id_full = [session_name_full(1:8), session_suffix, cell_id{i}];
                idx_in_table = find(strcmp(unitID_all, cell_id_full));
                if ~isempty(idx_in_table)
                    % calc value
                    % posi_I = 101, posi_N = 102, posi_IN=103
                    % nega_I = 110, nega_N = 120, nega_IN=130
                    cell_type(i) = idxUncI_posi(idx_in_table)+...
                        idxUncN_posi(idx_in_table)*2+...
                        idxUncI_nega(idx_in_table)*10+...
                        idxUncN_nega(idx_in_table)*20+100;
                end
            end
        end
        type_file = ['../GLM_data/MuscimolPre_full/celltype_MuscimolPre_full_', int2str(session), ...
            '.mat'];
        save(type_file, "cell_type");
    end
end