%% PDS_dataset.m - Register PDS datasets and save to meta folder.
% output: 
% /Data/Working/Meta/PDS_dataset_info.mat
% contains:
% dataset_num: number of datasets
% dataset_names: names of datasets
% session_nums: number of sessions in each dataset
% cortex_files: cell array of cortex file names for each dataset
% thalamus_files: cell array of thalamus file names for each dataset
% eyeID_files: cell array of eyeID file names for each dataset


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

%% Main
% dataset info
dataset_names = {'SlayerMus', 'SlayerSal', 'SlayerNoinj',...
                 'ZeppelinMus', 'ZeppelinSal', 'ZeppelinNoinj',...
                 'EmperorMus', 'EmperorSal', 'EmperorNoinj'};
                 
cortex_files = {...
    {'10272023-008', '11172023-008', '12012023-008',...
     '12082023-008', '12152023-008', '12292023-008', '01052024-008', '01122024-008'},... % Slayer Muscimol
    {'01302024-008', '02022024-008', '02092024-008', '02162024-008', '02292024-008'},... % Slayer Saline
    {'08112023-003', '08142023-001', '08152023-001', '08162023-001', '08172023-001'},... % Slayer No injection
    {'03072024-008'},... % Zeppelin Muscimol
    {'03122024-008'},... % Zeppelin Saline
    {'', '', '', '', '',...
     '', '', ''},... % Zeppelin No injection, these are thal only sessions
    {'01092026-008', '01162026-008', '01232026-008', '01302026-008'},... % Emperor Muscimol
    {'12122025-008', '12262025-008', '01022026-008'},... % Emperor Saline
    {},... % Emperor No injection
    };
thalamus_files = {...
    {'10272023-001', '11172023-001', '12012023-001',...
     '12082023-001', '12152023-001', '12292023-001', '01052024-001', '01122024-001'},... % Slayer Muscimol
    {'01302024-001', '02022024-001', '02092024-001', '02162024-001', '02292024-001'},... % Slayer Saline
    {'08112023-003', '08142023-001', '08152023-001', '08162023-001', '08172023-001'},... % Slayer No injection
    {'03072024-001'},... % Zeppelin Muscimol
    {'03122024-001'},... % Zeppelin Saline
    {'02122024-004', '02142024-002', '02212024-001', '02212024-002', '02222024-002',...
     '02282024-001', '02292024-002', '02292024-003'},... % Zeppelin No injection
    {'01092026-001', '01162026-001', '01232026-001', '01302026-001'},... % Emperor Muscimol
    {'12122025-001', '12262025-001', '01022026-001'},... % Emperor Saline
    {},... % Emperor No injection
    };
eyeID_files = {...
    {'10272023-008', '11172023-008', '12012023-008',...
     '12082023-008', '12152023-008', '12292023-008', '01052024-008', '01122024-008'},... % Slayer Muscimol
    {'01302024-008', '02022024-008', '02092024-008', '02162024-008', '02292024-008'},... % Slayer Saline
    {'', '', '', '', ''},... % Slayer No injection, no eye ID data.
    {'03072024-008'},... % Zeppelin Muscimol
    {'03122024-008'},... % Zeppelin Saline
    {'02122024-004', '02142024-002', '02212024-001', '02212024-002', '02222024-002',...
     '02282024-001', '02292024-002', '02292024-003'},... % Zeppelin No injection
    {'01092026-008', '01162026-008', '01232026-008', '01302026-008'},... % Emperor Muscimol
    {'12122025-008', '12262025-008', '01022026-008'},... % Emperor Saline
    {},... % Emperor No injection
    };


dataset_num = length(dataset_names);
session_nums = zeros(dataset_num, 1);
for i = 1:dataset_num
    session_nums(i) = length(cortex_files{i});
    assert(session_nums(i) == length(thalamus_files{i}), 'Number of cortex and thalamus files must match!');
    assert(session_nums(i) == length(eyeID_files{i}), 'Number of cortex and eyeID files must match!');
end

save_folder = fullfile(root, 'Data', 'Working', 'Meta');
check_path(save_folder);
save_path = fullfile(save_folder, 'PDS_dataset_info.mat');
save(save_path,'dataset_num', 'dataset_names', 'session_nums', 'cortex_files', 'thalamus_files', 'eyeID_files');