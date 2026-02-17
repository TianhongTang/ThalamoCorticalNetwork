%% Main script for preprocessing data

%% Get absolute folder
code_depth = 3;
script_path = mfilename('fullpath');
script_folder = fileparts(script_path);

root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
% include code folder and utils
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Main
fprintf('=====================\n');
fprintf('Preprocessing Phase 1/5: Registering PDS datasets...\n');
fprintf('=====================\n');
tic;
isoRun(fullfile(script_folder, 'PDS_dataset.m'));
toc;

fprintf('=====================\n');
fprintf('Preprocessing Phase 2/5: Matching resting ending time...\n');
fprintf('=====================\n');
tic;
isoRun(fullfile(script_folder, 'check_resting_dur.m'));
toc;

fprintf('=====================\n');
fprintf('Preprocessing Phase 3/5: Loading PDS files...\n');
fprintf('=====================\n');
tic;
isoRun(fullfile(script_folder, 'loadPDS.m'));
toc;

fprintf('=====================\n');
fprintf('Preprocessing Phase 4/5: Merging brain areas...\n');
fprintf('=====================\n');
tic;
isoRun(fullfile(script_folder, 'merge_session.m'));
toc;

fprintf('=====================\n');
fprintf('Preprocessing Phase 5/5: Aligning dataset to same length...\n');
fprintf('=====================\n');
tic;
isoRun(fullfile(script_folder, 'align_dataset.m'));
toc;