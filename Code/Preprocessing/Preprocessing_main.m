%% Main script for preprocessing data

%% Get absolute folder
script_path = mfilename('fullpath');
script_folder = fileparts(script_path);
addpath(script_folder);

%% Main
clear;
fprintf('=====================\n');
fprintf('Preprocessing Phase 1/3: Loading PDS files...\n');
fprintf('=====================\n');
tic;
run(fullfile(script_folder, 'loadPDS.m'));
toc;

fprintf('=====================\n');
fprintf('Preprocessing Phase 2/3: Merging brain areas...\n');
fprintf('=====================\n');
tic;
run(fullfile(script_folder, 'merge_session.m'));
toc;

fprintf('=====================\n');
fprintf('Preprocessing Phase 3/3: Aligning dataset to same length...\n');
fprintf('=====================\n');
tic;
run(fullfile(script_folder, 'align_dataset.m'));
toc;