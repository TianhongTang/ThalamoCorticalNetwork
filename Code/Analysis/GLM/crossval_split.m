function crossval_split(dataset_name, session, shuffle_id, fold_num, split_type)
% crossval_split - split shuffled raster data into training and testing sets for cross-validation
%
% Usage:
%   crossval_split(dataset_name, session, shuffle_id, fold_num, split_type)
%
% Inputs:
%   dataset_name - name of the dataset (string/char)
%   session      - session index (numeric)
%   shuffle_id   - shuffle identifier (numeric)
%   fold_num     - number of folds to split into (numeric)
%   split_type   - 'trial' or 'time' (string/char). Default: 'trial'
%                  'trial' - use whole trials as segments
%                  'time'  - split each trial into time segments (default segment length = 10 s)

% Notes:
%   Uses shuffled raster data saved by the shuffle.m function.

%% default parameters
if nargin < 5 || isempty(split_type)
    split_type = 'trial';
end
if isstring(split_type)
    split_type = char(split_type);
end

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
folder_name = fullfile(root, 'Data', 'Working', 'raster');
file_name = sprintf('shuffled_%s_%d_%d.mat', dataset_name, session, shuffle_id);
file_path = fullfile(folder_name, file_name);
if ~exist(file_path, 'file')
    error('crossval_split:FileNotFound', 'Shuffled raster file not found: %s', file_path);
end
load(file_path, "N", "n_trial", "rasters");

% Basic validation
if ~iscell(rasters) || numel(rasters) < n_trial
    error('crossval_split:InvalidRasters', 'Loaded ''rasters'' must be a cell array with at least n_trial elements.');
end

fold_rasters = cell(1, fold_num);
fold_total_lens = zeros(1, fold_num);
fold_trial_lens = cell(1, fold_num);
for i = 1:fold_num
    fold_rasters{i} = {};
    fold_trial_lens{i} = [];
end

assignments = {};
for i = 1:n_trial
    switch lower(split_type)
        case 'trial'
            % use the whole trial as one segment
            trial_rasters = rasters(i);
        case 'time'
            % split raster into fixed-length segments (in samples)
            % Default: 10000 samples per segment (~10 s at 1 kHz)
            segment_length_samples = 10000;
            [~, B] = size(rasters{i});
            segment_num = ceil(B / segment_length_samples);
            trial_rasters = cell(1, segment_num);
            for j = 1:segment_num
                start_idx = (j-1)*segment_length_samples + 1;
                end_idx = min(j*segment_length_samples, B);
                trial_rasters{j} = rasters{i}(:, start_idx:end_idx);
            end
        otherwise
            error('crossval_split:InvalidSplitType', 'Invalid split_type. Choose ''Trial'' or ''Time''.');
    end
    % assign segments to folds
    for j = 1:length(trial_rasters)
        shortest_fold = find(fold_total_lens == min(fold_total_lens), 1);
        fold_rasters{shortest_fold}{end+1} = trial_rasters{j};
        fold_total_lens(shortest_fold) = fold_total_lens(shortest_fold) + size(trial_rasters{j}, 2);
        fold_trial_lens{shortest_fold}{end+1} = size(trial_rasters{j}, 2);
        assignment = struct();
        assignment.trial_index = i;
        assignment.segment_index = j;
        assignment.fold = shortest_fold;
        assignment.length = size(trial_rasters{j}, 2);
        assignments{end+1} = assignment;
    end
end

save_folder = fullfile(root, 'Data', 'Working', 'crossval_split');
check_path(save_folder);
save_name = sprintf('crossval_%s_%d_%d.mat', dataset_name, session, shuffle_id);
save_path = fullfile(save_folder, save_name);

save(save_path, "N", "fold_num", "fold_rasters", "fold_trial_lens", "assignments", '-v7.3');