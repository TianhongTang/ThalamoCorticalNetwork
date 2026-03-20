function crossval_split(shuffled_meta, fold_num, split_type)
% crossval_split - split shuffled raster data into training and testing sets for cross-validation
%
% Usage:
%   crossval_split(shuffled_meta, fold_num, split_type)
%
% Inputs:
%   shuffled_meta - metadata for the shuffled raster data (struct)
%   fold_num     - number of folds to split into (numeric)
%   split_type   - 'trial' or 'time' (string/char). Default: 'trial'
%                  'trial' - use whole trials as segments
%                  'time'  - split each trial into time segments (default segment length = 10 s)

% Notes:
%   Uses shuffled raster data saved by the shuffle.m function.

%% default parameters
if nargin < 3 || isempty(split_type)
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
folder_name = fullfile(root, 'Data', 'Working', 'shuffled');
file_name = generate_filename('shuffled', shuffled_meta);
file_path = fullfile(folder_name, file_name);
if ~exist(file_path, 'file')
    error('crossval_split:FileNotFound', 'Shuffled raster file not found: %s', file_path);
end
shuffled_data = load(file_path, 'meta', 'data');
rasters   = shuffled_data.data.rasters;
N         = shuffled_data.meta.N;
trial_num = shuffled_data.meta.trial_num;


% Basic validation
if ~iscell(shuffled_data.data.rasters) || numel(shuffled_data.data.rasters) ~= shuffled_data.meta.trial_num
    error('crossval_split:InvalidRasters', 'Loaded ''rasters'' must be a cell array with trial_num elements.');
end

% data for each fold
fold_rasters    = cell(1, fold_num);
fold_total_lens = zeros(1, fold_num);
fold_trial_lens = cell(1, fold_num);
for i = 1:fold_num
    fold_rasters{i} = {};
    fold_trial_lens{i} = [];
end

assignments = {};
for i = 1:trial_num
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
        fold_trial_lens{shortest_fold}(end+1) = size(trial_rasters{j}, 2);
        assignment = struct();
        assignment.trial_index = i;
        assignment.segment_index = j;
        assignment.fold = shortest_fold;
        assignment.length = size(trial_rasters{j}, 2);
        assignments{end+1} = assignment; %#ok<AGROW>
    end
end

% construct meta and data for saving
meta = struct();
meta.animal_name     = shuffled_data.meta.animal_name;
meta.injection       = shuffled_data.meta.injection;
meta.prepost         = shuffled_data.meta.prepost;
meta.state           = shuffled_data.meta.state;
meta.area            = shuffled_data.meta.area;
meta.align           = shuffled_data.meta.align;
meta.session_idx     = shuffled_data.meta.session_idx;
meta.shuffle_idx     = shuffled_data.meta.shuffle_idx;
meta.file_name       = generate_filename('crossval', meta);
meta.N               = N;
meta.fold_num        = fold_num;
meta.assignment_num  = numel(assignments);
meta.fold_total_len  = fold_total_lens;
meta.fold_trial_lens = fold_trial_lens;

data = struct();
data.fold_rasters    = fold_rasters;
data.fold_trial_lens = fold_trial_lens;
data.assignments     = assignments;

save_folder = fullfile(root, 'Data', 'Working', 'crossval_split');
check_path(save_folder);
save_name = meta.file_name;
save_path = fullfile(save_folder, save_name);

save(save_path, "meta", "data", '-v7.3');
fprintf('Cross-validation split saved to: %s\n', save_path);
