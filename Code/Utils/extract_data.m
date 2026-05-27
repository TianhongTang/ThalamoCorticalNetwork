%% extract_data.m - Copy required data from original files to new folder for easy access

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

% set up output folder
output_folder = fullfile(root, 'Data', 'Extracted');
check_path(output_folder);
datetime_str = string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
output_folder = fullfile(output_folder, datetime_str);
check_path(output_folder);

% load metadata
load(fullfile(root, 'Data', 'Working', 'Meta', 'metadata.mat'), 'metadata');

% File filter
% data_type = 'raster';
% mt = struct2table(metadata.(data_type));
% mt = mt(strcmp(mt.area, "Cortex") & ...
%         strcmp(mt.align, 'Longest') & ...
%         cellfun(@(x) ~isempty(x) && x == 15, mt.resting_dur_threshold), :);

data_type = 'GLM';
mt = struct2table(metadata.(data_type));
mt = mt(mt.epoch == 3000, :);

% Copy files
for i = 1:height(mt)
    file_name = mt.file_name{i};

    src_path = fullfile(root, 'Data', 'Working', data_type, file_name);
    dst_path = fullfile(output_folder, file_name);
    copyfile(src_path, dst_path);
end