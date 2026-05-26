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
datetime_str = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
output_folder = fullfile(output_folder, datetime_str);
check_path(output_folder);

% load metadata
metadata = load(fullfile(root, 'Data', 'Meta', 'metadata.mat'));

% File filter
data_type = 'GLM';
mt = struct2table(metadata.(data_type));
mt = mt(mt.epoch == 3000, :);

% Copy files
for i = 1:height(mt)
    file_name = mt.file_name{i};

    src_path = fullfile(root, 'Data', data_type, file_name);
    dst_path = fullfile(output_folder, file_name);
    copyfile(src_path, dst_path);
end