%% update_meta.m - Scan current data folder and update metadata.mat.

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
metadata_folder = fullfile(root, 'Data', 'Meta');
backup_folder = fullfile(root, 'Data', 'Meta', 'backup');
data_folder = fullfile(root, 'Data', 'Working');
check_path(metadata_folder);
check_path(backup_folder);

data_folders = dir(fullfile(data_folder, '*'));
data_folders = data_folders([data_folders.isdir] & ~startsWith({data_folders.name}, '.'));

metadata = struct();
for i = 1:numel(data_folders)
    data_folder_name = data_folders(i).name;
    fprintf('Processing data folder: %s\n', data_folder_name);
    folder_path = fullfile(data_folder, data_folder_name);
    all_files = dir(fullfile(folder_path, '*.mat'));
    for file_idx = 1:numel(all_files)
        file_name = all_files(file_idx).name;
        file_path = fullfile(folder_path, file_name);
        metadata_entry = struct();
        metadata_entry.file_name = file_name;
        metadata_entry.folder_name = data_folder_name;
        metadata_entry.file_path = file_path;
        metadata_entry.processed_time = datetime('now');
        metadata.(sprintf('entry_%d', numel(fieldnames(metadata)) + 1)) = metadata_entry;
    end
end

current_time = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
metadata_backup_name = sprintf('metadata_backup_%s.mat', current_time);
metadata_backup_path = fullfile(backup_folder, metadata_backup_name);
if isfile(fullfile(metadata_folder, 'metadata.mat'))
    copyfile(fullfile(metadata_folder, 'metadata.mat'), metadata_backup_path);
    fprintf('Existing metadata backed up to: %s\n', metadata_backup_path);
else
    fprintf('No existing metadata found. No backup created.\n');
end