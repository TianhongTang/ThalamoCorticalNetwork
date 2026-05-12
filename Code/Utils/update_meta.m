%% update_meta.m - Scan current data folder and update metadata.mat.

% TODO: Compare new meta and current meta.

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
GET_DATA_FIELDS = false;
debug = false;

data_folder     = fullfile(root, 'Data', 'Working');
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');
backup_folder   = fullfile(root, 'Data', 'Working', 'Meta', 'backup');
check_path(metadata_folder);
check_path(backup_folder);
old_metadata_path = fullfile(metadata_folder, 'metadata.mat');
if isfile(old_metadata_path)
    load(old_metadata_path, 'metadata');
    fprintf('Existing metadata loaded from: %s\n', old_metadata_path);
    old_metadata = metadata;
else
    fprintf('No existing metadata found. Starting with empty metadata.\n');
end


data_folders = dir(fullfile(data_folder, '*'));
data_folders = data_folders([data_folders.isdir] & ~startsWith({data_folders.name}, '.'));

bad_files = {};

metadata = struct();
for i = 1:numel(data_folders)
    data_folder_name = data_folders(i).name;
    if strcmp(data_folder_name, 'Meta')
        continue; % skip metadata folder
    end

    fprintf('Processing data folder: %s\n', data_folder_name);
    folder_path = fullfile(data_folder, data_folder_name);
    all_files = dir(fullfile(folder_path, '*.mat'));

    metadata_entry = struct();
    for file_idx = 1:numel(all_files)
        fprintf('  Processing file %d/%d: %s\n', file_idx, numel(all_files), all_files(file_idx).name); 
        file_name = all_files(file_idx).name;
        file_path = fullfile(folder_path, file_name);

        try
            clear meta;
            load(file_path, "meta");
            % m = matfile(file_path);
            % vars = who(m);
            if ~exist('meta', 'var')
                warning('File %s does not contain meta variable. Skipping.', file_path);
                continue;
            end
            
            % meta = m.meta;
            if GET_DATA_FIELDS
                clear data;
                load(file_path, "data");

                if ~exist('data', 'var') %#ok<*UNRCH>
                    warning('File %s does not contain data variable. Skipping.', file_path);
                    continue;
                end
                meta.data_fields = fieldnames(data);
            end

            for field_name = fieldnames(meta)'
                fname = field_name{1};
                metadata_entry(file_idx).(fname) = meta.(fname);
            end

        catch ME
            fprintf("Failed: %s\n", ME.message);
            bad_files{end+1} = {file_name, ME.message}; %#ok<SAGROW> 
            if debug
                rethrow(ME); %#ok<*UNRCH>
            end
        end

    end
    metadata.(data_folder_name) = metadata_entry;
end

current_time = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
metadata_backup_name = sprintf('metadata_backup_%s.mat', current_time);
metadata_backup_path = fullfile(backup_folder, metadata_backup_name);
if isfile(old_metadata_path)
    copyfile(old_metadata_path, metadata_backup_path);
    fprintf('Existing metadata backed up to: %s\n', metadata_backup_path);
else
    fprintf('No existing metadata found. No backup created.\n');
end

save(fullfile(metadata_folder, 'metadata.mat'), 'metadata', '-v7.3');

bad_num = numel(bad_files);
fprintf('Bad file: %d\n', bad_num);
for i = 1:bad_num
    item = bad_files{i};
    fprintf('%d: %s\n%s', i, item{1}, item{2});
end