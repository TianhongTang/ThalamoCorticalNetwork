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
data_folder     = fullfile(root, 'Data', 'Working');
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');
backup_folder   = fullfile(root, 'Data', 'Working', 'Meta', 'backup');
check_path(metadata_folder);
check_path(backup_folder);
old_metadata_path = fullfile(metadata_folder, 'metadata.mat');
if isfile(old_metadata_path)
    load(old_metadata_path, 'metadata');
    fprintf('Existing metadata loaded from: %s\n', old_metadata_path);
else
    fprintf('No existing metadata found. Starting with empty metadata.\n');
    metadata = struct();
end

data_folders = dir(fullfile(data_folder, '*'));
data_folders = data_folders([data_folders.isdir] & ~startsWith({data_folders.name}, '.'));

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
        m = matfile(file_path);
        vars = who(m);
        if ~ismember('meta', vars)
            warning('File %s does not contain meta variable. Skipping.', file_path);
            continue;
        end
        meta = m.meta;
        meta.data_fields = fieldnames(m.data);

        % compare with current meta if exists
        if isfield(metadata, data_folder_name)
            % try to find matching filename
            existing_metas = metadata.(data_folder_name);
            match_idx = find(arrayfun(@(x) strcmp(x.file_name, meta.file_name), existing_metas), 1);
            if ~isempty(match_idx)
                existing_meta = existing_metas(match_idx);
                if isequal(existing_meta, meta)
                    fprintf('    Meta matches existing entry. No update needed.\n');
                    continue; % skip to next file
                else
                    fprintf('    Meta differs from existing entry. Updating metadata.\n');
                    % update the existing meta entry
                    metadata.(data_folder_name)(match_idx) = meta;
                end
            else
                fprintf('    No matching filename in existing metadata. Adding new entry.\n');
                % add new meta entry
                metadata.(data_folder_name)(end+1) = meta;
            end

        end

        for field_name = fieldnames(meta)'
            fname = field_name{1};
            metadata_entry(file_idx).(fname) = meta.(fname);
        end
    end
    metadata.(data_folder_name) = metadata_entry;
end

current_time = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
metadata_backup_name = sprintf('metadata_backup_%s.mat', current_time);
metadata_backup_path = fullfile(backup_folder, metadata_backup_name);
if isfile()
    copyfile(fullfile(metadata_folder, 'metadata.mat'), metadata_backup_path);
    fprintf('Existing metadata backed up to: %s\n', metadata_backup_path);
else
    fprintf('No existing metadata found. No backup created.\n');
end

save(fullfile(metadata_folder, 'metadata.mat'), 'metadata', '-v7.3');