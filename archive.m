%% archive.m - Archive old slurm out files

current_folder = pwd;
files = dir(fullfile(current_folder, 'slurm-*.out'));
archive_folder = fullfile(current_folder, 'archived');

if ~exist(archive_folder, 'dir')
    mkdir(archive_folder);
end

for i = 1:length(files)
    file_name = files(i).name;
    % if older than 1 days
    file_date = datetime(files(i).date);
    if (datetime('now') - file_date) > duration(24, 0, 0)
        movefile(fullfile(current_folder, file_name), fullfile(archive_folder, file_name));
        fprintf('Archived file: %s\n', file_name);
    else
        fprintf('Kept file: %s\n', file_name);
    end
end