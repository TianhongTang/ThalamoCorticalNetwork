function files_out = dir_mat_files(dir_in)

Files = dir(fullfile(dir_in, '*.mat'));
files_out = {Files.name}';

end