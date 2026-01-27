%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end

%% Main
folder = fullfile(root, 'Data', 'Working', 'border');
files = dir(fullfile(folder, '*.mat'));
% print file names
for i = 1:length(files)
    % match format: borders_Slayer*PreCortex.mat
    if startsWith(files(i).name, 'borders_Slayer') && contains(files(i).name, 'Pre')
        fprintf('%s\n', files(i).name);
        load(fullfile(folder, files(i).name), 'borders');
        new_filename = strrep(files(i).name, 'Pre', ''); % remove 'Pre' from filename
        save(fullfile(folder, new_filename), 'borders');
    end
    
end