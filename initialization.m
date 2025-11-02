% Setup folder structure
folders = {'Code', 'Data', 'Figures'};
subfolders = {{'Preprocessing', 'Analysis', 'Plotting'}, {'Experimental', 'Working', 'Plotting'}, {}};
for i = 1:length(folders)
    if ~isfolder(folders{i})
        mkdir(folders{i});
    end
    subfolder = subfolders{i};
    for j = 1:length(subfolder)
        if ~isfolder(fullfile(folders{i}, subfolder{j}))
            mkdir(fullfile(folders{i}, subfolder{j}));
        end
    end
end