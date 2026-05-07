function filenames = getDirFilenames(storageDir,label)

matFiles = dir(fullfile(storageDir,'*.mat'));
matFiles = {matFiles.name};
valid_files = ~cellfun('isempty',strfind(matFiles,label));
filenames = matFiles(valid_files);

end
