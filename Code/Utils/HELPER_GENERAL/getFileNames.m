function filenames = getFileNames(fullpath)
  
filenames = {};
tmp = dir(fullpath);
for i=1:length(tmp)  
  filenames{i} = tmp(i).name;
end

