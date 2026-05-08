function dgzdata = loadDGZFiles(dir, filelist)

tmp = size(filelist);
dgzdata = [];
for i=1:tmp(1)
   file = filelist{i};
   fprintf('Loading dgz file: %s.\n', file);
   fullpath 		= sprintf('%s%s',dir,file);
   tmpdgzdata = dg_read(fullpath);
   dgzdata = dgzAppend (dgzdata,tmpdgzdata);
end
