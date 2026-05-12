function fileroot = getFileRoot(filename)
% PURPOSE : To extract fileroot from fullpath
% USAGE :   fileroot = getFileRoot(filename)
% VERSION : 1.00 Apr-2000  YM

if nargin < 1
  fprintf('usage: fileroot = getFileRoot(filename)\n');
  return;
end

filename = getFileName(filename);
extpts = findstr(filename,'.');
if length(extpts)
   fileroot = filename(1:extpts(length(extpts))-1);
else 
   fileroot = filename;
end
