function fileext = getFileExt(filename)
% PURPOSE : To extract a file-extention from filename
% USAGE :   fileext = getFileExt(filename)
% VERSION : 1.00 Apr-2000  YM

if nargin < 1
  fprintf('usage: fileext = getFileExt(filename)\n');
end


fileext = '';
filename = getFileName(filename);
extpts = findstr(filename,'.');
if (length(extpts))
  % include '.'
  %fileext = filename(extpts(length(extpts)):length(filename));
  % not include '.'
  fileext = filename(extpts(length(extpts))+1:length(filename));
end
