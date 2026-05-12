function filename = getFileName(fullpath)
% PURPOSE : To extract a filename from fullpath
% USAGE :   filename = getFileName(fullpath)
% VERSION : 1.00 Apr-2000  YM


if nargin < 1
  fprintf('usage: filename = getFileName(fullpath)\n');
  return;
end

filename = fullpath;
token = findstr(strrep(fullpath,'\','/'), '/');
if length(token)
   filename = fullpath(token(length(token))+1:length(fullpath));
end
