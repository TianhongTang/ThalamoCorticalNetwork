function filepath = getFileDirectory(fullpath)
% PURPOSE : To extract directory from fullpath
% USAGE :   filepath = getFileDirectory(fullpath)
% VERSION : 1.00 Apr-2000  YM

if nargin < 1
  fprintf('usage: filepath = getFileDirectory(fullpath)\n');
  return;
end

filepath = '';
% set '\' to '/'
token = findstr(strrep(fullpath,'\','/'), '/');
if length(token)
  filepath = fullpath(1:token(length(token))-1);
end
