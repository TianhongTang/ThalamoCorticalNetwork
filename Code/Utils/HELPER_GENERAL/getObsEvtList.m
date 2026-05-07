function evtList = getObsEvtList(dgzFileList,pdm)
% PURPOSE : To make a list of event from a number of DGZ
% USAGE : evtObsList = getObsEvtList(dgzFileList,[pdm])
% VERSION : 1.00  14-Aug-2000  YM
%           1.01  07-Sep-2000  YM
  
global ENV

evtList = [];
if nargin < 1
  fprintf('Usage: evtObsList = getObsEvtList(dgzFileList,[pdm])\n');
  return;
end

if ~isa(ENV,'struct'), setEnvVars;  end

if ~isa(dgzFileList,'cell'),
  tmp = dgzFileList; dgzFileList = {};  dgzFileList{1} = tmp;
end

for i=1:length(dgzFileList)
  dgzdir = getFileDirectory(dgzFileList{i});
  dgzfile = getFileName(dgzFileList{i});
  if ~length(dgzdir), dgzdir = ENV.DGZPath;  end
  dgzfullpath = sprintf('%s/%s',dgzdir,dgzfile);
  dgz = dg_read(dgzfullpath);
  if ~exist('pdm')
    evtList{i} = getObsEvt(dgz);
  else
    evtList{i} = getObsEvt(dgz,pdm);
  end
end
