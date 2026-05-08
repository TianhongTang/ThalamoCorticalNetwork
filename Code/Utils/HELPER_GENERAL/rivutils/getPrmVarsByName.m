function [vars,vidx] = getPrmVarsByName(prmName,pdm)
% PURPOSE : To get parameter values from its name
% USAGE : [vars,vidx] = getPrmVarsByName(prmName,[pdm])
% VERSION : 1.00  10-Aug-2000  YM
  
global ENV FILES PDM

if nargin < 1
  fprintf('USAGE: [vars,vidx] = getPrmVarsByName(prmName,[pdm])\n');
  return;
end
if ~isa(ENV,'struct'), setEnvVars;  end;

if ~exist('pdm')
  if isa(PDM,'struct')
    pdm = PDM;
  else
    if isa(FILES,'struct')
      pdm = getPDMFileInfo(FILES.PDMFullPath);
    else
      pdm = getPDMFileInfo;
    end
  end
end

% main body of this function
vars = [];
vidx = 0;
for p=1:pdm.nParams
  if strcmp(pdm.prmNames{p},prmName)
    vidx = p;
    break;
  end
end

if vidx ~= 0, vars = pdm.prmVars{vidx};  end;
