function pdm = getPDMFileInfo(pdmfile)
% PURPOSE : This function retrieves PDM information.
% USAGE : pdm = getPDMFileInfo([pdmfile])
% VERSION : 1.00  May-2000     YM  modified from getPDMFileInfo.m by DAL
%           1.01  01-Sep-2000  YM  adds 'nprmVars'         
%           1.02  28-Mar-2001  YM  bug fix, use strmatch()

global PDM

% pickup file
if nargin < 1
  pdmfile = pickfile('Load PDM File',pwd,'*.pdm');
  if ~length(pdmfile), return; end
  fprintf('PDMFILE: %s\n',pdmfile);
end

pdmdir  = getFileDirectory(pdmfile);
pdmfile = getFileName(pdmfile);
if ~length(pdmdir), pdmdir = pwd; end

% init output
pdm = [];

% read text
pdmfullpath = sprintf('%s/%s',pdmdir,pdmfile);
pdm.fullpath = pdmfullpath;
pdm.filename = pdmfile;
pdmline = loadTclFile(pdmfullpath);

% process text
prmline = strmatch('newStmParam',pdmline);
nprm = 0;
ncmb = 1;
for i = 1:length(prmline)
  [t0,r0] = strtok(pdmline{prmline(i)});
  %fprintf(' %3d: %s   %s \n',i,pdmline{prmline(i)},t0);
  [t1,r1] = strtok(r0,'"');  [t1,r1] = strtok(r1,'"');
  [t2,r2] = strtok(r1,'"');  [t2,r2] = strtok(r2,'"');
  if length(t2)
    nprm = nprm + 1;
    pdm.prmNames{nprm} = t1;
    pdm.prmVars{nprm} = str2num(t2);
    pdm.nprmVars(nprm) = length(pdm.prmVars{nprm});
    ncmb = ncmb * length(pdm.prmVars{nprm});
  end
end

pdm.nParams = nprm;
pdm.nPattByPrms = ncmb;

if nargout == 0, PDM = pdm; end
