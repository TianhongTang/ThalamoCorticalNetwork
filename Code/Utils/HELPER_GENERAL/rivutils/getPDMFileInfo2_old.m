function pdm = getPDMFileInfo2(pdmfile)
% PURPOSE : This function retrieves PDM information.
% USAGE : pdm = getPDMFileInfo2([pdmfile])
% VERSION : 1.00  May-2000     YM  modified from getPDMFileInfo.m by DAL
%           1.01  01-Sep-2000  YM  adds 'nprmVars'         

global ENV PDM

if ~isa(ENV,'struct'), setEnvVars;  end;

% pickup file
if nargin < 1
  patt = strrep(sprintf('%s/*.%s',ENV.PDMPath,ENV.PDMExt),'/','\');
  [pdmfile, pdmdir]  = uigetfile(patt,'Load PDM File');
  if ~pdmfile, return; end
  pdmfile = strrep(sprintf('%s%s',pdmdir,pdmfile),'\','/');
  fprintf('PDMFILE: %s\n',pdmfile);
end

pdmdir  = getFileDirectory(pdmfile);
pdmfile = getFileName(pdmfile);
if ~length(pdmdir), pdmdir = ENV.PDMPath; end

% init output
pdm = [];

% read text
pdmfullpath = sprintf('%s/%s',pdmdir,pdmfile);
pdm.fullpath = pdmfullpath;
pdm.filename = pdmfile;
fid     = fopen(pdmfullpath,'r');
i = 1;
while 1
  line = fgets(fid,80);
  line = line(1:length(line)-1);
  if ~isstr(line),break,end;
  % remove comments following '#'
  ci = findstr(line,'#');
  if length(ci), line = line(1:ci(1)); end
  pdmline{i} = line;
  i=i+1;
end
fclose(fid);

% process text
nprm = 0;
ncmb = 1;
for i = 1:length(pdmline)
  [t0,r0] = strtok(pdmline{i});	
  if (findstr(t0,'newStmParam'))
    nprm = nprm+1;
    [t1,r1] = strtok(r0,'"');
    [t2,r2] = strtok(r1,'"');
    [t3,r3] = strtok(r2,'"');
    pdm.prmNames{nprm} = t2;
    r4 = r3(2:length(r3));
    r5 = r4(1:findstr(r4,'"')-1);
    tmpi = findstr(r5,' ');
    pdm.prmVars{nprm} = str2num(r5);
    pdm.nprmVars(nprm) = length(pdm.prmVars{nprm});
    ncmb = ncmb * length(pdm.prmVars{nprm});
  end
end

pdm.nParams = nprm;
pdm.nPattByPrms = ncmb;

if nargout == 0, PDM = pdm; end
