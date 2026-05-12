function obsevt = getObsEvt(dgz,pdm)
% PURPOSE : Picks up interesting events in observations
% USAGE :   obsevt = getObsEvt(dgz,[pdm])
% NOTE :    Use dg_read.dll modified by YM or dg_read2.m
%           Some float value may have a small offset like 1.2e-8.
%           So, I use 'round(x*10000)/10000' as you see below.
%
% ver 1.00  17-May-2000  Yusuke MURAYAMA, MPI
% ver 1.02  13-Oct-2000  Yusuke MURAYAMA, MPI

global ENV VERBOSE DATA DGZ EVT OBS PDM

if ~isa(ENV,'struct'), setEnvVars; end;
if ~isa(EVT,'struct'), essevents;  end;

if nargin < 1
  if class(DGZ) == 'struct', dgz = DGZ;
  elseif class(DATA.dgz) == 'struct', dgz = DATA.dgz;
  else
    error('Neither local or global dgzdata');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get names of stimulus parameter
if exist('pdm')
  stmPrmNames = pdm.prmNames;
elseif class(PDM) == 'struct'
  stmPrmNames = PDM.prmNames;
else
  fprintf(' getObsEvt: getting parameter names from dgzdata...\n');
  for i=1:length(dgz.e_params{1})
    if dgz.e_types{1}(i) == EVT.STRINGS_1
      tmpPrmNames = dgz.e_params{1}{i};
      break;
    end
  end
  % verify names to avoid duplicated name
  n=size(tmpPrmNames,1);
  rm_indx = zeros(1,n);
  for i=1:n
    for j=(i+1):n
      if strcmp(tmpPrmNames(i,:),tmpPrmNames(j,:)), rm_indx(j) = 1; end
    end
  end
  k = 0;
  for i=1:n
    if rm_indx(i) ~= 1
      k = k + 1;
      stmPrmNames{k} = deblank(tmpPrmNames(i,:));
    end
  end
end
obsevt.stmPrmNames = stmPrmNames;
obsevt.stmNPrm     = length(stmPrmNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get events
nobs = length(dgz.e_types);
obsevt.nobs = nobs;
for i=1:nobs
  % timings
  obsevt.obslen(i)   = getEVTObsTime(dgz,i,EVT.ENDOBS);
  obsevt.patofft{i}  = getEVTObsTime(dgz,i,EVT.PATTERN,0);
  obsevt.patont{i}   = getEVTObsTime(dgz,i,EVT.PATTERN,1);
  % stimulus parametes
  for j=1:length(stmPrmNames)
    p = getEVTObsParam(dgz,i,EVT.FLOATS_1,j);
    % fix bugs, some float value may have a small offset like
    % 1.2e-8 maybe,due to float->double conversion.
    p = round(p*10000)/10000;
    if j==1, tmpprms = zeros(length(p),length(stmPrmNames)); end
    tmpprms(:,j) = p;
    % following thing is to be compatibile with old version
    %str_field = sprintf('obsevt.%s{i}',stmPrmNames{j});
    str_field = sprintf('obsevt.prm%d{i}',j);
    eval(sprintf('%s = p;',str_field));
  end
  obsevt.prms{i}     = tmpprms;
  obsevt.riv{i}      = getEVTObsParam(dgz,i,EVT.PATTERN,1,1);
  obsevt.fade{i}     = getEVTObsParam(dgz,i,EVT.PATTERN,2,1);
  obsevt.trial{i}    = getEVTObsParam(dgz,i,EVT.PATTERN,3,1);
  obsevt.stimulus{i} = getEVTObsParam(dgz,i,EVT.PATTERN,4,1);
  obsevt.eye{i}      = getEVTObsParam(dgz,i,EVT.PATTERN,5,1);
  obsevt.duration{i} = getEVTObsParam(dgz,i,EVT.PATTERN,7,1);
  % before Aug-2000
  %obsevt.rivcfg{i}   = getEVTObsParam(dgz,i,EVT.FLOATS_1,3);
  % after Aug-2000
  %obsevt.rivcfg{i}   = getEVTObsParam(dgz,i,EVT.FLOATS_1,2);
end

obsevt.dgzfile = dgz.filename;

if nargout == 0, 
  OBS = catStruct(OBS,obsevt);
end
