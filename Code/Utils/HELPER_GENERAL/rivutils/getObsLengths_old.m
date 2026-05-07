function obslens = getObsLengths(dgz)
% PURPOSE : Gets observation lengths.
% USAGE :   obslens = getObsLengths(dgz)
% NOTE :    'dgz' can be either dgzdata or dgzfile.
% VERSION : 1.00  06-Dec-2000  Yusuke MURAYAMA, MPI

global ENV EVT DGZ DATA

if ~isa(ENV,'struct'), setEnvVars; end;
if ~isa(EVT,'struct'), essevents;  end;

if nargin == 1,
  if isa(dgz,'char')
    dgzfile = dgz;
    dgzdir = getFileDirectory(dgzfile);
    if ~length(dgzdir), dgzdir = ENV.DGZPath;  end
    dgzfile = sprintf('%s/%s.dgz',dgzdir,getFileRoot(dgzfile));
    dgz = dg_read(dgzfile);
  end
elseif nargin == 0
  if class(DGZ) == 'struct', dgz = DGZ;
  elseif class(DATA.dgz) == 'struct', dgz = DATA.dgz;
  else
    error('Neither local or global dgzdata');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get events
tmplens = getEVTTimes(dgz,EVT.ENDOBS);
nobs = length(tmplens);  obslens = zeros(1,nobs);
for i=1:nobs, obslens(i) = tmplens{i};  end
