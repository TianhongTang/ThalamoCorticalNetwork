function obslens = getObsLengths(dgz)
% PURPOSE : To get time length of all observation from dgz/adf.
% USAGE :   obslens = getObsLengths([dgz/dgzfile/adffile])
% NOTE :    argument can be either dgzdata or dgzfile/adffile.
% VERSION : 1.00  06-Dec-2000  Yusuke MURAYAMA, MPI
%         : 1.01  10-Oct-2002  Yusuke MURAYAMA, MPI

global EVT

if ~isa(EVT,'struct'), events;  end;

obslens = [];

if nargin == 0,
  dgz = pickfile('Select a dgz/adf file',pwd,'*.dgz;*.adf;*.adfw');
  if isempty(dgz),  return;  end
end

if isa(dgz,'char')
  if ~isempty(findstr(dgz,'.dgz')),
    % read the event file
    dgz = dg_read(dgz);
  else
    % assuming adf/adfw file, and get info
    [nchan,nobs,sampt,obslens] = adf_info(dgz);
    obslens = obslens*sampt;
    return;
  end
end

% get events from the dgzfile
tmplens = getEVTTimes(dgz,EVT.ENDOBS);
nobs = length(tmplens);  obslens = zeros(nobs,1);
for i=1:nobs, obslens(i) = tmplens{i};  end
