function retlist = getEVTTimes(dgzdata, evtype, evsubtype)
%
% Read event times from a dgz file
% DAL JAN-2000
% YM  05-Apr-2001  modified to use find()

if (nargin < 2)
  error('usage: evttimes = getEVTTimes(dgzdata,evtype,[evsubtype])')
end
if (nargin < 3),  evsubtype = -1;  end	

% select events
tmpsel = {};
if evsubtype == -1 
  %select based only on evtype
  for j=1:length(dgzdata.e_types)
    tmpsel{j} = find(dgzdata.e_types{j} == evtype);
  end	
else
  %select based on evtype and evsubtype
  for j=1:length(dgzdata.e_types)
    tmpsel{j} = find(dgzdata.e_types{j} == evtype & ...
		     dgzdata.e_subtypes{j} == evsubtype);
  end
end

% get times
retlist = {};
for j=1:length(tmpsel)
  retlist{j} = dgzdata.e_times{j}(tmpsel{j});
end

