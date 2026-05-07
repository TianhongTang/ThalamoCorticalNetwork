function retlist = getEVTObsSubtype(dgzdata, obs, evtype)
%
% Read event subtypes from a dgz file
% DAL JAN-2000
% 16-Aug-2000  YM modified from DAL version, using find()

if (nargin < 3)
  error('usage: evtobssubtype = getEVTObsSubtype(dgzdata,obs,evtype)')
end


retlist = [];
%select based only on evtype
sel = find(dgzdata.e_types{obs} == evtype);
if length(sel) ~= 0
  retlist = dgzdata.e_subtypes{obs}(sel);
end

