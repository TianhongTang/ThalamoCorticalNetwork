function retlist = getEVTObsTime(dgzdata, obs, evtype, evsubtype)
%
% 16-Aug-2000  YM modified from DAL version, using find()

if (nargin < 3)
  error('usage: evttimes = getEVTObsTime(dgzdata,obs,evtype,[evsubtype])')
end
if (nargin < 4),  evsubtype = -1;  end	


retlist = [];
if evsubtype == -1 
  %select based only on evtype
  sel = find(dgzdata.e_types{obs} == evtype);
else
  %select based on evtype and evsubtype
  %select based only on evtype
  sel = find(dgzdata.e_types{obs} == evtype &...
	     dgzdata.e_subtypes{obs} == evsubtype);
end

if length(sel) ~= 0, retlist = dgzdata.e_times{obs}(sel);  end
