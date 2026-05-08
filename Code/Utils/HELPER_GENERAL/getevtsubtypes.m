function retlist = getEVTSubtypes(dgzdata, evtype)
%
% Read event subtypes from a dgz file
% DAL JAN-2000
% YM  05-Apr-2001  modified to use find()

if (nargin < 2)
  error('usage: evtsubtypes = getEVTSubtypes(dgzdata,evtype)')
end

%
%  First get the main parameter field from the cell array
types    = getfield(dgzdata, 'e_types');
subtypes = getfield(dgzdata, 'e_subtypes');

retlist = {};
%select based only on evtype
for j=1:length(types)
  sel = find(dgzdata.e_types{j} == evtype);
  retlist{j} = dgzdata.e_subtypes{j}(sel);
end

