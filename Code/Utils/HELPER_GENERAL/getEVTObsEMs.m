function [hlist,vlist,emtimes] = getEVTObsEms(dgzdata,obs)
%
% Read one obsper's eye movement information from dgz structure
% DAL JAN-2000
%

if (nargin < 1)
	error('usage: [hlist,vlist,sample_rate] = getEVTObsEms(dgzdata,obs)')
end

if obs > length(dgzdata.ems)
   errordlg('Obs index (%d) exceeds total obspers (%d)\n', obs,length(ems))
end

hscales  = getEVTParams(dgzdata,26,1,0);
vscales  = getEVTParams(dgzdata,26,2,0);

ems     = dgzdata.ems;
if length(ems{obs})
   hlist   = ems{obs}{2}/hscales{obs};	
   vlist   = ems{obs}{3}/vscales{obs};	
   emtimes = ((1:length(vlist))*ems{obs}{1})';
end

