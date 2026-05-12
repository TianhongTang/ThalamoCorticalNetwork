function [hlist,vlist,emtimes] = getEVTEms(dgzdata)
%
% Read eye movement information from dgz structure
% DAL JAN-2000
%


if (nargin < 1)
	error('usage: [hlist,vlist,sample_rate] = getEVTEms(dgzdata)')
end

hscales  = getEVTParams(dgzdata,26,1,0);
vscales  = getEVTParams(dgzdata,26,2,0);

ems     = dgzdata.ems;
for j=1:length(ems)
   if length(ems{j})
	   hlist{j} = ems{j}{2}/hscales{j};	
 	   vlist{j} = ems{j}{3}/vscales{j};	
	   emtimes{j} = (1:length(vlist{j}))*ems{j}{1};
   end
end

