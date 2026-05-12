function [sort_data, areas, electrodes] = sortByArea(file,data)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% FUNCTION: 
% 	sortByArea(file,data);
%
% ARGUMENTS:
% 	file = name of the data file
%	data = must be sorted with the channels in dim 1
%
% RETURNS:
%	sort_data 	= resorted data, by area
% 	areas 		= areas in resorted order
%	electrodes	= electrodes in resorted order
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file
SI = getSessionInfo(file);
E  = SI.electrodes;
A  = SI.areas;
if isempty(E) | isempty(A)
    warndlg('Returning data without resorting');
    sort_data = data;	
    return;
end

[val,srtindx] 	= sort(A);
if ndims(data) == 3
   sort_data 	= data(:,:,srtindx);
elseif ndims(data) == 2
   sort_data 	= data(:,srtindx);
end
electrodes   	= E(srtindx);
areas	       	= val;



