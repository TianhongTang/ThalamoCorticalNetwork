function vname = varname(ad)
% ARCHDAT/VARNAME  Retrieve the variable in which archdat data is stored
% Usage:
%   var = varname(ad)
% For scalar archdat, var is a string containing the name of the variable
% in which the archdat data is stored. For non-scalar archdat, var is a
% cell array of the same size as ad, in which each element is the variable
% name of the corresponding element in ad.

if isempty(ad)
    vname = '';
elseif isscalar(ad)
    vname = ad.var;
else
    vname = reshape({ad.var},size(ad));
end