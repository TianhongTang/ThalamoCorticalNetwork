function tf = isready(ad)
% ARCHDAT/ISREADY  Determine whether an archdat object is "ready"
% Usage:
%   tf = isready(ad)
% ISREADY returns a logical array in which the value of each element
% denotes the readiness state of the corresponding element in ad.
% See also READY.

if isempty(ad)
    tf = logical([]);
elseif isscalar(ad)
    tf = ob.inmem;
else
    tf = reshape([ad.inmem],size(ad));
end