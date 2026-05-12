function tf = isnull(ad)
% ARCHDAT/ISNULL  Determine whether an archdat object is a null reference
% Usage:
%   tf = isnull(ad)
% ISNULL returns a logical array of the same size as ad, in which the value
% of each element denotes whether the corresponding element in ad is a null
% reference.

fname = filename(ad);
if isempty(fname)
    tf = logical([]);
elseif ~iscell(fname)
    tf = isempty(fname);
else
    % Call isempty on every element of fname and store the result in tf
    tf = cellfun(@isempty,fname);
end