function file = filename(ad)
% ARCHDAT/FILENAME  Retrieve the file name in which archdat data is stored
% Usage:
%   file = filename(ad)
% For scalar archdat, file is a string containing the name of the file in
% which the archdat data is stored. For non-scalar archdat, file is a cell
% array of the same size as ad, in which each element is the filename of
% the corresponding element in ad.

if isempty(ad)
    file = '';
elseif isscalar(ad)
    file = ad.file;
else
    file = reshape({ad.file},size(ad));
end