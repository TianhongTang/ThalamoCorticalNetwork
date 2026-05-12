function data = get(ad)
% ARCHDAT/GET  Retrieve data from an archdat object
% Usage:
%   data = get(ad)
% GET returns the data referred to by an archdat object. If the archdat is
% "ready", it returns the data stored in memory; otherwise, it loads the
% data from file and returns it.
%
% GET is not enabled for non-scalar archdat.

classname = classConst('classname');

if ~isscalar(ad)
    error([classname ':nonScalarGet'], ...
        'Can only get data on scalar archobj''s.');
end

if isnull(ad)
    % Null references contain no data
    data = [];
elseif ad.inmem
    % Return the data from memory
    data = ad.data;
else
    % Load the data from file and return it.
    checkArchive(ad.file,ad.var);
    s = load(ad.file,ad.var);
    data = s.(ad.var);
end
