function infoSt = info(ad)
% ARCHDAT/INFO  Get information on archdat objects
% 
% Usage:
%   infoSt = info(ad)
% infoSt is a structure array of the same size as the archdat array ad,
% with the following fields:
%    file: The filename in which the archive is stored
%     var: The variable name of the data in the archive
%   class: The data type of the archived data
%    size: The array size of the archived data
%   ready: True/false, whether this archdat is currently "ready"
%    null: True/false, whether this is a null archdat or not.
% The info structure for a null archdat object has empty values for the
% fields "file", "var", "class", and "size", the value false for "ready",
% and the value "true" for "null".
%
% Performance:
%   The information provided by INFO can also be retrieved with the methods
% FILENAME, VARNAME, DATACLASS, DATASIZE, ISREADY, and ISNULL. If the
% archdat is not "ready", DATACLASS and DATASIZE must access the archive on
% disk to retrieve the information. If both of these pieces of information
% are desired, it is faster to call INFO once rather than DATACLASS and
% DATASIZE separately. If info is called as in the following:
%   s = info(ad);
% then the following usages are equivalent for non-scalar ad:
%                  is equivalent to
%   filename(ad) ................... reshape({s.file},size(s))
%   varname(ad) .................... reshape({s.var},size(s))
%   dataclass(ad) .................. reshape({s.class},size(s))
%   datasize(ad) ................... reshape({s.size},size(s))
%   isready(ad) .................... reshape([s.ready],size(s))
%   isnull(ad) ..................... reshape([s.null],size(s))

% Define the struct that is returned if this is a null archdat:
nullInfo = struct('file',  '', ...
                  'var',   '', ...
                  'class', '', ...
                  'size',  [], ...
                  'ready', false, ...
                  'null',  true);

if isempty(ad)
    % If ad is an empty array, return an empty struct array:
    infoSt = struct('file',  {}, ...
                    'var',   {}, ...
                    'class', {}, ...
                    'size',  {}, ...
                    'ready', {}, ...
                    'null',  {});
elseif isscalar(ad)
    % If ad is scalar, our job is relatively easy:
    if isempty(ad.file)
        % Empty file string means this is a null archdat:
        infoSt = nullInfo;
    else
        if ad.inmem
            % The data is "ready" in memory, so just get its class and size
            dclass = class(ad.data);
            dsize  = size(ad.data);
        else
            % The data is not "ready", so get the class and size from file
            checkArchive(ad.file,ad.var);
            s = whos(ad.var,'-file',ad.file);
            dclass = s.class;
            dsize  = s.size;
        end
        infoSt.file  = ad.file;
        infoSt.var   = ad.var;
        infoSt.class = dclass;
        infoSt.size  = dsize;
        infoSt.ready = ad.inmem;
        infoSt.null  = false;
    end
else
    % This is a non-scalar archdat array, so we must get info for each
    % element
    infoSt = repmat(nullInfo,size(ad));
    dclass = cell(size(ad));
    dsize  = cell(size(ad));
    is_null = cellfun(@isempty,{ad.file});
    is_ready = [ad.inmem];
    % Get the class and size for data already in memory:
    dclass(is_ready) = cellfun(@class, {ad(is_ready).data}, ...
        'UniformOutput', false);
    dsize(is_ready)  = cellfun(@size, {ad(is_ready).data}, ...
        'UniformOutput', false);
    % Get the indices of the elements that are currently archived:
    archInds = find(~is_ready&~is_null);
    while ~isempty(archInds)
        % Determine which elements are stored in the same file:
        filename = ad(archInds(1)).file;
        sameFile = strcmp(filename,{ad(archInds).file});
        varList = {ad(archInds(sameFile)).var};
        checkArchive(filename,varList);
        % Get information on the archived variables in this file:
        s = whos(varList{:}, '-file', filename);
        fileVars = {s.name};
        % Reorder the results of WHOS according to varList:
        varOrder = cellfun(@(x)find(strcmp(x,fileVars)),varList);
        s = s(varOrder);
        dclass(archInds(sameFile)) = {s.class};
        dsize(archInds(sameFile)) = {s.size};
        archInds = archInds(~sameFile);
    end
    [infoSt(~is_null).file]  = deal(ad(~is_null).file);
    [infoSt(~is_null).var]   = deal(ad(~is_null).var);
    [infoSt(~is_null).class] = deal(dclass{~is_null});
    [infoSt(~is_null).size]  = deal(dsize{~is_null});
    [infoSt(is_ready).ready] = deal(true);
    [infoSt(~is_null).null]  = deal(false);
end