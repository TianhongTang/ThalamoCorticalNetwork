function sz = datasize(ob)
% ARCHDAT/DATASIZE  Get the size of archdat data
% Usage:
%   sz = datasize(ad)
% DATACLASS returns a vector specifying the size of the data referred to by
% the archdat object. For scalar archdat, it is essentially equivalent to
% calling SIZE on the data (e.g., "size(get(ad))"), but does not need to
% load the data from the archive to get the size information. For
% non-scalar archdat, DATASIZE returns a cell array of the same size as ad,
% with each element containing a string specifying the size of the data
% from the corresponding element of the archdat array.
%
% Performance:
%   If both class and size information are desired, it is faster to call
% INFO once rather than both DATACLASS and DATASIZE. See INFO.

nullSize = [0 0];

if isempty(ob)
    sz = nullSize;
elseif isscalar(ob)
    if isnull(ob)
        sz = nullSize;
    elseif ob.inmem
        sz = size(ob.data);
    else
        s = whos(ob.var,'-file',ob.file);
        sz = s.size;
    end
else
    sz = cell(size(ob));
    sz(isnull(ob)) = {nullSize};
    sz([ob.inmem]) = cellfun(@size,{ob([ob.inmem]).data}, ...
        'UniformOutput',false);
    archInds = find(cellfun(@isempty,sz));
    while ~isempty(archInds)
        filename = ob(archInds(1)).file;
        sameFile = strcmp(filename,{ob(archInds).file});
        varList = {ob(archInds(sameFile)).var};
        checkArchive(filename,varList);
        s = whos(varList{:}, '-file', filename);
        fileVars = {s.name};
        varOrder = cellfun(@(x)find(strcmp(x,fileVars)),varList);
        s = s(varOrder);
        sz(archInds(sameFile)) = {s.size};
        archInds = archInds(~sameFile);
    end
end
