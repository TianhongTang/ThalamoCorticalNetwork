function dc = dataclass(ob)
% ARCHDAT/DATACLASS  Get the class of archdat data
% Usage:
%   class = dataclass(ad)
% DATACLASS returns a string specifying the class of the archived data
% referred to by the archdat object. For scalar archdat, it is essentially
% equivalent to calling CLASS on the data (e.g., "class(get(ad))"), but
% does not need to load the data from the archive to get the class
% information. For non-scalar archdat, DATACLASS returns a cell array of
% the same size as ad, with each element containing a string specifying the
% class of the data from the corresponding element of the archdat array.
%
% Performance:
%   If both class and size information are desired, it is faster to call
% INFO once rather than both DATACLASS and DATASIZE. See INFO.

nullClass = '';

if isempty(ob)
    dc = nullClass;
elseif isscalar(ob)
    if isnull(ob)
        dc = nullClass;
    elseif ob.inmem
        dc = class(ob.data);
    else
        s = whos(ob.var,'-file',ob.file);
        dc = s.class;
    end
else
    dc = cell(size(ob));
    dc(isnull(ob)) = {nullClass};
    dc([ob.inmem]) = cellfun(@class,{ob([ob.inmem]).data}, ...
        'UniformOutput',false);
    archInds = find(cellfun(@isempty,dc));
    while ~isempty(archInds)
        filename = ob(archInds(1)).file;
        sameFile = strcmp(filename,{ob(archInds).file});
        varList = {ob(archInds(sameFile)).var};
        checkArchive(filename,varList);
        s = whos(varList{:}, '-file', filename);
        
        fileVars = {s.name};
        varOrder = cellfun(@(x)find(strcmp(x,fileVars)),varList);
        s = s(varOrder);
        dc(archInds(sameFile)) = {s.class};
        archInds = archInds(~sameFile);
    end
end