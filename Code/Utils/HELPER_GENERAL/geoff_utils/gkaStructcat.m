function s = gkaStructcat(varargin)

allStruct = all(cellfun(@isstruct,varargin));
if ~allStruct
    error('gkaStructcat:nonStruct','All inputs must be structs.');
end

dims = cellfun(@ndims,varargin);
if any(dims-dims(1))
    error('gkaStructcat:sizeMismatch','All inputs must be the same size.');
end
szs = cellfun(@size,varargin,'UniformOutput',false);
szs = vertcat(szs{:});
sz = szs(1,:);
sztest = bsxfun(@minus,szs,sz);
if any(sztest(:))
    error('gkaStructcat:sizeMismatch','All inputs must be the same size.');
end

allNames = cellfun(@fieldnames,varargin,'UniformOutput',false);
allNames = vertcat(allNames{:});
if numel(allNames)~=numel(unique(allNames))
    error('gkaStructcat:fieldRepeat', ...
        'All inputs must have different fields.');
end

%nEntries = numel(allNames);

allVals = cellfun(@struct2cell,varargin,'UniformOutput',false);
allVals = vertcat(allVals{:});

s = cell2struct(allVals,allNames,1);
