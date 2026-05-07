function klDendroClustChange(h,linkMat,idx,varargin)

% Set defaults

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-c','color'}
            colors = varargin{varStrInd(iv)+1};
    end
end


% Start by identifying the unique clusters
uClusts = unique(idx(~isnan(idx)));

% Set up handy constants
n = length(idx);
if ~exist('colors','var')
    colors = jet(length(uClusts));
end

% Loop through clusters to get an initial accounting of the idx indices
clustNums = cell(1,length(uClusts)+1);
for ic = 1:length(uClusts)
    clustNums{ic} = find(idx==uClusts(ic));
end
clustNums{ic+1} = find(isnan(idx));

% What we'll do is we'll set up a variable that tells us which cluster the
% stuff was in
hClust = nan(1,size(linkMat,1));

% Now, this may be time consuming but we'll have to loop through linkMat,
% at least until I think of a better solution...
for il = 1:size(linkMat,1)
    % Figure out which cluster these observations are in
    col1 = find(cellfun(@(x) any(ismember(x,linkMat(il,1))),clustNums));
    col2 = find(cellfun(@(x) any(ismember(x,linkMat(il,2))),clustNums));
    
    if col1 ~= col2
        break; %continue
    end
    
    hClust(il) = col1;
    clustNums{col1} = cat(1,clustNums{col1},il+n);
end
hClust(il:end) = length(clustNums);

% Now, we'll loop back through clusters and change the colors
for ic = 1:length(uClusts)
    set(h(hClust==uClusts(ic)),'color',colors(ic,:));
end
set(h(hClust==ic+1),'color','k');

% 
% 
% 
% 
% 
% for ic = 1:length(uClusts),
%     % Get this cluster index of idx
%     theseInds = find(idx==uClusts(ic));
%     
%     % Get the indices of linkMat corresponding to this cluster
%     hClustInds = find(ismember(linkMat(:,1),theseInds) | ismember(linkMat(:,2),theseInds));
%     
%     set(h(hClustInds),'color',colors(ic,:));
%     allHClusts = cat(2,allHClusts,hClustInds');
%     allHInds   = cat(2,allHInds,theseInds');
%     whichClust = cat(2,whichClust,ones(1,length(hClustInds)).*ic);
% end
% theseInds = find(isnan(idx));
% hClustInds = find(ismember(linkMat(:,1),theseInds) | ismember(linkMat(:,2),theseInds));
% set(h(hClustInds),'color','k');
% allHClusts = cat(2,allHClusts,hClustInds');
% allHInds   = cat(2,allHInds,theseInds');
% whichClust = cat(2,whichClust,ones(1,length(hClustInds)).*ic+1);
% 
% keyboard