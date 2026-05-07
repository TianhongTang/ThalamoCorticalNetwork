function [idx,idxDist,thisMin,clusMembers, linkMat, linkInds] = distMatAgglom(distMat,k,varargin)

nPrint = 100;
print = 1;
minN = 10;

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'min','-n'}
            minN = varargin{varStrInd(iv)+1};
        case {'-p'}
            print = varargin{varStrInd(iv)+1};
    end
end

if ~isnan(distMat(1,2)) && ~isnan(distMat(2,1))
    for ic = 1:size(distMat,2)
        for ir = ic:size(distMat,2)
            distMat(ir,ic) = nan;
        end
    end
end

% Save distance and observations in raw form in case I manipulate them in
% the next step...
distRaw = distMat;

% Initialize clusMembers
clusMembers = num2cell((1:size(distRaw,1))');
linkMat = [];
linkIDs(1:size(distRaw,1),1) = 1:size(distRaw,1);
nLoops = 0;
thisMin = nan(1,size(distMat,1));
while size(distMat,1) > 1
    if mod(size(distMat,1),nPrint) == 0 && print
        fprintf('%d members remaining...\n',size(distMat,1));
    end
    nLoops = nLoops+1;
    % Get the row and column corresponding to the minimum distance
    % remaining. If there are multiple, take just the first index
    % because other indices should be captured in the next loop
    % through.

    [minRow,minCol] = find(distMat == min(distMat(:)),1);
    if isnan(min(distMat(:))), break; end
    thisMin(nLoops+1) = min(distMat(:));
    
    % Get the indices of clusMembers to keep
    nonClust = 1:size(distMat,1); nonClust(ismember(nonClust,[minRow,minCol])) = [];
    newMems = clusMembers(nonClust,end);
    newMems{end+1} = [clusMembers{minRow,end},clusMembers{minCol,end}];
    theseMems = newMems{end};

    % Let's make a linkage matrix. Remember:
    % For example, suppose there are 30 initial nodes and at 
    % step 12 cluster 5 and cluster 7 are combined. Suppose their 
    % distance at that time is 1.5. Then Z(12,:) will be [5, 7, 1.5]. 
    % The newly formed cluster will have index 12 + 30 = 42. If 
    % cluster 42 appears in a later row, it means the cluster created 
    % at step 12 is being combined into some larger cluster.
    linkMat = cat(1,linkMat,[linkIDs(minRow),linkIDs(minCol),min(distMat(:))]);
    newLink = nLoops + size(distRaw,1);
    linkIDs = linkIDs(nonClust);
    linkIDs(end+1) = newLink;

    % Debugging to be sure we're not overcounting, should be worked
    % out...
    if any(ismember(clusMembers{minRow,end},clusMembers{minCol,end}))
        keyboard
    end

    % Save combined members in clusMembers
    clusMembers(1:length(newMems),size(clusMembers,2)+1) = newMems;

    % Cut the cluster members from distMat
    distMat = distMat(nonClust,nonClust);

    % Calculate distance between other clusters and the most recent
    % combination
    distMat = cat(1,distMat,nan(1,size(distMat,2)));
    distMat = cat(2,distMat,nan(size(distMat,1),1));
    for im = 1:(size(distMat,1)-1)
        tstMems = clusMembers{im,end};
        allCorrs = distRaw(tstMems,theseMems);
        allCorrsTrans = distRaw(theseMems,tstMems)';
        noNanCorrs = nan(size(allCorrs));
        noNanCorrs(~isnan(allCorrs)) = allCorrs(~isnan(allCorrs));
        noNanCorrs(~isnan(allCorrsTrans)) = allCorrsTrans(~isnan(allCorrsTrans));

        thisDist = nanmean(noNanCorrs(:));
        distMat(end,im) = thisDist;
    end    
end

% Get IDX from clustMembs
idx = nan(size(distRaw,1),max(k));
clustLens = cellfun(@length,clusMembers);
nkMeetsCrit = sum(clustLens >= minN,1);
sortK = sort(k);
idxDist = nan(1,max(k));
linkInds = nan(1,length(sortK));
for ik = 1:length(sortK)
    thisK = sortK(ik);
    thisCol = find(nkMeetsCrit == thisK,1,'last');
    if isempty(thisCol)
        fprintf('Criteria not met for %d clusters... Exiting loop\n',thisK);
        break;
    end
    linkInds(ik) = thisCol;
    idxDist(thisK) = thisMin(thisCol);
    clustRows = find(clustLens(:,thisCol) >= minN);
    for ic = 1:length(clustRows)
        idx(clusMembers{clustRows(ic),thisCol},thisK) = ic;
    end
end


    
    