%% klDoMetaClusteringv2 clusters SDFs according to desired conditions and analysis windows
%
%  Inputs
%     goodSDF: a cell array that contains SDFs. Each element of
%       the cell is a matrix that share a common alignment point (e.g.,
%       aligned on stimulus onset) and in the same condition (e.g., in RF).
%       Each row corresponds to a different unit. goodSDF can have as many 
%       elements as conditions desired for analysis.
%
%     allTimeCell: a cell array that contains time vectors that correspond
%       to the elements of goodSDF. For example, if goodSDF{1} is a matrix
%       containing SDFs from 200 ms prior to stimulus onset to 300 ms after
%       stimulus onset, allTimeCell{1} = -200:300.
%
%       Important/necessary options   
%           -e: A cell array of times that are used to define "epochs" within
%               the SDFs for different measurements (myEpocs is the
%               variable name)
%           -ei: Index into goodSDF/allTimeCell for which each epoch should be
%               selected. For example, if two epochs are used and myEpocInds = [1,2],
%               then the first epoch will be taken from goodSDF{1} and the
%               second will be taken from goodSDF{2}. Length of myEpocInds
%               and myEpocs must be identical.
%           -er: myEpocWinds is used to determine the overall time window that
%               should be used for calculating whole trial mean/SD for scaling
%
%
%       Tweaking options:
%           -k: Which range of K (number of clusters) should be tested for
%               in individual analysis pipelines
%           -n: The minimum number of units required per group in the
%               individual analysis pipelines
%           -mn: The minimum number of units required per group in the
%               consensus clustering
%           -r: The number of randomized bootstrapped iterations to calculate
%               Tibshirani gap
%           -in: If the individual pipeline structure has been calculated
%               elsewhere or has been saved, it can be passed in here
%           -u: "Uniform K" forces a given K for the individual pipelines
%           -m: Can be 'median' or 'mean' which refers to whether the
%               consensus distance matrix is formed by the mean or median across
%               individual pipelines
%           -sd: Which elements of goodSDF should be used to calculate
%               standard deviations for scaling. Defaults to all elements
%           -do: Distance only. This flag skips the actual clustering for each
%               individual pipeline and simply returns the appropriate distance
%               matrices. Defaults to 1 for time; individual pipeline
%               clusters are still supported but come from an earlier version
%               of the algorithm which compared individual results.
%           -c: Unclustered Criterion - the maximum percentage of units
%               that are acceptably unclustered. Used to determine myK.
%
%  Outputs
%     sortIDs: A matrix of IDs where column j corresponds to the results
%       for k = j
%     idxDist: A vector for which element i corresponds to the distance at
%       which clusters with k = i were formed.
%     raw: The composite distance matrix where the diagonal and lower
%       triangle are NaN.
%     respSumStruct: Summary structure for each individual pipeline. Useful
%       for accessing individual distance matrices or individual clustering
%       results if -do = 0.
%     rawLink: Linkage matrix for the agglomerative clustering procedure
%       that can be used to create dendrograms of the consensus clustering.
%     myK: Number of identified clusters as set by minimum membership and
%       maximum unclustered criteria



function [sortIDs, idxDist, raw, respSumStruct,rawLink, myK] = consensusCluster(goodSDF,allTimeCell,varargin)

%% Set options and defaults
% restrict = 0;
setK = 5;
minN = 15;
metaN = 10;
randReps = 5;
meanMedian = 'median';
distOnly = 1;
unClustCrit = .1;
normSD = 1:length(goodSDF);

%% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'-k','k'}
            setK = varargin{varStrInd(iv)+1};
        case {'-n'}
            minN = varargin{varStrInd(iv)+1};
        case {'-mn'}
            metaN = varargin{varStrInd(iv)+1};
        case {'-r'}
            randReps = varargin{varStrInd(iv)+1};
        case {'-in'}
            respSumStruct = varargin{varStrInd(iv)+1};
        case {'-m'}
            meanMedian = varargin{varStrInd(iv)+1};
        case {'-e'}
            myEpocs = varargin{varStrInd(iv)+1};
        case {'-ei'}
            myEpocInds = varargin{varStrInd(iv)+1};
        case {'-er'}
            myEpocWinds = varargin{varStrInd(iv)+1};
        case {'-sd'}
            normSD = varargin{varStrInd(iv)+1};
        case {'-do'}
            distOnly = varargin{varStrInd(iv)+1};
        case {'-c'}
            unClustCrit = varargin{varStrInd(iv)+1};
            
    end
end

% Here are the sets of options/combinations to use: They include scaling
% (respNorms), similarity measurements (respSims), how to calculate
% randomized spaces for individual clusters (respRand), which
% measurements of the SDFs to use (respInclude), and whether to z-score
% columns of each measurement (zResp).
respNorms = {'ztr','ztrbl','max','none','bl','zbl'};
respSims = {'corr','euc'};
respRand = {'pca'};
respInclude = {'mean','combo','slope','whole'};
zResp = 0;

% The above section has all possibilities, here are the same variables but
% subsets can be taken for the final clustering if desired.
whichNorms = respNorms;
whichSims = respSims;
whichRand = respRand;
whichInclude = respInclude;
whichZ = 0;


if ~exist('respSumStruct','var')
    %% Initialize output and printing variables
    isSameGroup = zeros(size(goodSDF{1},1),size(goodSDF{1},1),length(respNorms),length(respSims),length(respRand),length(respInclude),length(zResp));
    nNeeded = prod([length(respNorms),length(respSims),length(respRand),length(respInclude),length(zResp)]);
    nDone = 0;
    clear respSumStruct
    respSumStruct(1:nNeeded) = struct('ids',nan,'gap',nan,'gapErr',nan,'normType',nan,'simType',nan,'randType',nan,'include',nan,'k',nan);
    printStr = [];
    %% Loop on through and try the clustering procedures
    for in = 1:length(respNorms)
        for is = 1:length(respSims)
            for ir = 1:length(respRand)
                for ii = 1:length(respInclude)
                    for iz = 1:length(zResp)
                        nDone = nDone+1;
                        % Perform the clustering
                        fprintf(repmat('\b',1,length(printStr)));
                        printStr = sprintf('Clustering combination %d of %d...\n',nDone,nNeeded);
                        fprintf(printStr);
                        [respK,respClustIDs,respGap,respGapErr,linkMat,linkInd, distMat, respIn] = indivResponseCluster(goodSDF,allTimeCell,'norm',respNorms{in},'sim',respSims{is},'rand',respRand{ir},'resp',respInclude{ii},'z',zResp(iz),'-n',minN,'-r',randReps,'-e',myEpocs,'ei',myEpocInds,'rw',myEpocWinds,'sd',normSD,'do',distOnly);
                        % Save parameters and results
                        respSumStruct(nDone).ids        = respClustIDs;
                        respSumStruct(nDone).gap        = respGap;
                        respSumStruct(nDone).gapErr     = respGapErr;
                        respSumStruct(nDone).normType   = respNorms{in};
                        respSumStruct(nDone).simType    = respSims{is};
                        respSumStruct(nDone).randType   = respRand{ir};
                        respSumStruct(nDone).include    = respInclude{ii};
                        respSumStruct(nDone).k          = respK;
                        respSumStruct(nDone).z          = zResp(iz);
                        respSumStruct(nDone).linkMat    = linkMat;
                        respSumStruct(nDone).linkInd    = linkInd;
                        respSumStruct(nDone).distMat    = distMat;
                        respSumStruct(nDone).respIn     = respIn;
                        if size(respClustIDs,2) >= setK
                            for iii = 1:size(respClustIDs,1)
                                isSameGroup(iii,respClustIDs(:,setK)==respClustIDs(iii,setK),in,is,ir,ii,iz) = 1;
                            end
                        else
                            isSameGroup(:,:,in,is,ir,ii,iz) = nan;
                        end
                    end
                end
            end
        end
    end
end
fprintf(repmat('\b',1,length(printStr)));
fprintf('Finished pulling individual pipelines... Now clustering consensus matrix\n');

%% Pull out the relevant procedures
includes = {respSumStruct.include};
norms = {respSumStruct.normType};
rands = {respSumStruct.randType};
sims = {respSumStruct.simType};
zs  = [respSumStruct.z];
includeParams = find(ismember(includes,whichInclude) & ismember(norms,whichNorms) & ismember(rands,whichRand) & ismember(sims,whichSims) & ismember(zs,whichZ));

% Create and stack Z-scored distance matrix
allDist = nan(size(goodSDF{1},1),size(goodSDF{1},1),length(includeParams));
for ii = 1:length(includeParams)
    allDist(:,:,ii) = (respSumStruct(includeParams(ii)).distMat-nanmean(respSumStruct(includeParams(ii)).distMat(:)))./nanstd(respSumStruct(includeParams(ii)).distMat(:));
end

% Get the composite distance matrix my summarizing the stack
switch meanMedian
    case {'mean'}
        raw = nanmean(allDist,3);
    case {'median'}
        raw = nanmedian(allDist,3);
end

% Cluster the new distance matrix
[sortIDs,idxDist,~,~,rawLink] = distMatAgglom(raw,1:30,'-n',metaN);
% Cut off columns where k clusters didn't meet the clustering criteria
% (columns of all NaN)
while sum(isnan(sortIDs(:,end)))==size(sortIDs,1)
    sortIDs = sortIDs(:,1:(end-1));
end

% Adjust the linkage matrix in order to have no negative values (negative
% values throw off the "dendrogram" function
rawLink(:,3) = rawLink(:,3)-min(rawLink(:,3));

% Get the number of good and bad units where k = 1 (this eliminates any
% strage occurrences where a unit should never have been included in the
% first place)
nBad = sum(isnan(sortIDs(:,1)));
nGood = sum(~isnan(sortIDs(:,1)));

% Calculate percentages of NaN for each k
nNan = nan(1,size(sortIDs,2));
for ik = 1:size(sortIDs,2)
    nNan(ik) = (sum(isnan(sortIDs(:,ik)))-nBad)./nGood;
end

% Get a final K by finding the last K that meets the unclustered criterion
myK = find(nNan <= unClustCrit,1,'last');

