%% indivResponseCluster.m performs clustering according to specific parameters
%  set by consensusCluster.m



function [respK, respClustIDs, respGap, respGapErr, linkMat, linkInd, thisDistOut, respInOut] = indivResponseCluster(allSDF,allTimes,varargin)

%% Set Defaults
respNormType = 'none';
respSimType = 'corr';
respInclude = 'combo';
zResp = 1;
sdfK = 1:20;
blWind = -300:-100;
nRespMembers = 10;
randReps = 10; 
direction = 'first';
dirTypes = {'first','last'};
doGap = 1;
normSD = 1:length(allSDF);
distOnly = 0;

% Default analysis epochs
preVis = -100:0;
visTrans = 50:100;
visSust = 100:150;
preMov = -50:0;
postMov = 0:50;
nextVis = 50:100;
allEpocs = {preVis,visTrans,visSust,preMov,postMov,nextVis};
epocInds = [1,1,1,2,2,2];
respWinds = allTimes;

%% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'norm'}
            respNormType = varargin{varStrInd(iv)+1};
        case {'sim'}
            respSimType = varargin{varStrInd(iv)+1};
        case {'resp'}
            respInclude = varargin{varStrInd(iv)+1};
        case {'z'}
            zResp = varargin{varStrInd(iv)+1};
        case {'blwind'}
            blWind = varargin{varStrInd(iv)+1};
        case {'-n'}
            nRespMembers = varargin{varStrInd(iv)+1};
        case {'-d'}
            direction = varargin{varStrInd(iv)+1};
            if isnumeric(direction)
                direction = dirTypes{direction};
            end
        case {'-g','gap'}
            doGap = varargin{varStrInd(iv)+1};
        case {'-r','reps'}
            randReps = varargin{varStrInd(iv)+1};
        case {'-e','epocs'}
            allEpocs = varargin{varStrInd(iv)+1};
        case {'ei'}
            epocInds = varargin{varStrInd(iv)+1};
        case {'rw'}
            respWinds = varargin{varStrInd(iv)+1};
        case {'sd'}
            normSD = varargin{varStrInd(iv)+1};
        case {'do'}
            distOnly = varargin{varStrInd(iv)+1};
    end
end

%% Start the analysis

% Normalize/Scale the SDFs
normResp = scaleResp(allSDF,allTimes,respNormType,'-r',respWinds,'bl',blWind,'sd',normSD);

% Cut down the responses to relevant times and to get rid of nans
catNorms = [];
for i = 1:length(normResp)
    catNorms = cat(2,catNorms,normResp{i}(:,ismember(allTimes{i},respWinds{i})));
end
goodRows = sum(isfinite(catNorms),2) == size(catNorms,2);
for i = 1:length(normResp)
    normResp{i} = normResp{i}(goodRows,:);
end

% Set up epocs
epocSpks = cell(1,length(allEpocs));
epocTms = cell(1,length(allEpocs));
for i = 1:length(allEpocs)
    epocSpks{i} = normResp{epocInds(i)};
    epocTms{i} = allTimes{epocInds(i)};
end

% Smooth the SDFs a little bit in case we want to cluster the whole SDF
kern = makeGauss(10,'-x',-50:50);
convResp = [];
for i = 1:length(normResp)
    convResp = cat(2,convResp,conv2(normResp{i}(:,ismember(allTimes{i},respWinds{i})),kern,'same'));
end

%% Get the proper input ready
startTic = tic;
[respMat, respSlope] = parseEpochs(epocSpks,epocTms,allEpocs);

switch respInclude
    case {'mean','means'}
        respIn = respMat;
    case {'slope','slopes'}
        respIn = respSlope;
    case {'comb','combo'}
        respIn = [respMat,respSlope];
    case {'whole'}
        respIn = convResp;
end
if zResp
    respIn = (respIn-repmat(nanmean(respIn,1),size(respIn,1),1))./repmat(nanstd(respIn,[],1),size(respIn,1),1);
end

%% Make a distance matrix
thisDist = makeDistMat(respIn,respSimType,zeros(1,size(respIn,2)));
thisDistOut = nan(size(catNorms,1));
thisDistOut(goodRows,goodRows) = thisDist;

if distOnly
    [respK, respClustIDs, respGap, respGapErr, linkMat, linkInd] = deal(nan);
    respInOut = nan(size(catNorms,1),size(respIn,2));
    respInOut(goodRows,:) = respIn;
    return;
end



%% Cluster the distance matrix
[respClustIDs,~,~,~,linkMat, linkInd] = distMatAgglom(thisDist,sdfK,'-n',nRespMembers);
nBad = find(~goodRows);
for i = 1:length(nBad)
    check = linkMat(:,1:2);
    linkMat(check >= nBad(i)) = linkMat(check >= nBad(i))+1;
    linkMat = cat(1,linkMat,[nBad(i),max(check(:))+1,linkMat(end,3)]);
end
while sum(isnan(respClustIDs(:,end))) == size(respClustIDs,1)
    respClustIDs = respClustIDs(:,1:(end-1));
end

% Calculate W for each value of k (see Tibshirani et al 2001 for gap
% statistic calculations). 
W = nan(1,size(respClustIDs,2));
for ik = 1:size(respClustIDs,2)
    W(ik) = getW(respClustIDs(:,ik),thisDist);
end

if doGap
    %% Get the gap statistic (See Tibshirani et al 2001)
    wStar = nan(randReps,size(respClustIDs,2));
    for iRand = 1:randReps
        printStr = sprintf('Working on random repetition %d of %d...',iRand,randReps);
        fprintf(printStr);

        % Get a random set of data according to the mean/SD of the inputs
        randMat = makeRandomSet(size(respIn,1),nanmean(respIn,1),nanstd(respIn,1));
        randDist = makeDistMat(randMat,respSimType,zeros(1,size(randMat,2)));
        randIDs = distMatAgglom(randDist,sdfK,'-n',nRespMembers,'-p',0);
        while sum(isnan(randIDs(:,end))) == size(randIDs,1)
            randIDs = randIDs(:,1:(end-1));
        end
        % Calculate randomized W
        for ik = 1:min([size(wStar,2),size(randIDs,2)])
            wStar(iRand,ik) = getW(randIDs(:,ik),randDist);
        end
        for ib = 1:length(printStr)
            fprintf('\b');
        end
    end
    fprintf('Calculating on gap...\n');

    % First, calculate the gap itself
    respGap = nanmean(log(wStar)-repmat(log(W),size(wStar,1),1),1);

    % Now we need lBar
    lBar = nanmean(log(wStar),1);

    % Now s(respGapErr)/sd...
    sd = sqrt(nanmean((log(wStar)-repmat(lBar,size(wStar,1),1)).^2,1));
    respGapErr = sd.*sqrt(1+(1/randReps));

    % Now get the "official" k
    gapComp = [respGap(2:end)-respGapErr(2:end),nan];
    respK = find(respGap > gapComp,1,direction);
    if isempty(respK)
        respK = nan;
    end
else
    respGap = nan;
    respGapErr = nan;
    respK = nan;
end

% Set the variables aside...
respClustTmp = respClustIDs;
normRespTmp = cell(size(normResp));
for i = 1:length(normResp)
    normRespTmp{i} = normResp{i};
end
clear respClustIDs normResp;

% Initialize the output variables
respClustIDs = nan(size(catNorms,1),size(respClustTmp,2));
for i = 1:length(normRespTmp)
    normResp{i} = nan(size(catNorms,1),size(normRespTmp{i},2));
end
respInOut = nan(size(catNorms,1),size(respIn,2));

% Place values in the output variables
respClustIDs(goodRows,:) = respClustTmp;
for i = 1:length(normRespTmp)
    normResp{i}(goodRows,:) = normRespTmp{i};
end
respInOut(goodRows,:) = respIn;

%% And we're done!

function distMat = makeDistMat(obs,type,catVars,print)
    if nargin < 4
        print = 0;
    end
    
    % Start by computing pair-wise distances
    if print
        fprintf('Initializing pairwise distances...')
    end
    if ismember(type,{'corr'})
        % Clean up obs, if "corr"
        goodCols = zeros(1,size(obs,2));
        for ic = 1:size(obs,2)
            goodCols(ic) = sum(~isnan(obs(:,ic))) == size(obs,1);
        end
        obs = obs(:,logical(goodCols));
        distMat = 1-corr(obs',obs');
        for ir = 1:size(obs,1)
            for ic = ir:size(obs,1)
                distMat(ic,ir) = nan;
            end
        end
    elseif ismember(type,{'euc'})
        distMat = getEucDist(obs(:,~catVars));
        for ir = 1:size(obs,1)
            for ic = ir:size(obs,1)
                distMat(ic,ir) = nan;
            end
        end
    else
        error('Unsupported distance %s requested...\n',type);
    end
    if print
        fprintf('Done!\n');
    end
end

function timeStr = printTiming(ticVal)
    thisRawTime = toc(ticVal);
    thisTimeMin = floor(thisRawTime/60);
    if thisTimeMin > 60, thisTimeHour = floor(thisTimeMin/60); thisTimeMin = mod(thisTimeMin,60); else thisTimeHour = 0; end
    thisTimeSec = mod(thisRawTime,60);

    thisHrStr = num2str(thisTimeHour); while length(thisHrStr) < 2, thisHrStr = ['0',thisHrStr]; end
    thisMinStr = num2str(thisTimeMin); while length(thisMinStr) < 2, thisMinStr = ['0',thisMinStr]; end
    thisSecStr = sprintf('%.2f',thisTimeSec);  while length(thisSecStr) < 5, thisSecStr = ['0',thisSecStr]; end
    timeStr = sprintf('%s:%s:%s',thisHrStr,thisMinStr,thisSecStr);
end

% The below function grabs the sum of the peirwise distances amongst
% members of a cluster
function newD = dFromDist(dist, indices)
    % Nice, fast version by making a submatrix...
    newMat = dist(indices,indices);
    newD = nansum(newMat(:));
end

function w = getW(clusts,dist)
    uClusts = unique(clusts(~isnan(clusts)));
    allD = nan(1,length(uClusts));
    allN = nan(1,length(uClusts));
    for ic = 1:length(uClusts)
        clustInds = find(clusts==uClusts(ic));
        allD(ic) = dFromDist(dist,clustInds);
        allN(ic) = length(clustInds);
    end
    w = sum(allD./(2.*allN));
end

function out = getEucDist(x,y)
    if nargin < 2
        y = x;
    end
    if size(x,2) ~= size(y,2)
        error('Dimension mismatch');
    end
    runSum = zeros(size(x,1),size(y,1));
    for ic = 1:size(x,2)
        runSum = runSum + (abs(bsxfun(@minus,x(:,ic),y(:,ic)')).^2);
    end
    out = sqrt(runSum);
end

function [out, xvals] = makeGauss(sd,varargin)
    % Set defaults
    amp = 1;
    mu = 0;
    xvals = -100:100;
    scale = 1;
    const = 0;
    % Decode varargin
    varStrings = find(cellfun(@ischar,varargin));
    for is = 1:length(varStrings)
        switch varargin{varStrings(is)}
            case '-x'
                xvals = varargin{varStrings(is)+1};
            case '-m'
                mu = varargin{varStrings(is)+1};
            case '-a'
                amp = varargin{varStrings(is)+1};
            case '-s'
                scale = varargin{varStrings(is)+1};
            case '-c'
                const = varargin{varStrings(is)+1};
        end
    end
    normFun = ((1/sqrt(2*pi*sd.^2))*exp(-(xvals-mu).^2/(2*sd.^2)));
    if scale
        out = (normFun.*(amp/max(normFun)))+const;
    else
        out = (normFun.*amp)+const;
    end
end

function [epochMat, epochSlopes] = parseEpochs(spks,times,epochs,varargin)
epochMat    = cell2mat(cellfun(@getSpksFun,spks,times,epochs,'UniformOutput',0));
epochSlopes = cell2mat(cellfun(@getSlopesFun,spks,times,epochs,'UniformOutput',0));

    function outMat=getSpksFun(spks,times,epochs)
        outMat=nanmean(spks(:,ismember(times,epochs)),2);
    end

    function outSlope = getSlopesFun(spks,times,epochs)
        [~,outSlope,~] = regression(repmat(epochs,size(spks,1),1),spks(:,ismember(times,epochs)));
    end
end

function randVals = makeRandomSet(n,params1,params2,varargin)
    % Set defaults
    dists = ones(1,length(params1));
    % Decode varargin
    varStrings = find(cellfun(@ischar,varargin));
    for is = 1:length(varStrings)
        switch varargin{varStrings(is)}
            case {'-d','dists','dist'}
                dists = varargin{varStrings(is)+1};
        end
    end
    % Initialize output
    randVals = nan(n,length(params1));
    % Get gaussian distributions
    isGauss = dists==1;
    randTmp = randn(n,sum(isGauss));
    % Here, params2 = SD and params1 = mu
    randTmp = (randTmp.*repmat(params2(isGauss),n,1))+repmat(params1(isGauss),n,1);
    randVals(:,isGauss) = randTmp;
    clear randTmp;
    % Get uniform distributions
    isUnif = dists==2;
    if any(isUnif)
        randTmp = rand(n,sum(isUnif));
        % Here, params1 = min, params2=max
        randTmp = (randTmp.*repmat(params2(isUnif)-params1(isUnif),n,1))+repmat(params1(isUnif),n,1);
        randVals(:,isUnif) = randTmp;
    end
end


end