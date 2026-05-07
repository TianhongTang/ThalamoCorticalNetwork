function normAlign = scaleResp(respAlign,respTimes,respNormType,varargin)

blWind = -300:-100;
respWinds = cell(1,length(respAlign));
[respWinds{1:end}] = deal(-200:200);
calcSD = 1:length(respAlign);

% Decode varargin
varStrInd = find(cellfun(@ischar,varargin));
for iv = 1:length(varStrInd)
    switch varargin{varStrInd(iv)}
        case {'bl','-bl'}
            blWind = varargin{varStrInd(iv)+1};
        case {'r','-r'}
            respWinds = varargin{varStrInd(iv)+1};
        case {'sd'}
            calcSD = varargin{varStrInd(iv)+1};
    end
end

normAlign = cell(size(respAlign));
switch respNormType
    case {'zbl','z'}
        blMeans = nanmean(respAlign{1}(:,ismember(respTimes{1},blWind)),2);
        blStds = nanstd(respAlign{1}(:,ismember(respTimes{1},blWind)),[],2);
        for i = 1:length(respAlign)
            if ismember(i,calcSD)
                normAlign{i} = (respAlign{i}-repmat(blMeans,1,size(respAlign{i},2)))./repmat(blStds,1,size(respAlign{i},2));
            else
                normAlign{i} = respAlign{i}./repmat(blStds,1,size(respAlign{i},2));
            end
            
        end
    case {'ztr'}
        catVals = [];
        for i = 1:length(calcSD)
            catVals = cat(2,catVals,respAlign{calcSD(i)}(:,ismember(respTimes{calcSD(i)},respWinds{calcSD(i)})));
        end
        blMeans = nanmean(catVals,2);
        blStds = nanstd(catVals,[],2);
        for i = 1:length(respAlign)
            if ismember(i,calcSD)
                normAlign{i} = (respAlign{i}-repmat(blMeans,1,size(respAlign{i},2)))./repmat(blStds,1,size(respAlign{i},2));
            else
                normAlign{i} = respAlign{i}./repmat(blStds,1,size(respAlign{i},2));
            end
        end
    case {'ztrbl'}
        catVals = [];
        for i = 1:length(calcSD)
            catVals = cat(2,catVals,respAlign{calcSD(i)}(:,ismember(respTimes{calcSD(i)},respWinds{calcSD(i)})));
        end
        blMeans = nanmean(respAlign{1}(:,ismember(respTimes{1},blWind)),2);
        blStds = nanstd(catVals,[],2);
        for i = 1:length(respAlign)
            if ismember(i,calcSD)
                normAlign{i} = (respAlign{i}-repmat(blMeans,1,size(respAlign{i},2)))./repmat(blStds,1,size(respAlign{i},2));
            else
                normAlign{i} = respAlign{i}./repmat(blStds,1,size(respAlign{i},2));
            end
            
        end
    case {'max'}
        blMeans = nanmean(respAlign{1}(:,ismember(respTimes{1},blWind)),2);
        catVals = [];
        for i = 1:length(calcSD)
            catVals = cat(2,catVals,respAlign{calcSD(i)}(:,ismember(respTimes{calcSD(i)},respWinds{calcSD(i)})));
        end
        absMax = max(abs(catVals-repmat(blMeans,1,size(catVals,2))),[],2);
        for i = 1:length(respAlign)
            if ismember(i,calcSD)
                normAlign{i} = (respAlign{i}-repmat(blMeans,1,size(respAlign{i},2)))./repmat(absMax,1,size(respAlign{i},2));
            else
                normAlign{i} = respAlign{i}./repmat(absMax,1,size(respAlign{i},2));
            end
        end
    case {'scale'}
        catVals = [];
        for i = 1:length(calcSD)
            catVals = cat(2,catVals,respAlign{calcSD(i)}(:,ismember(respTimes{calcSD(i)},respWinds{calcSD(i)})));
        end
        minVals = min(catVals,[],2);
        maxVals = max(catVals,[],2);
        for i = 1:length(respAlign)
            if ismember(i,calcSD)
                normAlign{i} = (respAlign{i}-repmat(minVals,1,size(respAlign{i},2)))./repmat(maxVals-minVals,1,size(respAlign{i},2));
            else
                normAlign{i} = respAlign{i}./repmat(maxVals-minVals,1,size(respAlign{i},2));
            end
        end
    case {'base','bl'}
        blMeans = nanmean(respAlign{1}(:,ismember(respTimes{1},blWind)),2);
        for i = 1:length(respAlign)
            if ismember(i,calcSD)
                normAlign{i} = respAlign{i}-repmat(blMeans,1,size(respAlign{i},2));
            else
                normAlign{i} = respAlign{i};
            end
        end
    case {'sampmax'}
        overallVals = [];
        for i = 1:length(calcSD)
            overallVals = cat(3,overallVals,respAlign{calcSD(i)});
        end
        overallMn = nanmean(overallVals,3);
        maxVal = max(nanmean(overallMn,1));
        for i = 1:length(respAlign)
            normAlign{i} = respAlign{i}./maxVal;
        end
        
    case {'none'}
        for i = 1:length(respAlign)
            normAlign{i} = respAlign{i};
        end
end
