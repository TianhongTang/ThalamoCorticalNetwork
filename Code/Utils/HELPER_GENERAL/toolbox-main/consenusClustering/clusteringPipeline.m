function [clusterNeurons,sortIDs,normResp] = clusteringPipeline_SfN(inputSDF,sdfTimes,sdfEpoch,lineColors,colorMapping,nClusters)

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping);
normResp = scaleResp(inputSDF,sdfTimes,'max');

nClusters_manual = nClusters; clusterNeurons = [];
for i = 1:nClusters_manual
    clusterNeurons{i} = find(sortIDs(:,nClusters_manual) == i );
end

%% Plot clustering output
% Dendrogram
figure('Renderer', 'painters', 'Position', [100 100 500 400]);

subplot(1,5,5);
for ir = 1:size(raw,1)
    for ic = (ir+1):size(raw,2)
        raw(ic,ir) = raw(ir,ic);
    end
end
[h,~,outPerm] = dendrogram(rawLink,0,'Orientation','right');
set(gca,'YDir','Reverse');
klDendroClustChange(h,rawLink,sortIDs(:,nClusters_manual))
set(gca,'YTick',[]); xlabel('Similarity')
subplot(1,5,[1:4]);

imagesc(raw(outPerm,outPerm));
colormap(gray);
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
xticks([50:50:500]); yticks([50:50:500])

% SDFs
for i = 1:nClusters_manual
    figure('Renderer', 'painters', 'Position', [100 100 500 400]);hold on
    
    for ii = 1:size(inputSDF,2)
        plot(sdfTimes{ii},nanmean(normResp{ii}(clusterNeurons{i},:),1), 'color', lineColors{ii});
    end
    vline(0, 'k--'); xlim([sdfEpoch{1}(1) sdfEpoch{1}(end)])
    title(['Cluster ' int2str(i) ' - n: ' int2str(length(clusterNeurons{i}))])
    
end





