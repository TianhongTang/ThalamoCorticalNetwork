clear all; clc; close all;
addpath('/Users/Ilya/Dropbox/HELPER/HELPER_GENERAL/');
addpath('C:\Users\Ilya Monosov\Dropbox\HELPER\HELPER_GENERAL');
load('/Users/Ilya/Desktop/BFINAC_heatmap_datastruct.mat')


% % row#1 - Trial onset(after onset vs before onset)
% % row#2 - Go onset(6101vs6102)
% % row#3 - Go onset(6103vs6104)
% % row#4 - Flip onset(6101vs6102)
% % row#5 - Flip onset(6103vs6104)
% % row#6 - BEFORE flip onset - 6101(Left) vs 6102(Left)
% % row#7 - BEFORE flip onset - 6103(Left) vs 6104(Left)
% % row#8 - BEFORE flip onset - 6101(Right) vs 6102(Right)
% % row#9 - BEFORE flip onset - 6103(Right) vs 6104(Right)
% % row#10 - FreeReward(after onset vs before onset)
% % row#11 - Sound flash(after onset vs before onset)
% % row#12 - FreeReward vs Sound flash
% % row#13 - Go onset(Left vs Right)
% % row#14 - Flip onset(Left vs Right)


rng('default');  %this sets the "fixed" seed for clustering algorithms; basically so clustering is the same eachtime you do it

%find groups
for x=1:length(datastruct)
    temp=datastruct(x).CELL;
    c{x}=temp(1:3);       %%%%%%REMOVE TAKAYA's designation of "significant or not"
    datastruct(x).CELL=temp(1:3);
    clear temp
end
[Dd, id] = findgroups(c)
counts = histcounts(Dd)

datastruct=datastruct(find(Dd~=3))

%%
Firingrates=[datastruct(:).BaselineFiringrate];
Firingrates=(Firingrates-min(Firingrates));
Firingrates=(Firingrates./max(Firingrates));

savecellvectors=[];
for x=1:length(datastruct)
    
    test=contains(id,datastruct(x).CELL);
    GroupAssignment=find(test==1); clear test;
    
    conditions=[2 3 6 7 8 9 12 13];
    
    onecelldata=datastruct(x).ROC;
    onecelldata=onecelldata(conditions,:);
    
    onecelldataP=datastruct(x).Pvalue;
    onecelldataP=onecelldataP(conditions,:);
    %%%  lengthofeacharray=length(onecelldataP(2,:));
    onecelldata(find(onecelldataP>=0.05))=.5;
    %%
    
    
    onecelldataP=datastruct(x).Pvalue_Short;
    onecelldataP=onecelldataP(conditions,:);
    
    % lengthofeacharray=length(onecelldataP(2,:));
    %Threshold
    onecelldata(find(onecelldataP>=0.01),:)=.5;
    %%
    
    
    for xx=1:length(onecelldata(:,1))
        
        
        
        
        
        if length(find(onecelldata(xx,:)~=0.5))<50
            onecelldata(xx,1:length(onecelldata))=0.5;
        else
            onecelldata(xx,:)=smooth(onecelldata(xx,:),50);
        end
        
        
        temp=findseq(onecelldata(xx,:));
        if size(temp,1)>1
            if isempty(find(temp(:,1)~=0.5 & temp(:,4)>10))==1
                onecelldata(xx,:)=0.5;
                
            else
                %    onecelldata(xx,:)=smooth(onecelldata(xx,:),100);
            end
        end
        
    end
    onecelldata=onecelldata';
    onecelldata=onecelldata(:);
    
    % %     %add array of firing rates (same length as others)
    % %     temp(1:lengthofeacharray)=Firingrates(x);
    % %     temp=temp';
    % %     onecelldata=[onecelldata; temp]; clear temp;
    % %
    onecelldata(end+1)=GroupAssignment;
    savecellvectors=[savecellvectors; onecelldata'];
    clear onecelldata
    
end


%totalsize= round((size(savecellvectors(:,:),2)-1) ./ 400)



savecellvectors=savecellvectors(find(mean(savecellvectors(:,1:end-1)')~=0.5),:);

savecellvectors=savecellvectors(find(isnan(mean(savecellvectors(:,1:end-1)'))~=1),:);


G=savecellvectors(:,1:end-1); %remove the group assignment

[pc, zscores, pcvars] = pca(G,'Rows','pairwise'); %lets calculate Principal Comp

figuren;

%lets see which principal components explain most of the variance in your ndata
%nsubplot(20,30,1:4,1:4);set(gca,'ticklength',4*get(gca,'ticklength'));


VarExp=pcvars./sum(pcvars) * 100 %variance explained by each PC (PC1-PC6 i think)
%%%cumsum(pcvars./sum(pcvars) * 100) %cum sum of variance ...
plot(VarExp,'k','LineWidth',2);
scatter([1:length(VarExp)]', VarExp, 20, 'r', 'filled')
xlabel('PCs')
ylabel('% Variance Explained')
xlim([0 10])
ylim([0 100])

% nsubplot(20,30,11:14,1:4);set(gca,'ticklength',4*get(gca,'ticklength'));
% %%lets plot the PC weights // using matlab we can rotate this and explore it
% %%in 3d space; we can then use a simple setting to save a good orientation
% %%for plotting; however from the first plot we can see really we need more
% %%PCs to explain most of the data
% scatter3(zscores(:,1),zscores(:,2),zscores(:,3),20,'k','filled');
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% title('first 3 PCs. you can rotate using matlab ')

%note to Takaya: for the paper we may add a "shuffling" method to see which
%PCs are "significant" in terms of variance explained; we can then use that to choose number of PCs for clustering


X1=zscores(:,1:3)';
% VB fitting
[y1, model, L] = mixGaussVb(X1,1000);
%figure;
%plotClass(X1,y1);
%figure;
%plot(L)
% Predict testing data
[y2, R] = mixGaussVbPred(model,X1);
%figure;
%plotClass(X1,y2);


NUM =  max(y2)
D = pdist(zscores(:,1:10), 'euclidean'); %I chose the first 10 PCs because the plot above indicates that they explain most of the variance
T = linkage(D, 'ward');


figuren;

subplot(1,5,1); set(gca,'ticklength',4*get(gca,'ticklength'));
[H,T,outperm] = dendrogram(T, 0, 'colorthreshold',mean(T(end-NUM+1:end-NUM+2,3)),'Orientation','left');
set(H, 'LineWidth',1)
set(gca, 'XTickLabel',[], 'TickLength',[0 0])


subplot(1,5,2); set(gca,'ticklength',4*get(gca,'ticklength'));
savecormat_half=zscores(outperm,1:10);
imagesc(flipud(savecormat_half))
colormap('bone');colorbar; caxis([min(savecormat_half(:)) max(savecormat_half(:))]);
set(gca,'YTickLabel',[]);
xlabel('PCs')
ylabel('neurons')


subplot(1,5,3); set(gca,'ticklength',4*get(gca,'ticklength'));
savecormat_half=savecellvectors(outperm,1:end-1);
imagesc(flipud(savecormat_half))
colormap('bone');colorbar; caxis([0.2 0.8]);
ylabel('neurons')
xlabel('1000ms ROC vectors')
% names = {'pun val'; 'pun unc';  'rew val'; 'rew unc'; 'pun info'; 'rew info'; };
% set(gca,'xtick',[1:6],'xticklabel',names)
% xtickangle(45);
set(gca,'YTickLabel',[]);

%text(10,-4,['n=' mat2str(size(G1,1))  ])


subplot(1,5,4); set(gca,'ticklength',4*get(gca,'ticklength'));
savecormat_half=savecellvectors(outperm,end);
imagesc(flipud(savecormat_half))
caxis([min(Dd) max(Dd)]);

colormap('bone');colorbar;
ylabel('neurons')
xlabel(cell2mat(id))
set(gca,'YTickLabel',[]);














