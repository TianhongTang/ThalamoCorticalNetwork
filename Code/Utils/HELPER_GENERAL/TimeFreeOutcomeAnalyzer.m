%IEM HELPER
clc; clear all; close all; delete(gcf); clf;
addpath('HELPER_GENERAL'); %keep hepers functions here
CENTER=8001; %this is event alignment. this means that at 11,000 the event that we are studying (aligning spikes too) will be in the SDF *you will see in ploting
gauswindow_ms=50;
rastersplot=1;
limitsarray=[CENTER-1000:CENTER+1750];
xlimvar=length(limitsarray);
ylimvar=100;
load('Y:\MONKEYDATA\Robin_ongoing\TimingProcedureBF\TimingProcedureV3_19_12_2016_10_23.mat','PDS')


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%first get the data in the context where there are single rewards
% % %%%%%%%%and flash/sound events
% % 
% load('Y:\Kael\Batman\ProbAmtIsoLuminantV3\ProbAmtIsoLum_V3_07_07_2015_11_13.mat')
% 
% 
% close all;
% PDS.timesoffreeoutcomes_first(find(PDS.timesoffreeoutcomes_first>10))=NaN; %very few trials have a timing bug. just remove them
% free1=intersect(find(PDS.timesoffreeoutcomes_first>0),find(PDS.freeoutcometype==1)); %reward/juice
% free34=intersect(find(PDS.timesoffreeoutcomes_first>0),find(PDS.freeoutcometype==34)); %flash.sound combination
% %let's remove the reward noise
% cleanfreeoutcome=[free1 free34];
% CLEAN=cleanfreeoutcome;
% for x=1:length(CLEAN)
%     x
%     x=CLEAN(x);
%     spk=PDS(1).sptimes{x};
%     spk1=spk(find(spk<PDS.timesoffreeoutcomes_first(x)));
%     spk2=spk(find(spk>PDS.timesoffreeoutcomes_first(x)+10/1000));
%     spk=[spk1 spk2]; clear spk1 spk2
%     try
%         no_noise=find(spk>(PDS.timesoffreeoutcomes_first(x)+(10/1000)) & spk<(PDS.timesoffreeoutcomes_first(x)+(20/1000)));
%         no_noise=length(no_noise); %this many spikes in the window  so insert that many spikes in the deleted time range
%         grabnumbers=PDS.timesoffreeoutcomes_first(x)+randperm(10)/1000;
%         grabnumbers=grabnumbers(no_noise);
%         spk=[spk grabnumbers];
%     end
%     PDS(1).sptimes{x}=sort(spk);
%     %
%     clear y spk spk1 spk2 x
%     clc
% end
% %MAKE A RASTER AND A SPIKE DENSITY FUNCTION FOR EVERY TRIAL for
% %freeoutcomes
% Rasters=[];
% for x=1:length(PDS.fractals)
%     x
%     spk=PDS(1).sptimes{x}-PDS.timesoffreeoutcomes_first(x);
% 
%     spk=(spk*1000)+CENTER-1;
%     spk=fix(spk);
%     %
%     spk=spk(find(spk<CENTER*2));
%     %
%     temp(1:CENTER*2)=0;
%     try
%     temp(spk)=1;
%     catch
%     temp(1:length(temp))=NaN;    
%     end
%     Rasters=[Rasters; temp];
%     %
%     clear temp spk x
% end
% SDFfree=plot_mean_psth({Rasters},gauswindow_ms,1,(CENTER*2)-2,1); close all;
% Rastersfree=Rasters;
% clear PDS CLEAN cleanfreeoutcome
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %now load the timing trace for the same neuron -- in order to study how the
% %neuron responds to double versus single events (and airpuffs)
% close all
% load('Y:\Kael\Batman\TimingTrace\TimingTrace_V7_21_08_2015_13_31.mat')
% 
% PDS.timesoffreeoutcomes_first(find(PDS.timesoffreeoutcomes_first>10))=NaN; %very few trials have a timing bug. just remove them
% free1_=intersect(find(PDS.timesoffreeoutcomes_first>0),find(PDS.freeoutcometype==1)); %reward/juice
% free2_=intersect(find(PDS.timesoffreeoutcomes_first>0),find(PDS.freeoutcometype==2)); %aifpuff
% free34_=intersect(find(PDS.timesoffreeoutcomes_first>0),find(PDS.freeoutcometype==34)); %flash.sound combination
% durationsuntilreward=PDS(1).timeoutcome-PDS(1).timetargeton;
% durationsuntilreward=round(durationsuntilreward*10)./10;
% completedtrial=find(durationsuntilreward>0);
% deliv=find(PDS(1).deliveredornot==1)
% ndeliv=find(PDS(1).deliveredornot==0)
% trials6253=intersect(find(PDS(1).fractals==6253),completedtrial);
% trials6254=intersect(find(PDS(1).fractals==6254),completedtrial);
% trials6259=intersect(find(PDS(1).fractals==6259),completedtrial);
% trials6260=intersect(find(PDS(1).fractals==6260),completedtrial);
% trials6300=intersect(find(PDS(1).fractals==6300),completedtrial);
% trials6301=intersect(find(PDS(1).fractals==6301),completedtrial);
% trials6254d=intersect(trials6254,deliv);
% trials6260d=intersect(trials6260,deliv);
% trials6301d=intersect(trials6301,deliv);
% trials6254nd=intersect(trials6254,ndeliv);
% trials6260nd=intersect(trials6260,ndeliv);
% trials6301nd=intersect(trials6301,ndeliv);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cleanfreeoutcome=[free1_ free2_ free34_];
% CLEAN=cleanfreeoutcome;
% for x=1:length(CLEAN)
%     x=CLEAN(x);
%     spk=PDS(1).sptimes{x};
%     spk1=spk(find(spk<PDS.timesoffreeoutcomes_first(x)));
%     spk2=spk(find(spk>PDS.timesoffreeoutcomes_first(x)+10/1000));
%     spk=[spk1 spk2]; clear spk1 spk2
%     try
%         no_noise=find(spk>(PDS.timesoffreeoutcomes_first(x)+(10/1000)) & spk<(PDS.timesoffreeoutcomes_first(x)+(20/1000)));
%         no_noise=length(no_noise); %this many spikes in the window  so insert that many spikes in the deleted time range
%         grabnumbers=PDS.timesoffreeoutcomes_first(x)+randperm(10)/1000;
%         grabnumbers=grabnumbers(no_noise);
%         spk=[spk grabnumbers];
%     end
%     PDS(1).sptimes{x}=sort(spk);
%     %
%     clear y spk spk1 spk2 x
%     clc
% end
% for x=1:length(CLEAN)
%     x=CLEAN(x);
%     spk=PDS(1).sptimes{x};
%     spk1=spk(find(spk<PDS.timesoffreeoutcomes_second(x)));
%     spk2=spk(find(spk>PDS.timesoffreeoutcomes_second(x)+10/1000));
%     spk=[spk1 spk2]; clear spk1 spk2
%     try
%         no_noise=find(spk>(PDS.timesoffreeoutcomes_second(x)+(10/1000)) & spk<(PDS.timesoffreeoutcomes_second(x)+(20/1000)));
%         no_noise=length(no_noise); %this many spikes in the window  so insert that many spikes in the deleted time range
%         grabnumbers=PDS.timesoffreeoutcomes_second(x)+randperm(10)/1000;
%         grabnumbers=grabnumbers(no_noise);
%         spk=[spk grabnumbers];
%     end
%     PDS(1).sptimes{x}=sort(spk);
%     %
%     clear y spk spk1 spk2 x
%     clc
% end
% CLEAN=[trials6253 trials6259 trials6300];
% for x=1:length(CLEAN)
%     x=CLEAN(x);
%     spk=PDS(1).sptimes{x};
%     spk1=spk(find(spk<PDS.timeoutcome(x)));
%     spk2=spk(find(spk>PDS.timeoutcome(x)+10/1000));
%     spk=[spk1 spk2]; clear spk1 spk2
%     %
%     try
%         no_noise=find(spk>(PDS.timeoutcome(x)+(10/1000)) & spk<(PDS.timeoutcome(x)+(20/1000)));
%         no_noise=length(no_noise); %this many spikes in the window  so insert that many spikes in the deleted time range
%         grabnumbers=PDS.timeoutcome(x)+randperm(10)/1000;
%         grabnumbers=grabnumbers(no_noise);
%         spk=[spk grabnumbers];
%     end
%     PDS(1).sptimes{x}=sort(spk);
%     %
%     clear y spk spk1 spk2 x 
% end
% cleanreward50n=[trials6254nd trials6260nd trials6301nd];
% cleanreward50=[trials6254d trials6260d trials6301d];
% for x=1:length(cleanreward50)
%     x=cleanreward50(x);
%     spk=PDS(1).sptimes{x};
%     spk1=spk(find(spk<PDS.timeoutcome(x)));
%     spk2=spk(find(spk>PDS.timeoutcome(x)+10/1000));
%     spk=[spk1 spk2]; clear spk1 spk2
%     try
%         y=cleanreward50n(randperm(length(cleanreward50n))); y=y(1);
%         y=PDS(1).sptimes{y};
%         
%         y=y(find(y>PDS.timeoutcome(x) & y<PDS.timeoutcome(x)+10/1000));
%         spk=[spk; y];
%     end
%     PDS(1).sptimes{x}=sort(spk);
%     %
%     clear y spk spk1 spk2
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %MAKE A RASTER AND A SPIKE DENSITY FUNCTION FOR EVERY TRIAL for
% %freeoutcomes
% Rasters=[];
% for x=1:length(PDS.fractals)
%     spk=PDS(1).sptimes{x}-PDS.timesoffreeoutcomes_first(x);
% 
%     spk=(spk*1000)+CENTER-1;
%     spk=fix(spk);
%     %
%     spk=spk(find(spk<CENTER*2));
%     %
%     temp(1:CENTER*2)=0;
%     temp(spk)=1;
%     Rasters=[Rasters; temp];
%     %
%     clear temp spk x
% end
% SDFfree_=plot_mean_psth({Rasters},gauswindow_ms,1,(CENTER*2)-2,1); 
% close all;
% Rastersfree_=Rasters;
% 
% Rasters=[];
% for x=1:length(durationsuntilreward)
%     %
%     spk=PDS(1).sptimes{x}-PDS(1).timetargeton(x);
%     spk=(spk*1000)+CENTER-1; 
%     spk=fix(spk);
%     %
%     spk=spk(find(spk<CENTER*2));
%     %
%     temp(1:CENTER*2)=0;
%     temp(spk)=1; 
%     Rasters=[Rasters; temp];
%     %
%     clear temp spk x
% end
% SDFcstrace=plot_mean_psth({Rasters},gauswindow_ms,1,(CENTER*2)-2,1); 
% Rasters_cstrace=Rasters; clear Rasters
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% nsubplot(100,100,1:25,1:25)
% trialstoplot=free1;
% xuncert=[Rastersfree(trialstoplot,limitsarray)];
% SD=nanmean(SDFfree(trialstoplot,limitsarray));
% plot(SD,'r','LineWidth', 2); hold on;
% xt=[];
% rasts=[];
% for tq=1:size(xuncert,1)
%     Z=xuncert(tq,:);
%     Z(find(Z==1))=(find(Z==1));
%     Z(find(Z==0))=NaN;
%     xt_=length(find(Z>0));
%     if isempty(xt_)==1
%         xt_=NaN;
%     end
%     xt=[xt; xt_];
%     rasts=[rasts; Z];
%     clear Z tq xt_
% end
% MatPlot(1:size(xuncert,1),1:max(xt))=NaN
% for tq=1:size(xuncert,1)
%     temptq=find(rasts(tq,:)>0);
%     MatPlot(tq,1:length(temptq))=temptq;
% end
% rastList=MatPlot;
% rasIntv=1;
% LWidth=1;
% LColor='r';
% maxY_rast=tq+10;
% if rastersplot==1
%     for line = 1:size(rastList,1)
%         hold on
%         curY_rast = maxY_rast-rasIntv*line;
%         plot([rastList(line,:); rastList(line,:)],...
%             [(curY_rast+rasIntv/2)*ones(1,size(rastList,2)); ...
%             (curY_rast-rasIntv/2)*ones(1,size(rastList,2))],'-','linewidth',LWidth,'color',LColor );
%         hold on
%     end
% end
% clear rasts xt rastList MatPlot
% 
% trialstoplot=free34;
% xuncert=[Rastersfree(trialstoplot,limitsarray)];
% SD=nanmean(SDFfree(trialstoplot,limitsarray));
% plot(SD,'k','LineWidth', 2); hold on;
% xt=[];
% rasts=[];
% for tq=1:size(xuncert,1)
%     Z=xuncert(tq,:);
%     Z(find(Z==1))=(find(Z==1));
%     Z(find(Z==0))=NaN;
%     xt_=length(find(Z>0));
%     if isempty(xt_)==1
%         xt_=NaN;
%     end
%     xt=[xt; xt_];
%     rasts=[rasts; Z];
%     clear Z tq xt_
% end
% MatPlot(1:size(xuncert,1),1:max(xt))=NaN
% for tq=1:size(xuncert,1)
%     temptq=find(rasts(tq,:)>0);
%     MatPlot(tq,1:length(temptq))=temptq;
% end
% rastList=MatPlot;
% rasIntv=1;
% LWidth=1;
% LColor='k';
% maxY_rast=tq+40;
% if rastersplot==1
%     for line = 1:size(rastList,1)
%         hold on
%         curY_rast = maxY_rast-rasIntv*line;
%         plot([rastList(line,:); rastList(line,:)],...
%             [(curY_rast+rasIntv/2)*ones(1,size(rastList,2)); ...
%             (curY_rast-rasIntv/2)*ones(1,size(rastList,2))],'-','linewidth',LWidth,'color',LColor );
%         hold on
%     end
% end
% clear tq xt line rastList MatPlot curY_rast temptq rasts SD SDplotto xuncert
% x=[1000,1000]; y=[0,150]; plot(x,y,'k'); hold on;
% %xlim([0 2000])
% ylim([0 ylimvar])
% xlabel('free outcomes')
% axis square;
% clear rasts xt rastList MatPlot
% 
% 
% 
% 
% 
% nsubplot(100,100,1:25,31:55)
% trialstoplot=free1_;
% xuncert=[Rastersfree_(trialstoplot,limitsarray)];
% SD=nanmean(SDFfree_(trialstoplot,limitsarray));
% plot(SD,'r','LineWidth', 2); hold on;
% xt=[];
% rasts=[];
% for tq=1:size(xuncert,1)
%     Z=xuncert(tq,:);
%     Z(find(Z==1))=(find(Z==1));
%     Z(find(Z==0))=NaN;
%     xt_=length(find(Z>0));
%     if isempty(xt_)==1
%         xt_=NaN;
%     end
%     xt=[xt; xt_];
%     rasts=[rasts; Z];
%     clear Z tq xt_
% end
% MatPlot(1:size(xuncert,1),1:max(xt))=NaN
% for tq=1:size(xuncert,1)
%     temptq=find(rasts(tq,:)>0);
%     MatPlot(tq,1:length(temptq))=temptq;
% end
% rastList=MatPlot;
% rasIntv=1;
% LWidth=1;
% LColor='r';
% maxY_rast=tq+10;
% if rastersplot==1
%     for line = 1:size(rastList,1)
%         hold on
%         curY_rast = maxY_rast-rasIntv*line;
%         plot([rastList(line,:); rastList(line,:)],...
%             [(curY_rast+rasIntv/2)*ones(1,size(rastList,2)); ...
%             (curY_rast-rasIntv/2)*ones(1,size(rastList,2))],'-','linewidth',LWidth,'color',LColor );
%         hold on
%     end
% end
% clear rasts xt rastList MatPlot
% 
% trialstoplot=free34_;
% xuncert=[Rastersfree_(trialstoplot,limitsarray)];
% SD=nanmean(SDFfree_(trialstoplot,limitsarray));
% plot(SD,'k','LineWidth', 2); hold on;
% xt=[];
% rasts=[];
% for tq=1:size(xuncert,1)
%     Z=xuncert(tq,:);
%     Z(find(Z==1))=(find(Z==1));
%     Z(find(Z==0))=NaN;
%     xt_=length(find(Z>0));
%     if isempty(xt_)==1
%         xt_=NaN;
%     end
%     xt=[xt; xt_];
%     rasts=[rasts; Z];
%     clear Z tq xt_
% end
% MatPlot(1:size(xuncert,1),1:max(xt))=NaN
% for tq=1:size(xuncert,1)
%     temptq=find(rasts(tq,:)>0);
%     MatPlot(tq,1:length(temptq))=temptq;
% end
% rastList=MatPlot;
% rasIntv=1;
% LWidth=1;
% LColor='k';
% maxY_rast=tq+40;
% if rastersplot==1
%     for line = 1:size(rastList,1)
%         hold on
%         curY_rast = maxY_rast-rasIntv*line;
%         plot([rastList(line,:); rastList(line,:)],...
%             [(curY_rast+rasIntv/2)*ones(1,size(rastList,2)); ...
%             (curY_rast-rasIntv/2)*ones(1,size(rastList,2))],'-','linewidth',LWidth,'color',LColor );
%         hold on
%     end
% end
% clear tq xt line rastList MatPlot curY_rast temptq rasts SD SDplotto xuncert
% clear rasts xt rastList MatPlot
% 
% trialstoplot=free2_;
% xuncert=[Rastersfree_(trialstoplot,limitsarray)];
% SD=nanmean(SDFfree_(trialstoplot,limitsarray));
% plot(SD,'b','LineWidth', 2); hold on;
% xt=[];
% rasts=[];
% for tq=1:size(xuncert,1)
%     Z=xuncert(tq,:);
%     Z(find(Z==1))=(find(Z==1));
%     Z(find(Z==0))=NaN;
%     xt_=length(find(Z>0));
%     if isempty(xt_)==1
%         xt_=NaN;
%     end
%     xt=[xt; xt_];
%     rasts=[rasts; Z];
%     clear Z tq xt_
% end
% MatPlot(1:size(xuncert,1),1:max(xt))=NaN
% for tq=1:size(xuncert,1)
%     temptq=find(rasts(tq,:)>0);
%     MatPlot(tq,1:length(temptq))=temptq;
% end
% rastList=MatPlot;
% rasIntv=1;
% LWidth=1;
% LColor='b';
% maxY_rast=tq+70;
% if rastersplot==1
%     for line = 1:size(rastList,1)
%         hold on
%         curY_rast = maxY_rast-rasIntv*line;
%         plot([rastList(line,:); rastList(line,:)],...
%             [(curY_rast+rasIntv/2)*ones(1,size(rastList,2)); ...
%             (curY_rast-rasIntv/2)*ones(1,size(rastList,2))],'-','linewidth',LWidth,'color',LColor );
%         hold on
%     end
% end
% clear tq xt line rastList MatPlot curY_rast temptq rasts SD SDplotto xuncert
% x=[1000,1000]; y=[0,150]; plot(x,y,'k'); hold on;
% x=[1750,1750]; y=[0,150]; plot(x,y,'k'); hold on;
% %xlim([0 2000])
% ylim([0 ylimvar])
% xlabel('free outcomes')
% axis square;
% clear rasts xt rastList MatPlot
% 
% 
% 
% set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figure and orientation
% print('-dpdf', 'timing4' ); 
% 
% sdsdf
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


PDS.fractals(48:58)=NaN;
PDS.fractals(end-6:end)=NaN;

close all;
 

CENTER=6001; 
gauswindow_ms=50;
LINESIZE=3;
plotrange=[5500:12000];
durationsuntilreward=PDS(1).timeoutcome-PDS(1).timetargeton;
durationsuntilreward=round(durationsuntilreward*10)./10;
completedtrial=find(PDS.timetargeton>0);
deliv=find(PDS.rewardDuration>0)
ndeliv=find(PDS.rewardDuration==0)


%organize trial types into indices
trials6201=(find(PDS(1).fractals==6201)); 
trials6201=intersect(trials6201,completedtrial);
trials6201n=intersect(trials6201,deliv); 
trials6201nd=intersect(trials6201,ndeliv); 
%
trials6102=intersect(find(PDS(1).fractals==6102),completedtrial); %50% LONG 50 SHORT
trials6102_50s=intersect(find(durationsuntilreward==1.5),trials6102);
trials6102_50l=intersect(find(durationsuntilreward==4.5),trials6102);
%
trials6101=intersect(find(PDS(1).fractals==6101),completedtrial); %25SHORT
trials6101_25s=intersect(find(durationsuntilreward==1.5),trials6101);
trials6101_75l=intersect(find(durationsuntilreward==4.5),trials6101);
%
trials6103=intersect(find(PDS(1).fractals==6103),completedtrial);%25LONG
trials6103_75s=intersect(find(durationsuntilreward==1.5),trials6103);
trials6103_25l=intersect(find(durationsuntilreward==4.5),trials6103);
%
trials6104=intersect(find(PDS(1).fractals==6104),completedtrial); %SHORT
%
trials6105=intersect(find(PDS(1).fractals==6105),completedtrial);
trials6105_25s=intersect(find(durationsuntilreward==1.5),trials6105);
trials6105_25ms=intersect(find(durationsuntilreward==2.5),trials6105);
trials6105_25ml=intersect(find(durationsuntilreward==3.5),trials6105);
trials6105_25l=intersect(find(durationsuntilreward==4.5),trials6105);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAKE A RASTER AND A SPIKE DENSITY FUNCTION FOR EVERY TRIAL

CLEAN=[trials6101 trials6102 trials6103 trials6104 trials6105];
for x=1:length(CLEAN)
    x=CLEAN(x);
    spk=PDS(1).sptimes{x};
    spk1=spk(find(spk<PDS.timeoutcome(x)));
    spk2=spk(find(spk>PDS.timeoutcome(x)+10/1000));
    spk=[spk1 spk2]; clear spk1 spk2
    %
    try
        no_noise=find(spk>(PDS.timeoutcome(x)+(10/1000)) & spk<(PDS.timeoutcome(x)+(20/1000)));
        no_noise=length(no_noise); %this many spikes in the window  so insert that many spikes in the deleted time range
        grabnumbers=PDS.timeoutcome(x)+randperm(10)/1000;
        grabnumbers=grabnumbers(no_noise);
        spk=[spk grabnumbers];
    end
    PDS(1).sptimes{x}=sort(spk);
    %
    clear y spk spk1 spk2 x 
end
cleanreward50n=trials6201nd;
cleanreward50=trials6201n;
for x=1:length(cleanreward50)
    x=cleanreward50(x);
    spk=PDS(1).sptimes{x};
    spk1=spk(find(spk<PDS.timeoutcome(x)));
    spk2=spk(find(spk>PDS.timeoutcome(x)+10/1000));
    spk=[spk1 spk2]; clear spk1 spk2
    try
        y=cleanreward50n(randperm(length(cleanreward50n))); y=y(1);
        y=PDS(1).sptimes{y};
        
        y=y(find(y>PDS.timeoutcome(x) & y<PDS.timeoutcome(x)+10/1000));
        spk=[spk; y];
    end
    PDS(1).sptimes{x}=sort(spk);
    %
    clear y spk spk1 spk2
end

Rasters=[];
for x=1:length(durationsuntilreward)
    %
    spk=PDS(1).sptimes{x}-PDS(1).timetargeton(x);
    spk=(spk*1000)+CENTER-1; 
    spk=fix(spk);
    %
    spk=spk(find(spk<CENTER*2));
    %
    temp(1:CENTER*2)=0;
    temp(spk)=1; 
    Rasters=[Rasters; temp];
    %
    clear temp spk x
end
SDFtiming=plot_mean_psth({Rasters},gauswindow_ms,1,(CENTER*2)-2,1); 
Rasters_timing=Rasters; clear Rasters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

figure
nsubplot(100,100,1:100,1:100)
h = zeros(5,5);
h(:,1) =plot(nanmean(SDFtiming(trials6101_25s,6000:11000)),'r','LineWidth',LINESIZE); hold on
%h(:,2) =plot(nanmean(SDFtiming(trials6105_25s,6000:11000)),'m','LineWidth',LINESIZE); hold on
h(:,3) =plot(nanmean(SDFtiming(trials6102_50s,6000:11000)),'g','LineWidth',LINESIZE); hold on
%h(:,4) =plot(nanmean(SDFtiming(trials6104,6000:11000)),'y','LineWidth',LINESIZE); hold on
h(:,5) =plot(nanmean(SDFtiming(trials6103_75s,6000:11000)),'b','LineWidth',LINESIZE); hold on
x=[1500,1500]; y=[0,100]; plot(x,y,'k'); hold on;
%ylim([0 50])
xlim([1 2000])
ylim([-10 100])
axis square;

%for testing for buildup
%     y1=trialsregs(:);
%         x=trialsregs_nums(:);
%         [pval, r]=permutation_pair_test_fast(x,y1,permnumber,'rankcorr');
%         P = polyfit(x,y1,1);
%         yfit = P(1)*x+P(2);
%         linepar(zzzz).P34=P;
%         linepar(zzzz).x34=x;
%         linepar(zzzz).yfit34=yfit;
%         linepar(zzzz).pval34=pval;
%         linepar(zzzz).r34=r;

figure
nsubplot(100,100,1:100,1:100)
h = zeros(4,4);
h(:,1) =plot(nanmean(SDFtiming(trials6101_25s,6000:11000))-nanmean(SDFtiming(trials6101_75l,6000:11000)),'r','LineWidth',LINESIZE); hold on
%h(:,2) =plot(nanmean(SDFtiming(trials6105_25s,6000:11000))-nanmean(SDFtiming([trials6105_25ms trials6105_25ml trials6105_25l],6000:11000)),'m','LineWidth',LINESIZE); hold on
h(:,3) =plot(nanmean(SDFtiming(trials6102_50s,6000:11000))-nanmean(SDFtiming(trials6102_50l,6000:11000)),'g','LineWidth',LINESIZE); hold on
h(:,4) =plot(nanmean(SDFtiming(trials6103_75s,6000:11000))-nanmean(SDFtiming(trials6103_25l,6000:11000)),'b','LineWidth',LINESIZE); hold on
h(:,5) =plot(nanmean(SDFtiming(trials6201n,6000:11000))-nanmean(SDFtiming(trials6201nd,6000:11000)),'y','LineWidth',LINESIZE); hold on
x=[1500,1500]; y=[0,100]; plot(x,y,'k'); hold on;
%ylim([0 50])
xlim([0 3000])
ylim([-40 90])
xlim([1 2000])
axis square;




figure
nsubplot(100,100,1:100,1:100)
h = zeros(4,4);
h(:,1) =plot(nanmean(SDFtiming(trials6201n,6000:11000)),'r','LineWidth',LINESIZE); hold on
h(:,2) =plot(nanmean(SDFtiming(trials6201nd,6000:11000)),'b','LineWidth',LINESIZE); hold on
h(:,5) =plot(nanmean(SDFtiming(trials6104,6000:11000)),'k','LineWidth',LINESIZE); hold on
x=[1500,1500]; y=[0,100]; plot(x,y,'k'); hold on;
%ylim([0 50])
xlim([0 3000])
ylim([0 150])
xlim([1 2000])
axis square;












sadf



            figure
            plot(nanmean(SDFtiming(trials6101_25s,6000:11000)),'r'); hold on
            plot(nanmean(SDFtiming(trials6103_75s,6000:11000)),'m');  hold on
            plot(nanmean(SDFtiming(trials6102_50s,6000:11000)),'g');  hold on
            plot(nanmean(SDFtiming(trials6104,6000:11000)),'b'); hold on
            
            
            %%%%analyses on normalized activity (subtracted from non-delivery at short
            %%%%time)
            norm=nanmean(SDFtiming(trials6101_75l,6000:9000));
            A1=SDFtiming(trials6101_25s,6000:9000);
            ch4 = A1 - norm(ones(size(A1,1),1),:);
            %
            norm=nanmean(SDFtiming(trials6102_50l,6000:9000));
            A1=SDFtiming(trials6102_50s,6000:9000);
            ch3 = A1 - norm(ones(size(A1,1),1),:);
            %
            norm=nanmean(SDFtiming(trials6103_25l,6000:9000));
            A1=SDFtiming(trials6103_75s,6000:9000);
            ch2 = A1 - norm(ones(size(A1,1),1),:);
            %ch4=SDFtiming((trials6101_25s),6000:9000);
            %ch3=SDFtiming((trials6102_50s),6000:9000);
            %ch2=SDFtiming((trials6103_75s),6000:9000);
            ch1=SDFtiming((trials6104),6000:9000);
            
            figure
            plot(nanmean(ch2),'k')
            hold on
            plot(nanmean(ch3),'g')
            hold on
            plot(nanmean(ch4),'r')
            %
            %correlation done on a measure of variance in timing. i chose coeffecient
            %of variation as a reasonable one
            adval=1.5
            coeffofvariationCh4=std([0 0 0 3]+adval)./mean([0 0 0 3]+adval)
            coeffofvariationCh3=std([0 0 3 3]+adval)./mean([0 0 3 3]+adval)
            coeffofvariationCh2=std([0 3 3 3]+adval)./mean([0 3 3 3]+adval)
            %
            PvalueSave=[];
            PvalueSaveCor=[];
            SaveCor=[];
            for x = 1:size(ch1,2)
                %         t1=ch1(:,x);%-ch1mean(x);
                %         t1(1:length(t1),2)=coeffofvariationCh1; %1;
                t2=ch2(:,x);%-ch2mean(x);
                t2(1:length(t2),2)=coeffofvariationCh2; %2;
                t3=ch3(:,x);%-ch3mean(x);
                t3(1:length(t3),2)=coeffofvariationCh3; %3;
                t4=ch4(:,x);%-ch4mean(x);
                t4(1:length(t4),2)=coeffofvariationCh4; %4;
                temp=[t2; t3; t4;];
                P=kruskalwallis(temp(:,1),temp(:,2),'off');
                [pval, r]=permutation_pair_test_fast(temp(:,1),temp(:,2),1000,'corr'); %rankcorr for non parametric rank based; corr is linear
                PvalueSaveCor=[PvalueSaveCor; pval];
                PvalueSave=[PvalueSave; P];
                SaveCor=[SaveCor; r];
                clear P t1 t2 t3 t4 t5 t6 temp PlotPvalue pval r
            end
            
            S=SaveCor;
            figure
            S(PvalueSaveCor>0.01)=NaN;
            plot(S);




sdfsdf



%analyses -- my idea would be to perform running stats after we get all the
%sessions find the window in which there is variability (kruskal) and
%perform a correlation across the average activities during the sessions
ch1=SDFtiming((trials6101_25s),6000:9000);
ch2=SDFtiming((trials6105_25s),6000:9000);
ch3=SDFtiming((trials6102_50s),6000:9000);
ch4=SDFtiming((trials6103_75s),6000:9000);
ch5=SDFtiming((trials6104),6000:9000); 
%
ch1mean=nanmean(SDFtiming(trials6101_75l,6000:9000));
ch2mean=nanmean(SDFtiming([trials6105_25ms trials6105_25ml trials6105_25l],6000:9000));
ch3mean=nanmean(SDFtiming(trials6102_50l,6000:9000));
ch4mean=nanmean(SDFtiming(trials6103_25l,6000:9000));
%

PvalueSave=[];
PvalueSaveCor=[];
for x = 1:size(ch1,2)
    t1=ch1(:,x)%-ch1mean(x);
    t1(1:length(t1),2)=1;
    t2=ch2(:,x)%-ch2mean(x);
    t2(1:length(t2),2)=2;
    t3=ch3(:,x)%-ch3mean(x);
    t3(1:length(t3),2)=3;
    t4=ch4(:,x)%-ch4mean(x);
    t4(1:length(t4),2)=4;
    t5=ch5(:,x);
    t5(1:length(t5),2)=5;
    temp=[t1; t2; t3; t4; t5;];
    P=kruskalwallis(temp(:,1),temp(:,2),'off');
    [pval, r]=permutation_pair_test_fast(temp(:,1),temp(:,2),100,'rankcorr');
    PvalueSaveCor=[PvalueSaveCor; pval];
    PvalueSave=[PvalueSave; P]; 
    clear P t1 t2 t3 t4 t5 t6 temp PlotPvalue pval r
end
PlotPvalue(1:3001)=NaN;
PlotPvalue(find(PvalueSaveCor<0.01))=20; 
plot([1:3001],PlotPvalue,'r','LineWidth',2); hold on; 

warning off

% nsubplot(100,100,71:95,71:95)
% h = zeros(4,4);
% h(:,1) =plot(nanmean(SDFtiming(trials6101_25s,6000:11000))-nanmean(SDFtiming(trials6104,6000:11000)),'r','LineWidth',LINESIZE); hold on
% h(:,2) =plot(nanmean(SDFtiming(trials6105_25s,6000:11000))-nanmean(SDFtiming([trials6104],6000:11000)),'m','LineWidth',LINESIZE); hold on
% h(:,3) =plot(nanmean(SDFtiming(trials6102_50s,6000:11000))-nanmean(SDFtiming(trials6104,6000:11000)),'g','LineWidth',LINESIZE); hold on
% h(:,4) =plot(nanmean(SDFtiming(trials6103_75s,6000:11000))-nanmean(SDFtiming(trials6104,6000:11000)),'b','LineWidth',LINESIZE); hold on
% x=[1500,1500]; y=[0,100]; plot(x,y,'k'); hold on;
% %ylim([0 50])
% xlim([0 3000])
% ylim([-20 40])
% axis square;


nsubplot(100,100,41:65,1:35)
plot(nanmean(SDFcstrace(trials6253,8000:8000+3000)),'g','LineWidth',3); hold on
plot(nanmean(SDFcstrace(trials6259,8000:8000+3000)),'r','LineWidth',3); hold on
plot(nanmean(SDFcstrace(trials6300,8000:8000+3000)),'k','LineWidth',3); hold on
x=[2500,2500]; y=[0,100]; plot(x,y,'k'); hold on;
x=[1000,1000]; y=[0,100]; plot(x,y,'k'); hold on;
xlim([0 3000])
axis square;


nsubplot(100,100,71:95,1:35)
plot(nanmean(SDFcstrace(trials6254,8000:8000+2500)),'g','LineWidth',3); hold on
plot(nanmean(SDFcstrace(trials6260,8000:8000+2500)),'r','LineWidth',3); hold on
plot(nanmean(SDFcstrace(trials6301,8000:8000+2500)),'k','LineWidth',3); hold on
x=[2500,2500]; y=[0,100]; plot(x,y,'k'); hold on;
x=[1000,1000]; y=[0,100]; plot(x,y,'k'); hold on;
%plot outcome delivery
plot([2500:3000],nanmean(SDFcstrace(trials6254d,8000+2500:8000+2500+500)),'g','LineWidth',3); hold on
plot([2500:3000],nanmean(SDFcstrace(trials6260d,8000+2500:8000+2500+500)),'r','LineWidth',3); hold on
plot([2500:3000],nanmean(SDFcstrace(trials6301d,8000+2500:8000+2500+500)),'k','LineWidth',3); hold on
%plot outcome not deliver
plot([2500:3000],nanmean(SDFcstrace(trials6254nd,8000+2500:8000+2500+500)),'g','LineWidth',1); hold on
plot([2500:3000],nanmean(SDFcstrace(trials6260nd,8000+2500:8000+2500+500)),'r','LineWidth',1); hold on
plot([2500:3000],nanmean(SDFcstrace(trials6301nd,8000+2500:8000+2500+500)),'k','LineWidth',1); hold on
xlim([0 3000])
axis square;



%this is comparing 1.5 CS buildup with 2.5 CS normalized for reward by
%subtracing

%we may want to do this for trace too and compare teh time courses.. we may
%also need to do this in a separate figure for POST-OUTCOME responses to
%look at that. let's add that to this figure!

nsubplot(100,100,71:95,41:65)
plot(nanmean(SDFcstrace(trials6301,8000:8000+2500))-nanmean(SDFcstrace(trials6300,8000:8000+2500)),'k','LineWidth',2); hold on
%plot(nanmean(SDFcstrace(trials6260,8000:8000+2500))-nanmean(SDFcstrace(trials6259,8000:8000+2500)),'Color',[0.5 0.5 0.5],'LineWidth',1); hold on
%plot(nanmean(SDFcstrace(trials6254,8000:8000+2500))-nanmean(SDFcstrace(tri
%als6253,8000:8000+2500)),'Color',[0.8 0.8 0.8],'LineWidth',1); hold on
plot(nanmean(SDFtiming(trials6201,6000:6000+1500))-nanmean(SDFtiming(trials6104,6000:6000+1500)),'r','LineWidth',2); hold on
plot(nanmean(SDFtiming(trials6201,6000:6000+1500))-nanmean(SDFtiming(trials6104,6000:6000+1500)),'r','LineWidth',2); hold on
xlim([0 2500])
axis square;
x=[1500,1500]; y=[-100,100]; plot(x,y,'r'); hold on;
x=[2500,2500]; y=[-100,100]; plot(x,y,'k'); hold on;
ylim([-20 40])



set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figure and orientation
print('-dpdf', 'timing3' ); 


warning on


January 2015
January 2017

January 2019











