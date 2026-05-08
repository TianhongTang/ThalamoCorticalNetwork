clc; clear all;
addpath('HELPER_GENERAL'); %keep helpers functions here
close all; clear all;
CENTER=6001;
gauswindow_ms=50;
LINESIZE=3;
plotrange=[5500:12000];

%sort by AP / ML

[raw,string,xls1]=xlsread('Y:\Noah\WolverineData\Wolverine_Filelist_ActiveSheet.xls');
type2ornot= xls1(:,9);
unitNumber= xls1(:,16); unitNumber{1,1}=NaN;
%pull all timing procedures
indexTimingProcedure=[];
indexTimingTrace=[];
indexType2=[];
neuronID=[];
for x=1:length(string(:,1));
    tempname=string{x,1};
    tempneurontype=type2ornot{x,1}
    tempunit=unitNumber{x,1}
    
    if tempunit>0
        neuronID=[neuronID; tempunit];
    else
        neuronID=[neuronID; NaN];
    end
    
    
    
    try
        if strcmp(tempname(1:7),'TimingP')==1;
            indexTimingProcedure=[indexTimingProcedure; x];
        end
    end
    
    try
        if strcmp(tempname(1:7),'TimingT')==1;
            indexTimingTrace=[indexTimingTrace; x];
        end
    end
    
    try
        if tempneurontype==1;
            indexType2=[indexType2; x];
        end
    end
end


timeproceduretype2=intersect(indexType2,indexTimingProcedure)
unitID=[unitNumber{timeproceduretype2}]'

for xzv=1:length(timeproceduretype2)
   unitsave=unitID(xzv);
    try
        filename=string{timeproceduretype2(xzv),1}
        directory=string{timeproceduretype2(xzv),2}
        
        load(['y:\' directory '\' filename])
        
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
        PDSt=PDS;
        
        if isempty(trials6201)==1
            thisneuronID=unitNumber{timeproceduretype2(xzv),1};
            traceline=intersect(find(neuronID==14),indexTimingTrace);
            traceline=traceline(end);
            if ~isempty(traceline)==1
            %
            filenameT=string{(traceline),1}
            directoryT=string{(traceline),2}
            PDS1=load(['y:\' directory '\' filename]); clear c s
            PDSt=PDS1.PDS
            %
            trials6201=(find(PDSt(1).fractals==6201));
            trials6201=intersect(trials6201,completedtrial);
            trials6201n=intersect(trials6201,deliv);
            trials6201nd=intersect(trials6201,ndeliv);
            else
               pause; 
            end
        end
        
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
        %CLEAR REWARD NOISE
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        try
            checkthis=[length(trials6101_25s) length(trials6101_75l) length(trials6103_75s) length(trials6103_25l) length(trials6105_25s) length(trials6105_25ms) length(trials6105_25ml) length(trials6105_25l)]
            checkthis1=[isempty(trials6101_25s) isempty(trials6101_75l) isempty(trials6103_75s) isempty(trials6103_25l) isempty(trials6105_25s) isempty(trials6105_25ms) isempty(trials6105_25ml) isempty(trials6105_25l)]
            checkthis=[length(trials6101_25s) length(trials6101_75l) length(trials6103_75s) length(trials6103_25l)]
            checkthis1=[isempty(trials6101_25s) isempty(trials6101_75l) isempty(trials6103_75s) isempty(trials6103_25l)]
            
            if isempty(find(checkthis<2))==1 && isempty(find(checkthis1==1))==1
                dontsave=1;
            else
                dontsave=0;
            end
        catch
            dontsave=0;
        end
        
        if dontsave==1
         
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
            
            Rasters=[];
            for x=1:length(durationsuntilreward)
                %
                spk=PDS(1).sptimes{x}-PDSt(1).timetargeton(x);
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
            SDFtrace=plot_mean_psth({Rasters},gauswindow_ms,1,(CENTER*2)-2,1);
            Rasters_trace=Rasters; clear Rasters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
               
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
            S(PvalueSaveCor>0.05)=NaN;
            plot(S);
            close all;
            
            
%             %DOES THE INCREASED OVERALL VARIANCE MODIFY ACTIVITY AT THE FIRST 25%
%             %CHANCE OF REWARD?
%             norm=nanmean(SDFtiming(trials6101_75l,6000:9000));
%             A1=SDFtiming(trials6101_25s,6000:9000);
%             ch1 = A1 - norm(ones(size(A1,1),1),:);
%             
%             norm=nanmean(SDFtiming([trials6105_25ms trials6105_25ml trials6105_25l],6000:9000));
%             A2=SDFtiming(trials6105_25s,6000:9000);
%             ch2 = A2 - norm(ones(size(A2,1),1),:);
%             
%             [roc, proc]=rocarea3(A1,A2);
            
            
            if dontsave==1
                figure; 
                plot(nanmean(SDFtiming(trials6104,6000:11000)))
                hold on
                plot(nanmean(SDFtiming([trials6101_25s trials6102_50s trials6103_75s],6000:11000)),'r')
                hold on; plot(nanmean(SDFtiming([trials6101_75l trials6102_50l trials6103_25l],6000:11000)),'g')
                x=[1500,1500]; y=[0,60]; plot(x,y,'k'); hold on;
                close all;
                
                %%%SAVE THE DATA  -- all the average spike density functions
                savestruct(xzv).s6104=nanmean(SDFtiming(trials6104,6000:11000));
                
                savestruct(xzv).s6201d=nanmean(SDFtrace(trials6201n,6000:11000));
                savestruct(xzv).s6201nd=nanmean(SDFtrace(trials6201nd,6000:11000));
                
                savestruct(xzv).sTimingD=nanmean(SDFtiming([trials6101_25s trials6102_50s trials6103_75s],6000:11000));
                savestruct(xzv).sTimingND=nanmean(SDFtiming([trials6101_75l trials6102_50l trials6103_25l],6000:11000));
                
                savestruct(xzv).s6101_25s=nanmean(SDFtiming(trials6101_25s,6000:11000));
                savestruct(xzv).s6102_50s=nanmean(SDFtiming(trials6102_50s,6000:11000));
                savestruct(xzv).s6103_75s=nanmean(SDFtiming(trials6103_75s,6000:11000));
                
                savestruct(xzv).s6101_75l=nanmean(SDFtiming(trials6101_75l,6000:11000));
                savestruct(xzv).s6102_50l=nanmean(SDFtiming(trials6102_50l,6000:11000));
                savestruct(xzv).s6103_25l=nanmean(SDFtiming(trials6103_25l,6000:11000));
                
%                 savestruct(xzv).s6105_25s=nanmean(SDFtiming(trials6105_25s,6000:11000));
%                 savestruct(xzv).s6105_25ms=nanmean(SDFtiming(trials6105_25ms,6000:11000));
%                 savestruct(xzv).s6105_25ml=nanmean(SDFtiming(trials6105_25ml,6000:11000));
%                 savestruct(xzv).s6105_25l=nanmean(SDFtiming(trials6105_25l,6000:11000));
%                 
                
   %             savestruct(xzv).ROC_6105VS6101=roc;
   %             savestruct(xzv).P_6105VS6101=proc;
                savestruct(xzv).PvalueSaveCor=PvalueSaveCor;
                savestruct(xzv).PvalueSave=PvalueSave;
                savestruct(xzv).SaveCor=SaveCor;
                
                savestruct(xzv).filename=filename;
                savestruct(xzv).directory=directory;
                savestruct(xzv).unitsave=unitsave;
                
                
                
            end
        end
        
    end
    clc; clear PDS roc proc A1 A2 ch1 ch2 ch3 ch3
end

save timingproceduresummary1.mat savestruct
clc; clear all;
 

clc; clear all; close all; 
%load('giftforilya.mat')
load('timingproceduresummary1.mat')
linew=1.5;
linew1=0.5;


s6102_50s=[];
s6102_50l=[];
s6104=[];
s6103_75s=[];
s6103_25l=[];
s6101_25s=[];
s6101_75l=[];
PvalueSaveCor=[];
PvalueSave=[];
SaveCor=[];
ROC_6105VS6101=[];
P_6105VS6101=[];
s6105_25ms=[];
s6105_25ml=[];
s6105_25l=[];
s6105_25s=[];

%              s6104: [1x5001 double]
%          s6101_25s: [1x5001 double]
%          s6102_50s: [1x5001 double]
%          s6103_75s: [1x5001 double]
%          s6101_75l: [1x5001 double]
%          s6102_50l: [1x5001 double]
%          s6103_25l: [1x5001 double]
%          s6105_25s: [1x5001 double]
%         s6105_25ms: [1x5001 double]
%         s6105_25ml: [1x5001 double]
%          s6105_25l: [1x5001 double]
%     ROC_6105VS6101: [1x3001 double]
%       P_6105VS6101: [1x3001 double]
%      PvalueSaveCor: [3001x1 double]
%         PvalueSave: [3001x1 double]
%            SaveCor: [3001x1 double]
%           filename: 'TimingProcedure_12_04_2015_11_13.mat'
%          directory: 'Noah\WolverineData\TimingProcedure\oldTimingProcedure'


for x=1:21
    try

        s6102_50s=[s6102_50s;savestruct(x).s6102_50s];
        s6102_50l=[s6102_50l;savestruct(x).s6102_50l];
        s6104=[s6104;savestruct(x).s6104];
        s6103_75s=[s6103_75s;savestruct(x).s6103_75s];
        s6103_25l=[s6103_25l;savestruct(x).s6103_25l];
        s6101_25s=[s6101_25s;savestruct(x).s6101_25s];
        s6101_75l=[s6101_75l;savestruct(x).s6101_75l];
        SaveCor=[SaveCor;savestruct(x).SaveCor'];
        PvalueSave=[PvalueSave;savestruct(x).PvalueSave'];
        PvalueSaveCor=[PvalueSaveCor;savestruct(x).PvalueSaveCor'];
%         s6105_25ms=[s6105_25ms;savestruct(x).s6105_25ms];
%         s6105_25ml=[s6105_25ml;savestruct(x).s6105_25ml];
%         s6105_25l=[s6105_25l;savestruct(x).s6105_25l];
%         s6105_25s=[s6105_25s;savestruct(x).s6105_25s];

    end
end

figure
nsubplot(160,160,1:25,1:25)
plot(nanmean(s6101_75l),'r','LineWidth',linew1);
hold on
plot(nanmean(s6101_25s),'r','LineWidth',linew);
hold on
x=[1500,1500]; y=[0,60]; plot(x,y,'k'); hold on;
axis([0 2250 20 80])
axis square

nsubplot(160,160,1:25,31:55)
plot(nanmean(s6102_50s),'g','LineWidth',linew);
hold on
plot(nanmean(s6102_50l),'g','LineWidth',linew1);
hold on
x=[1500,1500]; y=[0,60]; plot(x,y,'k'); hold on;
axis([0 2250 20 80])
axis square

nsubplot(160,160,1:25,61:85)
plot(nanmean(s6103_75s),'b','LineWidth',linew);
hold on
plot(nanmean(s6103_25l),'b','LineWidth',linew1);
hold on
x=[1500,1500]; y=[0,60]; plot(x,y,'k'); hold on;
axis([0 2250 20 80])
axis square

nsubplot(160,160,1:25,91:115)
plot(nanmean(s6104),'k','LineWidth',linew);
x=[1500,1500]; y=[0,60]; plot(x,y,'k'); hold on;
axis([0 2250 20 80])
axis square


nsubplot(160,160,31:55,1:25)
plot(nanmean(s6101_25s)-nanmean(s6101_75l),'r','LineWidth',linew); hold on
plot(nanmean(s6102_50s)-nanmean(s6102_50l),'g','LineWidth',linew); hold on
plot(nanmean(s6103_75s)-nanmean(s6103_25l),'b','LineWidth',linew); hold on
x=[1500,1500]; y=[-20,60]; plot(x,y,'k'); hold on;
axis([0 2250 -20 30])
axis square

nsubplot(160,160,61:85,1:25)
plot(nanmean(SaveCor),'k','LineWidth',linew);
hold on; x=[1500,1500]; y=[0,1]; plot(x,y,'k'); hold on;
axis([0 2250 -0.2 0.4])
axis square
title('average r')

nsubplot(160,160,91:115,1:25)
PvalueSaveCor(find(PvalueSaveCor>0.01))=9; PvalueSaveCor(find(PvalueSaveCor~=9))=1;  PvalueSaveCor(find(PvalueSaveCor==9))=0;
plot(sum(PvalueSaveCor)./size(PvalueSaveCor,1),'k','LineWidth',linew);
hold on; x=[1500,1500]; y=[0,1]; plot(x,y,'k'); hold on;
axis([0 2250 0 0.7])
axis square
title('% cor significant')

nsubplot(160,160,121:145,1:25)
S=SaveCor;
S(PvalueSaveCor==0)=0;
image(colormapify(S,[-0.7 0.7],'b','w','r')); hold on
hold on; x=[1500,1500]; y=[0,20]; plot(x,y,'k','LineWidth',linew); hold on;
xlim([0 2250])
axis square
title('sig r for each neuron in time')



set(gcf,'Position',[1 41 2560 1484],'Paperposition',[0 0 26.6667 15.4583], 'Paperpositionmode','auto','Papersize',[26.6667 15.4583]);  % sets the size of the figure and orientation
print('-dpdf', 'timing1' ); 


ggg
figure
PvalueSave(find(PvalueSave>0.01))=9; PvalueSave(find(PvalueSave~=9))=1;  PvalueSave(find(PvalueSave==9))=0;
plot(sum(PvalueSave)./size(PvalueSave,1),'k','LineWidth',linew);
hold on; x=[1500,1500]; y=[0,1]; plot(x,y,'k'); hold on;
axis([0 2250 0 0.7])
axis square



jhgjhg
figure
temporary=([.3 .31 .32 .33 .34 .35 .36 .37 .38 .39 .4 .41 .42 .43 .44 .45 .46 .47 .48 .49 .5 .51 .52 .53 .54 .55 .56 .57 .58 .59 .6 .61 .62 .63 .64 .65 .66 .67 .68 .69 .7]);
image(colormapify(temporary,[.3 .7],'b','w','r')); clear temporary
axis square;





