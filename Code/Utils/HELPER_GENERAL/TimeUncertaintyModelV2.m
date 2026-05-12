%DEMO OF UNCERTAINTY AND TIME
%Ilya E Monosov 2015
clear all; clc; close all;

ShannonEntropy=-1*(.5*log(.5)); %entropy measure (nats) / but you can also set this to 1 and do it in bits (maximum uncertainty)
DIST=std([0.5 0.5 0.5 1])%./mean([0 0 0 1]);
DIST1=std([0.5 0.5 1 1])%./mean([0 0 1 1]);
DIST2=std([0.5 1 1 1])%./mean([0 0 1 1]);
VAL=mean([0.5 0.5 0.5 1])%./mean([0 0 0 1]);
VAL1=mean([0.5 0.5 1 1])%./mean([0 0 1 1]);
VAL2=mean([0.5 1 1 1])%./mean([0 0 1 1]);




iti=1; %we set this to 1 because the animal can't predict the time of ITI.. and so its baseline (this may change depending our task)
ts=1; %trial start fixation (let's say he does not have to fixate.. this is wrong if he does..because fixation takes away temporal uncertainty)
cs=1.5; %seconds CS stimulus 

%calculate the temporal parametsrs of the entire task in milliseconds
totalUncertain=(iti+ts+cs)*1000; %total for uncertain CS until reinforcement (outcome or no outcome)
totalCertain=(iti+cs+iti+ts+cs)*1000; %this is total for certain CS until the next reinforcement 
totalCertain=fliplr([1:totalCertain]);
%
total=totalUncertain;
total=fliplr([1:total]);
%

temp(1:length(total))=DIST;
UncertainCS=temp./total; clear temp
%temp(1:length(total))=VAL;
%plotValueCS=temp./(1:totalUncertain); clear temp
plotValueCS=VAL;
plot((UncertainCS-plotValueCS),'b', 'LineWidth', 1); hold on; clear temp

temp(1:length(total))=DIST1;
UncertainCS=temp./total; clear temp
%temp(1:length(total))=VAL1;
%plotValueCS=temp./(1:totalUncertain); clear temp
plotValueCS=VAL1;
plot((UncertainCS-plotValueCS),'m', 'LineWidth', 1); hold on; clear temp

sdf
ylim([-0.01 0.01])
xlim([0 3500])

sdfds

plot(plotUncertainCS,'b', 'LineWidth', 1); hold on; clear temp
%
temp(1:length(total))=DIST1;
plotUncertainCS=temp./total; clear temp
plot(plotUncertainCS,'m', 'LineWidth', 1); hold on; clear temp
%
temp(1:length(total))=DIST2;
plotUncertainCS=temp./total; clear temp
plot(plotUncertainCS,'r', 'LineWidth', 1); hold on; clear temp

ylim([0 0.01])
xlim([0 3500])

sdfdsfsd



temp(1:length(totalCertain))=DIST1;
plotCertainCS=temp./totalCertain; clear temp
plot([2000:4000],plotCertainCS(1:2001),'k-.','LineWidth', 1); hold on
hold on
plot([1:1999],plotUncertainCS(1:1999),'k-.','LineWidth', 1); hold on
ylim([0 0.01])
xlim([0 3500])

sdffsd
figure
plot(plotUncertainCS,'r', 'LineWidth', 3)
hold on
plot([2000:4000],plotCertainCS(1:2001),'k-.','LineWidth', 1)
hold on
plot([1:1999],plotUncertainCS(1:1999),'k-.','LineWidth', 1)
ylim([0 0.01])
xlim([0 3500])






total=(iti+ts+cs)*1000; %total for uncertain CS until reinforcement (outcome or no outcome)
totalCertain=(cs+iti+ts+cs)*1000; %this is total for certain CS until the next reinforcement 
total=fliplr([1:total]);
totalCertain=fliplr([1:totalCertain]);
temp(1:length(total))=DIST2;
plotUncertainCS=temp./total; clear temp
temp(1:length(totalCertain))=DIST2;
plotCertainCS=temp./totalCertain; clear temp
hold on
plot(plotUncertainCS,'g', 'LineWidth', 3)
hold on
plot([2000:4000],plotCertainCS(1:2001),'k-.','LineWidth', 1)
hold on
plot([1:1999],plotUncertainCS(1:1999),'k-.','LineWidth', 1)
ylim([0 0.01])
xlim([0 3500])





% %%%%%%%%%
% %this would be the case if we did not see the buildup...e.g. context
% uncertainty without time

% ITIVAL=ShannonEntropy/(iti+ts+cs);
% TS=ShannonEntropy/(ts+cs);
% CS50=ShannonEntropy/(cs);
% CS0_100=ShannonEntropy/(cs+iti+ts+cs); %this is because after 0 or 100 CS the time to the next reward is longer by the CS duration + iti_ts_cs period of the next trial
% figure
% plot([1:1500],ITIVAL)
% hold on
% plot([1500:2500],TS)
% hold on
% plot([2500:4000],CS50)
% hold on
% plot([2500:4000],CS0_100,'r')
% hold on
% plot([4000:5500],ITIVAL)
%%%%%%%%
