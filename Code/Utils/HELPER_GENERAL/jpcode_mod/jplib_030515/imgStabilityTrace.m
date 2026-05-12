%--------------------------------------------------------------
function [F, F134] = imgStabilityTrace(dir, filenum, indSignal, thresImg)
%
%  image series: stability testing
%
%  JP July 2002
%
%global STDPATH
%STDPATH.pv = 'M:/mridata/'     % wks4
%imgStabilityTrace('tgntest.d21', 5, [64 64], 2e4);   % #

global imgdat    %%% assumed to be 3D: (imgx,imgy,imgnum)
global acqp

flagReadDat = 1;
refImgNum = 1;
%%thresImg = 1e4;
   % parameters for SNR0
sizeSNR0 = 20;
indBackground = [1 1];  % start indices
%%indSignal = [64 64];    % center indices of image signal
DEV = 0;     % arbitrary cutoff to eliminate voxels with mean+/-DEV
ThresholdBad = 1.0;    % [%] specs for STD

if (flagReadDat | thresImg <= 0)
    img = PVrd2dseq(dir, filenum);
    s_d = size(img);
    imgdat = reshape(double(img(:,:,1,:)), s_d(1), s_d(2), s_d(4));
    refImg = imgdat(:,:,refImgNum);
    imgStat = imgStatistics(refImg);
else
    subplot(1,1,1)
    refImg = imgdat(:,:,refImgNum);
end

figure(1)
SNR0 = calcSNR0(indSignal, indBackground, sizeSNR0);
Fluctuation = FluctuationROI(indSignal, 4, 1);

ImgThresholded = (refImg > thresImg);
s_d = size(ImgThresholded);
figure(3)
imagesc(ImgThresholded);
axis image

fprintf('... calculation fluctuations ROI dependent ... ');
flagVerbose = 0;
nROI = 3;
F134 = zeros(s_d(1), s_d(2), nROI);
F = zeros(nROI,3);    % Ftheoretical, mean, std: statistics over image
for imgy=1:s_d(2)
 for imgx=1:s_d(1)
%for imgy=40:70
% for imgx=40:70
   if (ImgThresholded(imgx,imgy) == 1)
       F134(imgx,imgy,1) = FluctuationROI([imgx imgy], 1, flagVerbose);
       F134(imgx,imgy,2) = FluctuationROI([imgx imgy], 3, flagVerbose);
       F134(imgx,imgy,3) = FluctuationROI([imgx imgy], 4, flagVerbose);
   else  % do nothing
   end
 end
end
fprintf('done\n');

datF134 = reshape( F134, s_d(1)*s_d(2), nROI);
datInd = find(reshape( ImgThresholded, s_d(1)*s_d(2), 1 ) == 1);
datF134 = datF134(datInd,:);
mean0 = mean(datF134,1);
std0 = std(datF134,1);

   % eliminate voxels with > DEV
if (DEV > 0)
    datInd = find(datF134(:,1) < mean0(1,1)+DEV & datF134(:,1) > mean0(1,1)-DEV);
    datF134 = datF134(datInd,:);
    mean1 = mean(datF134,1);
    std1 = std(datF134,1);
else
    mean1 = mean0;
    std1 = std0;
end
nVoxelThres = length(datInd);
minval = min(datF134)
maxval = max(datF134)
nVoxelBad(1) = sum(datF134(:,1) > ThresholdBad);
nVoxelBad(2) = sum(datF134(:,2) > ThresholdBad);
nVoxelBad(3) = sum(datF134(:,3) > ThresholdBad)

F(:,1) = 100/SNR0*[1 1/3 1/4]';
F(:,2) = mean1';
F(:,3) = std1';

figure(1)
plotimg = F134(:,:,1);
subplot(3,1,1)
plot(plotimg(s_d(1)/2,:))
ylim([0 1.5])
xlabel('profile X');
str = sprintf('stabTest Statistics (n=%d):  SNR = %.0f,   Ft(1,3,4) = %.3f%%, %.3f%%, %.3f%%\nF(1) = %.2f%% (%.2f), F(3) = %.2f%% (%.2f), F(4) = %.2f%% (%.2f)', ...
    nVoxelThres, SNR0, F(1,1),F(2,1),F(3,1), F(1,2), F(1,3), F(2,2), F(2,3), F(3,2), F(3,3) );
title(str);
subplot(3,1,2)
plot(plotimg(:,s_d(2)/2))
ylim([0 1])
xlabel('profile Y');
subplot(3,1,3)
plot(datF134(:,1))
xlabel(sprintf('all fluctuations\n<%s> exp#%d  TR=%.1f s', dir, filenum, acqp.IMND_rep_time));

return
%--------------------------------------------------------------
function Weisskoff(dir, filenum, flagReadDat)

global imgdat    %%% assumed to be 3D: (imgx,imgy,imgnum)

refImgNum = 1;
%%thresImg = 1e4;
   % parameters for SNR0
sizeSNR0 = 20;
indBackground = [1 1];  % start indices
indSignal = [64 64];    % center indices of image signal

if (flagReadDat)
    img = PVrd2dseq(dir, filenum);
    s_d = size(img);
    imgdat = reshape(double(img(:,:,1,:)), s_d(1), s_d(2), s_d(4));
    refImg = imgdat(:,:,refImgNum);
    imgStat = imgStatistics(refImg);
    return
end

%SNR0 = calcSNR0(indSignal, indBackground, sizeSNR0);
%Fluctuation = FluctuationROI(indSignal, 4, 1);
subplot(1,1,1)
indROI = [1:10];
[Fn,Fnt,Fnt0,n] = WeisskoffTest(indROI, indSignal, indBackground, sizeSNR0);

return
%--------------------------------------------------------------
function [Fn,Fnt,Fnt0,n] = WeisskoffTest(indROI, indSignal, indBackground, sizeSNR0)
%
% Fn calculation for a series of ROI according Weisskoff MRM 36:643(1996)

global imgdat    %%% assumed to be 3D: (imgx,imgy,imgnum)

n = fix(indROI(fix(indROI)>=1));  % >=1
Fn = zeros(length(n),1);
Fnt = zeros(length(n),1);
Fnt0 = zeros(length(n),1);

SNR0 = calcSNR0(indSignal, indBackground, sizeSNR0);

for in=1:length(n)
    Fn(in) = FluctuationROI(indSignal, n(in), 1);
    Fnt(in) = 1/n(in)*Fn(1);
    Fnt0(in) = 100/(n(in)*SNR0);   %% rel. deviation (%)
end

loglog(n,Fnt)
ylim([0.01 1])
hold on
loglog(n,Fn,'--s')
loglog(n,Fnt0,'.')
hold off

return
%--------------------------------------------------------------
function Fluctuation = FluctuationROI(indSignal, sizeF, flagVerbose)
%
% Fn calculation according Weisskoff MRM 36:643(1996)
%     relative deviation (%)

global imgdat    %%% assumed to be 3D: (imgx,imgy,imgnum)

s_d = size(imgdat);
imgNum = s_d(3);
indSg = fix([indSignal(1)-sizeF/2+1 indSignal(2)-sizeF/2+1]);
indSg = [indSg(1) indSg(2) indSg(1)+sizeF-1 indSg(2)+sizeF-1];
indSg(1) = max([indSg(1) 1]);
indSg(2) = max([indSg(2) 1]);
indSg(3) = min([indSg(3) s_d(1)]);
indSg(4) = min([indSg(4) s_d(2)]);
sizeNew = (indSg(3)-indSg(1)+1) * (indSg(4)-indSg(2)+1);
arrSg = reshape( imgdat(indSg(1):indSg(3),indSg(2):indSg(4),:), sizeNew, imgNum);

meanI = mean(arrSg,1);                    % average over n.n voxels in single image
meanN = 1/(imgNum)*sum(meanI,2);   % average over N images
Fluctuation =  100*sqrt(1/(imgNum-1)*sum((meanI-meanN).^2))/meanN;

if flagVerbose
    fprintf('Fn(%2d) = %.2f%%:  Sg[%d %d %d %d]\n', sizeF, Fluctuation, indSg );
end

return
%--------------------------------------------------------------
function SNR0 = calcSNR0(indSignal, indBackground, sizeSNR0);
%
% SNR0 calculation according Weisskoff MRM 36:643(1996)

global imgdat    %%% assumed to be 3D: (imgx,imgy,imgnum)

flagVerbose = 1;

   % parameters for SNR0
%%sizeSNR0 = 20;
%%indBackground = [1 1];  % startIndices
%%indSignal = [64 64];    % center Indices

s_d = size(imgdat);
imgNum = s_d(3);
indBg = [indBackground(1) indBackground(2) indBackground(1)+sizeSNR0-1 indBackground(2)+sizeSNR0-1];
indSg = fix([indSignal(1)-sizeSNR0/2+1 indSignal(2)-sizeSNR0/2+1]);
indSg = [indSg(1) indSg(2) indSg(1)+sizeSNR0-1 indSg(2)+sizeSNR0-1];
arrBg = reshape( imgdat(indBg(1):indBg(3),indBg(2):indBg(4),:), sizeSNR0*sizeSNR0, imgNum);
arrSg = reshape( imgdat(indSg(1):indSg(3),indSg(2):indSg(4),:), sizeSNR0*sizeSNR0, imgNum);

SNR =  sum( mean(arrSg,1) ) / sum( mean(arrBg,1) );
SNR0 = sum( mean(arrSg,1) ) / (1.53* sum( std(arrBg,0,1) ));

if flagVerbose
    fprintf('SNR0 = %.1f, SNR = %.1f:  Sg[%d %d %d %d], Bg[%d %d %d %d], N=%d \n',...
        SNR0, SNR, indSg, indBg, sizeSNR0 );
end

return
%--------------------------------------------------------------
function imgStat = imgStatistics(img)

s_d = size(img);
img1D = reshape(img,prod(size(img)),1);
Imin = min(img1D);
Imax = max(img1D);

figure(3)
imagesc(img);
axis square

figure(2)
subplot(4,1,1)
plot(img(s_d(1)/2,:))
xlabel('profile X');
title(sprintf('refImg Statistics: min/max [%f,%f]',Imin, Imax))
subplot(4,1,2)
plot(img(:,s_d(2)/2))
xlabel('profile Y');
subplot(4,1,3)
plot(img(1,:))
xlabel('first line X');
subplot(4,1,4)
plot(img(:,1))
xlabel('first line Y');

imgStat = 0;

return
%--------------------------------------------------------------
function readSteData

dir = 'E:/STEtest/'

cutIndex = 10050;
cutIndex2 = cutIndex+30000;
dx0 = load([dir 'stex0.csv']);
dx1 = load([dir 'stex1.csv']);
dy0 = load([dir 'stey0.csv']);
dy1 = load([dir 'stey1.csv']);
dz0 = load([dir 'stez0.csv']);
dz1 = load([dir 'stez1.csv']);

dt = (dx0(2,1)-dx0(1,1))*1e3;    % ms
dx0 =     dx0(cutIndex+1:cutIndex2,2) - mean(dx0(1:cutIndex,2));
dy0 = -1*(dy0(cutIndex+1:cutIndex2,2) - mean(dy0(1:cutIndex,2)));
dz0 =     dz0(cutIndex+1:cutIndex2,2) - mean(dz0(1:cutIndex,2));
dx1 =     dx1(cutIndex+1:cutIndex2,2) - mean(dx1(1:cutIndex,2));
dy1 = -1*(dy1(cutIndex+1:cutIndex2,2) - mean(dy1(1:cutIndex,2)));
dz1 =     dz1(cutIndex+1:cutIndex2,2) - mean(dz1(1:cutIndex,2));
dt = ((1:length(dx0))*dt)';
subplot(3,1,1)
plot(dt,dx0)
subplot(3,1,2)
plot(dt,dy0)
subplot(3,1,3)
plot(dt,dz0)

%save('data/steTest.mat','dt','dx0','dy0','dz0','dx1','dy1','dz1');
stop
return
%--------------------------------------------------------------
function steTest

load('data/steTest.mat');

plot(dx0)
hold on
plot(dx1(14:end),'.')   % X:13us late 
hold off

plot(dy0)
hold on
plot(dy1(14:end),'.')   % Y:13us late 
hold off

plot(dz0)
hold on
plot(dz1(14:end),'.')   % Z:us late 
hold off

whos
return
%--------------------------------------------------------------