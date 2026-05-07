%% startup.m
%%
%%   JP Sep 2001

   %% mount lab::\\win49\E            on T:
   %% mount josef:nikos02::\\win58\Y  on Y:
nklDir = 'Y:/mri/'
%nklDir = 'D:/matlab/nikos/mri/'      %%% use, if Y::josef@win58 is down
addpath( strcat(nklDir, 'AdfSrc/adf_read'))
addpath( strcat(nklDir, 'dg_read'))
addpath( strcat(nklDir, 'MatLab/utils'))
addpath( strcat(nklDir, 'MatLab/sml'))
addpath( strcat(nklDir, 'MatLab/scat'))
addpath( strcat(nklDir, 'MatLab/pub'))
addpath( strcat(nklDir, 'MatLab/neu'))
addpath( strcat(nklDir, 'MatLab/erfmri'))
addpath( strcat(nklDir, 'MatLab/mri'))
addpath( strcat(nklDir, 'MatLab/figs'))
addpath( strcat(nklDir, 'MatLab/docs'))
addpath( strcat(nklDir, 'MatLab/anesth'))
addpath( strcat(nklDir, 'MatLab/ana'))
%%addpath D:/matlab/matlab.pfvdm/newfun
%%addpath D:/matlab/packages/nmrlab
%%addpath D:/matlab/packages/matpulse-2.4
%%addpath D:/matlab/packages/mseq
%%addpath D:/matlab/packages/spm99
addpath D:/matlab/libnkl    %% nkl routines modified by JP
addpath D:/matlab/lib
addpath D:/matlab

%Old CMRR Paths
%path(path,'/home/raid77/ku_grp/pfeuffer/matlab/WW_filter/');
%path(path,'/home/raid77/ku_grp/pfeuffer/matlab/spm99/');
%path(path,'/home/raid77/ku_grp/pfeuffer/matlab/lyngby/');
%path(path,'/home/raid77/ku_grp/pfeuffer/matlab/physiofix/');
%path(path,'/home/raid77/ku_grp/pfeuffer/matlab/pbanova/');
%path(path,'/home/raid77/ku_grp/pfeuffer/matlab/Functions_000110/');
%path(path,'/home/raid77/ku_grp/pfeuffer/matlab/lib/');
%path(path,'/usr/local/cmrr/m_src/');
%path(path,'/home/raid44/hu_grp/pfvdm/matlab/newfun/');
%path(path,'/home/raid133/hu_grp/shantanu/SPIRALS/4T/');

%%path(path,'/export/raid53/bill/Prog/')
%%path(path,'/export/raid53/bill/Prog/Functions/')
%%path(path,'/export/raid53/bill/Prog/Pbanova/Pb_new/')

%WavePath;      % WaveLab paths
%%spiralsort;   % test for spiralsetup
stdpath;       %
cd('data');
curdir = pwd
figure(1)
set(1,'Position',[58 371   470  260]);
colormap(gray);
figure(3)
set(3,'Position',[58  35   470  260]);
colormap(gray);

