function alignPeekachoooScans()
%procedure for AFNI for aligning MRI volumes for MN++ tracing -Ilya Monosov
%
%
clear all; close all; warning off; clc; clear all;

%setup path so we can run AFNI with matlab
PATH = getenv('PATH');
setenv('PATH', [PATH ':/Users/monosovi/Desktop/AFNI/']);
setenv DYLD_FALLBACK_LIBRARY_PATH ~/abin
unix('#!/bin/tcsh');
unix('# align_times.csh');



%%%%%%%%%%PREPARE BASELINE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%take BASELINE DAY (no MN) SCANS AND AVERAGE THEM - STEP1
cd '/Users/monosovi/Desktop/P/BRICKS/'
path = '/Users/monosovi/Desktop/P/BRICKS/';
unix('to3d -prefix daybscan1_peekads //Users/monosovi/Desktop/P/DICOM/MBF_PEEK_BASELINE/1/M*')
unix('to3d -prefix daybscan2_peekads //Users/monosovi/Desktop/P/DICOM/MBF_PEEK_BASELINE/2/M*')
unix('to3d -prefix daybscan3_peekads //Users/monosovi/Desktop/P/DICOM/MBF_PEEK_BASELINE/3/M*')
cd '/Users/monosovi/Desktop/P/MN _SCRIPT/'
[err, Day2_1, Infop, ErrMessage]= doBrikLoad(path,'daybscan1_peekads+orig');
[err, Day2_2, Infop, ErrMessage]= doBrikLoad(path,'daybscan2_peekads+orig');
[err, Day2_3, Infop, ErrMessage]= doBrikLoad(path,'daybscan3_peekads+orig');
AverageScan= [Day2_1+Day2_2+Day2_3]./3;
doBrikDump(AverageScan,path,'daybAVG_peekmbf');

%take baseline and uniformize  - STEP2
cd '/Users/monosovi/Desktop/P/BRICKS/'
unix('3dUniformize -anat daybAVG_peekmbf+orig -prefix daybAVG_peekmbf_unif -quiet');

%take resample baseline - STEP 3
cd '/Users/monosovi/Desktop/P/BRICKS/'
unix('3dresample -dxy 0.25 0.25 0.25 -rmode ''Cu'' -input daybAVG_peekmbf_unif+orig.HEAD -prefix -overwrite daybAVG_peekmbf_unif+orig');

%%%%%%%%%%PREPARE DAY 1 (24 hours after MN++ injection scan)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%take AVERAGE DAY - STEP1
cd '/Users/monosovi/Desktop/P/BRICKS/'
path = '/Users/monosovi/Desktop/P/BRICKS/';
unix('to3d -prefix day1scan1_peekads //Users/monosovi/Desktop/P/DICOM/MBF_PEEK_DAY_1/1/M*')
unix('to3d -prefix day1scan2_peekads //Users/monosovi/Desktop/P/DICOM/MBF_PEEK_DAY_1/2/M*')
unix('to3d -prefix day1scan3_peekads //Users/monosovi/Desktop/P/DICOM/MBF_PEEK_DAY_1/3/M*')
unix('to3d -prefix day1scan4_peekads //Users/monosovi/Desktop/P/DICOM/MBF_PEEK_DAY_1/4/M*')
unix('to3d -prefix day1scan5_peekads //Users/monosovi/Desktop/P/DICOM/MBF_PEEK_DAY_1/5/M*')
cd '/Users/monosovi/Desktop/P/MN _SCRIPT/'
[err, Day2_1, Infop, ErrMessage]= doBrikLoad(path,'day1scan1_peekads+orig');
[err, Day2_2, Infop, ErrMessage]= doBrikLoad(path,'day1scan2_peekads+orig');
[err, Day2_3, Infop, ErrMessage]= doBrikLoad(path,'day1scan3_peekads+orig');
[err, Day2_4, Infop, ErrMessage]= doBrikLoad(path,'day1scan4_peekads+orig');
[err, Day2_5, Infop, ErrMessage]= doBrikLoad(path,'day1scan5_peekads+orig');
AverageScan= [Day2_1+Day2_2+Day2_3+Day2_4+Day2_5]./5;
doBrikDump(AverageScan,path,'day1AVG_peekmbf');

%STEP 2 - same as for baseline
cd '/Users/monosovi/Desktop/P/BRICKS/'
unix('3dUniformize -anat day1AVG_peekmbf+orig -prefix day1AVG_peekmbf_unif -quiet');
%STEP 3 - same as for baseline
cd '/Users/monosovi/Desktop/P/BRICKS/'
unix('3dresample -dxy 0.25 0.25 0.25 -rmode ''Cu'' -input day1AVG_peekmbf_unif+orig.HEAD -prefix -overwrite day1AVG_peekmbf_unif+orig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%NOW ALIGN THE TWO PROCESSED AVERAGES IN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SEVERAL STEPS

%USE ALINEATE - STEP 1
cd '/Users/monosovi/Desktop/P/BRICKS/';
unix('3dAllineate -overwrite -maxscl 1.000001 -base daybAVG_peekmbf_unif+orig -input day1AVG_peekmbf_unif+orig.HEAD -prefix day1AVG_peekmbf_unif_align+orig -1Dmatrix_save day1matalign');

%USE UNIFIZE - STEP 2
cd '.../MBF';
unix('3dUnifize -prefix day0AVG_peekmbf_unif_align_uni+orig day0AVG_peekmbf_unif_align+orig*')
unix('3dUnifize -prefix day1AVG_peekmbf_unif_align_uni+orig day1AVG_peekmbf_unif_align+orig*')

%%%USE 3dQWARP - STEP 3 _...THIS TAKES FOREVER
unix('3dQwarp -prefix a0_aligned -blur 0 3 -workhard -base daybAVG_peekmbf_unif_uni+orig -source day0AVG_peekmbf_unif_align_uni+orig');
unix('3dQwarp -prefix a1_aligned -blur 0 3 -workhard -base daybAVG_peekmbf_unif_uni+orig -source day1AVG_peekmbf_unif_align_uni+orig');


%%%FINALLY APPLY THE MATRIX ACQUIRED USING THE 3dQWARP.... THIS MATRIX IS USED TO ALIGN/ROTATE THE non-unifized scans
 cd '/Users/monosovi/Desktop/P/BRICKS/';
 unix('3dNwarpApply -nwarp a2_aligned_WARP+orig -source day2AVG_peekmbf_unif_align+orig -prefix day2AVG_peekmbf_unif_aligned_warped');
 unix('3dNwarpApply -nwarp a0_aligned_WARP+orig -source day0AVG_peekmbf_unif_align+orig -prefix day0AVG_peekmbf_unif_aligned_warped');

 
function dumpNewBrik(brikdir,prefix,fmri_tc)

global HDR DATA

outname = sprintf('%s/%s',brikdir,prefix);

if (exist(outname))
    return;
end

M = fmri_tc;

[xd,yd,zd,t] = size(M);
Info.DATASET_DIMENSIONS=[xd yd zd];
Info.DATASET_RANK=[3,size(M,4)];
Info.TYPESTRING = '3DIM_HEAD_ANAT';  % functional
Info.SCENE_DATA = [0 2 0];       % normal
Info.ORIENT_SPECIFIC = [1 3 4];  % Right to Left, Ant to Post, Inf Sup
Info.ORIGIN =  [0 0 0];  %???
Info.BRICK_TYPES = 1*ones(1,t);
Info.TAXIS_NUMS = [Info.DATASET_RANK(2) 0 77002];
Info.TAXIS_FLOATS = [0 HDR.tr 0 0 0];

x = HDR.readres;
y = HDR.phaseres;
z = HDR.sliceres;

Info.DELTA =  [x y z];

Opt.Prefix = prefix;
Opt.NoCheck = 0;
disp('Writing volume as BRIK and HEAD files');
curdir = cd(brikdir);
WriteBrik(M,Info,Opt);  % needs to have a modified version to overwrite
% an existing file with "w+"
cd(curdir);
return
