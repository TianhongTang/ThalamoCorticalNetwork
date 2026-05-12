%------------------------------------------------------------
function procT1map(optin)
% %          procT1map(optin)
% %
% %  process T1/R1 maps with different ROIs
% %      uses global t1map
% %      loads roiArr
% %      save ASCII data
% %
% %  default options:
% %     opt(  
% %           'VERBOSE','1')
% %
% %  return: 
% %
% %  Functions called: 
% %  Tested for: 
% %
% %  Jul 2003 -  Josef Pfeuffer
% %
FCTNAME = 'procT1map';

global t1 t2 t2s t2epi t1map t2map t2smap t2epimap
        
savename = 'procT1map.dat';
roiname = 'roiArr.mat';

medfiltArr = [0 3 3];  % [3 3]: median filtering;  Arr(1) <=0: No filtering


      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.MODE    = 2;     % 1: mean/std,  2: gauss fit
dopt.THRES   = 0.10;  % >0: only pixel with T1Error < threshold;  <=0: thresholding off 
dopt.MAXITER = 20;
dopt.STDFAC  = 5;
dopt.NBINS   = 50;
dopt.NFIT    = 3;
dopt.PLOT    = 0;
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 0;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_mode    = dopt.MODE;   
thresT1Error = dopt.THRES; 
MAXITER   = dopt.MAXITER;
STDFAC    = dopt.STDFAC;
NBINS     = dopt.NBINS;
NFIT      = dopt.NFIT;
f_plot    = dopt.PLOT;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

load(roiname);
s_roi = [size(roiArr) 1 1 1];
numRoi = s_roi(3);
numDat = 2;

for idat=1:numDat
    if idat== 1
        t1dat = t1map(:,:,1,2);   % T1 values
        result = zeros(2,2,numRoi);  % mean/std : T1/R1 : #roi
    else
        t1dat = t1map(:,:,1,3);   % R1 values
    end
    t1Error = t1map(:,:,1,6);   % T1 single voxel fit errors
  for iRoi=1:numRoi
    if iRoi== 1
        datTmp = t1dat;
        if medfiltArr(1) > 0
            datTmp = medfilt2(datTmp,medfiltArr);
            fprintf('--- median filter [%d %d]\n',medfiltArr);
        end
    end
    
    mask = roiArr(:,:,iRoi);
    if thresT1Error > 0
        maskind = find(mask > 0 & datTmp > 0 & t1Error < thresT1Error);
    else
        maskind = find(mask > 0 & datTmp > 0);
    end
    dat = datTmp(maskind);

    [meanVal, stdVal, iiter] = stdThres(dat,opt('MODE',f_mode,'STDFAC',STDFAC,'NBINS',NBINS,'NFIT',NFIT,...
        'VERBOSE',f_verbose,'PLOT',f_plot));
    
    result(1,idat,iRoi) = meanVal;
    result(2,idat,iRoi) = stdVal;    
  end
end

key('SAVE as ascii? ')
%%fid = fopen(savename,'w');
%%fprintf(fid,'T1 T1std R1 R1std\n');
%%fclose(fid);
dat2save = reshape(result,2*numDat,numRoi)'
save(savename,'dat2save','-ascii');

%------------------------------------------------------------
%------------------------------------------------------------

