function [roiArr, refImg] = mgetMultRoi(img, numRoi, optin);
% %  [roiArr, refImg] = mgetMultRoi(img, numRoi, optin)
% %
% %  define interactively multiple ROIs and save 
% %
% %  default options:
% %     opt(  'SAVENAME','roiArr.mat',
% %           'REFIMG',1,
% %           'VERBOSE',1)
% %
% %  return: 
% %
% %  Functions called:  mgetroi
% %  Tested for: 
% %
% %  Jul 2003 -  Josef Pfeuffer
% %
FCTNAME = 'mgetMultRoi';
         
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.SAVENAME = 'roiArr.mat';
dopt.REFIMG  = 1;
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

savename  = dopt.SAVENAME;
refNum    = dopt.REFIMG;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

s_dat = [size(img) 1 1 1 1];
refNum = min([refNum s_dat(3)]);
refImg = img(:,:,refNum,1);
for iRoi=1:numRoi
    if iRoi == 1
        roiArr = zeros(s_dat(1),s_dat(2),numRoi);
    end
    str = sprintf('----- select ROI #%d',iRoi);
    fprintf('%s\n',str)
    DEF.IROINAME = str;
    roi = mgetroi(refImg, DEF);
    maskSignal = roi.mask';   % rotate mask !!

    roiArr(:,:,iRoi) = maskSignal;
    
    figure(1)
    imagesc(refImg')
    figure(3)
    imagesc(maskSignal')
end
key('SAVE ? ')
save(savename,'roiArr');
fprintf('\n<%s> saved\n', savename)


%------------------------------------------------------------
