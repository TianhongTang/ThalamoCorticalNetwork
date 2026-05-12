function procT1T2Clean(optin);
% %  procT1T2Clean(optin)
% %
% %  description
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
% %  May 2003 -  Josef Pfeuffer
% %
FCTNAME = 'procT1T2Clean';

global t1 t2 t2s t1map t2map t2smap

savename = 't1t2.mat';
maskBrainStr = 'maskBrain';
maskCsfStr   = 'maskCSF';
mask1Str   = 'mask1';
mask2Str   = 'mask2';
mask3Str   = 'mask3';

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.THRES1   = [9.99e10 6]';      
dopt.THRES2   = [9.99e10 6]';      
dopt.TYPE    = {'t1' 't2' 't2s'};
dopt.WRITE   = 1;           % 0: no write,  1: clean and writeSDT,  2: No fit, writeSDT only
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 0;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

threshold1 = dopt.THRES1;
threshold2 = dopt.THRES2;
fittype   = dopt.TYPE;
f_write   = dopt.WRITE;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

for itype = 1:length(fittype)
    
    type = char( fittype(itype) );
    savenameFit = sprintf('%smap.mat', type);
    save(savenameFit, sprintf('%smap', type));
    fprintf('\n<%s> saved\n', savenameFit)

        %% read 'map'
    filename = sprintf('%s%s/%d_%smap', ...
            eval(sprintf('%s.Ses.DataMri', type)),...
            eval(sprintf('%s.Ses.dirname', type)),...
            eval(sprintf('%s.Ses.expp(%s.Ses.grp.%s.exps(1)).scanreco(1)', type, type, type)),...
            type);
    datFit = STIMrdSdt( filename);
    eval(sprintf('%smap = datFit;', type));
    
    %% read 'mask'
    filename = sprintf('%s%s/%s', ...
            eval(sprintf('%s.Ses.DataMri', type)),...
            eval(sprintf('%s.Ses.dirname', type)),...
            maskBrainStr);
    maskBrain = STIMrdSdt( filename);
    filename = sprintf('%s%s/%s', ...
            eval(sprintf('%s.Ses.DataMri', type)),...
            eval(sprintf('%s.Ses.dirname', type)),...
            maskCsfStr);
    maskCsf = STIMrdSdt( filename);
    
    if f_write < 2
        fprintf('--- cleaning <%s> map ...\n', type);
        [datFit1, datFit2, mask] = t1t2Clean(threshold1, threshold2, maskBrain, maskCsf, type );
        eval(sprintf('%smap = datFit;', type));
        
        datFit1( find(datFit1 ==  Inf) ) = 0;
        datFit1( find(datFit1 == -Inf) ) = 0;
        datFit1( find(datFit1 >  2^31) ) = 0;
        datFit1( find(datFit1 < 0) ) = 0;

        datFit2( find(datFit2 ==  Inf) ) = 0;
        datFit2( find(datFit2 == -Inf) ) = 0;
        datFit2( find(datFit2 >  2^31) ) = 0;
        datFit2( find(datFit2 < 0) ) = 0;
    end
    if f_write > 0
        savenameFit = sprintf('%smap.mat', type);
        save(savenameFit, sprintf('%smap', type));
        fprintf('\n<%s> saved\n', savenameFit)

        filename = sprintf('%s%s/%d_%smap1', ...
            eval(sprintf('%s.Ses.DataMri', type)),...
            eval(sprintf('%s.Ses.dirname', type)),...
            eval(sprintf('%s.Ses.expp(%s.Ses.grp.%s.exps(1)).scanreco(1)', type, type, type)),...
            type);
        STIMwrSdt( datFit1, filename, opt('ROTATE',0,'OVERWRITE',0));

        filename = sprintf('%s%s/%d_%smap2', ...
            eval(sprintf('%s.Ses.DataMri', type)),...
            eval(sprintf('%s.Ses.dirname', type)),...
            eval(sprintf('%s.Ses.expp(%s.Ses.grp.%s.exps(1)).scanreco(1)', type, type, type)),...
            type);
        STIMwrSdt( datFit2, filename, opt('ROTATE',0,'OVERWRITE',0));
        
        filename = sprintf('%s%s/%s', ...
            eval(sprintf('%s.Ses.DataMri', type)),...
            eval(sprintf('%s.Ses.dirname', type)),...
            mask1Str);
        STIMwrSdt( mask(:,:,:,1), filename, opt('ROTATE',0,'OVERWRITE',0));
    
        filename = sprintf('%s%s/%s', ...
            eval(sprintf('%s.Ses.DataMri', type)),...
            eval(sprintf('%s.Ses.dirname', type)),...
            mask2Str);
        STIMwrSdt( mask(:,:,:,2), filename, opt('ROTATE',0,'OVERWRITE',0));
    
        filename = sprintf('%s%s/%s', ...
            eval(sprintf('%s.Ses.DataMri', type)),...
            eval(sprintf('%s.Ses.dirname', type)),...
            mask3Str);
        STIMwrSdt( mask(:,:,:,3), filename, opt('ROTATE',0,'OVERWRITE',0));
    end
end
    
return
if f_write > 0
    key('SAVE ? ')
    save(savename,'t1','t2','t2s','t1map','t2map','t2smap');
    fprintf('\n<%s> saved\n', savename)
end


%------------------------------------------------------------
%------------------------------------------------------------
function  [datFit, datFit2, mask] = t1t2Clean(threshold1, threshold2, maskBrain, maskCsf, type, optin)
% %  [datFit, datFit2, mask] = t1t2Clean(threshold, thresImageNum, maskBrain, maskCsf, type, optin)
% %
% %  description
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
% %  May 2003 -  Josef Pfeuffer
% %
FCTNAME = 't1t2Clean';

%%%global idat acqp datFit STDPATH
global t1 t2 t2s t1map t2map t2smap

plotSliceNum = 1;
thresImageVal1 = threshold1(:,1);  
thresImageNum1 = threshold1(:,2)+1;    % count from zero on (STIMULATE convention)
thresImageVal2 = threshold2(:,1);   
thresImageNum2 = threshold2(:,2)+1;    % count from zero on (STIMULATE convention)

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.VERBOSE = 0;

      % --- arg handling
nargVars = 5;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

eval(sprintf('datFit = %smap;', type));
eval(sprintf('datFit2 = %smap;', type));
s_dat = size(datFit);

thresMask1 = ones(s_dat(1:3));
for iThres = 1:size(thresImageVal1)
    thresMaskTmp = abs( datFit(:,:,:,thresImageNum1(iThres))) >= thresImageVal1(iThres);   
    thresMask1( find(thresMaskTmp) ) = 0;

    figure(3)
    imagesc( ~thresMask1(:,:,plotSliceNum,1) );
    title(sprintf('thresmask1 [%d,%.2f]',thresImageNum1(iThres),thresImageVal1(iThres)))
    key
end
thresMask2 = ones(s_dat(1:3));
for iThres = 1:size(thresImageVal2)
    thresMaskTmp = abs( datFit(:,:,:,thresImageNum2(iThres)) ) >= thresImageVal2(iThres);   
    thresMask2( find(thresMaskTmp) ) = 0;

    figure(3)
    imagesc( ~thresMask2(:,:,plotSliceNum,1) );
    title(sprintf('thresmask2 [%d,%.2f]',thresImageNum2(iThres),thresImageVal2(iThres)))
    key
end
imagesc( maskBrain(:,:,plotSliceNum) );
title('maskBrain')
key
imagesc( maskCsf(:,:,plotSliceNum) );
title('maskCsf')
key

thresMask2 = thresMask2 & maskBrain & ~maskCsf;

for i=1:s_dat(4)
    datTmp = datFit(:,:,:,i);
    datTmp( find(~thresMask1) ) = 0;
    datFit(:,:,:,i) = datTmp;
    
    datTmp2 = datFit2(:,:,:,i);
    datTmp2( find(~thresMask2) ) = 0;
%    datTmp2( find(~maskBrain) ) = 0;
%    datTmp2( find(maskCsf) ) = 0;
    datFit2(:,:,:,i) = datTmp2;
end
mask(:,:,:,1) = thresMask1;
mask(:,:,:,2) = thresMask2;
mask(:,:,:,3) = datFit2(:,:,:,3) > 0;   % R2 mask
%------------------------------------------------------------
