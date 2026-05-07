function t1t2Fit(optin);
% %  x = x(dir, filenum, optin)
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
% %  Jan 2003 -  Josef Pfeuffer
% %
FCTNAME = 't1t2Fit';

global t1 t2 t2s t1map t2map t2smap

savename = 't1t2.mat';

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.TYPE    = {'t1' 't2' 't2s'};
dopt.WRITE   = 1;           % 0: no write,  1: fit and writeSDT,  2: No fit, writeSDT only
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 0;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

fittype   = dopt.TYPE;
f_write   = dopt.WRITE;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

for itype = 1:length(fittype)
    
    type = char( fittype(itype) );
    if f_write < 2
        fprintf('--- fitting <%s> map ...\n', type);
        [datFit, xdat] = t1t2FitSingleType( type );
        eval(sprintf('%smap = datFit;', type));
    end
    if f_write > 0
        savenameFit = sprintf('%smap.mat', type);
        save(savenameFit, sprintf('%smap', type));
        fprintf('\n<%s> saved\n', savenameFit)

        filename = sprintf('%s%s/%d_%smap', ...
            eval(sprintf('%s.Ses.DataMri', type)),...
            eval(sprintf('%s.Ses.dirname', type)),...
            eval(sprintf('%s.Ses.expp(%s.Ses.grp.%s.exps(1)).scanreco(1)', type, type, type)),...
            type);
        STIMwrSdt( eval(sprintf('%smap', type)), filename, opt('ROTATE',0,'OVERWRITE',1));
    end
end
      
if f_write > 0
    key('SAVE ? ')
    save(savename,'t1','t2','t2s','t1map','t2map','t2smap');
    fprintf('\n<%s> saved\n', savename)
end


%------------------------------------------------------------
%------------------------------------------------------------
function  [datFit, xdat] = t1t2FitSingleType(type, optin)
% %  [datFit, xdat] = t1t2FitSingleType(type, optin)
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
% %  Jan 2003 -  Josef Pfeuffer
% %
FCTNAME = 't1t2FitSingleType';

%%%global idat acqp datFit STDPATH
global t1 t2 t2s t1map t2map t2smap

FITTHRES = 10;    %% SD of noise
noiseInd = [1 1 10 10];
f_plot = 0 ;
f_avg = 0;
f_write = 0;
f_crop = 0;
cropI = [125 127 125 127];  % 
cropI = [73 93 73 93];  % 
 
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.VERBOSE = 0;

      % --- arg handling
nargVars = 1;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

type = char(type);
if strcmp(type, 't1')
    refImage = 0;    %% 0: last image in series (T1 recovery)
elseif strcmp(type, 't2')
    refImage = 1;    %% 1: first image in series (T2 decay)
elseif strcmp(type, 't2s')
    refImage = 1;    %% 1: first image in series (T2 decay)
else
    error(sprintf('%s: unknown TYPE <%s>\n', FCTNAME, type));
end
eval(sprintf('typestruct = %s;', type));

    % format image data
idat = typestruct.dat;  % [nx,ny,nrecovtr,nslices,nr]
s_dat = size(idat);
idat = reshape( idat, s_dat(1),s_dat(2),s_dat(3),prod(s_dat)/s_dat(1)/s_dat(2)/s_dat(3) );
s_dat = [size(idat) 1 1 1 1];
s_dat = s_dat(1:4);

    % prepare fit data
xdat = typestruct.time;   %%
datFit = zeros(s_dat(1), s_dat(2), s_dat(4), 8);    %% x,y,slices,[A, tau, 1/tau, M0, resnorm, Aerr, tauErr, M0err]

if refImage <= 0
    refImage = length(xdat);
end
noiseDat = reshape( double(idat(noiseInd(1):noiseInd(3),noiseInd(2):noiseInd(4),1,1)), ...
    (noiseInd(3)-noiseInd(1)+1)*(noiseInd(4)-noiseInd(2)+1),1 );
noisemean = mean(noiseDat);
noisestd  = std(noiseDat);
imageThres = noisemean + FITTHRES*noisestd;

if f_plot
    figure(1)
    imagesc(double(idat(:,:,refImage,1)) > imageThres)
    if f_verbose
        key
    end
    figure(3)
    warning on
else
    warning off
end

if f_crop
    idat = idat(cropI(1):cropI(2), cropI(3):cropI(4),:,:);
    s_dat = [size(idat) 1 1 1 1];
    s_dat = s_dat(1:4);
    imagesc(double(idat(:,:,refImage,1)) > imageThres)
end
if f_avg
    idatind = find(double(idat(:,:,refImage,1)) > imageThres);
    idat = reshape(idat,s_dat(1)*s_dat(2), s_dat(3), s_dat(4));
    idat = idat(idatind,:,:);
    idat = reshape( avg(idat,1), 1, 1, s_dat(3), s_dat(4));
    s_dat = [size(idat) 1 1 1 1];
    s_dat = s_dat(1:4);
    plot(xdat,squeeze(idat(:,:,:,1)))
end
imageMax = max( squeeze(reshape(double(idat(:,:,refImage,:)),s_dat(1)*s_dat(2)*s_dat(4),1) ))
fprintf('noise: %g  +/-  %g\n',noisemean,noisestd);
fprintf('imageMax = %g, thres = %g (%.2f%%)\n',imageMax,imageThres,imageThres/imageMax*100);

fprintf('slice #');
for isl = 1:s_dat(4)
   fprintf(' %d',isl);
for iy = 1:s_dat(2)
for ix = 1:s_dat(1)
    
    idatVoxel = double( squeeze( idat(ix,iy,:,isl) ));
    if idatVoxel(refImage) > imageThres
		if strcmp(type, 't1')
%%%         RESULTS = fitLsq('funExpRecov',P,xdat,idatVoxel,opt('DISPLAY','off','MAXITER',50,'MAXFUNEVALS',50,'PLOT',f_plot));
            P = [idatVoxel(refImage) xdat(ceil(length(xdat)/5)) 1.0];    %% initial estimate
            RESULTS = fitNonLin('funExpRecov',P,xdat,idatVoxel,opt('PLOT',f_plot));
        elseif strcmp(type, 't2')
%%%         RESULTS = fitLsq('funExpMono',P,xdat,idatVoxel,opt('DISPLAY','off','MAXITER',50,'MAXFUNEVALS',50,'PLOT',f_plot));
            P = [idatVoxel(refImage) xdat(ceil(length(xdat)/5)) 0.0];    %% initial estimate
            RESULTS = fitNonLin('funExpMono',P,xdat,idatVoxel,opt('PLOT',f_plot));
		elseif strcmp(type, 't2s')
%%%         RESULTS = fitLsq('funExpMono',P,xdat,idatVoxel,opt('DISPLAY','off','MAXITER',50,'MAXFUNEVALS',50,'PLOT',f_plot));
            P = [idatVoxel(refImage) xdat(ceil(length(xdat)/5)) 0.0];    %% initial estimate
            RESULTS = fitNonLin('funExpMono',P,xdat,idatVoxel,opt('PLOT',f_plot));
		else
            error(sprintf('%s: unknown TYPE <%s>\n', FCTNAME, type));
		end
        Pfit   = [RESULTS.Pfit' 0.0]';   % extend in case of only 2-parameter-fit
        PfitErr = [(RESULTS.Pfit-RESULTS.PfitCI(:,1))' 0.0]';   % extend in case of only 2-parameter-fit
        if RESULTS.exitflag >= 0
            datFit(ix,iy,isl,1) = Pfit(1);
            datFit(ix,iy,isl,2) = Pfit(2);
            datFit(ix,iy,isl,3) = 1/Pfit(2);
            datFit(ix,iy,isl,4) = Pfit(3);
            datFit(ix,iy,isl,5) = RESULTS.resnorm;
            datFit(ix,iy,isl,6) = PfitErr(1)/Pfit(1);
            datFit(ix,iy,isl,7) = PfitErr(2)/Pfit(2);
            datFit(ix,iy,isl,8) = PfitErr(3)/Pfit(3);
            if f_verbose
                fprintf('%d %d %d: %8.0f (%.0f)   %8.4f (%.4f)  %8.1f (%.1f)\n',ix,iy,isl,Pfit(1),PfitErr(1),Pfit(2),PfitErr(2),Pfit(3),PfitErr(3))
                key
            end
        end
    end
end
end
end
warning on
fprintf('\n');


%------------------------------------------------------------
