function procT1T2histo(optin);
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
FCTNAME = 'procT1T2histo';

global t1 t2 t2s t1map t2map t2smap

savename = 't1t2.mat';

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.MAPINDEX = 0;      
dopt.MAPNAME = 'map2';
dopt.TYPE    = {'t1' 't2' 't2s'};
dopt.WRITE   = 0;           % 0: no write,  1: clean and writeSDT,  2: No fit, writeSDT only
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 0;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

MapIndex  = dopt.MAPINDEX;   % count from zero on (STIMULATE convention)
MapName   = dopt.MAPNAME;
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
    filename = sprintf('%s%s/%d_%s%s', ...
            eval(sprintf('%s.Ses.DataMri', type)),...
            eval(sprintf('%s.Ses.dirname', type)),...
            eval(sprintf('%s.Ses.expp(%s.Ses.grp.%s.exps(1)).scanreco(1)', type, type, type)),...
            type, MapName);
    datFit = STIMrdSdt( filename);
    eval(sprintf('%smap = datFit;', type));
       
    [datFit1] = t1t2Histogram(type, MapIndex);
    
    if f_write > 0
%         savenameFit = sprintf('%smap.mat', type);
%         save(savenameFit, sprintf('%smap', type));
%         fprintf('\n<%s> saved\n', savenameFit)
    end
end
    
return

%------------------------------------------------------------
%------------------------------------------------------------
function  [datFit1] = t1t2Histogram( type, plotMapIndex )
% %  
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
FCTNAME = 't1t2Histogram';

%%%global idat acqp datFit STDPATH
global t1 t2 t2s t1map t2map t2smap

%%%plotMapIndex = 0;  % count from zero on (STIMULATE convention)
 
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.NBINS   = 200;
dopt.VERBOSE = 0;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

nbins     = dopt.NBINS;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

eval(sprintf('datFit = %smap;', type));
s_dat = size(datFit);

maxIntensity = max(max(max(datFit(:,:,:,plotMapIndex+1))))

%for islice = 1:1
for islice = 1:s_dat(3)
    img = datFit(:,:,islice,plotMapIndex+1);    % count from zero on (STIMULATE convention)
    
    figure(3)
    imagesc( img );
    figure(1)
    
    imgdat = reshape(img, prod(size(img)), 1);
    [n,xout] = hist(imgdat,nbins);
    n(1) = 0;   % assume Magnitude images, don't count zeros
    bar(xout,n,'s');
    %xlim([0 maxIntensity]);
    title(sprintf('Histogram slice<%d> MapID0<%d>',islice,plotMapIndex))  
    key
end

datFit1 = 0;
return
%------------------------------------------------------------
