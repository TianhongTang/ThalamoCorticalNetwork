function datProc = erProc(dir, filenum, optin);
% %  dat = erAvg(dir, filenum, optin)
% %
% %  Event-Related Processing: extracting epochs and processing
% %  read 2dseq file, write SDT file
% %
% %  default options:
% %     opt(  'NS',[0 0],       % slice range
% %           'NR',[0 0],       % range in time series
% %           'XY',[0 0],       % single xy point
% %           'TYPE',[control,(stim,control),reptimes,control,TR]
% %           'RECO',1,         % number of reco dir
% %           'AVGEPOCHS',1     %
% %           'ORIENT','ax',    % for writing SDT
% %           'VERBOSE','1')
% %
% %  return: processed data
% %          global acqp reco
% %
% %
% %  Apr 2002 -  Josef Pfeuffer
% %
% %  changes: 
% %  x 2002: 
% %
FCTNAME = 'erProc';

global STDPATH
global acqp reco

APPAVG = '_av';

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.NS = [0 0];
dopt.NR = [0 0]; 
dopt.XY = [0 0];
dopt.TYPE = [0];
dopt.RECO = 1;
dopt.AVGEPOCHS = 1;
dopt.ORIENT = 'ax';
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

nsrange = dopt.NS;
nrrange = dopt.NR; 
xy = dopt.XY;
reconum = dopt.RECO;      
f_avgEpochs = 1;        % average epochs
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

if (filenum > 0)    % stdpath handling
   file = sprintf('%s%s/%s', STDPATH.pv,num2str(dir),num2str(filenum));
else              % give full path as arg1
   file = sprintf('%s', num2str(filenum));
end
savename = file;

    %%% convert TYPE code into full TYPE array
if ( length(dopt.TYPE) ~= 6 )
   typecode = dopt.TYPE(1);
        %%% define your OWN type codes HERE:
   if (typecode == 0) type = [0 0 0 1 0 1]; end
else
    type = dopt.TYPE;
end
typeCont1 = type(1);
typeStimOn = type(2);
typeStimOff = type(3);
typeStimRep = type(4);
typeCont2 = type(5);
typeTR = type(6);

     %% retrieve info struct, acqp, reco pars
imgInfo = PVrd2dseq(dir, filenum, opt('GETINFO',1,'RECO',reconum) );
nx = imgInfo.nx;
ny = imgInfo.ny;
ns = imgInfo.nslices;
nr = imgInfo.nr;

    %% check consistency of TYPE with data
if ( (typeCont1+(typeStimOn+typeStimOff)*typeStimRep+typeCont2) ~= nr )
    fprintf('%s: TYPE is [%d (%d, %d) x %d, %d, %f]\n', FCTNAME, type);
    error(sprintf('%s: TYPE not consistent with data length <%d>', FCTNAME, nr));
end

if ( f_avgEpochs )
    dataAvgEpochs = zeros(nx,ny,ns,(typeStimOn+typeStimOff));    
    readEpoch1 = 1;
    readEpoch2 = typeStimRep;
end

for iEpoch=readEpoch1:readEpoch2
        %%%% read each epoch separately and process it
    nrStart = typeCont1 +(iEpoch-1)*(typeStimOn+typeStimOff) + 1;
    nrEnd   = typeCont1 +(iEpoch)*(typeStimOn+typeStimOff);
    fprintf('---- reading epoch NR [%d %d]\n', nrStart, nrEnd);
    dataEpoch = PVrd2dseq(dir, filenum, opt('NR',[nrStart nrEnd],'RECO',reconum,'VERBOSE',0) );
    dataEpoch = double(dataEpoch);
    
    if ( f_avgEpochs )
        dataAvgEpochs = dataAvgEpochs + dataEpoch;
        dataEpoch = 0;   %% release memory
    end
end

if ( f_avgEpochs )
    datProc = dataAvgEpochs/(readEpoch1+readEpoch2-1);   % normalize
    savename = strcat(savename, APPAVG);
end

if ( ~isempty(savename) )
    STIMwrSdt( datProc, savename, opt('ORIENT',dopt.ORIENT,'ROTATE',0) );
end

return
%--------------------------------------------------------------------
