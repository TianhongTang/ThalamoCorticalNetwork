%--------------------------------------------------------------
function dorkTest

     % Dork SETUP
dataPath.stdpath = 'M:/mpidata/';
dataPath.dir     = 'N02.gh1';
dataPath.filenum = 8;
optstruct        = opt('NR',[4 0],'VERBOSE',0);

     % call Dork routine
DATA = getDorkData(dataPath, optstruct)

     % use DATA further on ...
figure(10)
islice=1
subplot(2,1,1)
plot(squeeze(DATA.deltaAmp(islice,1,:)))
ylabel('D amplitude [%]')
subplot(2,1,2)
plot(squeeze(DATA.deltaFreq(islice,1,:)))
ylabel('D offres [Hz]')

end
%--------------------------------------------------------------
function DorkData = getDorkData(dataPath, optstruct)
% %
% % Analysis of Dynamic Off-Resonance Effects in K-space (DORK)
% % Interface to function dorkTrace()
% % 
% % input:
% % dataPath.stdpath = 'M:/mpidata/';    !! MUST have '/' at the end !!
% % dataPath.dir     = 'N02.gh1';
% % dataPath.filenum = 8;
% % default options:
% %     opt(  'NS',[0 0],       % slice range
% %           'NR',[0 0],       % range in time series
% %           'REFNUM',1,       % reference FID #
% %           'FIDFILE','fid',  % fid/fid.orig/ser
% %           'DETREND',0,
% %           'FUDGE' = [0 0 0 0];
% %           'VERBOSE','1')
% %
% %  return: dat.deltaAmp  (nslices, nsegments, nrepetitions);
% %          dat.deltaFreq (nslices, nsegments, nrepetitions);
% %          dat.acqp      (ACQP struct)
% %          global acqp/reco parameters
% %
% %  Tested for: onepulse, EPI
% %
% %  Oct 2002 - JP
% %
FCTNAME = 'getDorkData';

global acqp navFidDat navDat
global STDPATH
STDPATH.pv = dataPath.stdpath;
dir        = dataPath.dir;
filenum    = dataPath.filenum;

DorkData = dorkTrace(dir, filenum, optstruct);


end
%--------------------------------------------------------------