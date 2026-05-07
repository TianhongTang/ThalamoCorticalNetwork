function SesSystemTest(SESSION, optin);
% %  SesSystemTest(SESSION);
% %  get imaging parameters directly from the acqp/reco
% %  create data array/files
% %
% %  default options:
% %     opt(  
% %           'VERBOSE','1')
% %
% %  return: global t1 t2 t2s
% %
% %  Functions called: 
% %  Tested for: 
% %
% %  Apr 2003 -  Josef Pfeuffer
% %
FCTNAME = 'SesSystemTest';

global STDPATH acqp reco
global noise spike

savename = 'SesSystemTest.mat';

%%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.READMAT = 0;
dopt.WRITE   = 0;
dopt.VERBOSE = 0;

      % --- arg handling
nargVars = 1;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_readmat = dopt.READMAT;
f_write   = dopt.WRITE;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

Ses = goto(SESSION);
%%infogen = getpvpars(SESSION,'t1',1)

if f_readmat
    fprintf('\n<%s> loading\n', savename)
    load(savename);
    f_t1 = 0;
    f_t2 = 0;
    f_t2s = 0;
    f_write = 0;
else
    noiseind = Ses.grp.noise.exps;
    if noiseind(1) <= 0
        f_noise = 0;
    else
        f_noise = 1;
    end
    
    spikeind = Ses.grp.spike.exps;
    if spikeind(1) <= 0
        f_spike = 0;
    else
        f_spike = 1;
    end
%     t1 = [];       %% loading new study: clear all data
%     t2 = [];
%     t2s = [];
%     t1map = [];    %% loading new study: clear all maps
%     t2map = [];
%     t2smap = [];
end

    %%%%% read noise data
if (f_noise == 1)
	for ind=1:length(noiseind)
		ExpNo = noiseind(ind);
		dataPath.stdpath	= Ses.DataMri;
		dataPath.dir		= Ses.dirname;
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);
		STDPATH.pv = dataPath.stdpath;
		
		if (ind == 1)
            info = PVrdFid(dataPath.dir, dataPath.filenum, opt('GETINFO',1,'VERBOSE',f_verbose) );
            noise.dat = zeros(info.nx/2, info.ny, length(noiseind), info.nslices, info.nr);
            noise.stat = zeros(3, length(noiseind));  % in percent
            noise.Ses = Ses;
            noise.acqp = acqp;
            noise.reco = reco;
            normfac = (2^16/2-1);
        end
        tdat = PVrdFid(dataPath.dir, dataPath.filenum, opt('VERBOSE',f_verbose) );
		noise.dat(:,:,ind,:,:) = real(double(tdat(:,1)));
        noise.stat(1,ind) = median(real(double(tdat(:,1)))) / normfac * 100;
        noise.stat(2,ind) = std(real(double(tdat(:,1)))) / normfac * 100;
        noise.stat(3,ind) = (max(real(double(tdat(:,1)))) - min(real(double(tdat(:,1))))) / 2^16 * 100;
        plot(real(double(tdat(:,1))),'k')
        ylim([-2^16/2 2^16/2-1])
        title( sprintf('median = %f %%,  std = %f %%,  p-p = %f %%', ...
             noise.stat(1,ind), noise.stat(2,ind), noise.stat(3,ind)) )
        xlabel( sprintf('<%s/%d> ExpNo<%d>', dataPath.dir, dataPath.filenum, ExpNo ))
        fprintf('%s %s, RG=%d, SP=<%.1f,%.1f>, <%.4f,%.2f,%.2f>%%\n', acqp.PULPROG, acqp.DIGMOD, acqp.RG, acqp.SP(1), acqp.SP(2),...
            noise.stat(1,ind), noise.stat(2,ind), noise.stat(3,ind))
        key
	end
end  % f_t1

    %%%%% read noise data
if (f_spike == 1)
	for ind=1:length(spikeind)
		ExpNo = spikeind(ind);
		dataPath.stdpath	= Ses.DataMri;
		dataPath.dir		= Ses.dirname;
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);
		STDPATH.pv = dataPath.stdpath;
		
		if (ind == 1)
            info = PVrdFid(dataPath.dir, dataPath.filenum, opt('GETINFO',1,'VERBOSE',f_verbose) );
            spike.dat = zeros(info.nx/2*info.ny*info.nslices*info.nr, length(spikeind) );
            spike.stat = zeros(3, length(spikeind));  % in percent
            spike.Ses = Ses;
            spike.acqp = acqp;
            spike.reco = reco;
            normfac = (2^16/2-1);
        end
        tdat = PVrdFid(dataPath.dir, dataPath.filenum, opt('VERBOSE',f_verbose) );
        tdat = reshape(tdat, prod(size(tdat)), 1);
        spike.dat(:,ind) = real(double(tdat(:,1)));
        spike.stat(1,ind) = median(real(double(tdat(:,1)))) / normfac * 100;
        spike.stat(2,ind) = std(real(double(tdat(:,1)))) / normfac * 100;
        spike.stat(3,ind) = (max(real(double(tdat(:,1)))) - min(real(double(tdat(:,1))))) / 2^16 * 100;
        plot(real(double(tdat(:,1))),'k')
        ylim([-2^16/2 2^16/2-1])
        title( sprintf('median = %f %%,  std = %f %%,  p-p = %f %%', ...
             spike.stat(1,ind), spike.stat(2,ind), spike.stat(3,ind)) )
        xlabel( sprintf('<%s/%d> ExpNo<%d>', dataPath.dir, dataPath.filenum, ExpNo ))
        fprintf('%s %s, RG=%d, SP=<%.1f,%.1f>, <%.4f,%.2f,%.2f>%%\n', acqp.PULPROG, acqp.DIGMOD, acqp.RG, acqp.SP(1), acqp.SP(2),...
            spike.stat(1,ind), spike.stat(2,ind), spike.stat(3,ind))
        key
	end
end  % f_t1

if f_write
    key('SAVE ? ')
    save(savename,'t1','t2','t2s','t1map','t2map','t2smap');
    fprintf('\n<%s> saved\n', savename)
    
    %load(savename);
    %%STIMwrSdt( t1map,sprintf('%s%s/%d_t1map', t1.Ses.DataMri,t1.Ses.dirname, t1.Ses.expp(t1.Ses.grp.t1.exps(1)).scanreco(1)), opt('ROTATE',0));

    if ~isempty(t1)
        s_dat = size(t1.dat);
        for inr=1:s_dat(3)
            t1dat(:,:,:,inr) = t1.dat(:,:,inr,:);
        end
        STIMwrSdt( t1dat, sprintf('%s%s/%d_t1dat', t1.Ses.DataMri,t1.Ses.dirname,t1.Ses.expp(t1.Ses.grp.t1.exps(1)).scanreco(1)), opt('ROTATE',0));
    end
end


%------------------------------------------------------------
