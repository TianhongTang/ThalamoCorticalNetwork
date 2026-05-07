function t1t2Ses(SESSION, optin);
% %  t1t2Ses(SESSION);
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
% %  Jan 2003 -  Josef Pfeuffer
% %
FCTNAME = 't1t2Ses';

global STDPATH acqp reco
global t1 t2 t2s t1map t2map t2smap

savename = 't1t2.mat';

%%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.READMAT = 0;
dopt.WRITE   = 1;
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
    t1ind = Ses.grp.t1.exps;
    t2ind = Ses.grp.t2.exps(1);
    t2sind = Ses.grp.t2s.exps;
    if t1ind(1) <= 0
        f_t1 = 0;
    else
        f_t1 = 1;
    end
    if t2ind(1) <= 0
        f_t2 = 0;
    else
        f_t2 = 1;
    end
    if t2sind(1) <= 0
        f_t2s = 0;
    else
        f_t2s = 1;
    end
    t1 = [];       %% loading new study: clear all data
    t2 = [];
    t2s = [];
    t1map = [];    %% loading new study: clear all maps
    t2map = [];
    t2smap = [];
end

    %%%%% read T1 images/pars
if (f_t1 == 1)
	for ind=1:length(t1ind)
		ExpNo = t1ind(ind);
		dataPath.stdpath	= Ses.DataMri;
		dataPath.dir		= Ses.dirname;
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);
		STDPATH.pv = dataPath.stdpath;
		
		if (ind == 1)
            info = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'GETINFO',1,'VERBOSE',f_verbose) );
            t1.dat = zeros(info.nx, info.ny, length(t1ind), info.nslices, info.nr);
            t1.time = zeros(length(t1ind),1);
            t1.Ses = Ses;
            t1.acqp = acqp;
            t1.reco = reco;
		end
		dat = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'VERBOSE',f_verbose) );
        t1.dat(:,:,ind,:,:) = dat;
        if acqp.IMND_recov_time > 0
            t1.time(ind) = acqp.IMND_recov_time/1000;   % [s]
        else
            t1.time(ind) = acqp.PVM_RepetitionTime/1000;   % [s]
        end
	end
end  % f_t1

    %%%%% read T2s images/pars
if (f_t2 == 1)
	%for ind=1:length(t2ind)
	for ind=1:1
		ExpNo = t2ind(ind);
		%fn					= getfilenames(Ses,ExpNo);
		dataPath.stdpath	= Ses.DataMri;
		dataPath.dir		= Ses.dirname;
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);
		STDPATH.pv = dataPath.stdpath;
		
		if (ind == 1)
            info = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'GETINFO',1,'VERBOSE',f_verbose) );
            t2.dat = zeros(info.nx, info.ny, length(t2ind), info.nechoes, info.nslices, info.nr);
            t2.time = zeros(length(t2ind),1);
            t2.Ses = Ses;
            t2.acqp = acqp;
            t2.reco = reco;
		end
		dat = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'VERBOSE',f_verbose) );
        t2.dat = reshape(dat, info.nx, info.ny, info.nechoes, info.nslices,info.nr);
        if acqp.IMND_recov_time > 0
            t2.time(ind) = acqp.IMND_echo_time/1000;   % [s]
            t2.time = t2.time(1)*(1:info.nechoes);
       else
            t2.time = acqp.EffectiveTE/1000.0;   % [s]
        end
	end
end  % f_t2

    % read T2s images/pars
if (f_t2s == 1)
	for ind=1:length(t2sind)
		ExpNo = t2sind(ind);
		dataPath.stdpath	= Ses.DataMri;
		dataPath.dir		= Ses.dirname;
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);
		STDPATH.pv = dataPath.stdpath;
		
		if (ind == 1)
            info = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'GETINFO',1,'VERBOSE',f_verbose) );
            t2s.dat = zeros(info.nx, info.ny, length(t2sind), info.nslices, info.nr);
            t2s.time = zeros(length(t2sind),1);
            t2s.Ses = Ses;
            t2s.acqp = acqp;
            t2s.reco = reco;
		end
		dat = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'VERBOSE',f_verbose) );
        t2s.dat(:,:,ind,:,:) = dat;
        if acqp.IMND_recov_time > 0
            t2s.time(ind) = acqp.IMND_echo_time/1000;   % [s]
        else
            t2s.time(ind) = acqp.PVM_EchoTime/1000;   % [s]
        end
	end
end  % f_t2s

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
    
    if ~isempty(t2)
        s_dat = size(t2.dat);
        for inr=1:s_dat(3)
            t2dat(:,:,:,inr) = t2.dat(:,:,inr,:);
        end
        STIMwrSdt( t2dat, sprintf('%s%s/%d_t2dat', t2.Ses.DataMri,t2.Ses.dirname,t2.Ses.expp(t2.Ses.grp.t2.exps(1)).scanreco(1)), opt('ROTATE',0));
    end
    
    if ~isempty(t2s)
        s_dat = size(t2s.dat);
        for inr=1:s_dat(3)
            t2sdat(:,:,:,inr) = t2s.dat(:,:,inr,:);
        end
        STIMwrSdt( t2sdat, sprintf('%s%s/%d_t2sdat', t2s.Ses.DataMri,t2s.Ses.dirname,t2s.Ses.expp(t2s.Ses.grp.t2s.exps(1)).scanreco(1)), opt('ROTATE',0));
    end
end


%------------------------------------------------------------
