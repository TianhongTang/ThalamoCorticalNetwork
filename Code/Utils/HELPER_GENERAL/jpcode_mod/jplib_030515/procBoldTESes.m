function procBoldTESes(SESSION, optin);
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
FCTNAME = 'procBoldTESes';

global STDPATH acqp reco
global t2s t2smap

savename = 'procBoldTE.mat';
avgNR = 2:10;   % average over time NR

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
    f_t2s = 0;
    f_write = 0;
else
    t2sind = Ses.grp.t2s.exps;
    if t2sind(1) <= 0
        f_t2s = 0;
    else
        f_t2s = 1;
    end
    t2s = [];
    t2smap = [];
end

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
            t2s.dat = zeros(info.nx, info.ny, length(t2sind), info.nslices, 1);
            t2s.time = zeros(length(t2sind),1);
            t2s.Ses = Ses;
            t2s.acqp = acqp;
            t2s.reco = reco;
		end
		dat = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'VERBOSE',f_verbose) );

            %--- averaging images
        dat = avg( dat(:,:,:,avgNR), 4);

        t2s.dat(:,:,ind,:,:) = dat;
        if acqp.EPI_TE_eff > 0
            t2s.time(ind) = acqp.EPI_TE_eff/1000;   % [s]
        else
            t2s.time(ind) = acqp.PVM_EchoTime/1000;   % [s]
        end
        
	end
end  % f_t2s

if f_write
    key('SAVE ? ')
    save(savename,'t2s','t2smap');
    fprintf('\n<%s> saved\n', savename)
    
    %load(savename);
    %%STIMwrSdt( t1map,sprintf('%s%s/%d_t1map', t1.Ses.DataMri,t1.Ses.dirname, t1.Ses.expp(t1.Ses.grp.t1.exps(1)).scanreco(1)), opt('ROTATE',0));

    if ~isempty(t2s)
        s_dat = size(t2s.dat);
        for inr=1:s_dat(3)
            t2sdat(:,:,:,inr) = t2s.dat(:,:,inr,:);
        end
        STIMwrSdt( t2sdat, sprintf('%s%s/%d_t2sdat', t2s.Ses.DataMri,t2s.Ses.dirname,t2s.Ses.expp(t2s.Ses.grp.t2s.exps(1)).scanreco(1)), opt('ROTATE',0));
    end
end


%------------------------------------------------------------
