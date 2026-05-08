function procT1T2Ses(SESSION, optin);
% %  procT1T2Ses(SESSION);
% %  get imaging parameters directly from the acqp/reco
% %  create data array/files
% %      sorting images by TE/TR ascending
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
% %  Jan/Jul 2003 -  Josef Pfeuffer
% %
FCTNAME = 'procT1T2Ses';

global STDPATH acqp reco
global t1 t2 t2s t2epi t1map t2map t2smap t2epimap

savename = 't1t2.mat';
f_avgt2epi = 1;

%%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.READMAT = 0;
dopt.WRITE   = 1;
dopt.SLICE   = 0;
dopt.AVGECHO = [0 0];
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
f_slicenum = dopt.SLICE;
f_avgecho = dopt.AVGECHO;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

Ses = goto(SESSION);
%%infogen = getpvpars(SESSION,'t1',1)
Ses.slicenum = f_slicenum;
dataPath.stdpath	= Ses.DataMri;
dataPath.dir		= Ses.dirname;
STDPATH.pv = dataPath.stdpath;

f_t1 = 0;
f_t2 = 0;
f_t2s = 0;
f_t2epi = 0;

if f_readmat
    fprintf('\n<%s> loading\n', savename)
    load(savename);
    f_write = 0;
else
    % check existence of groups
    grpnames = getgrpnames(SESSION);
    t1ind = 0;
    t2ind = 0;
    t2sind = 0;
    t2epiind = 0;
    if checkgrp(grpnames,'t1')
        t1ind = Ses.grp.t1.exps;
        if t1ind(1) > 0
            f_t1 = 1;
        end
    else
        fprintf('--- grp<%s> is not defined in description file <%s>!\n', 't1', SESSION);
    end
    if checkgrp(grpnames,'t2')
        t2ind = Ses.grp.t2.exps;
        if t2ind(1) > 0
            f_t2 = 1;
        end
    else
        fprintf('--- grp<%s> is not defined in description file <%s>!\n', 't2', SESSION);
    end
    if checkgrp(grpnames,'t2s')
        t2sind = Ses.grp.t2s.exps;
        if t2sind(1) > 0
            f_t2s = 1;
        end
    else
        fprintf('--- grp<%s> is not defined in description file <%s>!\n', 't2s', SESSION);
    end
    if checkgrp(grpnames,'t2epi')
        t2epiind = Ses.grp.t2epi.exps;
        if t2epiind(1) > 0
            f_t2epi = 1;
        end
    else
        fprintf('--- grp<%s> is not defined in description file <%s>!\n', 't2epi', SESSION);
    end
    t1 = [];       %% loading new study: clear all data
    t2 = [];
    t2s = [];
    t2epi = [];
    t1map = [];    %% loading new study: clear all maps
    t2map = [];
    t2smap = [];
    t2epimap = [];
end

    %%%%% read T1 images/pars
if (f_t1 == 1)
                % get sorting info 
	for ind=1:length(t1ind)
		ExpNo = t1ind(ind);
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);
        info = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'GETINFO',1,'VERBOSE',f_verbose) );
        if acqp.IMND_recov_time > 0
            time(ind) = acqp.IMND_recov_time/1000;   % [s]
        else
            time(ind) = acqp.PVM_RepetitionTime/1000;   % [s]
        end
	end
    timeSorted = sort(time);
	for ind=1:length(t1ind)
		ExpNo = t1ind(ind);
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);	
		if (ind == 1)
            info = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'GETINFO',1,'VERBOSE',f_verbose) );
            t1.dat = zeros(info.nx, info.ny, length(t1ind), info.nslices, info.nr);
            t1.time = zeros(length(t1ind),1);
            t1.Ses = Ses;
            t1.acqp = acqp;
            t1.reco = reco;
        end
		dat = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'VERBOSE',f_verbose) );
        findIndex = find(timeSorted == time(ind));
        sortIndex = findIndex(1);
        t1.sortindex(ind) = sortIndex;
        t1.dat(:,:,sortIndex,:,:) = dat;
        if acqp.IMND_recov_time > 0
            t1.time(sortIndex) = acqp.IMND_recov_time/1000;   % [s]
        else
            t1.time(sortIndex) = acqp.PVM_RepetitionTime/1000;   % [s]
        end
	end
end  % f_t1

    %%%%% read T2s images/pars
if (f_t2 == 1)
	for ind=1:length(t2ind)
		ExpNo = t2ind(ind);
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);		
		if (ind == 1)
            info = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'GETINFO',1,'VERBOSE',f_verbose) );
            if Ses.slicenum > 0 
                t2.dat = zeros(info.nx, info.ny, length(t2ind), info.nechoes, 1, info.nr);
            else
                t2.dat = zeros(info.nx, info.ny, length(t2ind), info.nechoes, info.nslices, info.nr);
            end
            t2.time = zeros(length(t2ind),1);
            t2.Ses = Ses;
            t2.acqp = acqp;
            t2.reco = reco;
            t2.sortindex = 0;
		end
        if Ses.slicenum > 0 
            dat = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('NS',[f_slicenum f_slicenum],'RECO',reconum,'VERBOSE',f_verbose) );
            dat = reshape(dat, info.nx, info.ny, info.nechoes, 1,info.nr);
        else
            dat = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'VERBOSE',f_verbose) );
            dat = reshape(dat, info.nx, info.ny, info.nechoes, info.nslices,info.nr);
        end
        t2.dat(:,:,ind,:,:) = dat;
        if acqp.IMND_recov_time > 0
            t2.time(ind) = acqp.IMND_echo_time/1000;   % [s]
            t2.time = t2.time(1)*(1:info.nechoes);
        else
            t2.time = acqp.EffectiveTE/1000.0;   % [s]
        end
	end
    if length(t2ind) > 1
        fprintf('t2.dat averaged\n');
        t2.dat = avg(t2.dat,3);
    end
    if Ses.slicenum > 0 
        t2.dat = reshape(t2.dat, info.nx, info.ny, info.nechoes, 1,info.nr);
    else
        t2.dat = reshape(t2.dat, info.nx, info.ny, info.nechoes, info.nslices,info.nr);
    end
end  % f_t2

    % read T2s images/pars
if (f_t2s == 1)
          % get sorting info
	for ind=1:length(t2sind)
		ExpNo = t2sind(ind);
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);		
        info = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'GETINFO',1,'VERBOSE',f_verbose) );
        if acqp.IMND_recov_time > 0
            time(ind) = acqp.IMND_echo_time/1000;   % [s]
        else
            time(ind) = acqp.PVM_EchoTime/1000;   % [s]
        end
	end
    timeSorted = sort(time);
	for ind=1:length(t2sind)
		ExpNo = t2sind(ind);
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);		
		if (ind == 1)
            info = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'GETINFO',1,'VERBOSE',f_verbose) );
            t2s.dat = zeros(info.nx, info.ny, length(t2sind), info.nslices, info.nr);
            t2s.time = zeros(length(t2sind),1);
            t2s.Ses = Ses;
            t2s.acqp = acqp;
            t2s.reco = reco;
		end
		dat = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'VERBOSE',f_verbose) );
        findIndex = find(timeSorted == time(ind));
        sortIndex = findIndex(1);
        t2s.sortindex(ind) = sortIndex;
        t2s.dat(:,:,sortIndex,:,:) = dat;
        if acqp.IMND_recov_time > 0
            t2s.time(sortIndex) = acqp.IMND_echo_time/1000;   % [s]
        else
            t2s.time(sortIndex) = acqp.PVM_EchoTime/1000;   % [s]
        end
	end
end  % f_t2s

    % read T2-EPI images/pars
if (f_t2epi == 1)
              % get sorting info
    for ind=1:length(t2epiind)
		ExpNo = t2epiind(ind);
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);		
        info = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'GETINFO',1,'VERBOSE',f_verbose) );
        if acqp.IMND_echo_time > 0
            time(ind) = acqp.IMND_echo_time/1000;   % [s]  for SE-EPI
        elseif acqp.EPI_TE_eff > 0
            time(ind) = acqp.EPI_TE_eff/1000;   % [s]      for GE-EPI
        else
           time(ind) = acqp.PVM_EchoTime/1000;   % [s]
        end
	end
    timeSorted = sort(time);
	for ind=1:length(t2epiind)
		ExpNo = t2epiind(ind);
		dataPath.filenum	= Ses.expp(ExpNo).scanreco(1);
		reconum				= Ses.expp(ExpNo).scanreco(2);		
		if (ind == 1)
            info = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'GETINFO',1,'VERBOSE',f_verbose) );
            if f_avgt2epi
                t2epi.dat = zeros(info.nx, info.ny, length(t2epiind), info.nslices, 1);
            else
                t2epi.dat = zeros(info.nx, info.ny, length(t2epiind), info.nslices, info.nr);
            end
            t2epi.time = zeros(length(t2epiind),1);
            t2epi.Ses = Ses;
            t2epi.acqp = acqp;
            t2epi.reco = reco;
		end
		dat = PVrd2dseq(dataPath.dir, dataPath.filenum, opt('RECO',reconum,'VERBOSE',f_verbose) );
        if f_avgt2epi
            s_dat = size(dat);
            if length(s_dat) >= 4
                fprintf('   averaging <t2epi.dat> over NR = %d\n',s_dat(4) );
                dat = avg(dat, 4);
            end
        end
        
        findIndex = find(timeSorted == time(ind));
        sortIndex = findIndex(1);
        t2epi.sortindex(ind) = sortIndex;
        t2epi.dat(:,:,sortIndex,:,:) = dat;
        if acqp.IMND_echo_time > 0
            t2epi.time(sortIndex) = acqp.IMND_echo_time/1000;   % [s]  for SE-EPI
        elseif acqp.EPI_TE_eff > 0
            t2epi.time(sortIndex) = acqp.EPI_TE_eff/1000;   % [s]      for GE-EPI
        else
            t2epi.time(sortIndex) = acqp.PVM_EchoTime/1000;   % [s]
        end
	end
end  % f_t2epi

if f_avgecho(1) > 0
    t2datavg = avg( t2.dat(:,:,f_avgecho(1):f_avgecho(2),:), 3);
    s_d = size(t2datavg);
    t2datavg = reshape(t2datavg, [s_d(1:2) s_d(4)]);
    imagesc(t2datavg(:,:,1,1))
    STIMwrSdt( t2datavg, sprintf('%s%s/%d_t2avg%d_%d', t2.Ses.DataMri,t2.Ses.dirname,t2.Ses.expp(t2.Ses.grp.t2.exps(1)).scanreco(1), f_avgecho(1), f_avgecho(2)), opt('ROTATE',0));
end    

if f_write
    key('SAVE ? ')
    save(savename,'t1','t2','t2s','t2epi','t1map','t2map','t2smap','t2epimap');
    fprintf('\n<%s> saved\n', savename)
    
    %load(savename);
    %%STIMwrSdt( t1map,sprintf('%s%s/%d_t1map', t1.Ses.DataMri,t1.Ses.dirname, t1.Ses.expp(t1.Ses.grp.t1.exps(1)).scanreco(1)), opt('ROTATE',0));

    if ~isempty(t1)
        s_dat = size(t1.dat);
        for inr=1:s_dat(3)
            t1dat(:,:,:,inr) = t1.dat(:,:,inr,:);
        end
        STIMwrSdt( t1dat, sprintf('%s%s/%d_t1dat', t1.Ses.DataMri,t1.Ses.dirname,t1.Ses.expp(t1.Ses.grp.t1.exps(1)).scanreco(1)), opt('ROTATE',0));
        t1dat=0;
    end
    
    if ~isempty(t2)
        s_dat = size(t2.dat);
        for inr=1:s_dat(3)
            t2dat(:,:,:,inr) = t2.dat(:,:,inr,:);
        end
        STIMwrSdt( t2dat, sprintf('%s%s/%d_t2dat', t2.Ses.DataMri,t2.Ses.dirname,t2.Ses.expp(t2.Ses.grp.t2.exps(1)).scanreco(1)), opt('ROTATE',0));
        t2dat=0;
    end
    
    if ~isempty(t2s)
        s_dat = size(t2s.dat);
        for inr=1:s_dat(3)
            t2sdat(:,:,:,inr) = t2s.dat(:,:,inr,:);
        end
        STIMwrSdt( t2sdat, sprintf('%s%s/%d_t2sdat', t2s.Ses.DataMri,t2s.Ses.dirname,t2s.Ses.expp(t2s.Ses.grp.t2s.exps(1)).scanreco(1)), opt('ROTATE',0));
        t2sdat=0;
    end
    
    if ~isempty(t2epi)
        s_dat = size(t2epi.dat);
        for inr=1:s_dat(3)
            t2epidat(:,:,:,inr) = t2epi.dat(:,:,inr,:);
        end
        STIMwrSdt( t2epidat, sprintf('%s%s/%d_t2epidat', t2epi.Ses.DataMri,t2epi.Ses.dirname,t2epi.Ses.expp(t2epi.Ses.grp.t2epi.exps(1)).scanreco(1)), opt('ROTATE',0));
        t2epidat=0
    end
end


%------------------------------------------------------------
function grpExists = checkgrp(grpnames, grpStr)
%
grpExists = 0;
for iGrp=1:length(grpnames)
    if strfind(char(grpnames(iGrp)),grpStr)
        grpExists = 1;
        return
    end
end
