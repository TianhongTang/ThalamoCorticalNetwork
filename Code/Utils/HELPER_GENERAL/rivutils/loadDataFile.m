function file = loadDataFile(dgzfilename,dgzdir,spontflag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loadDataFile 									% 
% This function finds all the relevant files, including:			%
%          DGZ File  (all events and eye movements)				%
%          LFP File  (local field potential file, sampled 500 Hz)		%
% 	       SPK File  (spike file extracted from the adf file)		%
% 	       SPC File  (files containing spectral density information)	%
% 	       SAC File  (files containing extracted EM information)		%
% 	       STM File  (stimulus files)					%
% 	       PDM File  (paradigm files)					%
% 	       MUL File  (preprocessed multiunit data)				%
% While only the DGZ file is loaded automatically, the user is informed		%
% as to the availability of the other files					%
%										%
% DAL  APR-2000									%
%										%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global GH FILES DATA CUR STM PDM

setUserFilePaths;
events


if nargin < 3
   % for spontaneous files	
   spontflag = 0;
end
%
% If there is no argument, it gives a typical 'windows explorer' type
% selection box
%
if ~nargin
  [dgzfilename, dgzdir]  = pickfile('Load DGZ File',FILES.DGZPath,'*.dgz');
elseif nargin == 1
  dgzdir = sprintf('%s',FILES.DGZPath);
end
if ~length(dgzfilename)
  fprintf(' loadDataFile: no dgzfile\n');
  return;
end

if spontflag
  % If it's a spontaneous file analysis, then select the 'collect' dgz file.
  frr = dgzfilename(1:7);
  si  = getSessionInfo(frr);
  dgzfilename = sprintf('%s_collect_%03d.dgz',frr,si.spontfile);
end	

fileRoot = getFileRoot(dgzfilename);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILES is a global structure containing  information about filenames and paths  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FILES.curfiles 		= {};
FILES.LFPFile     	= {};
FILES.ADFFile     	= [];
FILES.SPKFile     	= [];
FILES.DGZFile     	= dgzfilename;		
FILES.MULTIFiles  	= [];
FILES.SPECTFiles  	= [];
FILES.SACFile  	  	= [];
FILES.SPONTFile     	= [];
FILES.FILERoot	   	= fileRoot;
FILES.FILERootRoot 	= fileRoot(1:15);
FILES.FILERootRootRoot 	= fileRoot(1:7)
fileRootRoot 		= FILES.FILERootRoot;

dgzfullpath      = sprintf('%s/%s',dgzdir,dgzfilename);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CUR is a global structure with information relevant to the current analysis    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CUR.lfp_samptime 	= 0.2240;                         % in milliseconds
CUR.passbands    	= {[1 8],[8 15],[15 25],[25 35],[35 45], [45 55], [55 75],...
				   [75 100], [100 150]};
CUR.bandcols    	= {[1 0 0],[1 0.2 0],[1 0.4 0],[0.6 0.6 0],[0.4 0.8 0],...
				   [0.0 1.0 0.2], [0.0 0.8 0.4], [0.0 0.6 0.6],...
				   [0.3 0.5 0.8]};
CUR.sdfsampletime  	= 5;
CUR.coltable       	= [0.5 0.5 1.0 ; 0.5 1.0 0.5 ; 0.8 0.8 0.3 ; 1.0 0.5 0.5 ; ...
		      	0.3 0.8 0.8 ; 0.8 0.3 0.8 ; 0.0 0.5 1.0 ; 1.0 0.0 0.5 ; ...
		      	0.5 1.0 0.0 ; 1.0 0.5 0.0 ; 0.0 1.0 0.5 ; 0.5 0.0 1.0 ;...
		      	0.4 0.8 0.2 ; 0.8 0.4 0.2 ; 0.2 0.8 0.4 ; 0.2 0.4 0.8];
CUR.colstep        	= [0.1 0.1 -0.2 ; 0.1 -0.2 0.1 ; -0.15 -0.15 0.15 ; -0.2 0.1 0.1 ;...
		      	0.15 -0.15 -0.15 ; -0.15 0.15 -0.15 ; 0.2 0.1 -0.2 ; -0.2 0.2 0.1 ; ...
		      	0.1 -0.2 0.2 ; -0.2 0.1 0.2 ; 0.2 -0.2 -0.1 ; -0.1 0.2 -0.2 ;...
		      	0.1 -0.15 0.15 ; -0.15 0.1 0.15 ; 0.15 -0.15 0.1 ; 0.15 0.1 -0.15];
CUR.fadecoltable   	= [0.5 0.8 0.2 ; 0.8 0.5 0.2 ; 0.75 0.75 0.2; 0.6 0.6 0.6];
CUR.respcoltable 	= [0.6 1.0 0.6 ; 1.0 0.6 0.6 ; 1.0 1.0 0.0];
CUR.stimcoltable 	= [0.2 1.0 0.2 ; 1.0 0.2 0.2 ; 0.85 0.85 0.2; 0.6 0.6 0.6];
CUR.maxstims      = 25;
CUR.maxcells      = 5;
CUR.maxchan	  = 16;
CUR.s_pretime       = 1000;
CUR.s_posttime      = 1000;
CUR.s_binwidth      = 50;
CUR.r_pretime       = 1500;
CUR.r_posttime      = 1500;
CUR.e_pretime       = 200;
CUR.e_posttime      = 500;
CUR.e_binwidth      = 10;
CUR.ce_pretime      = 500;
CUR.ce_posttime     = 800;
CUR.r_binwidth      = 50;
CUR.o_binwidth      = 100;
CUR.minemwin        = -1.0;
CUR.maxemwin        = 1.0;
CUR.sp_curobs       = 1;
CUR.spk_visible     = ones(CUR.maxchan,CUR.maxcells);
CUR.nchan 	    = 0;
CUR.minsaccamp      = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the main file								 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = dir(dgzfullpath);
bytes = tmp.bytes;
totfiles = 0;
if bytes < 5000
  FILES.datafiletype = 'CONCATFILE';
  indx = 0;
  fid = fopen(dgzfullpath);
  while 1
   indx = indx + 1;
   l = fgetl(fid);
   if l == -1 
     break; 
   else
    totfiles = totfiles+1;	
    full = sprintf('%s/%s',dgzdir,l);
    dgz(indx) = dg_read(full);
    FILES.curfiles{totfiles} = l(1:19);		
   end
  end
  fclose(fid);
  DATA.dgz = dgzCat(dgz);
else
 FILES.datafiletype = 'SINGLEFILE';
 DATA.dgz    = dg_read(dgzfullpath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the correct names and paths						 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if ~strcmp(FILES.datafiletype,'CONCATFILE')
  eegfilename     = sprintf('%s.%s', fileRoot, FILES.EEGExt);
  adffilename     = sprintf('%s.%s', fileRoot, FILES.ADFExt);
  spkfilename     = sprintf('%s.adfspk', fileRoot);
  stmfilename     = sprintf('%s.stm',    fileRoot);
  pdmfilename     = sprintf('%s.pdm',    fileRoot);
  lfpfilename     = sprintf('%s.%s', 	 fileRoot, FILES.LFPExt);
  multifilename   = sprintf('%s.%s',     fileRoot, FILES.MULTIExt);
  spectfilename   = sprintf('%s.spect',  fileRoot);
  sacfilename     = sprintf('%s.sacc',   fileRoot);

  eegfullpath      = sprintf('%s/%s',	FILES.EEGPath,lfpfilename);
  adffullpath      = sprintf('%s/%s',	FILES.ADFPath,adffilename);
  stmfullpath      = sprintf('%s/%s',	FILES.STMPath,stmfilename);
  pdmfullpath      = sprintf('%s/%s',	FILES.PDMPath,pdmfilename);

  spectfullpath    = sprintf('%s/%s',	FILES.PROC_MULTISPECTPath,spectfilename);
  lfpfullpath      = sprintf('%s/%s',	FILES.PROC_LFPPath,lfpfilename);
  multifullpath    = sprintf('%s/%s',	FILES.PROC_MULTIPath,multifilename);
  sacfullpath      = sprintf('%s/%s',	FILES.PROC_SACPath,sacfilename);
  spkfullpath      = sprintf('%s/%s',	FILES.PROC_SPKPath,spkfilename);    % all chans

  lfpfile     = dir(lfpfullpath);
  multifile   = dir(multifullpath);
  sacfile     = dir(sacfullpath);
  adffile     = dir(adffullpath);
  spkfile     = dir(spkfullpath);
  stmfile     = dir(stmfullpath);
  pdmfile     = dir(pdmfullpath);
  spectfile   = dir(spectfullpath);
%end

if ~isfield(FILES,'notebook'), FILES.notebook = 'wally.txt';  end
nbfullpath  = sprintf('%s/%s', FILES.NBPath, FILES.notebook);
nbfile      = dir(nbfullpath);



%################################################################################
%## Here we start loading all the files that are available  			#
%################################################################################


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD Preprocessed Saccade File (s)						%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(FILES.datafiletype,'CONCATFILE'))
  for i=1:length(FILES.curfiles)		
     sacfilenames{i} = sprintf('%s.sacc', FILES.curfiles{i});
  end	
else
  if length(sacfile)
    sacfilenames   = cellstr(getFileName(sacfile.name))
  else
    sacfilenames = [];
  end
end
nsacfiles = length(sacfilenames);
if nsacfiles > 0	
  sacfullpath = [];	
  for i=1:nsacfiles
    sacfullpath{i}  = sprintf('%s/%s',FILES.PROC_SACPath,sacfilenames{i});
  end	

  DATA.sac  = loadSaccFiles(sacfullpath);
  FILES.SACFile = sacfilenames;
  if exist('GH') & isfield(GH,'sacflag')	
    set(GH.sacflag,'ForegroundColor',[0.1 0.4 0.4]);
  end
else 
  FILES.SACFile = [];
  DATA.sac     = [];
  if exist('GH') & isfield(GH,'sacflag')	
    set(GH.sacflag,'ForegroundColor',[0.65 0.65 0.65]);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD ADFSPK File(s)								%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(FILES.datafiletype,'CONCATFILE'))
  for i=1:length(FILES.curfiles)		
     for j=1:CUR.maxchan
        spkfilenames{i,j} = sprintf('%s.ch%02d.adfspk', FILES.curfiles{i},j);
     end
  end	
else
  spkfilenames   = cellstr(char(spkfile.name));
end

nspkfiles = length(spkfilenames);
if nspkfiles > 0	
  spkfullpath = [];	
  for i=1:nspkfiles
%    spkfullpath{i}  = sprintf('%s/%s',FILES.PROC_SPKPath,spkfilenames{i});
    spkfullpath{i}  = sprintf('%s/%s',FILES.PROC_SPKPath,getFileName(spkfilenames{i}));
  end	
  DATA.spk = [];	
  FILES.SPKFile = spkfullpath;
  if exist('GH') & isfield(GH,'mulflag')	
     set(GH.mulflag,'ForegroundColor',[0.1 0.4 0.4]);
  end
else 
  FILES.SPKFile    = [];
  DATA.adfspk      = [];
  if exist('GH') & isfield(GH,'mulflag')	
     set(GH.spkflag,'ForegroundColor',[0.65 0.65 0.65]);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD Local Field Potential File(s)						%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(FILES.datafiletype,'CONCATFILE'))
  for i=1:length(FILES.curfiles)		
     lfpfilenames{i} = sprintf('%s.lfp', FILES.curfiles{i});
  end	
else
  lfpfilenames   = {lfpfilename};
end
nlfpfiles = length(lfpfilenames);
if nlfpfiles > 0	
  lfpfullpath = [];	
  for i=1:nlfpfiles
    [pt,nm,ext] = fileparts(lfpfilenames{i});
    lfpfullpath{i}  = sprintf('%s/%s%s',FILES.PROC_LFPPath,nm,ext);
  end	
  DATA.lfp = [];	
  FILES.LFPFile = lfpfullpath;
  if exist('GH') & isfield(GH,'lfpflag')	
     set(GH.lfpflag,'ForegroundColor',[0.1 0.4 0.4]);
  end
else 
  FILES.LFPFile = [];
  DATA.lfp      = [];
  if exist('GH') & isfield(GH,'lfpflag')	
     set(GH.lfpflag,'ForegroundColor',[0.65 0.65 0.65]);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD Preprocessed Multiunit File(s)						%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(FILES.datafiletype,'CONCATFILE'))
  for i=1:length(FILES.curfiles)		
     multifilenames{i} = sprintf('%s.multi', FILES.curfiles{i});
  end	
else
  multifilenames   = {multifilename};
end
nmultifiles = length(multifilenames);
if nmultifiles > 0	
  multifullpath = [];	
  for i=1:nmultifiles
    multifullpath{i}  = sprintf('%s/%s',FILES.PROC_MULTIPath,getFileName(multifilenames{i}));
  end	
  DATA.multi = [];	
  FILES.MULTIFile = multifullpath;
  if exist('GH') & isfield(GH,'mulflag')	
     set(GH.mulflag,'ForegroundColor',[0.1 0.4 0.4]);
  end
else 
  FILES.MULTIFile = [];
  DATA.multi      = [];
  if exist('GH') & isfield(GH,'mulflag')	
      set(GH.mulflag,'ForegroundColor',[0.65 0.65 0.65]);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD Preprocessed Spect File(s)						%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(FILES.datafiletype,'CONCATFILE')
  if size(spectfile,1)
    FILES.SPECTFiles = getFileNames(spectfullpath);   % gets the names only
    for i=1:length(FILES.SPECTFiles)
      DATA.spect{i}     = [];
    end
    if exist('GH') & isfield(GH,'spectflag')	
      set(GH.spectflag,'ForegroundColor',[0.1 0.1 0.7]);
    end	
  else 
    FILES.SPECTFiles = [];
    for i=1:16
      DATA.spect{i}     = [];
    end
    if exist('GH') & isfield(GH,'spectflag')	
       set(GH.spectflag,'ForegroundColor',[0.65 0.65 0.65]);
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD Stimulus File(s)								%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(FILES.datafiletype,'CONCATFILE'))
  stmfilename   = sprintf('%s_ALL.stm',FILES.FILERootRoot);
else
  if length(stmfile)
    stmfilename   = stmfile(1).name;
  else
    stmfilename = [];
  end
end
stmfullpath      = sprintf('%s/%s',	FILES.STMPath,stmfilename);
if length(stmfilename)
  FILES.STMFile = stmfilename;
  STM = getSTMFILEInfo2(stmfullpath);
  if exist('GH') & isfield(GH,'stmflag')	
     set(GH.stmflag,'ForegroundColor',[0.4 0.4 0.1]);
  end
else 
  FILES.STMFile = [];
  STM     = [];
  if exist('GH') & isfield(GH,'stmflag')	
     set(GH.stmflag,'ForegroundColor',[0.65 0.65 0.65]);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD Paradigm File(s)								%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(FILES.datafiletype,'CONCATFILE'))
  pdmfilename   = sprintf('%s_ALL.pdm',FILES.FILERootRoot);
else
  if length(pdmfile)
    pdmfilename   = pdmfile(1).name;
  else
    pdmfilename = [];
  end
end
pdmfullpath      = sprintf('%s/%s',	FILES.PDMPath,pdmfilename);
if length(pdmfilename)
  FILES.PDMFile = pdmfilename;
  PDM = getPDMFILEInfo2(pdmfullpath);
  if exist('GH') & isfield(GH,'pdmflag')	
     set(GH.pdmflag,'ForegroundColor',[0.4 0.1 0.1]);
  end	
else 
  FILES.PDMFile = [];
  PDM     = [];
  if exist('GH') & isfield(GH,'pdmflag')	
     set(GH.pdmflag,'ForegroundColor',[0.65 0.65 0.65]);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD Notebook File(s)								%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(nbfile,1)
  [gNB, magic_line] = getNBFileInfo(nbfullpath, FILES.FILERootRoot);
  if exist('GH') & isfield(GH,'nb')	
     set(GH.nb,'String',gNB','ListBoxTop',magic_line-10);
  end
else 
  gNB     = [];
end


										
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not actually load the adf file, but be ready to.				%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('adffile') & size(adffile,1)
  FILES.ADFFile = adffilename;
else 
  FILES.ADFFile = [];
end

if exist('GH') & isfield(GH,'filename')	
   set(GH.filename,'String',FILES.DGZFile);
end
%
% read in the eye movement data into a more convenient form
%
[DATA.emh,DATA.emv,DATA.emt] = getEVTEms(DATA.dgz);



%
% These are important initializations
%
DATA.spk_riv = [];
DATA.tr_sac  = [];
DATA.spikesobs = 0;
DATA.multi_sac = [];
DATA.multi_sac_corr = [];

CUR.multi.single    = [1 1];
CUR.multi.condition = 1;
CUR.multi.epoch     = 1;
CUR.multi.lag       = 512;
CUR.multi.lag_incr  = 64;
CUR.shuffle_correct = 1;
CUR.multi.allags    = [16 32 64 128 192 256 384 512];

CUR.totobs       	= length(DATA.emh);
SI = getSessionInfo(FILES.FILERootRootRoot);
CUR.electrodes 	= SI.electrodes;
CUR.areas	= SI.areas;
CUR.angle	= SI.angle;
CUR.rf		= SI.rf;
if ~isempty(SI.spontfile)
  FILES.SPONTFile = sprintf('%s_collect_%03d',...
	FILES.FILERootRootRoot,SI.spontfile);
  DATA.spontdgz    = dg_read(sprintf('%s/%s.dgz',...
	FILES.DGZPath,FILES.SPONTFile));
end
CUR.nchan = length(CUR.electrodes);
DATA.chan = [1:CUR.nchan];
