function PV2lcmTest(dir, filenum, optin);


if nargin < 1
    %PV2lcm('jptest.en1',4, opt('BASIS','4T', 'TRAMP',1,'VOLUME',1,...
    %    'DEGZER',0,'DEGPPM',0,'PPMST',4.2,'PPMEND',0.1,'VERBOSE',2) );
    PV2lcm('juchem_tes.dv1', 14, opt('BASIS','4T_MAC', 'TRAMP',1,'VOLUME',1,...
        'LSHIFT',63,'DEGZER',0,'DEGPPM',0,'PPMST',4.2,'PPMEND',0.1,'VERBOSE',2) );
end
    %juchem_tes.dv1/14
    
%--------------------------------------------------------------
function PV2lcm(dir, filenum, optin);
% %  PV2lcm(dir, filenum, optin)
% %
% %  Bruker ParaVision
% %  converts FID file to LCModel format
% %
% %  default options:
% %     opt(   'DIRAPP', '',
% %            'TWODIM', 0,
% %            'LSHIFT', 0,
% %            'CT', 0,
% %            'WDWEM', 0,
% %            'BASIS', '',
% %            'SHOWRAW', 0,
% %            'MACHINE', 0,
% %            'XSIZE', 2048,     % number of complex data points
% %            'NODISPLAY', 0,
% %            'VERBOSE','1')
% %                     %%--- LCModel input parameter: (see manual for use)
% %                     %%  default    example
% %            'TRAMP',  1.0,  
% %            'VOLUME', 1.0,
% %            'DEGZER', 0,
% %            'DEGPPM', 0,
% %            'PPMST',  999.,
% %            'PPMEND', -999,
% %            'FWHMBA', 0.009,
% %            'RFWHM',  2.5,
% %            'SHIFMN', [-0.2,-0.1],
% %            'PPMSHF', 0,
% %            'NEACH',  999,
% %            'CHOMIT', '',
% %            'CHKEEP', '',
% %            'VITRO',  0,
% %            'DKNTMN', 0,
% %            'SDDEGP', 0,
% %            'NOXXT2', 0  )
% %
% %  Tested for: PC (not Linux)
% %
% %  June 2002 -  Josef Pfeuffer
% %
FCTNAME = 'PV2lcm';

global STDPATH
global acqp
global reco

FIDFILE = 'fid';        
SDTEXT = '.sdt';
SPREXT = '.spr';
STDPATHLOCAL = STDPATH.lcm; 
PATHBASISSET = '/usr/local/mpi/LCModel/bin/basis-sets/';
VERSION = ' PV2lcm() - JP June 2002';

PLOTRAW = 'plotraw';  %commands
LCMODEL = 'lcm'; 

STDEXT = '.fid';      %VNMR  data directory
LCMEXT  = '.lcm';     %of LCM data directory
LCMEXT1  = '.RAW';
LCMEXT2  = '.PLOTIN';
LCMEXT3  = '.CONTROL';
FNAMEDEF = 'lcm';	% of LCModel raw data filename
SLICEAPP = 'S';		% slice appendix
MACHINEDEF = ['wks4' 'wks4'];	% strarray of possible machines
LCMBATFILE = 'lcmodel.bat';

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.DIRAPP = '';
dopt.TWODIM = 0;
dopt.LSHIFT = 0;
dopt.CT     = 0;
dopt.WDWEM  = 0;
dopt.BASIS  = 'unknown';
dopt.SHOWRAW = 0;
dopt.MACHINE = MACHINEDEF;
dopt.XSIZE  = 2048;
dopt.NODISPLAY = 0;
dopt.VERBOSE = 1;
dopt.TRAMP  = 1.0;    %%--- LCModel input parameter: (see manual for use)
dopt.VOLUME = 1.0;
dopt.DEGZER = 0;
dopt.DEGPPM = 0;
dopt.PPMST  = 999.;
dopt.PPMEND = -999;
dopt.FWHMBA = 0.009;
dopt.RFWHM  = 2.5;
dopt.SHIFMN = [-0.2,-0.1];
dopt.PPMSHF = 0;
dopt.NEACH  = 999;
dopt.CHOMIT = '';
dopt.CHKEEP = '';
dopt.VITRO  = 0;
dopt.DKNTMN = 0;
dopt.SDDEGP = 0;
dopt.NOXXT2 = 0;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

dirapp    = dopt.DIRAPP;
f_slices  = dopt.TWODIM;
flagLS    = dopt.LSHIFT;
f_ct      = dopt.CT;
wdwem     = dopt.WDWEM;
basisfile = dopt.BASIS;
flagShowRaw = dopt.SHOWRAW;
machines  = dopt.MACHINE;
xsize     = dopt.XSIZE;				% complex data points
f_nodisplay = dopt.NODISPLAY;
f_verbose = dopt.VERBOSE;
tramp     = dopt.TRAMP;        %%--- LCModel input parameter: (see manual for use)
volume    = dopt.VOLUME;
degzer    = dopt.DEGZER;
degppm    = dopt.DEGPPM;
ppmst     = dopt.PPMST;
ppmend    = dopt.PPMEND;
fwhmba    = dopt.FWHMBA;
rfwhm     = dopt.RFWHM;
shifmn    = dopt.SHIFMN;
if length(shifmn) ~= 2
    fprintf,'! n_elements(SHIFMN) ~= 2 -> default setting used !\n';
    shifmn = [-0.2,-0.1]
end
ppmshf    = dopt.PPMSHF;
neach     = dopt.NEACH;
chomit    = dopt.CHOMIT;
chkeep    = dopt.CHKEEP;
f_vitro   = dopt.VITRO;
dkntmn    = dopt.DKNTMN;
sddegp    = dopt.SDDEGP;
f_noxxt2  = dopt.NOXXT2;
      %%%%  end: handling options  %%%%%
      
if (filenum > 0)     % stdpath handling
   file = sprintf('%s%s/%s', STDPATHLOCAL, num2str(dir), ...
      num2str(filenum));
else              % give full path as arg1
   file = sprintf('%s', num2str(dir));
end


%IF(strpos(dir, STDEXT) GT -1) THEN $
%   fdir = strmid(dir, 0, strpos(dir, STDEXT)) $
%ELSE $
fdir = dir;

	% ---- determine vndir
%IF (n_params() EQ 2) THEN BEGIN
%      vndir = STDPATH.vnp + fdir + "_" + mystring(filenum) $
%		+ dirapp + STDEXT $

	% ---- determine lcmdir
	%      whole path: +SLICEAPP+islice+LCMEXT+'/'
% IF (n_params() EQ 2) THEN begin
%    j = checkfile(STDPATH.lcm + fdir, /write, IS_DIR=is_dir)
%    if is_dir NE 1 then begin		; create directory
%       spawn,'mkdir '+STDPATH.lcm + fdir, result
%       if result(0) NE '' then print,'MKDIR :'+result
%    endif
%    j = checkfile(STDPATH.lcm + fdir, /write, IS_DIR=is_dir)
%    if is_dir EQ 1 then $
%       dirstr = "/" $
%    else $
%       dirstr = "_"
%    lcmdir = STDPATH.lcm + fdir + dirstr + mystring(filenum) + dirapp 
lcmdir = strcat(fdir, dirapp);

      % read Parameter file(s)
acqp = PVrdPar(dir, filenum);
%%reco = PVrdParReco(dir, filenum, opt('RECO',dopt.RECO));

%data = readvnmr( dir, DIRAPP=dirapp, ct=ctcount )
data = PVrdFid(dir, filenum, opt('FIDFILE','ser') );

%   wdwem = double(f_wdwem)/procpar.sw $
wdwem = 0;

%bc, data, LSHIFT = f_ls, WDWEM=wdwem, /BCFT	; conversion to float

% 		; ---- calibrate according to ctcount
% if f_ct NE 0 then begin
%    data = data / float(f_ct) 
%    print, string(f_ct, format='("--- CT = ",G0.0,"  Calibration")')
% endif else begin
%    j = where(ctcount LT 1, c_ct)
%    if c_ct EQ 0 then begin
%       for ict=0,n_elements(ctcount)-1 do $
%           data(*,ict) = data(*,ict) / float(ctcount(ict))
%       print,  "--- CT = " + strjoin(string(ctcount,format='(G0.0)'))  $
% 		+ "  Calibration"
%    endif else begin
%       print, 'CTCOUNT is zero !'
%       print, ctcount
%       print, '!!! use CT=ctcount to calibrate !'
%       stop
%    endelse
% endelse

if flagLS >= 1
    data = double(data(flagLS:end,:));
else
    data = double(data);
end

 		% ---- truncate/fill according XSIZE
s_td = size(data);
datalen = s_td(1);
datalenNew = xsize;
dataNew = zeros(datalenNew, prod(s_td(2:end)));
if datalen <= datalenNew
    dataNew(1:datalen,:) = data;
    data = dataNew;
else 
    data = data(1:datalenNew,:);
    if datalenNew ~= datalen 
       fprintf('-> data truncated/filled up to XSIZE = %d', xsize);
    end
end

h = defineLCMheader;
				% ---- RAW parameters
h.comm = {strcat('written by ', VERSION) ' ' ' '};	    % cell array of comments
h.id = sprintf('%s/%s', num2str(dir),num2str(filenum));
h.tramp = 1;  %tramp;
h.volume = 1; %volume;
				% ---- PLOTIN parameters
h.hzpppm = acqp.SFO1;
h.nunfil = datalenNew;
h.deltat = 1/acqp.SW_h;   %% maybe factor 2 wrong !!!
h.filraw = strcat(FNAMEDEF, '.RAW');
h.filps  = strcat(FNAMEDEF, '_RAW.PS');
h.ppmst  = ppmst; 
h.ppmend = ppmend;
h.degzer = degzer;
h.degppm = degppm;
				% ---- CONTROL parameters
str_title = strcat(acqp.PULPROG, ' '); 
h.title   = str_title;		% max length is 120 chars
h.filbas  = strcat(PATHBASISSET, basisfile, '.BASIS');
h.filps2  = strcat(FNAMEDEF, '.PS');
h.filcoo  = strcat(FNAMEDEF, '.COORD');
h.vitro   = f_vitro;
h.dkntmn  = dkntmn;
h.sddegp  = sddegp;
h.fwhmba  = fwhmba;
h.rfwhm   = rfwhm;
h.shifmn  = shifmn;
h.ppmshf  = ppmshf;
h.neach   = neach;
h.chomit  = chomit;
h.chkeep  = chkeep;
h.noxxt2  = f_noxxt2;

    %---- write LCModel files
LCMwrRaw( data, h, dir, filenum, opt('DIRAPP','','TWODIM',0,'VERBOSE',2) );

return
%--------------------------------------------------------------
function LCMheader = defineLCMheader

   % define LCM struct here:
   %     = 0   type: double
   %     = ''  type: string
h.comm   = '';	    % arbitrary length
h.id     = '';      % ID string
h.tramp  = 1.0;
h.volume = 1.0;

h.hzpppm = 0;       %---- PLOTIN parameters
h.nunfil = 0;
h.deltat = 0;
h.filraw = 0;
h.filps  = 0;
h.ppmst  = 0;
h.ppmend = 0;
h.degzer = 0;
h.degppm = 0;

h.title  = 0;       %---- CONTROL parameters
h.filbas = '';
h.filps2 = '';
h.filcoo = '';
h.vitro  = 0;
h.dkntmn = 0;
h.sddegp = 0;
h.fwhmba = 0;
h.rfwhm  = 2.5;
h.shifmn =[-0.2 -0.1];
h.ppmshf = 0;
h.neach  = 1;
h.chomit = '';
h.chkeep = '';
h.noxxt2 = 0;

LCMheader = h;

return
%--------------------------------------------------------------
