%--------------------------------------------------------------
function LCMwrRaw(data, header, dir, filenum, optin);
% %  --- write files in LCModel format  ---
% %      called by PV2lcm()
% % 
% %  INPUT:
% % 	data
% %  	header
% % 	dir
% % 	filenum 
% % 
% %  default options:
% %     opt(  'DIRAPP','',
% %           'TWODIM',0,
% %           'VERBOSE',2)
% %
% %  OUTPUT:
% % 	LCModel files are written
% % 
% %  Tested for: PC (not Linux)
% %
% %  June 2002 -  Josef Pfeuffer
% %
FCTNAME = 'LCMwrRaw';

global STDPATH

STDPATHLOCAL = STDPATH.lcm; 
STDEXT  = '.fid';	% of VN data directory
LCMEXT  = '.lcm';	% of LCM data directory
FNAMEDEF = 'lcm';	% of LCModel raw data filename
TWODIMAPP = 'S'; 	% TwoDim appendix

			% definitions for LCM files
LCMEXT1  = '.RAW';
NLBEGIN = '$NMID';
NLEND   = '$END';
FMTDAT  = '(8E14.5)';		% to define FORTRAN format for RAW data: '(8E13.5)'
NCOLS   = 8;

LCMEXT2  = '.PLOTIN';
NLBEGIN2 = '$PLTRAW';
NLEND2   = '$END';

LCMEXT3  = '.CONTROL';
NLBEGIN3 = '$LCMODL';
NLEND3   = '$END';

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.DIRAPP  = '';
dopt.TWODIM  = 0;
dopt.VERBOSE = 2;

      % --- arg handling
nargVars = 4;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

dirapp    = dopt.DIRAPP;
flagTwoDim  = dopt.TWODIM;      
flagVerbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

if (filenum > 0)     % stdpath handling
   file = sprintf('%s%s/%s', STDPATHLOCAL, num2str(dir), ...
      num2str(filenum));
else              % give full path as arg1
   file = sprintf('%s', num2str(dir));
end

% IF(strpos(dir, STDEXT) GT -1) THEN $
%    fdir = strmid(dir, 0, strpos(dir, STDEXT)) $
% ELSE IF(strpos(dir, LCMEXT1) GT -1) THEN $
%    fdir = strmid(dir, 0, strpos(dir, LCMEXT1)) $
% ELSE IF(strpos(dir, LCMEXT2) GT -1) THEN $
%    fdir = strmid(dir, 0, strpos(dir, LCMEXT2)) $
% ELSE IF(strpos(dir, LCMEXT3) GT -1) THEN $
%    fdir = strmid(dir, 0, strpos(dir, LCMEXT3)) $
% ELSE IF(strpos(dir, LCMEXT) GT -1) THEN $
%    fdir = strmid(dir, 0, strpos(dir, LCMEXT)) $
% ELSE $
fdir = dir;

		%% --- whole path: +SLICEAPP+islice+LCMEXT+'/'
if narg == 4 
%     j = checkfile(STDPATH.lcm + fdir, /write, IS_DIR=is_dir)
%     if is_dir EQ 1 then $
%        dirstr = "/" $
%    else $
%       dirstr = "_"
    dirstr = '/';
    fdir = strcat(STDPATH.lcm, fdir, dirstr, sprintf('%d',filenum), dirapp);
else
    fdir = strcat(fdir, dirapp);
end

% j = findfile(fdir + strsl + LCMEXT + '/' + FNAMEDEF + '*', count=c_fdir)
% if (c_fdir GT 0) then begin
%    ch = ''
%    read, '!! overwrite <'+fdir+strsl+LCMEXT+'> ?  (yes/no)[no] ', ch
%    if strmid(strupcase(strtrim(ch, 1)), 0, 3) NE 'YES' then begin
%       return
%    endif
% endif 
if flagVerbose > 1 
   header
end

s_data = size(data);
numTwoDim = prod(s_data(2:end));
ldata = reshape(data, s_data(1), numTwoDim );	% --- lcm data
if flagTwoDim == 0 
    numTwoDim = 1;
end
	
for itwodim=1:numTwoDim 		% ---- LOOP:  2D data set
    
    if flagTwoDim ~= 0
        strTwoDim = sprintf('%s%d', TWODIMAPP, itwodim);
    else
        strTwoDim = '';
    end
    fdirslice = strcat(fdir, strTwoDim, LCMEXT, '/');

%    j = findfile(fdirslice + FNAMEDEF + '*', count=c_fdir)
%    if (c_fdir LE 0) then begin
%       spawn, 'mkdir ' + fdirslice, result
%       if result NE '' then $
%          print, result
%    endif

    header.title = [fdirslice ' ' header.title];
    header.title = header.title(1:min(end,120));    % ----  ; max length is 120 chars
    isdata = ldata(:, itwodim);
         %% convert to alternating real/imag
    isdatatmp = zeros(2, length(isdata));
    isdatatmp(1,:) = real(isdata)';
    isdatatmp(2,:) = imag(isdata)';
    isdata = isdatatmp(:);

			% ---- writing RAW file
file = strcat(fdirslice, FNAMEDEF, LCMEXT1);
unit = fopen(file,'w');
if unit < 0
    error(sprintf('Can not open <%s>', file));
end
if flagVerbose ~= 0
    fprintf('writing <%s>\n', file);
end

for i=1:length(header.comm) 
    fprintf(unit, ' %s\n', header.comm{i});
end
fprintf(unit, ' %s\n', NLBEGIN);
str_id = deblank(upper(header.id));
fprintf(unit, ' ID=''%s''\n', str_id(1:min(16,end)));
fprintf(unit, ' FMTDAT=''%s''\n', FMTDAT);
fprintf(unit, ' TRAMP=%14.5E\n', header.tramp);
fprintf(unit, ' VOLUME=%14.5E\n', header.volume);
fprintf(unit, ' BRUKER=T\n');
fprintf(unit, ' %s\n', NLEND);

nrows = fix( length(isdata) / NCOLS );
for irow=1:nrows
    fprintf(unit, ' %+12.5E', isdata( ((irow-1)*NCOLS+1):irow*NCOLS));
    fprintf(unit, '\n');
end
fclose(unit);

ndatalost = mod(length(isdata), NCOLS);
if ndatalost ~= 0
   fprint('--- writelcmraw: %d data points lost !', ndatalost);
end

			% ---- writing PLOTIN file
file = strcat(fdirslice, FNAMEDEF, LCMEXT2);
unit = fopen(file,'w');
if unit < 0
    error(sprintf('Can not open <%s>', file));
end
if flagVerbose ~= 0
    fprintf('writing <%s>\n', file);
end

fprintf(unit, ' %s\n', NLBEGIN2);
fprintf(unit, ' HZPPPM=%g\n', header.hzpppm);
fprintf(unit, ' NUNFIL=%d\n', header.nunfil);
fprintf(unit, ' DELTAT=%g\n', header.deltat);
fprintf(unit, ' FILRAW=''%s''\n', header.filraw);
fprintf(unit, ' FILPS=''%s''\n', header.filps);
fprintf(unit, ' PPMST=%g\n', header.ppmst);
fprintf(unit, ' PPMEND=%g\n', header.ppmend);
fprintf(unit, ' DEGPPM=%g\n', header.degppm);
fprintf(unit, ' DEGZER=%g\n', header.degzer);
fprintf(unit, ' %s\n', NLEND2);

fclose(unit);
 
 			% ---- writing CONTROL file
file = strcat(fdirslice, FNAMEDEF, LCMEXT3);
unit = fopen(file,'w');
if unit < 0
    error(sprintf('Can not open <%s>', file));
end
if flagVerbose ~= 0
    fprintf('writing <%s>\n', file);
end
 
fprintf(unit, ' %s\n', NLBEGIN3);
fprintf(unit, ' TITLE=''%s''\n', header.title);
fprintf(unit, ' OWNER=''MR Imaging and Spectroscopy, Max Planck Institute for Biological Cybernetics''\n');
fprintf(unit, ' KEY=   316040592, 551241956\n');
fprintf(unit, ' PGNORM=''A4''\n');
fprintf(unit, ' FILPS=''%s''\n', header.filps2);
fprintf(unit, ' FILCOO=''%s''\n', header.filcoo);
if header.noxxt2 == 0
    fprintf(unit, ' ALSDT2(5)=0.01\n');
    fprintf(unit, ' ALSDT2(4)=0.01\n');
    fprintf(unit, ' ALSDT2(3)=0.01\n');
    fprintf(unit, ' ALSDT2(2)=0.01\n');
    fprintf(unit, ' ALSDT2(1)=0.01\n');
    fprintf(unit, ' ALEXT2(5)=0.5\n');
    fprintf(unit, ' ALEXT2(4)=0.5\n');
    fprintf(unit, ' ALEXT2(3)=0.5\n');
    fprintf(unit, ' ALEXT2(2)=0.5\n');
    fprintf(unit, ' ALEXT2(1)=0.5\n');
    fprintf(unit, ' CHSDT2(5)=''Ala''\n');
    fprintf(unit, ' CHSDT2(4)=''Gln''\n');
    fprintf(unit, ' CHSDT2(3)=''Lac''\n');
    fprintf(unit, ' CHSDT2(2)=''Tau''\n');
    fprintf(unit, ' CHSDT2(1)=''Ins''\n');
    fprintf(unit, ' CHEXT2(5)=''Ala''\n');
    fprintf(unit, ' CHEXT2(4)=''Gln''\n');
    fprintf(unit, ' CHEXT2(3)=''Lac''\n');
    fprintf(unit, ' CHEXT2(2)=''Tau''\n');
    fprintf(unit, ' CHEXT2(1)=''Ins''\n');
    fprintf(unit, ' NSDT2=5\n');
    fprintf(unit, ' NEXT2=5\n');
end
fprintf(unit, ' CHCOMB(3)=''Glu+Gln''\n');
fprintf(unit, ' CHCOMB(2)=''GPC+PCho''\n');
fprintf(unit, ' NCOMBI=8\n');
fprintf(unit, ' NAMREL=''Cr+PCr''\n');
fprintf(unit, ' LPRINT=6\n');
fprintf(unit, ' FILPRI=''lcm.PRI''\n');
fprintf(unit, ' LCOORD=9\n');
if ~isempty(header.chkeep)
    for ikeep=length(header.chkeep):-1:1
        fprintf(unit, ' CHKEEP(%d)=''%s''\n', ikeep, header.chkeep{ikeep});
    end
    fprintf(unit, ' NKEEP=%d\n', length(header.chkeep));
end
if ~isempty(header.chomit)
    for iomit=length(header.chomit):-1:1
        fprintf(unit, ' CHOMIT(%d)=''%s''\n', iomit, header.chomit{iomit});
    end
    fprintf(unit, ' NOMIT=%d\n', length(header.chomit));
end
fprintf(unit, ' NAMEAC(5)=''GABA''\n');
fprintf(unit, ' NAMEAC(4)=''GSH''\n');
fprintf(unit, ' NAMEAC(3)=''Ins''\n');
fprintf(unit, ' NAMEAC(2)=''Glc''\n');
fprintf(unit, ' NAMEAC(1)=''Mac''\n');
fprintf(unit, ' NEACH=%d\n', header.neach);
fprintf(unit, ' DEGPPM=%g\n', header.degppm);
if header.sddegp >= 0
    fprintf(unit, ' SDDEGP=%g\n', header.sddegp);
end
fprintf(unit, ' SHIFMN=%g,%g\n', header.shifmn(1),header.shifmn(2) );
if header.ppmshf ~= 0
    fprintf(unit, ' PPMSHF=%g\n', header.ppmshf);
end
fprintf(unit, ' FWHMBA=%g\n', header.fwhmba);
fprintf(unit, ' RFWHM=%g\n', header.rfwhm);
if header.dkntmn ~= 0
    fprintf(unit, ' DKNTMN=%g\n', header.dkntmn);
else
    fprintf(unit, ' ECCDON=.TRUE.\n');
    if header.vitro ~= 0
        fprintf(unit, ' VITRO=.TRUE.\n');
    end
end
fprintf(unit, ' PPMEND=%g\n', header.ppmend);
fprintf(unit, ' PPMST=%g\n', header.ppmst);
fprintf(unit, ' FILRAW=''%s''\n', header.filraw);
fprintf(unit, ' FILBAS=''%s''\n', header.filbas);
fprintf(unit, ' DELTAT=%g\n', header.deltat);
fprintf(unit, ' NUNFIL=%d\n', header.nunfil);
fprintf(unit, ' HZPPPM=%g\n', header.hzpppm);
fprintf(unit, ' %s\n', NLEND3);

fclose(unit);

end 	% for: --- itwodim

return
%--------------------------------------------------------------
