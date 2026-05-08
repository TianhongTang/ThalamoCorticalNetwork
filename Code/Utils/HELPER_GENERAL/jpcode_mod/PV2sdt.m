function PV2sdt(dir, filenum, optin);
% %  PV2sdt(dir, filenum, optin)
% %
% %  Bruker ParaVision
% %  converts 2dseq RECO file to stimulate SDT format
% %
% %  default options:
% %     opt(   'ORIENT','sag',
% %            'RECO',2,
% %            'VERBOSE','1')
% %
% %  Tested for: PC (not Linux)
% %
% %  Sep 2001 -  Josef Pfeuffer
% %
FCTNAME = 'PV2sdt';

global STDPATH
global acqp
global reco

IMGFILE = '2dseq';        
SDTEXT = '.sdt';
SPREXT = '.spr';

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.ORIENT = 'sag';
dopt.RECO = 1;
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

orient = dopt.ORIENT;
reconum = dopt.RECO;      
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

if (filenum > 0)     % stdpath handling
   file = sprintf('%s%s/%s/pdata/%s/%s', STDPATH.pv,num2str(dir), ...
      num2str(filenum),num2str(reconum),IMGFILE);
   if (reconum <=1)
      outfile = sprintf('%s%s/%s', STDPATH.pv,num2str(dir),num2str(filenum));
      outfileS = sprintf('%s/%s', num2str(dir),num2str(filenum));
   else
      outfile = sprintf('%s%s/%s_%s', STDPATH.pv,num2str(dir),num2str(filenum), ...
                           num2str(reconum) );
      outfileS = sprintf('%s/%s_%s', num2str(dir),num2str(filenum), ...
                           num2str(reconum) );
   end
else              % give full path as arg1
   file = sprintf('%s/%s', num2str(dir), IMGFILE);
   outfile = sprintf('%s', num2str(dir));
   error([FCTNAME ': !! dopt() not yet implemented !!'])
end
      % read Parameter file(s)
acqp = PVrdPar(dir, filenum);
reco = PVrdParReco(dir, filenum, opt('RECO',dopt.RECO));

nx = reco.RECO_size(1);
ny = reco.RECO_size(2);
if (length(reco.RECO_size) > 2)
   reco.RECO_size
   error(sprintf('%s: unexpected case of parameters !', FCTNAME ));
end   
nslices = acqp.NSLICES;
nr = acqp.NR;
if strcmp(reco.RECO_wordtype, '_16BIT_SGN_INT')
   NBYTE = 2;              % 16 bit 
   precind = 2;
elseif strcmp(reco.RECO_wordtype, '_32BIT_SGN_INT')
   NBYTE = 4;              % 32 bit 
   precind = 3;
else
   error(sprintf('%s: unknown type <%s>', FCTNAME, reco.RECO_wordtype))
end

      % consistency check
fsize = checkfilesize(file);
if ( fsize/NBYTE == nx*ny*nslices*nr)
   fprintf('reading <%s>\n', file);
   fprintf('%s %s\n', acqp.PULPROG, acqp. GRDPROG);
   fprintf('dim <%d><%d><%d><%d>  filesize %dk\n', ...
            nx, ny, nslices, nr, fsize/1024);
else
   error(sprintf('%s: nint <%d> ~= dimensions <%d><%d><%d><%d>\n', ...
            FCTNAME, fsize/NBYTE, nx, ny, nslices, nr));
end   

   % --- defs
%precind=    1       2        3       4          5       6
precStr = {'';   'int16'; 'int32'; 'float32'; '';      ''};
typeStr = {'BYTE'; 'WORD';'LWORD'; 'REAL'; 'COMPLEX'; 'TBD'};

      % ----- copy SDT (SignalDaTa) file (raw binary data)
fn = [outfile SDTEXT];
if ( exist(fn,'file') ) 
   ch = input(['!! overwrite <' outfileS SDTEXT '> ?  (yes/no)[no] '], 's');
   if ( strncmpi(strjust(ch,'left'), 'yes', 3) ~= 1)
      return
   else
      if isunix
         unix(sprintf('rm %s ',fn));    % to be able to link file
      else
         unix(sprintf('del %s ',fn));    % to be able to link file
      end
   end
end

if isunix
   unix(sprintf('ln -s %s %s',file,fn));
else
   unix(sprintf('ln %s %s',file,fn));
end
disp(['--- file <' fn '> written']);

      % ----- write SPR (SignalPaRameters) file (ASCII)
dim = [nx ny nslices nr];
fn = [outfile SPREXT];
fp = fopen( fn, 'w' );

fprintf( fp, 'numDim: %d\n', length(dim));
fprintf( fp, 'dim: ');
for idx=1:length(dim)
   fprintf( fp, '%d ', dim(idx));
end
fprintf( fp, '\n');
fprintf( fp, 'dataType: %s\n', char(typeStr(precind)) );
if ( ~isempty(orient) )
   fprintf( fp, 'sdtOrient:%s\n', orient);
end

fclose( fp );
disp(['--- file <' fn '> written']);
type(fn);


