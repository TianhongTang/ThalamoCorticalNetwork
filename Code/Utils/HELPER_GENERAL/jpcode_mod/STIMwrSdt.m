function STIMwrSdt( data, filename, optin)
% %  STIMwrSdt(data, filename, optin)
% %
% %  STIMULATE:
% %  --- writes data in SDT/SPR format  ---
% %
% %  default options:
% %     opt(  'ROTATE',1,        % slice range
% %           'ORIENT','sag',    % range in time series
% %           'PRECISION',4,     % 3:LWORD, 4:REAL
% %           'VERBOSE','1')
% %
% %  Sep 2001 -  Josef Pfeuffer
% %
FCTNAME = 'STIMwrSdt';
BYTEORDER = 'ieee-be';

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.ROTATE = 0;	% rotate from MATLAB to STIMULATE
dopt.ORIENT = 'sag';
dopt.PRECISION = 4;       % float32
dopt.OVERWRITE = 0;
dopt.VERBOSE = 1;

% --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_rotate = dopt.ROTATE;	
orient = dopt.ORIENT;
precind = dopt.PRECISION; 
f_overwrite = dopt.OVERWRITE;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

   % --- defs
%precind=    1       2        3       4          5       6
precStr = {'';   'int16'; 'int32'; 'float32'; '';      ''};
typeStr = {'BYTE'; 'WORD';'LWORD'; 'REAL'; 'COMPLEX'; 'TBD'};

tdat = data;
dim = size(tdat);
if ( ndims(tdat) < 2 )
   size(data)
   error('--- not possible to handle 1D data')
end
if ( ndims(tdat) < 4 )
   dimtmp = [1 1 1 1];
   dimtmp(1:ndims(tdat)) = dim;
   dim = dimtmp;
end
ni = prod(dim)/(dim(1)*dim(2));

if (f_rotate) 
    if dim(1) == dim(2)
        disp( sprintf('--- rotating (%d)', f_rotate) );
        tdat = reshape(tdat, dim(1), dim(2), ni);
        for ii=1:ni
            tdat(:,:,ii) = tdat(:,:,ii)';    
        end
    else
        disp( sprintf('--- !!! CAN NOT rotate (%d) matrix with dim[%d,%d]', f_rotate, dim(1), dim(2)) );
    end
end

      % ----- write SDT (SignalDaTa) file (raw binary data)
fn = [filename '.sdt'];
if ( exist(fn,'file') & f_overwrite == 0) 
   ch = input(['!! overwrite <' fn '> ?  (yes/no)[no] '], 's');
   if ( strncmpi(strjust(ch,'left'), 'yes', 3) ~= 1)
      return
   end
end

disp(['writing <' fn '> ']);
fp = fopen( fn, 'w', BYTEORDER);
fwrite( fp, tdat, char(precStr(precind)) );
fclose( fp );

      % ----- write SPR (SignalPaRameters) file (ASCII)
fn = [filename '.spr'];
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
disp(['writing <' fn '>']);
type(fn);


%--------------------------------------------------------------------
