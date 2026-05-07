function mrsWrite(mrsDat, optin);
% %  mrsWrite(mrsDat, optin)
% %
% %  write FID data for MRS
% %
% %  default options:
% %     opt(  
% %           'VERBOSE','1')
% %
% %  return: 
% %
% %  Functions called: 
% %  Tested for: 
% %
% %  Nov 2002 -  Josef Pfeuffer
% %
FCTNAME = 'mrsWrite';
       
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.WRITEAPP = '_mrs';
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 1;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

WRAPP = dopt.WRITEAPP;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

%  mrsDat.info
%          file: 'M:/mridata/B02.gL1/22/fid'
%         fsize: 2097152
%     precision: 'int32'
%  info.byteorder
%            nx: 2048
%            ny: 1
%       nslices: 1
%            nr: 256


s_tdat = size(mrsDat.tdat);   %check dim not larger than nx/nr and be 2D!
if prod(s_tdat) > mrsDat.info.nx*mrsDat.info.nr; stop; end
%%%tdat = zeros( 2, mrsDat.info.nx/2, mrsDat.info.ny*mrsDat.info.nslices*mrsDat.info.nr );   %% keep full length
tdat = zeros( 2, mrsDat.info.nx/2, s_tdat(2) );     %% keep length only in 1st DIM
tdat(1, 1:s_tdat(1), 1:s_tdat(2)) = real(mrsDat.tdat);
tdat(2, 1:s_tdat(1), 1:s_tdat(2)) = imag(mrsDat.tdat);

file = strcat(mrsDat.info.file, WRAPP);
BYTEORDER = mrsDat.info.byteorder;
PRECISION = mrsDat.info.precision;

if f_verbose
    fprintf('writing <%s>\n', file);
end
fid = fopen(file, 'w', BYTEORDER);
flen = prod(size(tdat));
fcount = fwrite(fid, tdat, PRECISION); 
if (fcount < flen)
    error(sprintf('%s: only %d / %d written', FCTNAME, fcount, flen))
end

return
%------------------------------------------------------------

