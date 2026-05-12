function mrsDat = mrsRead(dir, filenum, optin);
% %  mrsDat = mrsRead(dir, filenum, optin)
% %
% %  reads FID data for MRS
% %
% %  default options:
% %     opt(  
% %           'VERBOSE','1')
% %
% %  return: struct mrsDat
% %
% %  Functions called: mrsDatInit
% %  Tested for: 
% %
% %  Nov 2002 -  Josef Pfeuffer
% %
FCTNAME = 'mrsRead';

global acqp
      
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.FIDFILE = {'fid.orig' 'fid' 'ser'};    %% order of preference to read FID data;
dopt.CONVDTA = 0;
dopt.VERBOSE = 0;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

fidFileArr = dopt.FIDFILE;
f_convd2a  = dopt.CONVDTA;
f_verbose  = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

   %% get info about target data
   %% try different FID file names
fidFile = '';
iarr = 1;
while isempty(fidFile) & iarr <= length(fidFileArr)
    fidFileTmp = char( fidFileArr(iarr) );    
    fidInfoTmp = PVrdFid(dir, filenum, opt('FIDFILE',fidFileTmp,'GETINFO',1,'VERBOSE',f_verbose));
    if fidInfoTmp.fsize > 0
        fidFile = fidFileTmp;
        fidInfo = fidInfoTmp;
        iarr = length(fidFileArr);
    end
    iarr = iarr + 1;
end
if isempty(fidFile)
    error(sprintf('%s: file NOT found <%s>', FCTNAME, fidInfoTmp.file));
end

mrsDat = mrsDatInit;
mrsDat.dir     = dir;
mrsDat.filenum = filenum;
   %% read data
mrsDat.tdat = PVrdFid(dir, filenum, opt('FIDFILE',fidFile,'VERBOSE',f_verbose));
mrsDat.acqp = acqp;
mrsDat.info = fidInfo;

    %% MRS convention: 2 dimensions in double (complex) format
mrsDat.tdat = double(mrsDat.tdat);   % don't expect memory problems FOR NOW !
s_tdat = size(mrsDat.tdat);
mrsDat.tdat = reshape(mrsDat.tdat, s_tdat(1), prod(s_tdat)/s_tdat(1));

if f_convd2a
    mrsDat = mrsConvdta(mrsDat,opt('SHIFT',72.5));
end
if f_verbose
    mrsDat
end

return

%------------------------------------------------------------
function mrsDat = mrsDatInit;
% %  mrsDat = mrsDatInit
% %
% %  defines/inits mrsDat
% %
% %  return: struct mrsDat
% %
% %  Functions called: 
% %  Tested for: 
% %
% %  Nov 2002 -  Josef Pfeuffer
% %
FCTNAME = 'mrsDatInit';

%--- define struct mrsDat 
mrsDat.dir     = '';
mrsDat.filenum = 0;
mrsDat.tdat = 0;
mrsDat.acqp = 0;
mrsDat.info = 0;
mrsDat.eccDat  = 0;
mrsDat.freqDat = 0;
mrsDat.flagEccDone     = 0;
mrsDat.flagConvdtaDone = 0;
mrsDat.flagAvgDone     = 0;
mrsDat.flagFreqDone    = 0;
%--- end struct mrsDat

return
