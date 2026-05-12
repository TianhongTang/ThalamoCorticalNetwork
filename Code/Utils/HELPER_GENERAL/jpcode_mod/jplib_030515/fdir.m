function sdir = fdir(arg1, arg2, arg3)
%
%%   constructs path/filename according options and STDPATH
%%   sdir = fdir(dir, filenum, optin)
%%   sdir = fdir(dir, optin)
%%
%%   JP Apr 2000

global STDPATH
   
narg = nargin;
error(nargchk(2,3,narg));

    % default options: define whole struct
dopt.DIRAPP = '';
dopt.FILEEXT = '';
dopt.VNMR = 0;
dopt.VNP = 0;
dopt.SDT = 0;
dopt.ROI = 0;
dopt.PATH = 0;
dopt.VERBOSE = 0;

   % handle arguments
dir = arg1;
if (narg == 2)
   optin = arg2;
else
   filenum = arg2;
   optin = arg3;
end
dopt = setopt(dopt,optin);
narg = narg - 1;    % not including options

STDEXT = '.fid';    %of data directory
STDDEF = 'fid';	    %of VNMR filename
SDTEXT = 'P';       %sdt data extension
SDTDEF = '.sdt';    %sdt file extension
ROIEXT = 'P';       %roi data extension
ROIDEF = '.roi';    %roi file extension

dirapp = dopt.DIRAPP;
fileext = dopt.FILEEXT;
f_vnmr = dopt.VNMR;
f_vnp = dopt.VNP;
f_sdt = dopt.SDT;
f_roi = dopt.ROI;
f_path = dopt.PATH;

if (~f_vnmr & ~f_vnp & ~f_sdt & ~f_roi & ~f_path)
   f_vnmr = 1;
end
sdir = '';
f_set = 0;

if (f_set == 0 & f_vnmr) 
   if ~isempty( findstr(dir, STDEXT))
      sdir = dir(1:findstr(dir, STDEXT)-1);
   else
      sdir = dir;
   end

   if (narg == 2) 
      if (~isempty(dirapp))
         sdir = [STDPATH.vnp sdir '_' num2str(filenum) ...
                dirapp STDEXT]; 
      else 
         sdir = [STDPATH.vnmr sdir '/' num2str(filenum) STDEXT];
      end
   else
      sdir = [sdir dirapp STDEXT];
   end
   f_set = 1;
end

if (f_set == 0 & f_sdt) 
   if (narg == 2) 
      sdir = [STDPATH.vnmr dir SDTEXT];
      if ( ~exist(sdir,'dir') ) 
         unix(['mkdir ' sdir])
      end
      sdir = [sdir '/' num2str(filenum) fileext dirapp];
   else
      sdir = [dir fileext dirapp];
   end
   f_set = 1;
end

if (f_set == 0 & f_roi)
   if (narg == 2) 
      sdir = [STDPATH.vnmr dir ROIEXT];
      if ( ~exist(sdir,'dir') ) 
         unix(['mkdir ' sdir])
      end
      sdir = [sdir '/' num2str(filenum) dirapp ROIDEF];
   else
      sdir = [dir dirapp];
   end
   f_set = 1;
end

if (f_set == 0 & f_path) 
   if ~isempty( findstr(dir, STDEXT))
      sdir = dir(1:findstr(dir, STDEXT)-1);
   else
      sdir = dir;
   end

   if (narg == 2) 
      sdir = [STDPATH.vnmr sdir '/'];
   else
      sdir = sdir;
   end
   f_set = 1;
end

%--------------------------------------------------------------
