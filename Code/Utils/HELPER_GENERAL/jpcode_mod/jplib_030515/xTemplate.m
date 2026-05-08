function x = x(dir, filenum, optin);
% %  x = x(dir, filenum, optin)
% %
% %  description
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
FCTNAME = 'x';

f_x = 0;          % 
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.X = [0];
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

x = dopt.X;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

%------------------------------------------------------------
