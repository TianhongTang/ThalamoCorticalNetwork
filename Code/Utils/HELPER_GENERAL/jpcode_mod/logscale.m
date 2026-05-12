function x = logscale(logarr, optin);
% %  x = logscale(logarr, optin)
% %
% %  logscale = [xStart xEnd xNsteps]
% %
% %  create isodistant x values on a log scale
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
% %  Dec 2002 -  Josef Pfeuffer
% %
FCTNAME = 'logscale';
       
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 1;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%
      
xlogstart = log( logarr(1) );
xlogend   = log( logarr(2) );
xnum   = logarr(3);
xlogstep = (xlogend - xlogstart)/(xnum-1);

xlog = xlogstart:xlogstep:xlogend;
x = exp(xlog);

%------------------------------------------------------------
