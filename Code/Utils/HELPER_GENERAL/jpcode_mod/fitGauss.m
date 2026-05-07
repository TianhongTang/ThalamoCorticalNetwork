function results = fitGauss(xdata, ydata, optin)
% %        results = fitGauss(xdata, ydata, optin)
% %
% %  fit (x,y) to a gaussian model with nfit parameters
% %     [amplitude xmean sigma poly(a0,a1,a2) ]
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
% %  Jul 2003 -  Josef Pfeuffer
% %
FCTNAME = 'fitGauss';

fun = 'funGauss';
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.PSTART = zeros(1,7);
dopt.NFIT = 3;
dopt.PLOT = 1;
dopt.MAXITER = 10;
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

Pstart    = dopt.PSTART;
nfit      = dopt.NFIT;
f_plot    = dopt.PLOT;
MAXITER   = dopt.MAXITER;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%
      
if (length(Pstart) < nfit)
    Pstart = [Pstart 0 0 0 0 0 0 0]; 
end

     %% compute_estimates
if Pstart(1) == 0
   x = xdata;
   y = ydata;
   c = polyfit(x,y,1);		 %Fit a straight line
   yf = polyval(c,x);
   yd = y - yf;
   ymax = max(yd);               % x,y and subscript of extrema
   imax = find(yd == ymax);
   xmax	= x(imax);
   ymin = min(yd);          
   imin = find(yd == ymin);
   xmin	= x(imin);

   if (abs(ymax) > abs(ymin))    % emiss or absorp?
      i0=imax;
   else 
      i0=imin; 
   end
   i0 = min([i0 (length(y)-4)]);	 %never take edges
   i0 = max([i0 3]);	
   dy = yd(i0);			%diff between extreme and mean
   del = dy/exp(1);		%1/e value
   il = 0;
   while ( ((i0+il) < length(y)) &... 	% guess at 1/2 width.
      ((i0-il) > 0) &...
      (abs(yd(i0+il)) >= abs(del)) &...
      (abs(yd(i0-il)) >= abs(del)) )
      il = il + 1;
   end
   a = [yd(i0) x(i0) abs(x(i0)-x(i0+il))];
   if (nfit > 3) 	%estimates
      a = [a c(1)]; 
   end
   if (nfit > 4) 
      a = [a c(2)]; 
   end
   if (nfit > 5) 
      a = [a 0.]; 
   end
   Pstart = a;
end

if f_verbose >= 2
    f_display = 'iter';
else
    f_display = 'off';    
end

results = fitLsq(fun,Pstart,xdata,ydata,opt('DISPLAY',f_display,'MAXITER',MAXITER,'PLOT',f_plot));

if ( f_verbose > 0)
   fprintf('\nstartfit: x0 = [');
   fprintf('%10.3g',Pstart);
   fprintf(']\n');
   fprintf('  endfit: x0 = [');
   fprintf('%10.3g',results.Pfit);
   fprintf(']\n');
end
