function [x0fit, yfit] = fitGauss(arg1,arg2,arg3,arg4)
%
%%  fit to a gaussian model with nfit parameters
%%     [amplitude xmean sigma poly(a0,a1,a2) ]
%%  [x0fit, yfit] = fitgauss(x0start,xdata,ydata,nfit)
%%  [x0fit, yfit] = fitgauss(xdata,ydata,nfit)
%%
%%   JP Apr 2000

f_verbose = 1;
fun = 'funGauss';

narg = nargin;
error(nargchk(3,4,narg))
if (narg == 4)
   x0start = arg1;
   xdata = arg2;
   ydata = arg3; 
   nfit = arg4; 
   f_compute_estimates = 0;
   if (length(x0start) < nfit)
      x0start(nfit) = 0 
   end
else
   xdata = arg1;
   ydata = arg2;
   nfit = arg3; 
   f_compute_estimates = 1;
end

if (f_compute_estimates)
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
   x0start = a;
end

     % options used by curvefit: 1,2,3,5,7,9,14,16,17
% options(1) = 0;    % output: 0,1,-1
% options(2) = 1e-6; % termination for f
% options(3) = 1e-6; % termination for g
% options(5) = 0;    % 0:LM method  	1: GN method
% options(7) = 1;    % 0:cubic+quad 	1: only cubic
% options(14) = 1500;
lb = [];
up = [];

x0fit = lsqcurvefit(fun,x0start,xdata,ydata,lb,up);
yfit = feval(fun,x0fit,xdata);

if ( f_verbose )
   fprintf('startfit: x0 = [');
   fprintf('%10.3g',x0start);
   fprintf(']\n');
   fprintf('  endfit: x0 = [');
   fprintf('%10.3g',x0fit);
   fprintf(']\n');
end
