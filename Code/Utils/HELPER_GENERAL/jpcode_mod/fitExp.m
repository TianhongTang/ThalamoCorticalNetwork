function RESULTS = fitExp(xdata, ydata, Pstart, optin);
% %  RESULTS = fitExp(xdata, ydata, Pstart, optin)
% %
% %  curvefit of Exponential
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
% %  Jan 2003 -  Josef Pfeuffer
% %
FCTNAME = 'fitExp';

fun = 'funExpMono';  % function used for fitting
f_test = 0;          % 
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.PLOT    = 1;
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 3;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_plot    = dopt.PLOT;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%
      
if f_test
	Pstart = [10,1];
	xdata=[0,1,2,3,4,5,6,7];
	ydata=[10,7,5,4,3.5,3.3,3.5,1.7];
	plot(xdata,ydata,'.')
	
	lb = [];
	up = [];
	[Pfit, resnorm] = lsqcurvefit(fun,Pstart,xdata,ydata,lb,up)
	yfit = feval(fun,Pfit,xdata);
	hold on 
	plot(xdata,yfit)
	hold off
end

if f_plot
    plot(xdata,ydata,'.')
end
if f_plot
    plot(xdata,ydata,'.')
end
	
% ======================================================================
% REMINDER: HERE ARE THE DEFAULT MATLAB OPTIONS FOR lsqcurvefit
% ======================================================================
% defaultopt = optimset('display','final','LargeScale','on', ...
%   'TolX',1e-6,'TolFun',1e-6,'DerivativeCheck','off',...
%   'Diagnostics','off',...
%   'Jacobian','off','MaxFunEvals','100*numberOfVariables',...
%   'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
%   'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
%	'MaxPCGIter','max(1,floor(numberOfVariables/2))', ...
%   'TolPCG',0.1,'MaxIter',400,'JacobPattern',[], ...
%   'LineSearchType','quadcubic','LevenbergMarq','on'); 
% ======================================================================

MAXITER		= 20;
MAXFUNEVALS = 20*length(Pstart)*2;
TOLFUN		= 1.0000e-006;

options = optimset('MaxIter',MAXITER,'MaxFunEvals',MAXFUNEVALS,'TolFun',TOLFUN,'Display','off');

% LOWER END UPPER BOUNDS
lb = []; %[.1 1 1 1 0 -2];
ub = []; %[2 8 8 4 45 2];

% exitflag > 0 -- function converged to a solution x
% exitflag = 0 -- max number of evaluation/iterrations was reached
% exitflag < 0 -- function did not converge to a solution

[Pfit,resnorm,residual,exitflag,output,LAMBDA,JACOB] = ...
		lsqcurvefit(fun,Pstart,xdata,ydata,lb,ub,options);

yfit = feval(fun,Pfit,xdata);

if nargout,
	RESULTS.Pfit = Pfit;
	RESULTS.resnorm = resnorm;
	RESULTS.residual = residual;
	RESULTS.exitflag = exitflag;
	RESULTS.output = output;
	RESULTS.yfit = yfit;
end;

if f_plot
	hold on 
	plot(xdata,yfit)
	hold off
    pause(0.1)
end

