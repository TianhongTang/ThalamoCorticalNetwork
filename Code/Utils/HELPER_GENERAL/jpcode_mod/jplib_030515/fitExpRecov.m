function RESULTS = fitExpRecov(xdata, ydata, Pstart, optin);
% %  RESULTS = fitExpRecov(xdata, ydata, Pstart, optin)
% %
% %  curvefit of Exponential recovery (Saturation Recovery)
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
FCTNAME = 'fitExpRecov';

fun = 'funExpRecov';   % function used for fitting
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
	Pstart = [1,1,1];
	xdata=[   0.0925, 0.1615,    0.2765,    0.4705,    0.7985,    1.3495,    2.2870, 3.8445];
	ydata=[   0.2478, 0.2817,    0.3347,    0.4154,    0.5302,    0.6746,    0.8258, 0.9383];
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

MAXITER		= 7;
MAXFUNEVALS = 7*length(Pstart)*2;
TOLFUN		= 1.0000e-006;

options = optimset('MaxIter',MAXITER,'MaxFunEvals',MAXFUNEVALS,'TolFun',TOLFUN,'Display','iter');

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
%txt{1}	= 'No stimulus condition';
%txt{2}	= strcat('Monkeys...',sprintf('%s ',ir.monk{1}));
%txt{3}	= strcat('Variables Tested...',sprintf('%s ',ir.vars{:}));
%txt{4}	= strcat('Optimal Parameters...',sprintf('%3.2f, ',X));
% txt{4}	= strcat('Optimal Parameters...',sprintf('%3.2f, ',Pfit));
% txt{5}	= sprintf('ExitFlag = %d ',exitflag);
% txt{6}	='ExitFlag > 0 -- function converged to a solution x';
% txt{7}	='Exitflag = 0 -- max NO of evaluation/iterrations was reached';
% txt{8}	='Exitflag < 0 -- function did not converge to a solution';
% txt{9}	= sprintf('ResNorm = %f',resnorm);
% txt{10} = strcat('Algorithm: ',output.algorithm);
% txt{11}	= sprintf('FuncCount = %d ',output.funcCount);
% txt{12}	= sprintf('Iterations = %d ',output.iterations);
% txt{13}	= sprintf('PCG (Precond. Conj. Gradient) Iterations = %d ',output.iterations);
% txt{14}	= sprintf('Tolerance (TolFun) = %f',usroptions.TolFun);
% txt{15}	= sprintf('Tolerance (TolPCG) = %f',usroptions.TolPCG);
% txt{16} = sprintf('LevenbergMarquardt = %s',usroptions.LevenbergMarquardt);
% txt{17} = sprintf('LineSearchType = %s',usroptions.LineSearchType);


%------------------------------------------------------------
