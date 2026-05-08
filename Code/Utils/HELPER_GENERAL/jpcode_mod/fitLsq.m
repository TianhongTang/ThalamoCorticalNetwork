function RESULTS = fitLsq(fun, Pstart, xdata, ydata, optin);
% %  RESULTS = fitLsq(fun, Pstart, xdata, ydata, optin
% %
% %  least square curvefit (lsqcurvefit) interface
% %     ( from optimization toolbox )
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
FCTNAME = 'fitLsq';

%fun = 'funExpRecov';   % function used for fitting
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.LB          = [];
dopt.UB          = [];
dopt.PLOT    = 1;
dopt.VERBOSE = 1;
    % lsqoptions
dopt.MAXITER     = 10;
dopt.MAXFUNEVALS = 10;
dopt.TOLFUN		 = 1e-4;
dopt.DISPLAY     = 'off';

      % --- arg handling
nargVars = 4;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

lb          = dopt.LB;      % LOWER END UPPER BOUNDS
ub          = dopt.UB;
f_plot      = dopt.PLOT;
f_verbose   = dopt.VERBOSE;
    % options
MAXITER		= dopt.MAXITER;
MAXFUNEVALS = dopt.MAXFUNEVALS*length(Pstart)*2;
TOLFUN		= dopt.TOLFUN;
DISPLAY     = dopt.DISPLAY;
      %%%%  end: handling options  %%%%%
	
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
%        DerivativeCheck: [ on | {off} ]
%             Diagnostics: [ on | {off} ]
%           DiffMaxChange: [ positive scalar {1e-1} ]
%           DiffMinChange: [ positive scalar {1e-8} ]
%                 Display: [ off | iter | notify | final ]
%       GoalsExactAchieve: [ positive scalar | {0} ]
%              GradConstr: [ on | {off} ]
%                 GradObj: [ on | {off} ]
%                 Hessian: [ on | {off} ]
%                HessMult: [ function | {[]} ]
%             HessPattern: [ sparse matrix | {sparse(ones(NumberOfVariables))} ]
%              HessUpdate: [ dfp | gillmurray | steepdesc | {bfgs} ]
%                Jacobian: [ on | {off} ]
%               JacobMult: [ function | ([]) ]
%            JacobPattern: [ sparse matrix | {sparse(ones(Jrows,Jcols))} ]
%              LargeScale: [ {on} | off ]
%      LevenbergMarquardt: [ on | off ]
%          LineSearchType: [ cubicpoly | {quadcubic} ]
%             MaxFunEvals: [ positive scalar ]
%                 MaxIter: [ positive scalar ]
%              MaxPCGIter: [ positive scalar | {max(1,floor(numberOfVariables/2))}]
%           MeritFunction: [ singleobj | multiobj ]
%               MinAbsMax: [ positive scalar | {0} ]
%        PrecondBandWidth: [ positive scalar | {0} | Inf ]
%                  TolCon: [ positive scalar ]
%                  TolFun: [ positive scalar ]
%                  TolPCG: [ positive scalar | {0.1} ]
%                    TolX: [ positive scalar ]
%                TypicalX: [ vector | {ones(NumberOfVariables,1)} ]
% ======================================================================
% exitflag > 0 -- function converged to a solution x
% exitflag = 0 -- max number of evaluation/iterrations was reached
% exitflag < 0 -- function did not converge to a solution

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


if f_plot
    plot(xdata,ydata,'.')
end

options = optimset('MaxIter',MAXITER,'MaxFunEvals',MAXFUNEVALS,'TolFun',TOLFUN,'Display',DISPLAY);
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

%------------------------------------------------------------
