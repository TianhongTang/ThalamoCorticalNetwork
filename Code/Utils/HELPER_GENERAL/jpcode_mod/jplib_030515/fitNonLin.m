function RESULTS = fitNonLin(fun, Pstart, xdata, ydata, optin);
% %  RESULTS = fitNonLin(xdata, ydata, Pstart, optin)
% %
% %  Nonlinear least-squares data fitting by the Gauss-Newton method(nlinfit) interface
% %     ( from statistics toolbox )
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
% %  March 2003 -  Josef Pfeuffer
% %
FCTNAME = 'fitNonLin';

%fun = 'funExpRecov';   % function used for fitting
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.PLOT    = 1;
dopt.VERBOSE = 1;
    % lsqoptions
dopt.alpha   = 0.05;

      % --- arg handling
nargVars = 4;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_plot      = dopt.PLOT;
f_verbose   = dopt.VERBOSE;
    % lsqoptions
alpha       = dopt.alpha;   % confidence interval (1-alpha)
    %%%%  end: handling options  %%%%%
	
if f_plot
    plot(xdata,ydata,'.')
end

[Pfit,residual,JACOB] = nlinfitJP(xdata,ydata,fun,Pstart);
if length(JACOB(:,1)) >= length(Pfit(:,1))
       % confidence intervals (95%) on parameters
   PfitCI = nlparci(Pfit,residual,JACOB);
       % Confidence intervals on predictions of nonlinear models
       %   'simopt'    'on' for simultaneous intervals or 'off' (the default) for non-simultaneous intervals
       %   'predopt'   'curve' (the default) for confidence intervals for the function value at the inputs, or
       %                'observation' for confidence intervals for a new response value
   [yfit,yfitConf] = nlpredci(fun,xdata,Pfit,residual,JACOB,alpha,'off','curve');
else
       %% underdetermined: can not calculate
    PfitCI = zeros(length(Pfit(:,1)),2);  
    PfitCI(:,1) = Pfit;
    PfitCI(:,2) = Pfit;
    yfit = feval(fun,Pfit,xdata);
    yfitConf = zeros(size(yfit));
end
        
if nargout,
	RESULTS.Pfit = Pfit;
	RESULTS.PfitCI = PfitCI;
	RESULTS.residual = residual;
	RESULTS.resnorm =  sum(residual.^2);
	RESULTS.exitflag = 1;
	RESULTS.output = [];
	RESULTS.yfit = yfit;
	RESULTS.yfitConf = yfitConf;
end;

if f_plot
	hold on 
	plot(xdata,yfit,'-')
    plot(xdata,yfit-yfitConf,'r-')
	plot(xdata,yfit+yfitConf,'r-')
	hold off
    pause(0.1)
end

%------------------------------------------------------------
