function RESULTS = fitLsqTest(optin);
% %  RESULTS = fitLsqTest(xdata, ydata, Pstart, optin)
% %
% %  least square curvefit
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
FCTNAME = 'fitLsqTest';

f_test = 1;   % 'funExpRecov';   % Saturation recovery
f_test = 2;   % 'funExp';        % Exponetial decay
f_test = 3;   % 'funExp';        % Exponetial decay
          
      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.PLOT    = 1;
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 0;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_plot    = dopt.PLOT;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%
      
if f_test == 1
    fun = 'funExpRecov';
	Pstart = [1,1,1];
	xdata=[   0.0925, 0.1615,    0.2765,    0.4705,    0.7985,    1.3495,    2.2870, 3.8445];
	ydata=[   0.2478, 0.2817,    0.3347,    0.4154,    0.5302,    0.6746,    0.8258, 0.9383];
	plot(xdata,ydata,'.')
	
	results = fitLsq(fun,Pstart,xdata,ydata,opt('DISPLAY','iter','MAXITER',10))
end
if f_test == 2
    fun = 'funExpMono';
	Pstart = [10,1];
	xdata=[0,1,2,3,4,5,6,7];
	ydata=[10,7,5,4,3.5,3.3,3.5,1.7];
	plot(xdata,ydata,'.')

    results = fitLsq(fun,Pstart,xdata,ydata,opt('DISPLAY','iter','MAXITER',10))
end
if f_test == 3
    fun = 'funExpMono';
	Pstart = [10,1];
	xdata=[0,1,2,3,4,5,6,7];
	ydata=[10,7,5,4,3.5,3.3,3.5,1.7];
    alpha = 0.05;
    
    [beta,resids,J] = nlinfit(xdata,ydata,fun,Pstart);
    ci = nlparci(beta,resids,J);
    [ypred,delta] = nlpredci(fun,xdata,beta,resids,J,alpha,'off','curve')
    beta
    beta-ci(:,1)
    plot(xdata,ydata,'.')
    hold on
    plot(xdata,ypred-delta,'r-')
    plot(xdata,ypred,'-')
	plot(xdata,ypred+delta,'r-')
    
    %%%RESULTS = fitNonLin(fun, Pstart, xdata, ydata);

end

end
%------------------------------------------------------------


% load reaction
% betafit = nlinfit(reactants,rate,@hougen,beta)
% betafit =
%     1.2526
%     0.0628
%     0.0400
%     0.1124
%     1.1914
%     
% load reaction    
% [beta,resids,J] = nlinfit(reactants,rate,'hougen',beta);
% beta =
% 
%     1.2526
%     0.0628
%     0.0400
%     0.1124
%     1.1914
% ci = nlparci(beta,resids,J)
% ci =
%    -1.0798    3.3445
%    -0.0524    0.1689
%    -0.0437    0.1145
%    -0.0891    0.2941
%    -1.1719    3.7321
% 
% load reaction
% [beta,resids,J] = nlinfit(reactants,rate,@hougen,beta);
% [ypred,delta] = nlpredci(@hougen,[100 300 80],beta,resids,J)
% ypred =
%     13
% delta =
%        1.4277