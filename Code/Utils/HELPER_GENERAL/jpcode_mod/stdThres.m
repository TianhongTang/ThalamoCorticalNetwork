function [meanVal, stdVal, iiter] = stdThres(dat, optin);
% %  [mean, std, iiter] = stdThres(dat, optin)
% %
% %  calculate mean/std iteratively by thresholding outliers 
% %
% %  default options:
% %     opt(  
% %           'VERBOSE',1)
% %
% %  return: 
% %
% %  Functions called:  mean,std
% %  Tested for: 
% %
% %  Jul 2003 -  Josef Pfeuffer
% %
FCTNAME = 'stdThres';

PROB_FIT_END = 0.001;
STDVALMIN = 0.005;

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.MODE    = 2;     % 1: mean/std,  2: gauss fit
dopt.MAXITER = 20;
dopt.STDFAC  = 5;
dopt.NBINS   = 50;
dopt.NFIT    = 4;
dopt.PLOT    = 0;
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 1;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end

f_mode    = dopt.MODE;   
MAXITER   = dopt.MAXITER;
STDFAC    = dopt.STDFAC;
NBINS     = dopt.NBINS;
NFIT      = dopt.NFIT;
f_plot    = dopt.PLOT;
f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%
      
if f_plot
    f_verbose = f_verbose + 1;
end

dat = dat(:,1,1,1);   % only 1D data
s_dat = size(dat) ;

for iiter=1:MAXITER
    if iiter == 1
        datTmp = dat;
        thres = [0 0];
    else
        thresPrevious = thres;
        stdVal = max([stdVal STDVALMIN]);
        if meanVal <= 0
            meanVal = median(dat);    
        end
        thres = [meanVal-STDFAC*stdVal meanVal+STDFAC*stdVal]; 
        indThres = dat > thres(1) & dat < thres(2) & dat > 0;
        datTmp = dat(indThres);
        if isempty(datTmp)
            meanVal = median(dat);    
            stdVal = std(dat);
            thres = [meanVal-STDFAC*stdVal meanVal+STDFAC*stdVal]; 
            indThres = dat > thres(1) & dat < thres(2) & dat > 0;
            datTmp = dat(indThres);
        end
    end
    
    if f_mode == 2
        [ydata, xdata] = hist(datTmp, NBINS);
        if iiter == 1
            nfit = 3;
        else
            nfit = NFIT;    % 0th order baseline
        end
        if f_plot
            figure(3)
        end

        if f_verbose >=2
            results = fitGauss(xdata,ydata,opt('NFIT',nfit,'PLOT',f_plot,'VERBOSE',1));
        else
            results = fitGauss(xdata,ydata,opt('NFIT',nfit,'PLOT',f_plot,'VERBOSE',0));
        end

        if f_plot
            pause(0.1)
        end
        meanVal = results.Pfit(2);
        stdVal = results.Pfit(3);
    else
        meanVal = mean(datTmp);
        stdVal = std(datTmp);
    end
    if f_verbose 
        figure(1)
        plot(dat)
        if f_mode == 1
            [nhist, xout] = hist(datTmp, NBINS);
            figure(3)
            %bar(xout,nhist);
            plot(xout,nhist);
            title(sprintf('--- iter#%d: mean/std = %.4f +/- %.4f', iiter, meanVal, stdVal));
        end
    end
    if iiter > 1
        crit1 = abs( (thres(1) - thresPrevious(1))/thres(1) );
        crit2 = abs( (thres(2) - thresPrevious(2))/thres(2) );
%%%     if thres(1) == thresPrevious(1) & thres(2) == thresPrevious(2)
        if crit1 < PROB_FIT_END & crit2 < PROB_FIT_END     
            fprintf('--- iter#%d: mean/std = %.4f +/- %.4f\n', iiter, meanVal, stdVal);
            if f_verbose & f_mode == 2
                figure(3)
                plot(xdata,ydata,'.')
                hold on 
                plot(xdata,results.yfit)
                hold off
                title(sprintf('--- iter#%d: mean/std = %.4f +/- %.4f', iiter, meanVal, stdVal));
                pause(0.1)
            end
            return  
        end
    end

end
fprintf('--- iter MAX #%d: mean/std = %.4f +/- %.4f\n', iiter, meanVal, stdVal);

%------------------------------------------------------------
%------------------------------------------------------------
