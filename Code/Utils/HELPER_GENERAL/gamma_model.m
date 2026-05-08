function gamma_model(bins,vals)

  
  r0 = 5;
  lambda0 = 5;
  
  gamma_pars = nlinfit(PHASE_BINS,PHASE_HIST_VALS,'rivgamma',[r0 lambda0]);

  figure(100)
  bar(bins,vals,'b');
  hold on
  plot(rivgamma(bins,gamma_pars),'r');






%--------------------------------------------
function y = rivgamma(beta,x)
% y = rivgamma(beta,x)
% model for psychophysical probit function
% beta is the following vector: beta = [r lambda]
% DAL 
%
if (nargin < 2)
   error('usage: rivgamma([r lambda],x)');
end

r      = beta(1);
lambda = beta(2);


term1 = (lambda^r)/gamma(real(r));
term2 = (x.^(r-1));
term3 = exp(-lambda.*x);
y = term1.*term2.*term3;



%--------------------------------------------
function y = rivlognorm(beta,x)
% y = rivlognorm(beta,x)
% model for psychophysical probit function
% beta is the following vector: beta = [u std]
% DAL 
%
if (nargin < 2)
   error('usage: rivlognorm([u std],x)');
end

u      = beta(1);
std    = beta(2);
y = lognpdf(x,u,std);
