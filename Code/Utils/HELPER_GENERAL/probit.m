function y = probit(beta,x)
% y = probit(beta,x)
% model for psychophysical probit function
% beta is the following vector: beta = [MU SIGMA floor ceiling]
% DAL 
%
if (nargin < 2)
   error('usage: probit([mu sigma floor ceiling,x)');
end

mu      = beta(1);
sigma   = beta(2);
floor   = beta(3);
ceiling = beta(4);

y = floor+(ceiling-floor)*normcdf(x,mu,sigma);

