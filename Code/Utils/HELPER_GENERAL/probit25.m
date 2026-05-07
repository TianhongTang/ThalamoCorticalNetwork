function y = probit(beta,x)
%
% model for psychophysical probit function
% beta is the following vector: beta = [MU SIGMA]
% DAL 
%
if (nargin < 2)
   error('usage: probit(beta,x])');
end

floor = 0.25;
ceiling = 1.0;


y = floor+(ceiling-floor)*normcdf(x,beta(1),beta(2));
