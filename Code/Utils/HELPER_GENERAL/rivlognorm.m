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
