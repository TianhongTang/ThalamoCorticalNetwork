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
