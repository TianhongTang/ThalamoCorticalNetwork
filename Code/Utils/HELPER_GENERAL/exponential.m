function y = exponential(beta,x)
% y = exponential(beta,x)
% model for exponential function
% beta is the following vector: beta = [TAU floor]
% DAL 
%
if (nargin < 2)
   error('usage: exponential([TAU floor scale offset],x)');
end


tau     = beta(1);
floor   = beta(2);
scale   = beta(3);
offset  = beta(4);

y = floor+scale*exp(-(x-offset)/tau)';
