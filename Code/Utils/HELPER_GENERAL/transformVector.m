function out = transformVector(in,type)

if nargin < 2 
   errordlg('USAGE: out = transformVector(in, ''FISHERZ'' | ''FISHERZINV'' | ''ZSCORE'')');
   return
end

if (strcmp(type,'FISHERZ'))
   out = 0.5*log((1+in)./(1-in));
elseif (strcmp(type,'FISHERZINV'))
   out = (exp(2*in)-1)./(exp(2*in)+1);
end
