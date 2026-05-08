function l2 = repeat(v,n)
%
% FUNCTION
% l2 = repeat(v,n)
%
% ARGS
% v = input vector
% n  = number of times to repeat  i.e.  111222333444
%
% see also: replicate
  
l2 = [];
for i = 1:length(v)
  l2 = [l2 replicate(v(i),n)];
end

