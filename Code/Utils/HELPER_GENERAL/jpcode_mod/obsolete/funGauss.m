function f = funGauss(a,x)
% f = funGauss(a,x):
%	F = A(0)*EXP(-Z^2/2) + A(3) + A(4)*X + A(5)*X^2
%	Z = (X-A(1))/A(2)
%	Elements beyond A(2) are optional.
%
%   JP Apr 2000

if length(a) < 3
   error('funGauss: length(a) < 3 !')
end

if ( a(3) ~= 0.0 )
   z = (x-a(2))/a(3); 	%SET Z
   ez = exp(-z.^2/2.);	%GAUSSIAN PART
else 
   ez = 0.0;
end

switch length(a)
   case 3
 	f = a(1)*ez;
   case 4
 	f = a(1)*ez + a(4);
   case 5
 	f = a(1)*ez + a(4) + a(5)*x;
   otherwise
 	f = a(1)*ez + a(4) + a(5)*x + a(6)*x.^2; 
end
