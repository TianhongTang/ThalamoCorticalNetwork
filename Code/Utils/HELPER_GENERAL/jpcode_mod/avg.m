function avgdata = avg(data, dim);
%
%   averages subdims of a matrix
%
%   JP Sep 2001

narg = nargin;
error(nargchk(1,2,narg));

s_d = size(data);

if (narg == 2)
   avgdata = sum(data, dim) / s_d(dim);
else
   avgdata = sum(reshape( data, prod(s_d), 1), 1) / prod(s_d);
end

