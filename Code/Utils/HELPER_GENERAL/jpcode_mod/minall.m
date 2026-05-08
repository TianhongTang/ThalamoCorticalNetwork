function [minval,minind] = minall(data);
%
%   searches minimum in whole matrix
%
%   JP Sep 2001

s_d = size(data);
p_d = prod(s_d);
[minval,minind] = min( reshape( data, p_d, 1) );
