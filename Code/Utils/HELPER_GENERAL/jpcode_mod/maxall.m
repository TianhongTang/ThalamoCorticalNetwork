function [maxval,maxind,minval,minind] = maxall(data);
%
%   searches maximum/minimum in whole matrix
%
%   JP Sep 2001

s_d = size(data);
p_d = prod(s_d);
[maxval,maxind] = max( reshape( data, p_d, 1) );
[minval,minind] = min( reshape( data, p_d, 1) );
