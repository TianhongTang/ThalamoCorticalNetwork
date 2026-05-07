function [x_sem] = nansem(x,dim)
% [x_sem] = nansem(x,dim)
%
%calculate SEM of data along dimension 'dim', excluding nan values

if nargin < 2
    x_sem = nanstd(x,0) ./ sqrt(sum(~isnan(x)));
else
    x_sem = nanstd(x,0,dim) ./ sqrt(sum(~isnan(x),dim));
end;
