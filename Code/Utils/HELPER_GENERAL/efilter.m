function [y,f] = efilter(x,f,varargin)
%ethanfilter(x,f,filter_type)
%
% filters each row of x with a filter f
% for data points at the ends of x, uses only a portion of the
% filter, normalized to preserve its mean
% (thus you can filter a vector without 
%  'losing' data at the edges (as with convn(...,'same'))
%  or shrinking the data vector (as with convn(...,'valid'))
% 
%
%args:
% x ~ a matrix. Each row vector is filtered separately.
%      OR, a column vector. It is transposed, filtered, then transposed
%      again to be returned as a column vector. (for backward compatibility)
% f ~ a row vector. Must have an odd length.
%
% can use a built-in filter:
% specify the filter's length using 'f', and its type using 'filter_type'
% and additional arguments
%
% so far, implemented:
%  'average' - a simple running average (default: 1)
%    Args: window length.
%  'gauss' - gaussian filter. 
%    Args: standard deviation (default: 1), 
%          filter length (default: 1 + 2*(5*stdev))
%
%
%examples:
% figure; hold on;
% x = randn(2,50); plot(x','k');
% plot(ethanfilter(x,[.25 -.75 1 -.75 .25])','Color',[.7 .7 .7]);
% plot(ethanfilter(x,'gauss',1,25)','r');
% plot(ethanfilter(x,'average',5)','b');
%

[y,f] = ethanfilter(x,f,varargin{:});
