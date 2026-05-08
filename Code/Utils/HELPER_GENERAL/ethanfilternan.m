function [y,f] = ethanfilternan(x,f,varargin)
%ethanfilternan(x,f,filter_type)
%
% filters each row of x with a filter f
% for data points at the ends of x, uses only a portion of the
% filter, normalized to preserve its mean
% (thus you can filter a vector without 
%  'losing' data at the edges (as with convn(...,'same'))
%  or shrinking the data vector (as with convn(...,'valid'))
% does similar normalization for data points which are NaN
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
%  'gauss' - gaussian filter. 
%    Args: standard deviation (default: 1), 
%          filter length (default: 1 + 2*ceil(5*stdev))
%  'average' - a simple running average
%    Args: window length.  (default: 1)
%
%  'causalaverage' - causal running average. Same as 'average' but shifted
%          to only include data from the current column and later columns
%    Args: window length (default: 1)
%  'causalexpo' - causal exponential. Based on an exponential probability
%          distribution, cropped at a specified length
%    Args: mean of exponential distribution (default: 1)
%          filter length (default: until CDF of exponential distribution = 0.99)
%
%
%examples:
% figure; hold on;
% x = randn(2,50); plot(x','k');
% plot(ethanfilter(x,[.25 -.75 1 -.75 .25])','Color',[.7 .7 .7]);
% plot(ethanfilter(x,'gauss',1,25)','r');
% plot(ethanfilter(x,'average',5)','b');
%

if nargin < 2
    error('need at least two arguments');
elseif nargin >= 2
    %use f as the filter?
    if ~isstr(f)
        f = f;
    %use a built-in filter?
    elseif strcmp(f,'average') == 1
        if length(varargin) > 1 error('filter type ''average'' expects 0-1 args'); end;
        if length(varargin) >= 1 
            filter_length = varargin{1}; 
        else
            filter_length = 1; 
        end;
        
        f = ones(1,filter_length) ./ filter_length;
    elseif strcmp(f,'causalaverage') == 1
        if length(varargin) > 1 error('filter type ''causalaverage'' expects 0-1 args'); end;
        if length(varargin) >= 1 
            filter_length = varargin{1}; 
        else
            filter_length = 1; 
        end;
        
        f = [(ones(1,filter_length) ./ filter_length) zeros(1,filter_length-1)];
    elseif strcmp(f,'causalexpo') == 1
        if length(varargin) > 2 error('filter type ''causalaverage'' expects 0-2 args'); end;
        if length(varargin) >= 1
            filter_scale = varargin{1};
        else
            filter_scale = 1;
        end;
        if length(varargin) >= 2
            filter_length = varargin{2};
        else
            filter_length = ceil(icdf('exp',.99,filter_scale));
        end;
        f = pdf('exp',(filter_length-1):-1:0,filter_scale);
        f = normalize([f zeros(1,filter_length-1)],1);
    elseif strcmp(f,'gauss') == 1
        if length(varargin) > 2 error('filter type ''gauss'' expects 0-2 args'); end;
        if length(varargin) >= 1 
            filter_std = varargin{1};
        else 
            filter_std = 1;
        end;
        if length(varargin) >= 2 
            filter_length = varargin{2};
        else
            filter_length = 1 + 2*ceil(5 .* filter_std);
        end;
        
        f = diff(cdf('Normal',(-filter_length/2):1:(filter_length/2),0,filter_std));
        f = f ./ sum(f);
    else
        error(['unknown filter type: ' filter_type]);
    end;
end;

x_was_col_vector = false;
if isvector(x) && ~isrowvector(x)
    x = x';
    x_was_col_vector = true;
end;

if length(size(f)) ~= 2 || ~isrowvector(f) || mod(length(f),2) == 0
    error('filter must be a row vector with an odd length');
end;

if isempty(x)
    % no need to filter an empty input!
    y  = x;
    return;
end;

% NOTE: for some reason, convn applies f as if it was reversed...?
%  e.g. I would expect f = [1 1 1 0 0] to be a causal filter, but instead
%  it smears the data *forward* in time...
% so, here I manually reverse the filter 'f' so it works as I expect.
f = f(end:-1:1);

%run the filter over the data, neglecting NaNs
x_nan2zero = x; 
x_nan2zero(isnan(x_nan2zero)) = 0;
y_nan2zero = convn(x_nan2zero,f,'same');

%for each data point, find the sum of filter weights that were actually
% applied to it (e.g. points at the end of the data vector weren't 
% calculated using the full filter. Similar for points near NaNs)
x_validones = ones(size(x));
x_validones(isnan(x)) = 0;
y_validones = convn(x_validones,f,'same');

%y = y_nan2zero ./ y_validones;
% changed 2009-07-27 - now it *should* be mean-preserving...
y = y_nan2zero .* (sum(f) ./ y_validones);

if x_was_col_vector
    y = y';
    f = f';
end;
