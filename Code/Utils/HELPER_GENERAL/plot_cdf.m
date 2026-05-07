function [h] = plot_cdf(vals,varargin)
% [h] = plot_cdf(vals,...)
%  plots the CDF of a vector of values (excluding NaNs)
% 
% [h] = plot_vals(vals,'includenans',...)
%  same, but includes NaNs in determining the maximum level of the CDF
%  (e.g. if half the data is NaN, then the CDF will reach a maximum of 
%   0.5 instead of 1.0)

if ~isvector(vals)
    error('vals must be a vector');
end;

if ~isrowvector(vals)
    vals = vals';
end;

include_nans = 0;
if numel(varargin) >= 1 && ischar(varargin{1}) && strcmp(varargin{1},'includenans')
    include_nans = 1;
    varargin = varargin(2:end);
end;

num_nans = sum(isnan(vals));

vals = sort(vals(~isnan(vals)));
valquants = (1:length(vals)) ./ length(vals);

%replace vals([1 2 3 ...]) with vals([1 1 2 2 3 3 ...])
valids = [1:length(vals) ; 1:length(vals)];
vals = vals(valids(:));
valquants = valquants(valids);

%replace 1/N 1/N 2/N 2/N ... (N-1)/N (N-1)/N N/N N/N
% with 0 1/N 1/N ... (N-1)/N (N-1)/N N/N
valquants = [0 valquants(1:(end-1))];

% if included NaNs, scale down all y-variables proportionally
% (so that those data points are effectively 'missing', e.g. if half of the
% data in NaN then the y coords only go up to 0.5)
if include_nans
    valquants = valquants .* (numel(vals) ./ (numel(vals) + num_nans));
end;
if ~isempty(vals)
    h = plot(vals,valquants,varargin{:});
else
    h = [];
end;