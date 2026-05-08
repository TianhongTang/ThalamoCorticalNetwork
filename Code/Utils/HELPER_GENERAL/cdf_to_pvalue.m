function [p] = cdf_to_pvalue(y,tail,nsamples)
% [p] = cdf_to_pvalue(y,tail,proper_pvalues)
%  given the null distribution's cdf value at each observed test statistic,
%   convert them to p-values.
%
%   optional arguments:
%   tail - 'both' or 0: two-tailed test, 
%          'left' or '-' or -1: test the hypothesis stat < permstats
%          'right' or '+' or +1: test stat > permstats
%    (default: 'both')
%
%   nsamples - if the cdf was approximated using a set of samples 
%    from the null distribution, then we shouldn't report
%    p-values lower than 1/nsamples (or 2/nsamples for a two-tailed test).
%    (e.g. if do a two-tailed test where we approximate the null distribution 
%     using 100 samples, and *all* of the samples are larger than 
%     the observed statistic, then instead of reporting p == 0, 
%     we should report p < .02. Thus this function will report p = .02)
%    (default: inf (i.e. assuming cdf was computed analytically))
%
%   both tail and nsamples may be specified as matrices of size equal to y

if nargin < 2
    tail = 0;
end;
if nargin < 3
    nsamples = inf;
end;

if ischar(tail)
    if strcmp(tail,'left')
        tail = -1;
    elseif strcmp(tail,'right')
        tail = 1;
    elseif strcmp(tail,'both')
        tail = 0;
    else
        error('unknown tail');
    end;
end;
tail = sign(tail);


p = nans(size(y));

% for speed purposes (according to profiler), 
% split up function to allow special cases
% when tail/nsamples are scalar
if isscalar(tail)
    if tail < 0
        p = max(y, 1 ./ nsamples);
    elseif tail > 0
        p = max(1 - y, 1 ./ nsamples);
    else
        p = max(2*min(y,1-y), 2 ./ nsamples);
    end;
else
    ids_left = tail < 0;
    ids_right = tail > 0;
    ids_both = tail == 0;
    
    if isscalar(nsamples)
        if any(ids_left)
            p(ids_left) = max(y(ids_left), 1./nsamples);
        end;
        if any(ids_right)
            p(ids_right) = max(1 - y(ids_right), 1./nsamples);
        end;
        if any(ids_both)
            p(ids_both) = max(2*min(y(ids_both),1 - y(ids_both)), 2./nsamples);
        end;
    else
        if any(ids_left)
            p(ids_left) = max(y(ids_left), 1./nsamples(ids_left));
        end;
        if any(ids_right)
            p(ids_right) = max(1 - y(ids_right), 1./nsamples(ids_right));
        end;
        if any(ids_both)
            p(ids_both) = max(2*min(y(ids_both),1 - y(ids_both)), 2./nsamples(ids_both));
        end;
    end;
end;
