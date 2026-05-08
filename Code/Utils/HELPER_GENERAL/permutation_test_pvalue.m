function [p] = permutation_test_pvalue(stat,permstats,tail,proper_pvalues)
% [p] = permutation_test_pvalue(stat,permstats,tail,proper_pvalues)
%   compute p-value of a statistic 'stat', given a vector of samples
%    from the null distribution ('permstats').
%
% like calling
%  cdf_to_pvalue(empirical_cdf(stat,permstats),tail,nsamples)
%
%  if tail == [], defaults to 'both'
%  if proper_pvalues == [], defaults to true, i.e. computes nsamples based 
%   on the number of permutations. Otherwise, uses nsamples = inf.

if nargin < 3 || isempty(tail)
    tail = 'both';
end;
if nargin < 4 || isempty(proper_pvalues)
    proper_pvalues = true;
end;

y = empirical_cdf(stat,permstats);
nsamples = inf;
if proper_pvalues
    if iscell(permstats)
        nsamples = length(permstats);
    else
        nsamples = max(size(permstats) - size(stat) + 1);
    end;
end;

p = cdf_to_pvalue(y,tail,nsamples);
