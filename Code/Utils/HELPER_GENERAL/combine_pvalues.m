function [p] = combine_pvalues(pvals,dim)
%combine_pvalues(pvals)
%combine_pvalues(pvals,dim)
%computes a p-value from the product of the input p-values
%  along the dimension 'dim'
%  (if pvals is a vector, defaults to a produce along its length;
%   if pvals is a matrix, defaults to dimension = 1)
%
%e.g. say you have three datasets where condition A has a
% bigger mean than condition B, but only at p=.1 for each case.
% No single dataset is significant. But you don't think that
% A and B are the same, because if they were, it would be very unlikely to get 
% such low p-values three times in a row. You need a way to combine the
% three p-values together into an 'overall' p-value.
% 
%Apparently Fisher defined a test statistic for this, 
% -2*ln(p1*p2*p3*...*pn)
% if all null hypotheses were true, he says this'll have a chi-square 
% distribution with 2*n degrees of freedom. So, that's what this function
% does.
%
% See here:
% Fisher, R. A., Statistical Methods for Research Workers, London: Oliver and Boyd, 11th ed., 1950, 
% pages 103-105.
% http://www.haghish.com/resources/materials/Statistical_Methods_for_Research_Workers.pdf
% 
%This method has been used by neuroscientists: (Rolls ET, Xiang JZ. 2005)
%
%Apparently (Littell and Folks 1971) proved that this combined p-value is
% "Asymptotically Bahadour Optimal", which as far as I can tell, means that if the null
% hypothesis (that all p-vals are U(0,1)) is false, then as we gather more individual
% p-values, the combined p-value becomes small as fast as it possibly can
% (at least asymptotically, i.e. as the number of p-values we gather grows to
% infinity...is this meaningful for small sample sizes?)
% 
% (cited in "Truncated Product Method for Combining P-Values" by 
% D.V. Zaykin, Lev A. Zhivotovsky, P.H. Westfall, and B.S. Weir,
% Genetic Epidemiology 22:170–185 (2002)
% ...who end up inventing a different method more suited to when you're combining hundreds of
% p-values, most of which you think follow the null hypothesis of no effect; but for which you
% think a real effect is lurking in a small number of cases (e.g. whole-genome scans). Probably not suited
% to my case, where I'm combining ~2-6 p-values, or expect an effect in most cases.)

if nargin < 2
    if isrowvector(pvals)
        dim = 2;
    elseif iscolvector(pvals)
        dim = 1;
    else
        dim = 1;
    end;
end;

p = 1 - chi2cdf( -2*sum(log(pvals),dim), 2*size(pvals,dim));
