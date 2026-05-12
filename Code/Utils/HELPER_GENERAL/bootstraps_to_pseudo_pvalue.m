function [pseudop,p_floor] = bootstraps_to_pseudo_pvalue(bootval,nullval,enforce_floor)
% [pseudop,p_floor] = bootstraps_to_psuedo_pvalue(bootvals,nullval)
% given a vector of bootstrap values, and a null value (default: 0),
% return a pseudo p value, i.e. the widest bootstrap CI that excludes the
% null value
%
% if enforce_pvalue_floor is true (default: true) do not allow the pseudo p
% value to be below 2/numel(bootval). E.g. if you do 200 bootstraps, you
% cannot get a pseudo p value below 0.01. The second output returns the
% floor p value that was used
%fron ESBM:
%I think it is possible to get something but it is not a proper p-value. It is something I have sometimes called a %pseudo p-value’, 
%based on bootstrap CIs. Basically pseduop = 0.01 means "The widest bootstrap confidence interval excluding 0 is at 0.99% confidence”.  
%  I write something in figure caption/methods like “*, **, *** indicates the 95%, 99%, 99.9% boostrap CI excludes 0". 
% Since these are technically slightly diff from a p-value (in some conditions they can be the same, 
% but other conditions they can be diff, e.g. cases with asymmetrical CIs).

if nargin < 2 || isempty(nullval)
    nullval = 0;
end
if nargin < 3 || isempty(enforce_floor) || enforce_floor
    p_floor = min(1,2 ./ numel(bootval));
else
    p_floor = 0;
end

pseudop = min(1,2*min(mean(bootval <= nullval),mean(bootval >= nullval)));

pseudop = max(pseudop,p_floor);