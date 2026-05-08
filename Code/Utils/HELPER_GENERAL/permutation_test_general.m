function [pval,truestat,permstats,perms] = permutation_test_general(data,label,permset,statfxn,nperms,tail,cutoff_low_pvalues,perms)
% [pval,truestat,permstats,perms] = permutation_test_general(data,label,permset,statfxn,nperms,tail,cutoff_low_pvalues,perms)
%
% input:
%  data - vector of data values
%
%  label - vector of labels sorting the data into groups. The permuter
%   doesn't touch the actual data, it just shuffles these labels.
%    e.g. to compare means of data point that have diff labels, try the following:
%     statfxn = @(data,label,pset) mean(data(label==1)) - mean(data(label==2))
%
%  permset - vector. the permuter is only allowed to shuffle labels if they have
%   the same 'permset' value. e.g. if you want to compare the distribution
%   of firing rate differences between two conditions, you would assign one label
%   to each condition, then shuffle labels only *within neurons* so that
%   you can calculate within-neuron modulation indexes
%
% 'statfxn' is the function used to compute the test statistic. 
%   it must accept three inputs: (data, label, and permset)
%   and it must return a scalar output.
%   (default: returns the average of the pairwise differences between 
%     the "means" of the labeled groups, where the "mean" of a labeled group 
%     is not the true mean, but rather is weighted so that each permutation 
%     set contributes equally (as long as the permutation set contains 
%     at least one data point with that label).
%    In a way, this tests whether the labels are ordered so that 
%     small labels are assigned to small data values
%  
%    e.g., if data have labels in {1,2,3}, and permutation sets in {8,9},
%     and Mij denotes the mean of data with label i and permutation set j,
%     then this returns:
%    ((mean([M38,M39])-mean([M28,M29]))
%    +(mean([M38,M39])-mean([M18,M19]))
%    +(mean([M28,M29])-mean([M18,M19]))) / 3
%   )
%
% 'nperms' is the number of permutations to use (default: 500)
%
% 'tail' is whether the test is 'left', 'right', or 'both' -tailed
%  (default: 'both')
%
% 'cutoff_low_pvalues' if true, doesn't allow p-values to go below
%  the lowest possible value that can be calculated by this test
%  (1/nperms for a 1-tailed test, 2/nperms for a 2-tailed test).
%  After all, if the true statistic's value is smaller than 100 permuted statistics
%   for a 2-tailed test, then it is wrong to say p=0. Rather, p<.01, but
%   we don't know *how* much less than .01. So, the function will
%   just report p=.01
%  (default=true)
%
% 'perms' is a matrix of permuted labels
%  it has size [length(label) nperms]
%  perms(v,p) is the label of the v-th observed variable
%   in the p-th permutation.
%  (default: use a randomly generated permutation matrix)
%
%
% Written by ESBM
% edited 2018-01-12 to reduce dependencies on my custom code

if nargin < 2
    error('need at least 2 arguments');
end

is_default_arg = @(x) ischar(x) && strcmp(x,'auto');

if ~isvector(data) || ~isvector(label) error('inputs must be vectors'); end

data = data(:);
label = label(:);
ulabels = unique(label);

if nargin < 3 || is_default_arg(permset) permset = ones(size(data)); end

upermsets = unique(permset);

if nargin < 4 || is_default_arg(statfxn) statfxn = @(d,lab,p) default_statfxn(d,lab,p); end
if nargin < 5 || is_default_arg(nperms) nperms = 500; end
if nargin < 6 || is_default_arg(tail) tail = 'both'; end
if nargin < 7 || is_default_arg(cutoff_low_pvalues) cutoff_low_pvalues = true; end
if nargin < 8 || is_default_arg(perms)
    perms = nan*ones(length(data),nperms);
    % generate permutations by shuffling labels within each unique permset
    for u = 1:length(upermsets)
        upermids = permset == upermsets(u);
        upermlabels = label(upermids);
        for p = 1:nperms
            perms(upermids,p) = upermlabels(randperm(length(upermlabels))');
        end
    end
else
    %verify that the uesr-input permutations are valid
    if any(size(perms) ~= [length(data) nperms])
        error(['user-input permutations have the wrong size, should have size ' num2str([length(data) nperms])]);
    end
    for ups = 1:length(upermsets)
        upermids = permset == upermsets(ups);
        upermlabels = label(upermids);
        for ul = 1:length(ulabels)
            for p = 1:nperms
                if sum(perms(upermids,p) == ulabels(ul)) ~= sum(upermlabels == ulabels(ul))
                    error(['user-input permutation ' num2str(p) ' has incorrect number of label ' num2str(ulabels(ul)) ' in permset ' upermsets(ups)]);
                end
            end
        end
    end
end

% calculate statistic on true data 
truestat = statfxn(data,label,permset);

% calculate statistic after permuting labels
permstats = nan(1,nperms);
for p = 1:nperms
    permstats(p) = statfxn(data,perms(:,p),permset);
end

%calc p-value
left_pval = sum(permstats <= truestat) / nperms; 
right_pval = sum(permstats >= truestat) / nperms; 

if strcmp(tail,'left')
    pval = left_pval;
elseif strcmp(tail,'right')
    pval = right_pval;
else %strcmp(tail,'both')
    pval = min(1,2*min(left_pval,right_pval));
end

if cutoff_low_pvalues
    if strcmp(tail,'both')
        pval = max(pval,2/nperms);
    else
        pval = max(pval,1/nperms);
    end
end


%--------------------------------------------------------------------------
% function stat = default_statfxn(d,lab,p)
% %returns the average of the pairwise absolute differences between 
% % the means of the labeled groups
% 
% if length(ulabels) <= 1
%     stat = 0;
% else
%     umeans = nan*ones(size(ulabels));
%     for u = 1:length(ulabels)
%         umeans(u) = nanmean(d(lab == ulabels(u)));
%     end
%     sumdiffs = 0;
%     for u = 1:length(ulabels)
%         for u2 = 1:(u-1)
%             sumdiffs = sumdiffs + abs(umeans(u) - umeans(u2));
%         end
%     end
% 
%     stat = sumdiffs / (length(ulabels)*(length(ulabels)-1) / 2);
% end

%--------------------------------------------------------------------------
function stat = default_statfxn(d,lab,p)
%returns the average of the pairwise differences between 
% the "means" of the labeled groups,
% where the "mean" of a labeled group is not the true mean,
% but rather is normalized so that each permutation set contributes
% equally (as long as it contains at least one data point with that label)

ulabels = unique(lab);
upermsets = unique(p);
if length(ulabels) <= 1
    stat = 0;
else
    umeans = nan*ones(length(ulabels),length(upermsets));
    for ulab = 1:length(ulabels)
        for ups = 1:length(upermsets)
            umeans(ulab,ups) = nanmean(d(lab == ulabels(ulab) & p == upermsets(ups)));
        end
    end
    umeans = nanmean(umeans,2);
    sumdiffs = 0;
    for ulab = 1:length(ulabels)
        for ulab2 = 1:(ulab-1)
            sumdiffs = sumdiffs + (umeans(ulab) - umeans(ulab2));
        end
    end

    stat = sumdiffs / (length(ulabels)*(length(ulabels)-1) / 2);
end
