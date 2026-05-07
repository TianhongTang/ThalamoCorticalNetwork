function [p,index,permindexes] = roc_vs_roc_permutation_test(noise1,signal1,noise2,signal2,nperms)
% [p,index] = roc_vs_roc_permutation_test(noise1,signal1,noise2,signal2,nperms)
%
% calculate ROC area (signal vs noise) separately for two conditions (1 vs
% 2). Return an index of the differnece between ROC areas (ROC1 - ROC2)
% and a p-value against the null hypothesis that signal1 and signal2 are
% exchangeable and noise1 and noise2 are exchangeable.
%
% input:
%  noise1 - noise data from condition 1
%  signal1 - signal data from condition 1
%  noise2 - noise data from condition 2
%  signal2 - signal data from condition 2
%  nperms - number of permutation to use for shuffling
%
% output:
%  p - permutation p-value
%  index - true index (ROC area 1 - ROC area 2)
%  permindexes - vector of nperms indexes from permuted datasets

% calculate Index and its p-value
rate1sig = signal1;
rate1noi = noise1;

rate2sig = signal2;
rate2noi = noise2;

% construct proper input to permutation_test_general
% to allow a permutation-based test of this index
% (which shuffles the correct conditions to represent the null
% hypothesis that the uncertainty signal does not change
% between the two conditions+epochs being compared)
addset = @(rate,lab,permgroup,newrate,newlab,newpermgroup) deal([rate ; newrate],[lab ; newlab*ones(numel(newrate),1)],[permgroup ; newpermgroup*ones(numel(newrate),1)]);

rate = [];
lab = [];
permgroup = [];
[rate,lab,permgroup] = addset(rate,lab,permgroup,rate1noi,1,1);
[rate,lab,permgroup] = addset(rate,lab,permgroup,rate1sig,2,2);
[rate,lab,permgroup] = addset(rate,lab,permgroup,rate2noi,3,1);
[rate,lab,permgroup] = addset(rate,lab,permgroup,rate2sig,4,2);
[p,index,permindexes] = permutation_test_general(rate,lab,permgroup, ...
    @(rate,lab,permgroup) rocarea3(rate(lab==1),rate(lab==2))-rocarea3(rate(lab==3),rate(lab==4)), ...
    nperms);
