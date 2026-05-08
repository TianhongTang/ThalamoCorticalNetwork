function [h,b] = getSpikeHist(data,ntrials,pars)

%
% getSpikeHist 
% USAGE: getSpikeHist(data, [min max bw]);
%
% Generate PSTH
%
% DAL

bins = [pars(1)-pars(3):pars(3):pars(2)+pars(3)];
[h,b] = hist(data,bins);
h = h(2:length(h)-2);
b = b(2:length(b)-2);

h = (1000.*h)/(ntrials*pars(3));
