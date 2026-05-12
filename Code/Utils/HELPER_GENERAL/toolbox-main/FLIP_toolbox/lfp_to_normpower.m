function [relpow, nonnormpow] = lfp_to_normpower(lfp)

% lfp should be in format nchans x trialtime x ntrials
% assumes lfp data in 1 kHz resolution
% require fieldtrip toolbox

foi = 1:150; % frequencies of interest

nchannels = size(lfp, 1); ntrials = size(lfp, 3); ntrialtime =  size(lfp, 2);
% prep ft format
dataLFP = [];
for t = 1:ntrials
    dataLFP.time{t,1} = 0.001:0.001:ntrialtime*0.001; % in seconds
    dataLFP.trial{t,1} = zeros(nchannels, ntrialtime); %nchannels x ntrialtime
    dataLFP.trial{t,1} = lfp(:, :, t);

end
dataLFP.fsample = 1000;
dataLFP.label = {};
%make each channel unique
for c = 1:nchannels
    dataLFP.label{c} = ['ch' num2str(c)];
end
dataLFP.trialinfo = zeros(ntrials, 1);
dataLFP.trialinfo(:, 1) = 1:ntrials;


cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'pow';
cfg.keeptrials = 'yes';
cfg.foi = foi;
cfg.pad = 'nextpow2';
pow = ft_freqanalysis(cfg, dataLFP);
meanpow = squeeze(mean(pow.powspctrm));
relpow = meanpow ./ max(meanpow);
nonnormpow = meanpow;
end
