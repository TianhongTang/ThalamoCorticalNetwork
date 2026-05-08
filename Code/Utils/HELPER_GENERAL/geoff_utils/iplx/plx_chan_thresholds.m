function [n,thresholds] = plx_chan_thresholds(filename)
% PLX_CHAN_THRESHOLDS  Mimic the behavior of Plexon's file reading function
%
% [n,thresholds] = plx_chan_thresholds(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   thresholds - array of tresholds, expressed in raw A/D counts
%   n - number of channel

fullplx = iplx_mimichelper(filename);
n = numel(fullplx.Spikes);
thresholds = [fullplx.Spikes.Threshold]';