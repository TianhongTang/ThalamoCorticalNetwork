function [n,freqs] = plx_adchan_freqs(filename)
% PLX_ADCHAN_FREQS  Mimic the behavior of Plexon's file reading function
%
% [n,freqs] = plx_adchan_freqs(filename)
%
% OUTPUT:
%   freqs - array of frequencies
%   n - number of channels

fullplx = iplx_mimichelper(filename);
n = numel(fullplx.Slow);
freqs = [fullplx.Slow.SampleFreq]';