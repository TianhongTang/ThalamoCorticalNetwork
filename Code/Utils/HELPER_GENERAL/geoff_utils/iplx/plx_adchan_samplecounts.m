function [n,samplecounts] = plx_adchan_samplecounts(filename)
% PLX_ADCHAN_SAMPLECOUNTS  Mimic the behavior of Plexon's function
% Read the per-channel sample counts for analog channels from a .plx file
%
% [n,samplecounts] = plx_adchan_samplecounts(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   n - number of channels
%   samplecounts - array of sample counts

fullplx = iplx_mimichelper(filename);
n = numel(fullplx.Slow);
samplecounts = zeros(n,1);
for k=1:n
    samplecounts(k) = sum(cellfun(@numel,fullplx.Slow(k).Data));
end