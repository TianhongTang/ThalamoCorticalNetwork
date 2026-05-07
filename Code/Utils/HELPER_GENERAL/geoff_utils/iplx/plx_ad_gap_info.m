function [adfreq, n, ts, fn] = plx_ad_gap_info(filename, ch)
% PLX_AD_GAP_INFO  Mimic the behavior of Plexon's file reading function
%
% [adfreq, n, ts, fn] = plx_ad_gap_info(filename, ch)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 0-based channel number
%
%           a/d data come in fragments. Each fragment has a timestamp
%           and a number of a/d data points. The timestamp corresponds to
%           the time of recording of the first a/d value in this fragment.
%           All the data values stored in the vector ad.
% 
% OUTPUT:
%   adfreq - digitization frequency for this channel
%   n - total number of data points 
%   ts - array of fragment timestamps (one timestamp per fragment, in seconds)
%   fn - number of data points in each fragment

fullplx = iplx_mimichelper(filename);
chindx = find([fullplx.Spikes.Channel]==ch);
if isempty(chindx)
    error('No such channel in this data file.');
elseif numel(chindx)>1
    error('Data integrity failure; multiple channels with this number.');
end
adfreq = fullplx.Slow(chindx).SampleFreq;
ts = fullplx.Slow(chindx).Timestamps;
fn = cellfun(@numel,fullplx.Slow(chindx).Data);
n = sum(fn);