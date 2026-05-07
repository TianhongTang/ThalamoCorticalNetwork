function [n, ts] = plx_ts(filename, channel, unit)
% PLX_TS  Mimic the behavior of Plexon's .plx file reading function
% Read spike timestamps from a .plx file
%
% [n, ts] = plx_ts(filename, channel, unit)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 1-based channel number
%   unit  - unit number (0- unsorted, 1-4 units a-d)
% OUTPUT:
%   n - number of timestamps
%   ts - array of timestamps (in seconds)

fullplx = iplx_mimichelper(filename);
chindx = find([fullplx.Spikes.Channel]==channel);
if isempty(chindx)
    error('No such channel in this data file.');
elseif numel(chindx)>1
    error('Data integrity failure; multiple channels with this number.');
end
isunit = fullplx.Spikes(chindx).Units==unit;
ts = fullplx.Spikes(chindx).Timestamps(isunit);
n = numel(ts);