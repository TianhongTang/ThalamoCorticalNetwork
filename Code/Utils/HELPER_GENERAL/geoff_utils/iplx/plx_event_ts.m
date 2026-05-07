function [n, ts, sv] = plx_event_ts(filename, ch)
% PLX_EVENT_TS Mimic the behavior of Plexon's .plx file reading function
%
% [n, ts, sv] = plx_event_ts(filename, channel)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 1-based external channel number
%             strobed channel has channel number 257  
% OUTPUT:
%   n - number of timestamps
%   ts - array of timestamps (in seconds)
%   sv - array of strobed event values (filled only if channel is 257)

fullplx = iplx_mimichelper(filename);
chindx = find([fullplx.Events.Channel]==ch);
if isempty(chindx)
    error('No such channel in this data file.');
elseif numel(chindx)>1
    error('Data integrity failure; multiple channels with this number.');
end
ts = fullplx.Events(chindx).Timestamps;
n = numel(ts);
sv = fullplx.Events(chindx).Strobed;