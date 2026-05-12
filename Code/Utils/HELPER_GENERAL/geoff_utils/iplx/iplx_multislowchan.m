function [slow, ts] = iplx_multislowchan(fullplx)
% IPLX_MULTISLOWCHAN  Collect slow data from multiple channels
% Usage:
%   [slow, ts] = iplx_multislowchan(fullplx)
% When possible, IPLX_MULTISLOWCHAN collects all the data from multiple
% slow channels in an experiment and joins it together as a single matrix,
% (# of samples)-by-(# of channels), each column representing the full data
% collected on a single channel. If the experiment contains more than one
% block of slow data per channel, or if the timing is not aligned across
% channels properly, IPLX_MULTISLOWCHAN will generate an error. Otherwise,
%   slow    contains the multichannel data block
%   ts      is the timestamp of the beginning of the data block (make sure
%           to take this value into account when comparing with spike or
%           event data).

fullplx = iplx_tidychans(fullplx);
if any(cellfun(@numel,{fullplx.Slow.Timestamps})~=1)
    error('Multiple slow data blocks, cannot collect');
end

allts = [fullplx.Slow.Timestamps];
if range(allts)~=0
    error('Slow data channels start at different times, cannot collect');
end
ts = allts(1);

allslow = [fullplx.Slow.Data];
if range(cellfun(@numel,allslow))~=0
    error('Slow data channels end at different times, cannot collect');
end
slow = [allslow{:}];