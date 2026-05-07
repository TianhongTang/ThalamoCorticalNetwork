function [ts, waves] = iplx_getunit(fullplx,channel,unit)
% IPLX_GETUNIT  Extract spikes from a specific channel and unit
% Usage:
%   [ts, waves] = iplx_getunit(fullplx,chan,unit)
% 
%   ts     is a vector of spike timestamps
%   waves  is a (# samples/wave)-by-(# spikes) matrix, each column
%          containing the waveform (in mV) of a single spike from the
%          specified unit on the specified channel

chindx = find([fullplx.Spikes.Channel]==channel);
if isempty(chindx)
    error('IPLX:notChannel','No such channel in this data file.');
elseif numel(chindx)>1
    error('IPLX:tooManyChannels', ...
        'Data integrity failure; multiple channels with this number.');
end
isunit = fullplx.Spikes(chindx).Units==unit;
ts = fullplx.Spikes(chindx).Timestamps(isunit);
waves = fullplx.Spikes(chindx).Waveforms(:,isunit);