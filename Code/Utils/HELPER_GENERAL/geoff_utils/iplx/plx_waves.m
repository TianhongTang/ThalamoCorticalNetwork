function [n, npw, ts, wave] = plx_waves(filename, ch, u)
% PLX_WAVES  Mimic the behavior of Plexon's .plx file reading function
% Read waveform data from a .plx file
%
% [n, npw, ts, wave] = plx_waves(filename, channel, unit)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 1-based channel number
%   unit  - unit number (0- unsorted, 1-4 units a-d)
% OUTPUT:
%   n - number of waveforms
%   npw - number of points in each waveform
%   ts - array of timestamps (in seconds) 
%   wave - array of waveforms [npw, n], raw a/d values

fullplx = iplx_mimichelper(filename);
chindx = find([fullplx.Spikes.Channel]==ch);
if isempty(chindx)
    error('IPLX:notChannel','No such channel in this data file.');
elseif numel(chindx)>1
    error('IPLX:tooManyChannels', ...
        'Data integrity failure; multiple channels with this number.');
end
[ts, wave] = iplx_getunit(fullplx,ch,u);
wave = wave'/fullplx.Spikes(chindx).ADConversion;
[n,npw] = size(wave);