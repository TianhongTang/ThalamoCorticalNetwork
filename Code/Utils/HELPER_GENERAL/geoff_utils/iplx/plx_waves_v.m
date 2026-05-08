function [n, npw, ts, wave] = plx_waves_v(filename, ch, u)
% PLX_WAVES_V  Mimic the behavior of Plexon's .plx file reading function
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
%   wave - array of waveforms [npw, n] converted to mV

fullplx = iplx_mimichelper(filename);
[ts,wave] = iplx_getunit(fullplx,ch,u);
wave = wave';
[n,npw] = size(wave);