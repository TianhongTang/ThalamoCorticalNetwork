function [n,gains] = plx_adchan_gains(filename)
% PLX_ADCHAN_GAINS  Mimic the behavior of Plexon's file reading function
% Read analog channel gains from .plx file
%
% [n,gains] = plx_adchan_gains(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
% 
% OUTPUT:
%  gains - array of total gains
%  n - number of channels

if isempty(filename)
    [filename,pathname] = uigetfile('*.plx', 'Choose a PLX file');
    filename = fullfile(pathname, filename);
end

hdr = iplx_header(filename);
n = numel(hdr.SlowHeaders);
gains = ([hdr.SlowHeaders.Gain].*[hdr.SlowHeaders.PreAmpGain])';