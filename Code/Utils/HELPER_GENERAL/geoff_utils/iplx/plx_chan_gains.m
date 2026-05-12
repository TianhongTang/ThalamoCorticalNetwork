function [n,gains] = plx_chan_gains(filename)
% PLX_CHAN_GAINS  Mimic the behavior of Plexon's .plx file reading function
% Read channel gains from .plx file
%
% [gains] = plx_chan_gains(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
% 
% OUTPUT:
%  gains - array of total gains
%   n - number of channels

if isempty(filename)
    [filename,pathname] = uigetfile('*.plx', 'Choose a PLX file');
    filename = fullfile(pathname, filename);
end

hdr = iplx_header(filename);
n = numel(hdr.SpikeHeaders);
gains = [hdr.SpikeHeaders.Gain]';