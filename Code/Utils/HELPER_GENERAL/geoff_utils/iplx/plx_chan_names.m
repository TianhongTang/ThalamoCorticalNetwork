function [n,names] = plx_chan_names(filename)
% PLX_CHAN_NAMES  Mimic the behavior of Plexon's file reading function
% Read name for each spike channel from a .plx file
%
% [n,names] = plx_chan_names(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   names - array of channel name strings
%   n - number of channels

fullplx = iplx_mimichelper(filename);
names = char({fullplx.Spikes.Name});
n = size(names,1);