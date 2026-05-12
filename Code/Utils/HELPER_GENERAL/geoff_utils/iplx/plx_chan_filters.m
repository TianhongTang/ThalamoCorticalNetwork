function [n,filters] = plx_chan_filters(filename)
% PLX_CHAN_FILTERS  Mimic the behavior of Plexon's file reading function
%
% [n,filters] = plx_chan_filters(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
% 
% OUTPUT:
%   filter - array of filter values (0 or 1)
%   n - number of channels

fullplx = iplx_mimichelper(filename);
n = numel(fullplx.Spikes);
filters = double([fullplx.Spikes.Filter]');