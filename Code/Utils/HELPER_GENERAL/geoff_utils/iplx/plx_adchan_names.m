function [n,names] = plx_adchan_names(filename)
% PLX_ADCHAN_NAMES  Mimic the behavior of Plexon's file reading function
% Read name for each a/d channel from a .plx file
%
% [n,names] = plx_adchan_names(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   names - array of a/d channel name strings
%   n - number of channels

fullplx = iplx_mimichelper(filename);
names = char({fullplx.Slow.Name});
n = size(names,1);