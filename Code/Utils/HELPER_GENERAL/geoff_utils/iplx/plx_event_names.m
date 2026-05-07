function [n,names] = plx_event_names(filename)
% PLX_EVENT_NAMES  Mimic the behavior of Plexon's file reading function
% Read name for each event type from a .plx file
%
% [n,names] = plx_event_names(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
% 
% OUTPUT:
%   names - array of event name strings
%   n - number of channels

fullplx = iplx_mimichelper(filename);
names = char({fullplx.Events.Name});
n = size(names,1);