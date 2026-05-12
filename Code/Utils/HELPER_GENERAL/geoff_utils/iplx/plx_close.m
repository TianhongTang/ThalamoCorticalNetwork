function [n] = plx_close(filename)
% PLX_CLOSE  Clean up the mimic function persistent workspace
% 
% It is very important that you call this function after you are done using
% the Plexon mimic functions to read in your file! Otherwise, the entire
% file (which may be several hundred megabytes) remains in persistent
% memory.
%
% [n] = plx_close(filename)
%
% INPUT:
%   filename - if empty string, will close any open files
% OUTPUT:
%   n - always 0

n = 0;
iplx_mimichelper(filename,true);