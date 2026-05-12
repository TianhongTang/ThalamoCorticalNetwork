function [adfreq, n, ad] = plx_ad_span(filename, ch, startCount,endCount)
% PLX_AD_SPAN  Mimic the behavior of Plexon's .plx file reading function
% Read a span of a/d data from a .plx file
%
% [adfreq, n, ad] = plx_ad_span(filename, ch,startCount,endCount)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   startCount - index of first sample to fetch
%   endCount - index of last sample to fetch
%   channel - 0 - based channel number
%
% OUTPUT:
%   adfreq - digitization frequency for this channel
%   n - total number of data points 
%   ad - array of raw a/d values

[adfreq, n, ts, fn, ad] = plx_ad(filename, ch);
ad = ad(startCount:endCount);