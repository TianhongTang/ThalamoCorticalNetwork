function  [n,adchans] = plx_ad_chanmap(filename)
% PLX_AD_CHANMAP  Mimic the behavior of Plexon's .plx file reading function
% 
% [n,adchans] = plx_ad_chanmap(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
% OUTPUT:
%   n - number of continuous channels
%   adchans - 1 x n array of continuous channel numbers
%
% Normally, there is one channel entry in the .plx for for each raw continuous channel,
% so the mapping is trivial adchans[i] = i-1 (because continuous channels start at 0).
% However, for certain .plx files saved in some ways from OFS (notably after
% loading data files from other vendors), the mapping can be more complex.
% E.g. there may be only 2 non-empty channels in a .plx file, but those channels
% correspond to raw channel numbers 7 and 34. So in this case NChans = 2, 
% and adchans[1] = 7, adchans[2] = 34.
% The plx_ routines that return arrays always return arrays of size NChans. However,
% routines that take channels numbers as arguments always expect the raw  
% channel number.  So in the above example, to get the data from  
% the second channel, use
%   [adfreq, n, ts, fn, ad] = plx_ad(filename, adchans[2])

fullplx = iplx_mimichelper(filename);
n = numel(fullplx.Slow);
adchans = reshape([fullplx.Slow.Channel],1,[]);