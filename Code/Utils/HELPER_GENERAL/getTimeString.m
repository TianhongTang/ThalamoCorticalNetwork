function tstr = getTimeString(t)
% PURPOSE : To get a string representing time.
% USAGE : tstr = getTimeString([t])
% ARG :   't' is a return value of 'clock'.
% VERSION : 1.00  14-May-2000 YM

if nargin < 1, t = clock; end
t = fix(t);
if length(t) == 1
  h = fix(t/3600);
  m = mod(fix(t/60),60);
  s = mod(t,60);
else
  h = t(4);
  m = t(5);
  s = t(6);
end
tstr = sprintf('%02d:%02d:%02d',h,m,s);
