function [v] = iscolorspec(cs)
% [v] = colorspec(cs)
%  returns true if cs is a color specification
%   (a string like 'm' or 'magenta', or an RGB triplet)
%
%  note: if cannot interpret as a colorspec,
%  strips it down to the characters 'ymcrgbkw' and tries again
%  (this is to handle arguments to 'plot' correctly)
%  So ':go' counts as the colorspec 'g'

v = true;

try
    rgb = colorspec_to_rgb(cs);
catch
    v = false;
end;
