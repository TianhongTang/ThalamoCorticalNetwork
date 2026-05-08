function [icolor] = interpcolor(color1,color2,alpha)
% [icolor] = interpcolor(color1,color2,alpha)
% interpolate between two colors, simply using:
%  icolor = color1.*((1-alpha).*[1 1 1]) + color2.*((alpha).*[1 1 1])
% default: alpha = .5

if nargin < 2
    color2 = color1;
end;
if nargin < 3
    alpha = .5;
end;
if isequal(color1,'none') || isequal(color2,'none')
    icolor = 'none';
else
    icolor = colorspec_to_rgb(color1).*((1-alpha).*[1 1 1]) + colorspec_to_rgb(color2).*((alpha).*[1 1 1]);
end;
