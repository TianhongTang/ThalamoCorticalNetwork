function h = errorshade(x, y, errorx, errory, plot_option, shade_color)
% errorshade by Vincent S. Huang
% usage: errorshade(x, y, errorx, errory, plot_option, shade_color)
%
% like errorbar, but instead of plotting bars, fills the the area enclosed by
% [x-errorx x+errorx, y-yerror, y+yerror] with shade_color, and plots the 
% x,y as a line or as indicated by plot_option

default_shade_color = [0.8 0.8 0.8];
if nargin < 6
    shade_color = default_shade_color;
end;

x_low = x;
y_low = y -  errory;
x_high = x;
y_high =  y + errory;
fill_error_area(x,y,x_low, x_high, y_low, y_high, shade_color);

if nargin > 3
    x_low = x - errorx;
    y_low = y;
    x_high = x + errory;
    y_high =  y;
end;
fill_error_area(x,y,x_low, x_high, y_low, y_high, shade_color);

if nargin < 5
    plot_option = '';
end;
plot(x,y, plot_option)

figure(gcf); 
hold_stat = ishold;
if hold_stat == 0
    hold off;
else 
    hold on;
end;

function fill_error_area(x,y,x_low, x_high, y_low, y_high, shade_color)
x_high = x_high(end:-1:1);
y_high = y_high(end:-1:1);

if size(x,1) ~= 1 
x_vec = [x_low; x_high];
y_vec = [y_low; y_high];
else
x_vec = [x_low x_high];
y_vec = [y_low y_high];
end;

x_vec = x_vec(find(~isnan(x_vec)));
y_vec = y_vec(find(~isnan(y_vec)));

figure(gcf); 
hold_stat = ishold;
hold on;
if nargin > 5
h = fill(x_vec, y_vec, shade_color);
else
h = fill(x_vec, y_vec, default_shade_color);
end;
set(h, 'LineStyle', 'none')
