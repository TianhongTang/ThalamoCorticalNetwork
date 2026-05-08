function new_color = convertColorToL48D(color)
% Converts a color of format R to a
% [R G -] using the MSB of R and G

new_color = [0 0 0];
new_color(1) = bitand(color, 240);
new_color(2) = bitshift(bitand(color, 15), 4);

return
