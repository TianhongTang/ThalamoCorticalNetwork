function imageL48D = convertIndexedImageToL48DImage(image)
% Takes an images in index format (indexes to a CLUT), and return it
% as an indexed image using R and G MSB to avoid dithering.

image_size = size(image);
imageL48D = zeros(image_size(1), image_size(2), 3, 'uint8');
imageL48D(:, :, 1) = bitand(image(:, :), 240); % Red MSB
imageL48D(:, :, 2) = bitshift(bitand(image(:, :), 15), 4); % Green MSB
return 

