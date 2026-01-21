function [r, err] = calc_synchrony_corr_across(raster_area1, raster_area2)
% Calculate pairwise Pearson correlation across two areas

N1 = size(raster_area1, 1);
N2 = size(raster_area2, 1);
assert(size(raster_area1, 2) == size(raster_area2, 2), 'Rasters must have the same number of time points.');
T = size(raster_area1, 2);

corr_values = zeros(N1, N2);
for i = 1:N1
    for j = 1:N2
        neuron1 = raster_area1(i, :);
        neuron2 = raster_area2(j, :);
        if std(neuron1) == 0 || std(neuron2) == 0
            corr_values(i, j) = NaN; % undefined correlation
        else
            corr_values(i, j) = corr(neuron1', neuron2');
        end
    end
end
r = mean(corr_values(:), 'omitnan');
err = std(corr_values(:), 0, 'omitnan') / sqrt(sum(~isnan(corr_values(:))));

end