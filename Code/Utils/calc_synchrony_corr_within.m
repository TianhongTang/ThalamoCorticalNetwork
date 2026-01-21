function [r, err] = calc_synchrony_corr_within(raster_area)
% Calculate pairwise Pearson correlation within an area
N = size(raster_area, 1);
T = size(raster_area, 2);
corr_values = NaN(N, N);
for i = 1:N
    for j = i+1:N
        neuron1 = raster_area(i, :);
        neuron2 = raster_area(j, :);
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