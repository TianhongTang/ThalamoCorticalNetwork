function raster_visualization(data)
%%   Raster Visualization Function
%  This function takes in a 2D array of data and visualizes it as a raster plot.
%  Inputs:
%  - data: (N, T) array of raster. Values are binary (0 or 1).
%  Outputs:
%  - None. Displays a raster plot.

%% Main
[N, T] = size(data);
figure;
hold on;
for i = 1:N
    plot(data(i,:) + i - 0.5);
end
hold off;
xlabel("Time (ms)");
ylabel("Neuron No.")