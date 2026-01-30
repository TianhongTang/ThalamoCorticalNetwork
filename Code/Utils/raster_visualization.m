function f=raster_visualization(data, colors, trial_borders, visible)
%%   Raster Visualization Function
%  This function takes in a 2D array of data and visualizes it as a raster plot.
%  Inputs:
%  - data: (N, T) array of raster. Values are binary (0 or 1).
%  - colors: (N, 3) array of RGB colors for each neuron.
%  - trial_borders: (1, M) array of time points indicating trial borders (optional).
%  - visible: 'on' or 'off' to control figure visibility (optional, default 'on').
%  Outputs:
%  - Figure handle of the raster plot.

%% Input handling
if nargin < 4
    visible = 'on';
end
if nargin < 3
    trial_borders = [];
end

%% Main
[N, T] = size(data);
f = figure("Visible", visible);
hold on;
for i = 1:N
    plot(data(i,:) + i - 0.5, 'Color', colors(i, :));
end
if exist('trial_borders', 'var')
    for b = 1:length(trial_borders)-1
        xline(trial_borders(b), 'k--', 'LineWidth', 1);
    end
end
hold off;
xlabel("Time (ms)");
ylabel("Neuron No.");