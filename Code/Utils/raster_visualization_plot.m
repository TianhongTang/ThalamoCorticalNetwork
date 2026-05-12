function raster_visualization_plot(f, data, colors, trial_borders)
%%   Raster Visualization Function
%  This function takes in a 2D array of data and plot it as a raster plot on the given figure handle.
%  Inputs:
%  - data: (N, T) array of raster. Values are binary (0 or 1).
%  - colors: (N, 3) array of RGB colors for each neuron.
%  - trial_borders: (1, M) array of time points indicating trial borders (optional).
%  - visible: 'on' or 'off' to control figure visibility (optional, default 'on').
%  Outputs:
%  - Figure handle of the raster plot.

%% Input handling
% if nargin < 3
%     trial_borders = [];
% end
if nargin < 4
    trial_borders = [];
end

%% Main
[N, T] = size(data);
hold on;
for i = 1:N
    plot(f, data(i,:) + i - 0.5, 'Color', colors(i, :));
end
if exist('trial_borders', 'var')
    for b = 1:length(trial_borders)-1
        xline(f, trial_borders(b), 'k--', 'LineWidth', 1);
    end
end
hold off;
xlabel("Time (ms)");
ylabel("Neuron No.");