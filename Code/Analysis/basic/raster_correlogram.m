%% raster_correlogram.m - Compute and plot cross-correlogram between pairs of neurons from raster data.

clear;
%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
% include code folder and utils
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Main

% Example session: Slayer Mus 6. Selected t_range. Justification: Best sleep period for Pre-eyeclose.
meta = struct();
meta.animal_name = 'Slayer';
meta.injection = 'Muscimol';
meta.align = 'Last';
meta.session_idx = 6;
meta.resting_dur_threshold = 15;
t_range = 1:60000;


% areas    = {'Full',     'Full',      'Cortex',   'Cortex'};
areas    = {'Cortex',   'Cortex',    'Cortex',   'Cortex'};
preposts = {'Pre',      'Pre',       'Post',     'Post'};
states   = {'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};

%% Raster
f = figure();
tiles = tiledlayout(4, 2, "TileSpacing", "Compact", "Padding", "Compact");

% selected_neurons = 1:95;
% selected_neurons = [16:21, 39, 44, 61:67, 83:93];
% selected_neurons = [14:19, 33:40, 42, 44, 46, 49, 58:73, 83:95];
selected_neurons = [15, 73];
% selected_neurons = [16, 17];

for i = 1:length(areas)
    % Load data
    meta.area = areas{i};
    meta.prepost = preposts{i};
    meta.state = states{i};
    meta.filename = generate_filename('raster', meta);

    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    % Plot
    % N = raster_data.meta.N;

    cell_area = raster_data.data.cell_area;
    tile = nexttile(i*2-1);
    % tile = nexttile(i);
    raster = raster_data.data.rasters{1}(:, 1:60000);

    cell_area = cell_area(selected_neurons);
    raster = raster(selected_neurons, :);

    N = numel(selected_neurons);
    colors = zeros(N, 3);
    for j = 1:N
        switch cell_area{j}
            case 'Thalamus'
                colors(j, :) = [1, 0, 1];
            case 'ACC'
                colors(j, :) = [0, 0, 1]; % blue
            case 'VLPFC'
                colors(j, :) = [1, 0, 0]; % red
            otherwise
                colors(j, :) = [0, 0, 0]; % black
        end
    end
    raster_visualization_plot(tile, raster, colors)
    title(sprintf('%s-%s, %s', meta.prepost, meta.injection, meta.state));
    xlabel('Time (ms)');
    ylabel('Neuron No.');
    % ylim([-1.5, 4.5]);
    % ylim([-1.5, 4.5]);
end

%% Pairwise correlogram
% f = figure();
% tiles = tiledlayout(2, 2, "TileSpacing", "Compact", "Padding", "Compact");

% selected_neurons = [85, 86];
% selected_neurons = [15, 73];
selected_neurons = [35, 60];
smooth_window = 30; % ms

for i = 1:length(areas)
    % Load data
    meta.area = areas{i};
    meta.prepost = preposts{i};
    meta.state = states{i};
    meta.filename = generate_filename('raster', meta);

    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    raster = raster_data.data.rasters{1};

    r1 = raster(selected_neurons(1), 1:60000);
    r2 = raster(selected_neurons(2), 1:60000);

    % Compute correlogram
    [correlogram, lags] = norm_xcorr(r1, r2, 1000);
    smooth_correlogram = movmean(correlogram, smooth_window);

    % Compute shuffled correlogram for control
    shuffle_N = 10;
    shuffle_correlograms = zeros(shuffle_N, length(correlogram));
    for j = 1:shuffle_N
        fprintf('%d/%d shuffles started\n', j, shuffle_N);
        shuffled_r2 = r2(randperm(length(r2)));
        [shuffle_corr, ~] = norm_xcorr(r1, shuffled_r2, 1000);
        shuffle_correlograms(j, :) = movmean(shuffle_corr, smooth_window);
        fprintf('%d/%d shuffles finished\n', j, shuffle_N);
    end
    shuffle_mean = mean(shuffle_correlograms, 1);
    shuffle_std = std(shuffle_correlograms, [], 1);
    % shuffle_mean = movmean(shuffle_mean, smooth_window);
    % shuffle_std = movmean(shuffle_std, smooth_window);

    % Plot
    nexttile(i*2);
    std_multiplier = 2;
    shuffle_upper = shuffle_mean + std_multiplier * shuffle_std;
    shuffle_lower = shuffle_mean - std_multiplier * shuffle_std;
    % plot shuffled correlogram as a shaded area
    fill([lags, fliplr(lags)], [shuffle_upper, fliplr(shuffle_lower)], [0.8, 0.8, 0.8], 'EdgeColor', 'none');

    % plot correlogram
    hold on;
    xline(0, 'k--');
    yline(1, 'k--');
    plot(lags, smooth_correlogram);
    hold off;

    title(sprintf('Correlogram: %s-%s, %s', meta.prepost, meta.injection, meta.state));
    xlabel('Lag (ms)');
    ylabel('Normalized correlation');
    legend({'Shuffled 2SD', '', '', 'correlation'});
    ylim([-1, 5]);
end

function [correlogram, lags] = norm_xcorr(r1, r2, max_lag)
    % Compute normalized cross-correlogram, not only by the number of overlapping points, but also by the overall firing rates of the two neurons.
    % Compute raw cross-correlogram
    [raw_correlogram, lags] = xcorr(r1, r2, max_lag, 'unbiased');
    ave_rate = zeros(size(lags));
    for i = 1:length(lags)
        lag = lags(i);
        if lag < 0
            ave_rate(i) = mean(r1(1:end+lag)) * mean(r2(1-lag:end));
        else
            ave_rate(i) = mean(r1(1+lag:end)) * mean(r2(1:end-lag));
        end
    end
    valid_filter = ave_rate > 0;
    correlogram = zeros(size(raw_correlogram));
    correlogram(valid_filter) = raw_correlogram(valid_filter) ./ ave_rate(valid_filter);

    if any(~valid_filter)
        correlogram(~valid_filter) = NaN;
        warning('Some lags have zero average firing rate product, resulting in NaN in the normalized correlogram.');
    end
end