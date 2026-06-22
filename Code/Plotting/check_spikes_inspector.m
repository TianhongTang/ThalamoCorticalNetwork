%% check_spikes_inspector.m - targeted raster + area firing-rate inspector
%
% This inspector version generates one targeted figure only.
% It plots:
%   1. Raster plot for the selected concatenated time window.
%   2. Smoothed population firing rate for each detected area.
%
% Output folder:
%   root/Figures/Rasters_PDS_inspector/

clear;

%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end

addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Inspector parameters
% Required identifiers for generate_filename('raster', meta)
animal_name = 'Slayer';
injection = 'Muscimol';   % 'Muscimol', 'Saline', or 'No injection'
session_idx = 5;

prepost_str = 'Pre';      % 'Pre' or 'Post'
state = 'RestClose';       % 'RestOpen' or 'RestClose'
area_type = 'Full';       % 'Full' or 'Cortex'
align = 'None';        % 'None', 'Last', or 'Longest'

% Time window in the concatenated raster, in ms.
% Example: 1:10000 means the first 10 s after concatenating all trials.
plot_start_ms = 260001;
plot_end_ms = 290000;

% Sorting and smoothing
sort_criterion = 'channel';
smooth_sigma_ms = 50;
smooth_half_width_ms = 200;

% Output behavior
REPLOT = true;
figure_visible = 'off';
output_resolution = 300;

%% Build metadata and load raster
meta = struct();
meta.animal_name = animal_name;
meta.injection = injection;
meta.prepost = prepost_str;
meta.state = state;
meta.area = area_type;
meta.align = align;
meta.session_idx = session_idx;
meta.file_name = generate_filename('raster', meta);

data_folder = fullfile(root, 'Data', 'Working', 'raster');
data_path = fullfile(data_folder, meta.file_name);

if ~isfile(data_path)
    error('Raster data file not found: %s', data_path);
end

d = load(data_path);
rasters = d.data.rasters;
cell_area = normalize_cell_area(d.data.cell_area);
N = d.meta.N;

fprintf('Loaded raster: %s\n', data_path);
fprintf('Animal: %s | Injection: %s | Session: %d | %s %s | Area: %s | Align: %s\n', ...
    animal_name, injection, session_idx, prepost_str, state, area_type, align);
fprintf('N neurons: %d | Trials: %d\n', N, numel(rasters));

%% Load sort index if available
sortidx_folder = fullfile(root, 'Data', 'Working', 'sortidx');
meta.criterion = sort_criterion;
meta.file_name = generate_filename('sortidx', meta);
sortidx_path = fullfile(sortidx_folder, meta.file_name);

if isfile(sortidx_path)
    sort_idx = load(sortidx_path).data.sort_idx;
    fprintf('Loaded sort index: %s\n', sortidx_path);
else
    warning('Sort index file not found: %s. Using original neuron order.', sortidx_path);
    sort_idx = 1:numel(cell_area);
end

cell_area = cell_area(sort_idx);
for r_idx = 1:numel(rasters)
    rasters{r_idx} = rasters{r_idx}(sort_idx, :);
end

%% Concatenate rasters and select target window
concatenated_raster = cell2mat(rasters);
trial_borders = cumsum(cellfun(@(x) size(x, 2), rasters));

total_len = size(concatenated_raster, 2);
plot_start_ms = max(1, round(plot_start_ms));
plot_end_ms = min(total_len, round(plot_end_ms));

if plot_start_ms > plot_end_ms
    error('Invalid time window after clipping: start=%d, end=%d, total_len=%d.', ...
        plot_start_ms, plot_end_ms, total_len);
end

plot_range = plot_start_ms:plot_end_ms;
segment_raster = concatenated_raster(:, plot_range);

segment_trial_borders = trial_borders( ...
    trial_borders >= plot_start_ms & trial_borders <= plot_end_ms);

fprintf('Plot window: %d-%d ms of concatenated raster. Window length = %d ms.\n', ...
    plot_start_ms, plot_end_ms, numel(plot_range));

%% Detect areas and assign colors
area_names = get_area_order(cell_area);
fprintf('Detected areas:\n');
for area_i = 1:numel(area_names)
    area_name = area_names{area_i};
    fprintf('  %s: %d neurons\n', area_name, sum(strcmp(cell_area, area_name)));
end

neuron_colors = zeros(numel(cell_area), 3);
for neuron_i = 1:numel(cell_area)
    neuron_colors(neuron_i, :) = get_area_color(cell_area{neuron_i});
end

%% Compute smoothed area firing rates
smooth_kernel = make_gaussian_kernel(smooth_half_width_ms, smooth_sigma_ms);
area_fr = struct([]);

for area_i = 1:numel(area_names)
    area_name = area_names{area_i};
    area_filter = strcmp(cell_area, area_name);

    % Firing rate per neuron for this area:
    % mean spike value across neurons at each ms * 1000 ms/s.
    raw_fr = mean(segment_raster(area_filter, :), 1) * 1000;
    smooth_fr = conv(raw_fr, smooth_kernel, 'same');

    area_fr(area_i).name = area_name;
    area_fr(area_i).N = sum(area_filter);
    area_fr(area_i).raw_fr = raw_fr;
    area_fr(area_i).smooth_fr = smooth_fr;
    area_fr(area_i).color = get_area_color(area_name);
end

%% Plot raster + aligned firing rate
fig_width = 14;
fig_height = 9;

f = figure('Color', 'w', 'Visible', figure_visible, ...
    'Units', 'inches', 'Position', [1, 1, fig_width, fig_height]);

ax_raster = axes(f, 'Position', [0.08, 0.36, 0.80, 0.55]);
plot_raster_panel(ax_raster, segment_raster, cell_area, neuron_colors, ...
    plot_start_ms, plot_end_ms, segment_trial_borders, area_names);
title(ax_raster, sprintf('%s %s session %d | %s %s | %s | %s | %d-%d ms', ...
    animal_name, injection, session_idx, prepost_str, state, area_type, align, ...
    plot_start_ms, plot_end_ms), ...
    'Interpreter', 'none');

ax_fr = axes(f, 'Position', [0.08, 0.11, 0.80, 0.18]);
plot_area_fr_panel(ax_fr, area_fr, plot_range, plot_start_ms, plot_end_ms, segment_trial_borders);
xlabel(ax_fr, 'Concatenated time (ms)');
ylabel(ax_fr, 'Firing rate (Hz/neuron)');
title(ax_fr, sprintf('Smoothed area firing rate, Gaussian \\sigma = %d ms', smooth_sigma_ms), ...
    'Interpreter', 'tex');

linkaxes([ax_raster, ax_fr], 'x');
xlim(ax_raster, [plot_start_ms, plot_end_ms]);
xlim(ax_fr, [plot_start_ms, plot_end_ms]);

%% Export
figure_folder = fullfile(root, 'Figures', 'Rasters_PDS_inspector');
check_path(figure_folder);

figure_name = sprintf('inspector_raster_fr_%s_%s_s%d_%s_%s_%s_%s_%06d-%06dms.png', ...
    sanitize_filename(animal_name), sanitize_filename(injection), session_idx, ...
    sanitize_filename(prepost_str), sanitize_filename(state), ...
    sanitize_filename(area_type), sanitize_filename(align), ...
    plot_start_ms, plot_end_ms);

figure_path = fullfile(figure_folder, figure_name);

if isfile(figure_path) && ~REPLOT
    fprintf('Figure already exists: %s\n', figure_path);
else
    exportgraphics(f, figure_path, ...
        'ContentType', 'image', ...
        'BackgroundColor', 'white', ...
        'Resolution', output_resolution);
    fprintf('Saved figure: %s\n', figure_path);
end

close(f);

%% Helper functions

function cell_area = normalize_cell_area(cell_area)
    if isstring(cell_area)
        cell_area = cellstr(cell_area);
    elseif iscategorical(cell_area)
        cell_area = cellstr(cell_area);
    elseif ischar(cell_area)
        cell_area = {cell_area};
    end

    if ~iscell(cell_area)
        error('Unsupported cell_area format.');
    end

    cell_area = cellfun(@char, cell_area, 'UniformOutput', false);
end

function area_names = get_area_order(cell_area)
    canonical_order = {'Thalamus', 'ACC', 'VLPFC'};
    present = unique(cell_area, 'stable');
    area_names = {};

    for i = 1:numel(canonical_order)
        if any(strcmp(present, canonical_order{i}))
            area_names{end+1} = canonical_order{i}; %#ok<AGROW>
        end
    end

    for i = 1:numel(present)
        if ~any(strcmp(area_names, present{i}))
            area_names{end+1} = present{i}; %#ok<AGROW>
        end
    end
end

function color = get_area_color(area_name)
    switch char(string(area_name))
        case 'Thalamus'
            color = [1.0, 0.0, 1.0]; % magenta
        case 'ACC'
            color = [0.0, 0.0, 1.0]; % blue
        case 'VLPFC'
            color = [1.0, 0.0, 0.0]; % red
        otherwise
            color = [0.0, 0.0, 0.0]; % black
    end
end

function kernel = make_gaussian_kernel(half_width_ms, sigma_ms)
    x = -half_width_ms:half_width_ms;
    kernel = exp(-(x.^2) / (2 * sigma_ms^2));
    kernel = kernel / sum(kernel);
end

function plot_raster_panel(ax, segment_raster, cell_area, neuron_colors, ...
    plot_start_ms, plot_end_ms, segment_trial_borders, area_names) %#ok<INUSD>

    axes(ax); %#ok<LAXES>
    cla(ax);
    hold(ax, 'on');

    for area_i = 1:numel(area_names)
        area_name = area_names{area_i};
        neuron_idx = find(strcmp(cell_area, area_name));
        if isempty(neuron_idx)
            continue;
        end

        sub_raster = segment_raster(neuron_idx, :);
        [row_sub, col_sub] = find(sub_raster);

        if isempty(row_sub)
            continue;
        end

        x = plot_start_ms + col_sub - 1;
        y = neuron_idx(row_sub);
        this_color = get_area_color(area_name);

        scatter(ax, x, y, 3, ...
            'Marker', '.', ...
            'MarkerEdgeColor', this_color, ...
            'DisplayName', sprintf('%s, N=%d', area_name, numel(neuron_idx)));
    end

    draw_trial_borders(ax, segment_trial_borders);

    hold(ax, 'off');
    xlim(ax, [plot_start_ms, plot_end_ms]);
    ylim(ax, [0.5, size(segment_raster, 1) + 0.5]);
    % set(ax, 'YDir', 'reverse');
    ylabel(ax, 'Neuron #');
    xlabel(ax, '');
    set(ax, 'XTickLabel', []);
    box(ax, 'on');
    % legend(ax, 'Location', 'northoutside');

    hold(ax, 'on');
    x_text = plot_start_ms + 0.005 * (plot_end_ms - plot_start_ms);
    for area_i = 1:numel(area_names)
        area_name = area_names{area_i};
        neuron_idx = find(strcmp(cell_area, area_name));
        if isempty(neuron_idx)
            continue;
        end
        y_text = median(neuron_idx);
        text(ax, x_text, y_text, area_name, ...
            'Color', get_area_color(area_name), ...
            'FontWeight', 'bold', ...
            'Interpreter', 'none', ...
            'VerticalAlignment', 'middle');
    end
    hold(ax, 'off');
end

function plot_area_fr_panel(ax, area_fr, plot_range, plot_start_ms, plot_end_ms, segment_trial_borders)
    axes(ax); %#ok<LAXES>
    cla(ax);
    hold(ax, 'on');

    for area_i = 1:numel(area_fr)
        plot(ax, plot_range, area_fr(area_i).smooth_fr, ...
            'Color', area_fr(area_i).color, ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s, N=%d', area_fr(area_i).name, area_fr(area_i).N));
    end

    draw_trial_borders(ax, segment_trial_borders);

    hold(ax, 'off');
    xlim(ax, [plot_start_ms, plot_end_ms]);
    box(ax, 'on');
    legend(ax, 'Location', 'southoutside');
end

function draw_trial_borders(ax, segment_trial_borders)
    yl = ylim(ax);
    for border_i = 1:numel(segment_trial_borders)
        xline(ax, segment_trial_borders(border_i), ...
            'Color', [0.5, 0.5, 0.5], ...
            'LineStyle', '--', ...
            'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end
    ylim(ax, yl);
end

function safe_name = sanitize_filename(name)
    safe_name = regexprep(char(string(name)), '[^A-Za-z0-9_-]', '_');
end
