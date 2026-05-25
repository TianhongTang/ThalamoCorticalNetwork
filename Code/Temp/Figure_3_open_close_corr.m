%% Figure 3: Example session network structure and Open-Close correlation
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

% Example session: Slayer Mus 6.
meta = struct();
meta.animal_name = 'Slayer';
meta.injection = 'Muscimol';
meta.align = 'Last';
meta.session_idx = 6;
meta.resting_dur_threshold = 15;
meta.area = 'Cortex';
meta.shuffle_idx = 0;
meta.kernel_name = 'DeltaPure';
meta.reg_name = 'L2=0_2';
meta.epoch = 3000;
meta.fold_idx = 0;

preposts = {'Pre', 'Pre', 'Post', 'Post'};
states   = {'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};
n_state = numel(states);

% figure
f = figure('Color', 'w');
tiles = tiledlayout(2, 4, "TileSpacing", "Compact", "Padding", "Compact");

%% parameters
err_multi = 2; % threshold for significant J, in multiples of the error estimate from GLM.
density_nbin = 60;
scatter_marker_size = 8;
scatter_alpha = 0.25;

%% Load all four states first
state_data = struct();
for i = 1:n_state
    sd = load_state_connectivity(root, meta, preposts{i}, states{i}); %#ok<SAGROW>
    if i == 1
        state_data = sd;
    else
        state_data(i) = sd; % Store loaded state data in the array
    end
end

%% Row 1: network plot for each state
for i = 1:n_state
    tile = nexttile(i);
    call_plot_network(tile, ...
        state_data(i).J12, state_data(i).J21, ...
        state_data(i).err12, state_data(i).err21, err_multi, [], []);
    title(tile, sprintf('Network plot: %s, %s', state_data(i).prepost, state_data(i).state));
end

%% Row 2: Open-Close correlation within Pre and Post
% Pre = compare RestOpen vs RestClose within Pre
pre_open  = state_data(find(strcmp({state_data.prepost}, 'Pre')  & strcmp({state_data.state}, 'RestOpen'), 1));
pre_close = state_data(find(strcmp({state_data.prepost}, 'Pre')  & strcmp({state_data.state}, 'RestClose'), 1));
post_open  = state_data(find(strcmp({state_data.prepost}, 'Post') & strcmp({state_data.state}, 'RestOpen'), 1));
post_close = state_data(find(strcmp({state_data.prepost}, 'Post') & strcmp({state_data.state}, 'RestClose'), 1));

[x_pre, y_pre, label_pre] = make_open_close_vectors(pre_open, pre_close);
[x_post, y_post, label_post] = make_open_close_vectors(post_open, post_close);

[rho_pre, p_pre] = pearson_stats(x_pre, y_pre);
[rho_post, p_post] = pearson_stats(x_post, y_post);

% Row 2 col 1: Pre scatter
ax = nexttile(5);
plot_open_close_scatter(ax, x_pre, y_pre, rho_pre, p_pre, scatter_marker_size, scatter_alpha, ...
    sprintf('Pre scatter: %s', label_pre));

% Row 2 col 2: Pre density
ax = nexttile(6);
plot_open_close_density(ax, x_pre, y_pre, rho_pre, p_pre, density_nbin, ...
    sprintf('Pre density: %s', label_pre));

% Row 2 col 3: Post scatter
ax = nexttile(7);
plot_open_close_scatter(ax, x_post, y_post, rho_post, p_post, scatter_marker_size, scatter_alpha, ...
    sprintf('Post scatter: %s', label_post));

% Row 2 col 4: Post density
ax = nexttile(8);
plot_open_close_density(ax, x_post, y_post, rho_post, p_post, density_nbin, ...
    sprintf('Post density: %s', label_post));

%% Export to pdf
fig = gcf;

save_folder = fullfile(root, 'Figures', 'Paper');
check_path(save_folder);
filename = fullfile(save_folder, 'Figure3.pdf');

figWidth  = 16.0;   % inches
figHeight = 10.0;   % inches
resolution = 300;  % dpi; mainly affects rasterized components

set(fig, 'Units', 'inches');
fig.Position(3:4) = [figWidth, figHeight];

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [figWidth, figHeight]);
set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);
set(fig, 'Color', 'w');

exportgraphics(fig, filename, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);

filename = fullfile(save_folder, 'Figure3_preview.jpg');
exportgraphics(fig, filename, ...
    'ContentType', 'image', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);

close(fig);

%% functions
function state_struct = load_state_connectivity(root, meta, prepost, state)
    state_struct = struct();
    state_struct.prepost = prepost;
    state_struct.state = state;

    meta.prepost = prepost;
    meta.state = state;

    % Load raster data for neuron area info
    meta.filename = generate_filename('raster', meta);
    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    cell_area = raster_data.data.cell_area;
    filter1 = ismember(cell_area, {'ACC'});
    filter2 = ismember(cell_area, {'VLPFC'});

    % Load GLM data for connectivity info
    meta.filename = generate_filename('GLM', meta);
    GLM_data = load(fullfile(root, 'Data', 'Working', 'GLM', meta.filename));
    fprintf('Loaded GLM data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    N = GLM_data.meta.N;
    J = GLM_data.data.model_par(:, (2:N+1));
    err = GLM_data.data.model_err.total(:, (2:N+1));

    state_struct.cell_area = cell_area;
    state_struct.filter1 = filter1;
    state_struct.filter2 = filter2;
    state_struct.J = J;
    state_struct.err = err;
    state_struct.J12 = J(filter1, filter2);
    state_struct.J21 = J(filter2, filter1);
    state_struct.err12 = err(filter1, filter2);
    state_struct.err21 = err(filter2, filter1);
end

function [x, y, label_text] = make_open_close_vectors(open_state, close_state)
    validate_matching_filters(open_state, close_state);

    % Compare corresponding cross-area directed Jij values.
    x12 = open_state.J12(:);
    y12 = close_state.J12(:);
    x21 = open_state.J21(:);
    y21 = close_state.J21(:);

    x = [x12; x21];
    y = [y12; y21];

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    label_text = 'ACC↔VLPFC Jij';
end

function validate_matching_filters(state_a, state_b)
    if ~isequal(state_a.filter1(:), state_b.filter1(:)) || ~isequal(state_a.filter2(:), state_b.filter2(:))
        error('ACC/VLPFC neuron identities do not match between %s %s and %s %s.', ...
            state_a.prepost, state_a.state, state_b.prepost, state_b.state);
    end

    if ~isequal(size(state_a.J12), size(state_b.J12)) || ~isequal(size(state_a.J21), size(state_b.J21))
        error('Cross-area J matrix sizes do not match between %s %s and %s %s.', ...
            state_a.prepost, state_a.state, state_b.prepost, state_b.state);
    end
end

function [rho, pval] = pearson_stats(x, y)
    if numel(x) < 2 || numel(y) < 2
        rho = NaN;
        pval = NaN;
        return;
    end

    [R, P] = corrcoef(x, y, 'Rows', 'complete');
    rho = R(1, 2);
    pval = P(1, 2);
end

function plot_open_close_scatter(ax, x, y, rho, pval, marker_size, marker_alpha, title_text)
    scatter(ax, x, y, marker_size, 'filled', ...
        'MarkerFaceColor', [0.2, 0.2, 0.2], ...
        'MarkerFaceAlpha', marker_alpha, ...
        'MarkerEdgeAlpha', marker_alpha, ...
        'DisplayName', 'J_{ij}');
    hold(ax, 'on');
    plot_identity_line(ax, x, y);
    hold(ax, 'off');

    axis(ax, 'square');
    xlabel(ax, 'Open J_{ij}');
    ylabel(ax, 'Close J_{ij}');
    title(ax, title_text);
    add_stats_text(ax, rho, pval);
end

function plot_open_close_density(ax, x, y, rho, pval, nbin, title_text)
    [edges_x, edges_y] = make_density_edges(x, y, nbin);
    N = histcounts2(x, y, edges_x, edges_y);
    imagesc(ax, edges_x, edges_y, N.');
    set(ax, 'YDir', 'normal');
    hold(ax, 'on');
    plot_identity_line(ax, x, y);
    hold(ax, 'off');
    colorbar(ax);

    axis(ax, 'square');
    xlim(ax, [edges_x(1), edges_x(end)]);
    ylim(ax, [edges_y(1), edges_y(end)]);
    xlabel(ax, 'Open J_{ij}');
    ylabel(ax, 'Close J_{ij}');
    title(ax, title_text);
    add_stats_text(ax, rho, pval);
end

function [edges_x, edges_y] = make_density_edges(x, y, nbin)
    all_vals = [x(:); y(:)];
    all_vals = all_vals(isfinite(all_vals));
    if isempty(all_vals)
        all_vals = [-1; 1];
    end
    vmin = min(all_vals);
    vmax = max(all_vals);
    if vmin == vmax
        delta = max(abs(vmin) * 0.1, 1e-6);
        vmin = vmin - delta;
        vmax = vmax + delta;
    end
    edges_x = linspace(vmin, vmax, nbin + 1);
    edges_y = linspace(vmin, vmax, nbin + 1);
end

function plot_identity_line(ax, x, y)
    vals = [x(:); y(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        vals = [-1; 1];
    end
    vmin = min(vals);
    vmax = max(vals);
    if vmin == vmax
        delta = max(abs(vmin) * 0.1, 1e-6);
        vmin = vmin - delta;
        vmax = vmax + delta;
    end
    plot(ax, [vmin, vmax], [vmin, vmax], 'r--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlim(ax, [vmin, vmax]);
    ylim(ax, [vmin, vmax]);
end

function add_stats_text(ax, rho, pval)
    if isnan(rho) || isnan(pval)
        stats_text = 'rho = NaN\np = NaN';
    else
        if pval < 1e-3
            p_str = sprintf('%.2e', pval);
        else
            p_str = sprintf('%.3f', pval);
        end
        stats_text = sprintf('rho = %.3f\np = %s', rho, p_str);
    end
    text(ax, 0.04, 0.96, stats_text, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 4, 'FontSize', 9);
end

function call_plot_network(ax, J12, J21, err12, err21, err_multi, highlight_i, highlight_j)
    if ~isempty(highlight_i) && ~isempty(highlight_j)
        try
            plot_network(ax, J12, J21, err12, err21, err_multi, highlight_i, highlight_j);
            return;
        catch ME_highlight
            warning('plot_network with highlight failed: %s. Retrying without highlight.', ME_highlight.message);
        end
    end

    try
        plot_network(ax, J12, J21, err12, err21, err_multi);
    catch ME_no_highlight
        if ~isempty(J12)
            warning('plot_network without highlight failed: %s. Retrying with fallback highlight (1,1).', ME_no_highlight.message);
            plot_network(ax, J12, J21, err12, err21, err_multi, 1, 1);
        else
            rethrow(ME_no_highlight);
        end
    end
end
