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
tiles = tiledlayout(3, 4, "TileSpacing", "Compact", "Padding", "Compact");

%% parameters
err_multi = 1; % threshold for significant J, in multiples of the error estimate from GLM.
density_nbin = 60;
scatter_marker_size = 8;
scatter_alpha = 0.25;
density_clip_percentile = [0.5, 99.5]; % Robust density axis limits; reduces outlier domination.
density_use_log_count = true; % Plot log10(count+1) so sparse bins are visible.
category_labels = {'Negative', 'Non-sig', 'Positive'};

%% Load from metadata
mt = load_meta(root, 'table'); % metadata table

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

[data_pre, label_pre] = make_open_close_vectors(pre_open, pre_close);
[data_post, label_post] = make_open_close_vectors(post_open, post_close);

[rho_pre, p_pre, n_pre] = pearson_stats(abs(data_pre.x), abs(data_pre.y));
[rho_post, p_post, n_post] = pearson_stats(abs(data_post.x), abs(data_post.y));

[cat_counts_pre, agreement_pre, kappa_pre, cat_n_pre] = make_open_close_category_counts(pre_open, pre_close, err_multi);
[cat_counts_post, agreement_post, kappa_post, cat_n_post] = make_open_close_category_counts(post_open, post_close, err_multi);

fprintf('Pre Open-Close Pearson: rho = %.6f, p = %.3e, n = %d\n', rho_pre, p_pre, n_pre);
fprintf('Post Open-Close Pearson: rho = %.6f, p = %.3e, n = %d\n', rho_post, p_post, n_post);
fprintf('Pre categorical agreement: agreement = %.6f, kappa = %.6f, n = %d\n', agreement_pre, kappa_pre, cat_n_pre);
fprintf('Post categorical agreement: agreement = %.6f, kappa = %.6f, n = %d\n', agreement_post, kappa_post, cat_n_post);

% Row 2 col 1: Pre scatter
ax = nexttile(5);
plot_open_close_scatter(ax, data_pre, rho_pre, p_pre, n_pre, scatter_marker_size, scatter_alpha, ...
    sprintf('Pre scatter: %s', label_pre));

% Row 2 col 2: Pre density
ax = nexttile(6);
plot_open_close_density(ax, data_pre.x, data_pre.y, rho_pre, p_pre, n_pre, density_nbin, density_clip_percentile, density_use_log_count, ...
    sprintf('Pre density: %s', label_pre));

% Row 2 col 3: Post scatter
ax = nexttile(7);
plot_open_close_scatter(ax, data_post, rho_post, p_post, n_post, scatter_marker_size, scatter_alpha, ...
    sprintf('Post scatter: %s', label_post));

% Row 2 col 4: Post density
ax = nexttile(8);
plot_open_close_density(ax, data_post.x, data_post.y, rho_post, p_post, n_post, density_nbin, density_clip_percentile, density_use_log_count, ...
    sprintf('Post density: %s', label_post));

%% Row 3: categorical Open-Close transition tables
% Rows = Open category; columns = Close category.
ax = nexttile(9, [1, 2]);
plot_category_transition_table(ax, cat_counts_pre, category_labels, agreement_pre, kappa_pre, cat_n_pre, ...
    'Pre categorical transition');

ax = nexttile(11, [1, 2]);
plot_category_transition_table(ax, cat_counts_post, category_labels, agreement_post, kappa_post, cat_n_post, ...
    'Post categorical transition');

%% Export to pdf
fig = gcf;

save_folder = fullfile(root, 'Figures', 'Paper');
check_path(save_folder);
filename = fullfile(save_folder, 'Figure3.pdf');

figWidth  = 16.0;   % inches
figHeight = 15.0;   % inches
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

function [data, label_text] = make_open_close_vectors(open_state, close_state)
    validate_matching_filters(open_state, close_state);

    % Compare corresponding cross-area directed Jij values.
    x12 = open_state.J12(:);
    y12 = close_state.J12(:);
    x21 = open_state.J21(:);
    y21 = close_state.J21(:);
    xerr12 = open_state.err12(:);
    yerr12 = close_state.err12(:);
    xerr21 = open_state.err21(:);
    yerr21 = close_state.err21(:);

    x = [x12; x21];
    y = [y12; y21];
    xerr = [xerr12; xerr21];
    yerr = [yerr12; yerr21];

    % Finite filter
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    xerr = xerr(valid);
    yerr = yerr(valid);

    data = struct();

    % Significance filter
    % significant = abs(x) > xerr & abs(y) > yerr;
    pos = x > xerr & y > yerr;
    neg = x < -xerr & y < -yerr;
    
    data.xpos = x(pos);
    data.ypos = y(pos);
    data.xneg = x(neg);
    data.yneg = y(neg);
    data.x = [data.xpos; data.xneg];
    data.y = [data.ypos; data.yneg];

    label_text = 'Significant ACC↔VLPFC Jij';
end


function [counts, agreement, kappa, n_valid] = make_open_close_category_counts(open_state, close_state, err_multi)
    validate_matching_filters(open_state, close_state);

    open_cat12 = classify_connections(open_state.J12(:), open_state.err12(:), err_multi);
    close_cat12 = classify_connections(close_state.J12(:), close_state.err12(:), err_multi);

    open_cat21 = classify_connections(open_state.J21(:), open_state.err21(:), err_multi);
    close_cat21 = classify_connections(close_state.J21(:), close_state.err21(:), err_multi);

    open_cat = [open_cat12; open_cat21];
    close_cat = [close_cat12; close_cat21];

    valid = isfinite(open_cat) & isfinite(close_cat);
    open_cat = open_cat(valid);
    close_cat = close_cat(valid);
    n_valid = numel(open_cat);

    class_values = [-1, 0, 1];
    counts = zeros(3, 3);
    for i_class = 1:3
        for j_class = 1:3
            counts(i_class, j_class) = sum(open_cat == class_values(i_class) & ...
                                           close_cat == class_values(j_class));
        end
    end

    agreement = compute_raw_agreement(counts);
    kappa = compute_cohen_kappa(counts);
end

function cat = classify_connections(J, err, err_multi)
    % Categories:
    %   -1 = significantly negative
    %    0 = non-significant
    %    1 = significantly positive
    J = J(:);
    err = err(:);
    cat = nan(size(J));

    valid = isfinite(J) & isfinite(err) & err >= 0;
    cat(valid) = 0;
    cat(valid & J >  err_multi * err) = 1;
    cat(valid & J < -err_multi * err) = -1;
end

function agreement = compute_raw_agreement(counts)
    total_n = sum(counts(:));
    if total_n == 0
        agreement = NaN;
    else
        agreement = trace(counts) / total_n;
    end
end

function kappa = compute_cohen_kappa(counts)
    total_n = sum(counts(:));
    if total_n == 0
        kappa = NaN;
        return;
    end

    p_observed = trace(counts) / total_n;
    row_marginal = sum(counts, 2) / total_n;
    col_marginal = sum(counts, 1) / total_n;
    p_expected = row_marginal.' * col_marginal.';

    if abs(1 - p_expected) < eps
        kappa = NaN;
    else
        kappa = (p_observed - p_expected) / (1 - p_expected);
    end
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

function [rho, pval, n_valid] = pearson_stats(x, y)
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    n_valid = numel(x);

    if n_valid < 3
        rho = NaN;
        pval = NaN;
        return;
    end

    [R, P] = corrcoef(x, y);
    rho = R(1, 2);
    pval = P(1, 2);
end

function plot_open_close_scatter(ax, data, rho, pval, n_valid, marker_size, marker_alpha, title_text)
    scatter(ax, data.xpos, data.ypos, marker_size, 'filled', ...
        'MarkerFaceColor', [1, 0.2, 0.2], ...
        'MarkerFaceAlpha', marker_alpha, ...
        'MarkerEdgeAlpha', marker_alpha, ...
        'DisplayName', 'positive');
    hold(ax, 'on');
    scatter(ax, -data.xneg, -data.yneg, marker_size, 'filled', ...
        'MarkerFaceColor', [0.2, 0.2, 1], ...
        'MarkerFaceAlpha', marker_alpha, ...
        'MarkerEdgeAlpha', marker_alpha, ...
        'DisplayName', 'negative');
    plot_identity_line(ax, data.xpos, data.ypos);
    hold(ax, 'off');

    axis(ax, 'square');
    xlabel(ax, 'Open |J_{ij}|');
    ylabel(ax, 'Close |J_{ij}|');
    title(ax, title_text);
    legend(ax, 'Location', 'best');
    add_stats_text(ax, rho, pval, n_valid);
end

function plot_open_close_density(ax, x, y, rho, pval, n_valid, nbin, clip_percentile, use_log_count, title_text)
    [edges_x, edges_y, n_in_range] = make_density_edges(x, y, nbin, clip_percentile);
    count_mat = histcounts2(x, y, edges_x, edges_y);

    centers_x = 0.5 * (edges_x(1:end-1) + edges_x(2:end));
    centers_y = 0.5 * (edges_y(1:end-1) + edges_y(2:end));

    if use_log_count
        plot_mat = log10(count_mat.' + 1);
    else
        plot_mat = count_mat.';
    end

    imagesc(ax, centers_x, centers_y, plot_mat);
    set(ax, 'YDir', 'normal');
    hold(ax, 'on');
    plot_identity_line_with_limits(ax, edges_x(1), edges_x(end), edges_y(1), edges_y(end));
    hold(ax, 'off');
    cb = colorbar(ax);
    if use_log_count
        ylabel(cb, 'log_{10}(count + 1)');
    else
        ylabel(cb, 'count');
    end

    axis(ax, 'square');
    xlim(ax, [edges_x(1), edges_x(end)]);
    ylim(ax, [edges_y(1), edges_y(end)]);
    xlabel(ax, 'Open J_{ij}');
    ylabel(ax, 'Close J_{ij}');
    title(ax, title_text);
    add_stats_text(ax, rho, pval, n_valid);

    fprintf('%s density: total n = %d, in density range = %d, max bin count = %d\n', ...
        title_text, n_valid, n_in_range, max(count_mat(:)));
end
function [edges_x, edges_y, n_in_range] = make_density_edges(x, y, nbin, clip_percentile)
    valid = isfinite(x) & isfinite(y);
    x_valid = x(valid);
    y_valid = y(valid);
    all_vals = [x_valid(:); y_valid(:)];

    if isempty(all_vals)
        all_vals = [-1; 1];
        x_valid = [];
        y_valid = [];
    end

    if nargin < 4 || isempty(clip_percentile)
        clip_percentile = [0, 100];
    end

    vmin = local_percentile(all_vals, clip_percentile(1));
    vmax = local_percentile(all_vals, clip_percentile(2));

    if ~isfinite(vmin) || ~isfinite(vmax) || vmin == vmax
        vmin = min(all_vals);
        vmax = max(all_vals);
    end

    if vmin == vmax
        delta = max(abs(vmin) * 0.1, 1e-6);
        vmin = vmin - delta;
        vmax = vmax + delta;
    end

    edges_x = linspace(vmin, vmax, nbin + 1);
    edges_y = linspace(vmin, vmax, nbin + 1);

    in_range = x_valid >= edges_x(1) & x_valid <= edges_x(end) & ...
               y_valid >= edges_y(1) & y_valid <= edges_y(end);
    n_in_range = sum(in_range);
end

function q = local_percentile(vals, pct)
    vals = sort(vals(isfinite(vals)));
    if isempty(vals)
        q = NaN;
        return;
    end

    pct = max(0, min(100, pct));
    if numel(vals) == 1
        q = vals(1);
        return;
    end

    pos = 1 + (numel(vals) - 1) * pct / 100;
    lo = floor(pos);
    hi = ceil(pos);
    if lo == hi
        q = vals(lo);
    else
        q = vals(lo) + (pos - lo) * (vals(hi) - vals(lo));
    end
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

function plot_identity_line_with_limits(ax, xmin, xmax, ymin, ymax)
    vmin = max(xmin, ymin);
    vmax = min(xmax, ymax);
    if isfinite(vmin) && isfinite(vmax) && vmin < vmax
        plot(ax, [vmin, vmax], [vmin, vmax], 'r--', 'LineWidth', 1, 'HandleVisibility', 'off');
    else
        vals_min = min([xmin, ymin]);
        vals_max = max([xmax, ymax]);
        plot(ax, [vals_min, vals_max], [vals_min, vals_max], 'r--', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

function add_stats_text(ax, rho, pval, n_valid)
    if nargin < 4
        n_valid = NaN;
    end

    if isnan(rho) || isnan(pval)
        stats_text = sprintf('rho = NaN\np = NaN\nn = %g', n_valid);
    else
        if pval < 1e-3
            p_str = sprintf('%.2e', pval);
        else
            p_str = sprintf('%.3f', pval);
        end
        stats_text = sprintf('rho = %.3f\np = %s\nn = %d', rho, p_str, n_valid);
    end
    text(ax, 0.04, 0.96, stats_text, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 4, 'FontSize', 9);
end



function plot_category_transition_table(ax, counts, category_labels, agreement, kappa, n_valid, title_text)
    imagesc(ax, counts);
    axis(ax, 'square');

    xticks(ax, 1:3);
    yticks(ax, 1:3);
    xticklabels(ax, category_labels);
    yticklabels(ax, category_labels);
    xtickangle(ax, 30);

    xlabel(ax, 'Close category');
    ylabel(ax, 'Open category');

    cb = colorbar(ax);
    ylabel(cb, 'count');

    max_count = max(counts(:));
    for row_idx = 1:3
        for col_idx = 1:3
            this_count = counts(row_idx, col_idx);
            if max_count > 0 && this_count > 0.55 * max_count
                text_color = 'w';
            else
                text_color = 'k';
            end
            text(ax, col_idx, row_idx, sprintf('%d', this_count), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', ...
                'Color', text_color);
        end
    end

    if isnan(agreement)
        agreement_str = 'NaN';
    else
        agreement_str = sprintf('%.3f', agreement);
    end

    if isnan(kappa)
        kappa_str = 'NaN';
    else
        kappa_str = sprintf('%.3f', kappa);
    end

    title(ax, sprintf('%s\\nAgreement = %s, kappa = %s, n = %d', ...
        title_text, agreement_str, kappa_str, n_valid));
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
