%% Figure 3 supplement: by-session Open-Close scatter reference
% Each row is one selected session.
% Columns:
%   1. Pre abs scatter, same-sign significant connections only
%   2. Post abs scatter, same-sign significant connections only
%   3. Pre original signed J scatter, all connections, colored by Open/Close significance
%   4. Post original signed J scatter, all connections, colored by Open/Close significance

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

%% State definitions
preposts = {'Pre', 'Pre', 'Post', 'Post'};
states   = {'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};
n_state = numel(states);

%% Parameters
kernel_idx = 1;
err_multi = 1; % threshold for significant J, in multiples of the GLM error estimate.
scatter_marker_size = 8;
scatter_alpha = 0.25;
abs_axis_limit = [0, 3.5]; % set [] for dynamic axis limits.
signed_axis_limit = [];    % set [] for dynamic symmetric axis limits.
skip_failed_sessions = false;
max_sessions_to_include = inf; % set smaller while debugging.

%% Load and filter metadata
mt = load_meta(root, 'table');
mt = mt.GLM;

selected_rows = default_metadata_filter(mt);
selected_mt = mt(selected_rows, :);
meta_array = table2struct(selected_mt);

if isfinite(max_sessions_to_include)
    meta_array = meta_array(1:min(numel(meta_array), max_sessions_to_include));
end

if isempty(meta_array)
    error('No metadata rows selected.');
end

%% Load each selected session and store scatter data
session_results = struct([]);
valid_session_count = 0;
failed_session_count = 0;

for session_i = 1:numel(meta_array)
    meta = meta_array(session_i);
    session_label = make_session_label(meta);
    fprintf('\n===== Loading %s (%d/%d) =====\n', session_label, session_i, numel(meta_array));

    try
        state_data = struct([]);
        for state_i = 1:n_state
            sd = load_state_connectivity(root, meta, preposts{state_i}, states{state_i}, kernel_idx);
            if state_i == 1
                state_data = sd;
            else
                state_data = [state_data; sd]; %#ok<AGROW>
            end
        end

        pre_open  = state_data(find(strcmp({state_data.prepost}, 'Pre')  & strcmp({state_data.state}, 'RestOpen'), 1));
        pre_close = state_data(find(strcmp({state_data.prepost}, 'Pre')  & strcmp({state_data.state}, 'RestClose'), 1));
        post_open  = state_data(find(strcmp({state_data.prepost}, 'Post') & strcmp({state_data.state}, 'RestOpen'), 1));
        post_close = state_data(find(strcmp({state_data.prepost}, 'Post') & strcmp({state_data.state}, 'RestClose'), 1));

        [data_pre, ~] = make_open_close_vectors(pre_open, pre_close, err_multi);
        [data_post, ~] = make_open_close_vectors(post_open, post_close, err_multi);

        valid_session_count = valid_session_count + 1;
        session_results(valid_session_count).label = session_label; %#ok<SAGROW>
        session_results(valid_session_count).data_pre = data_pre;
        session_results(valid_session_count).data_post = data_post;

    catch ME
        failed_session_count = failed_session_count + 1;
        if skip_failed_sessions
            warning('Skipping %s because loading/processing failed: %s', session_label, ME.message);
            continue;
        else
            rethrow(ME);
        end
    end
end

if valid_session_count == 0
    error('No valid sessions were loaded.');
end

fprintf('\nValid sessions: %d. Failed sessions: %d.\n', valid_session_count, failed_session_count);

%% Figure
f = figure('Color', 'w', 'Visible', 'off');
tiles = tiledlayout(valid_session_count, 4, "TileSpacing", "Compact", "Padding", "Compact");

for session_i = 1:valid_session_count
    session_label = session_results(session_i).label;
    data_pre = session_results(session_i).data_pre;
    data_post = session_results(session_i).data_post;

    [rho_pre_abs, p_pre_abs, n_pre_abs] = pearson_stats(data_pre.x_abs, data_pre.y_abs);
    [rho_post_abs, p_post_abs, n_post_abs] = pearson_stats(data_post.x_abs, data_post.y_abs);
    [rho_pre_orig, p_pre_orig, n_pre_orig] = pearson_stats(data_pre.x_orig, data_pre.y_orig);
    [rho_post_orig, p_post_orig, n_post_orig] = pearson_stats(data_post.x_orig, data_post.y_orig);

    ax = nexttile((session_i - 1) * 4 + 1);
    plot_abs_scatter(ax, data_pre, rho_pre_abs, p_pre_abs, n_pre_abs, ...
        scatter_marker_size, scatter_alpha, abs_axis_limit, ...
        sprintf('%s\nPre abs same-sign sig', session_label));

    ax = nexttile((session_i - 1) * 4 + 2);
    plot_abs_scatter(ax, data_post, rho_post_abs, p_post_abs, n_post_abs, ...
        scatter_marker_size, scatter_alpha, abs_axis_limit, ...
        sprintf('%s\nPost abs same-sign sig', session_label));

    ax = nexttile((session_i - 1) * 4 + 3);
    plot_signed_scatter(ax, data_pre, rho_pre_orig, p_pre_orig, n_pre_orig, ...
        scatter_marker_size, scatter_alpha, signed_axis_limit, ...
        sprintf('%s\nPre original signed J', session_label));

    ax = nexttile((session_i - 1) * 4 + 4);
    plot_signed_scatter(ax, data_post, rho_post_orig, p_post_orig, n_post_orig, ...
        scatter_marker_size, scatter_alpha, signed_axis_limit, ...
        sprintf('%s\nPost original signed J', session_label));
end

%% Export to pdf
fig = gcf;

save_folder = fullfile(root, 'Figures', 'Paper');
check_path(save_folder);

figWidth  = 16;   % inches
figHeight = max(4.0, 3.2 * valid_session_count);
resolution = 300;  % dpi; mainly affects rasterized components

set(fig, 'Units', 'inches');
fig.Position(3:4) = [figWidth, figHeight];

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [figWidth, figHeight]);
set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);
set(fig, 'Color', 'w');

filename = fullfile(save_folder, sprintf('Figure3_by_session_scatter_k%d_preview.jpg', kernel_idx));
exportgraphics(fig, filename, ...
    'ContentType', 'image', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);

filename = fullfile(save_folder, sprintf('Figure3_by_session_scatter_k%d.pdf', kernel_idx));
exportgraphics(fig, filename, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);


close(fig);

%% functions
function selected_rows = default_metadata_filter(mt)
    base_filter = strcmp(mt.kernel_name, "DeltaPure") & ...
                  strcmp(mt.align, 'Last') & ...
                  strcmp(mt.area, "Cortex") & ...
                  strcmp(mt.injection, 'Muscimol') & ...
                  (mt.epoch == 3000) & ...
                  (mt.fold_idx == 0) & ...
                  (mt.shuffle_idx == 0) & ...
                  cellfun(@(x) ~isempty(x) && x == 15, mt.resting_dur_threshold);

    anchor_filter = base_filter & ...
                    strcmp(mt.prepost, 'Pre') & ...
                    strcmp(mt.state, 'RestOpen');

    selected_rows = false(height(mt), 1);
    anchor_idx = find(anchor_filter);

    required_preposts = {'Pre', 'Pre', 'Post', 'Post'};
    required_states   = {'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};

    match_fields = {'animal_name', 'injection', 'align', 'session_idx', ...
                    'resting_dur_threshold', 'area', 'kernel_name', ...
                    'reg_name', 'epoch', 'fold_idx', 'shuffle_idx'};

    for k = 1:numel(anchor_idx)
        idx = anchor_idx(k);
        same_session = base_filter;

        for f = 1:numel(match_fields)
            field = match_fields{f};
            if ~ismember(field, mt.Properties.VariableNames)
                continue;
            end

            anchor_value = mt.(field)(idx);

            if iscell(mt.(field))
                anchor_value = anchor_value{1};
                same_session = same_session & cellfun(@(x) isequal(x, anchor_value), mt.(field));
            else
                same_session = same_session & arrayfun(@(x) isequal(x, anchor_value), mt.(field));
            end
        end

        has_all_required = true;
        for q = 1:numel(required_preposts)
            has_this = any(same_session & ...
                           strcmp(mt.prepost, required_preposts{q}) & ...
                           strcmp(mt.state, required_states{q}));

            if ~has_this
                has_all_required = false;
                if ismember('file_name', mt.Properties.VariableNames)
                    anchor_name = mt.file_name{idx};
                else
                    anchor_name = sprintf('row %d', idx);
                end
                warning('Session %s is missing required Pre/Post x Open/Close combination: %s %s. Skipping this session.', ...
                    anchor_name, required_preposts{q}, required_states{q});
                break;
            end
        end

        if has_all_required
            selected_rows(idx) = true;
        end
    end

    fprintf('Default metadata filter selected %d/%d rows.\n', sum(selected_rows), sum(base_filter));
    for k = 1:height(mt)
        if selected_rows(k)
            if ismember('file_name', mt.Properties.VariableNames)
                fprintf('  %s\n', mt.file_name{k});
            else
                fprintf('  row %d\n', k);
            end
        end
    end

    if ~any(selected_rows)
        warning('Default metadata filter selected no complete Pre/Post x Open/Close sets. Falling back to Pre RestOpen anchors.');
        selected_rows = anchor_filter;
    end
end

function label = make_session_label(meta)
    animal_str = char(string(meta.animal_name));
    injection_str = char(string(meta.injection));
    session_str = char(string(meta.session_idx));
    label = sprintf('%s %s session %s', animal_str, injection_str, session_str);
end

function state_struct = load_state_connectivity(root, meta, prepost, state, kernel_idx)
    state_struct = struct();
    state_struct.prepost = prepost;
    state_struct.state = state;

    meta.prepost = prepost;
    meta.state = state;

    meta.filename = generate_filename('raster', meta);
    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    cell_area = raster_data.data.cell_area;
    filter1 = ismember(cell_area, {'ACC'});
    filter2 = ismember(cell_area, {'VLPFC'});

    meta.filename = generate_filename('GLM', meta);
    GLM_data = load(fullfile(root, 'Data', 'Working', 'GLM', meta.filename));
    fprintf('Loaded GLM data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    N = GLM_data.meta.N;
    J = GLM_data.data.model_par(:, ((2+N*(kernel_idx-1)) : (1+N*kernel_idx)));
    err = GLM_data.data.model_err.total(:, ((2+N*(kernel_idx-1)) : (1+N*kernel_idx)));

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

function [data, label_text] = make_open_close_vectors(open_state, close_state, err_multi)
    validate_matching_filters(open_state, close_state);

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

    valid = isfinite(x) & isfinite(y) & isfinite(xerr) & isfinite(yerr);
    x = x(valid);
    y = y(valid);
    xerr = xerr(valid);
    yerr = yerr(valid);

    open_sig = abs(x) > err_multi * xerr;
    close_sig = abs(y) > err_multi * yerr;

    % For the abs scatter, keep same-sign significant connections only.
    pos = x > err_multi * xerr & y > err_multi * yerr;
    neg = x < -err_multi * xerr & y < -err_multi * yerr;

    data = struct();
    data.x_orig = x;
    data.y_orig = y;
    data.open_sig = open_sig;
    data.close_sig = close_sig;

    data.xpos = x(pos);
    data.ypos = y(pos);
    data.xneg = x(neg);
    data.yneg = y(neg);

    data.xpos_abs = data.xpos;
    data.ypos_abs = data.ypos;
    data.xneg_abs = -data.xneg;
    data.yneg_abs = -data.yneg;

    data.x_abs = [data.xpos_abs; data.xneg_abs];
    data.y_abs = [data.ypos_abs; data.yneg_abs];

    label_text = 'same-sign significant ACC↔VLPFC |J_{ij}|';
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

function plot_abs_scatter(ax, data, rho, pval, n_valid, marker_size, marker_alpha, axis_limit, title_text)
    cla(ax);
    scatter(ax, data.xpos_abs, data.ypos_abs, marker_size, 'filled', ...
        'MarkerFaceColor', [1, 0.2, 0.2], ...
        'MarkerFaceAlpha', marker_alpha, ...
        'MarkerEdgeAlpha', marker_alpha, ...
        'DisplayName', 'positive');
    hold(ax, 'on');
    scatter(ax, data.xneg_abs, data.yneg_abs, marker_size, 'filled', ...
        'MarkerFaceColor', [0.2, 0.2, 1], ...
        'MarkerFaceAlpha', marker_alpha, ...
        'MarkerEdgeAlpha', marker_alpha, ...
        'DisplayName', 'negative');

    if isempty(axis_limit)
        axis_limit = dynamic_positive_axis_limit(data.x_abs, data.y_abs);
    end

    plot(ax, axis_limit, axis_limit, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlim(ax, axis_limit);
    ylim(ax, axis_limit);
    hold(ax, 'off');

    axis(ax, 'square');
    xlabel(ax, 'Open |J_{ij}|');
    ylabel(ax, 'Close |J_{ij}|');
    title(ax, sprintf('%s\nr = %.3f, p = %.4f, n = %d', title_text, rho, pval, n_valid), ...
        'Interpreter', 'none');
    legend(ax, 'Location', 'northwest');
end

function plot_signed_scatter(ax, data, rho, pval, n_valid, marker_size, marker_alpha, axis_limit, title_text)
    cla(ax);

    x = data.x_orig;
    y = data.y_orig;
    open_sig = data.open_sig;
    close_sig = data.close_sig;

    none_sig = ~open_sig & ~close_sig;
    open_only = open_sig & ~close_sig;
    close_only = ~open_sig & close_sig;
    both_sig = open_sig & close_sig;

    hold(ax, 'on');
    scatter(ax, x(none_sig), y(none_sig), marker_size, 'filled', ...
        'MarkerFaceColor', [0.65, 0.65, 0.65], ...
        'MarkerFaceAlpha', marker_alpha, ...
        'MarkerEdgeAlpha', marker_alpha, ...
        'DisplayName', 'not sig');
    scatter(ax, x(open_only), y(open_only), marker_size, 'filled', ...
        'MarkerFaceColor', [0.1, 0.45, 0.95], ...
        'MarkerFaceAlpha', marker_alpha, ...
        'MarkerEdgeAlpha', marker_alpha, ...
        'DisplayName', 'open only');
    scatter(ax, x(close_only), y(close_only), marker_size, 'filled', ...
        'MarkerFaceColor', [0.95, 0.55, 0.1], ...
        'MarkerFaceAlpha', marker_alpha, ...
        'MarkerEdgeAlpha', marker_alpha, ...
        'DisplayName', 'close only');
    scatter(ax, x(both_sig), y(both_sig), marker_size, 'filled', ...
        'MarkerFaceColor', [0.45, 0.1, 0.75], ...
        'MarkerFaceAlpha', marker_alpha, ...
        'MarkerEdgeAlpha', marker_alpha, ...
        'DisplayName', 'both sig');

    if isempty(axis_limit)
        axis_limit = dynamic_symmetric_axis_limit(x, y);
    end

    plot(ax, axis_limit, axis_limit, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
    xline(ax, 0, 'k:', 'HandleVisibility', 'off');
    yline(ax, 0, 'k:', 'HandleVisibility', 'off');
    xlim(ax, axis_limit);
    ylim(ax, axis_limit);
    hold(ax, 'off');

    axis(ax, 'square');
    xlabel(ax, 'Open J_{ij}');
    ylabel(ax, 'Close J_{ij}');
    title(ax, sprintf('%s\nr = %.3f, p = %.4f, n = %d', title_text, rho, pval, n_valid), ...
        'Interpreter', 'none');
    % legend(ax, 'Location', 'southeast');
end

function axis_limit = dynamic_positive_axis_limit(x, y)
    vals = [x(:); y(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        axis_limit = [0, 1];
        return;
    end

    vmax = max(vals);
    if vmax <= 0
        vmax = 1;
    end
    axis_limit = [0, 1.05 * vmax];
end

function axis_limit = dynamic_symmetric_axis_limit(x, y)
    vals = [x(:); y(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        axis_limit = [-1, 1];
        return;
    end

    vmax = max(abs(vals));
    if vmax <= 0
        vmax = 1;
    end
    axis_limit = [-1.05 * vmax, 1.05 * vmax];
end
