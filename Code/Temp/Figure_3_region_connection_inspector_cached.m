%% Figure 3 helper: inspect selected scatter-region connections
% Load all session data once, then rerun the plot-parameter section to inspect different regions.

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

%% Load data parameters

data_preposts = {'Pre', 'Post'};
data_states = {'RestOpen', 'RestClose'};
data_kernel_indices = [1, 2];

skip_failed_sessions = false;
max_sessions_to_include = inf;

%% Load data

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

loaded_session_data = struct([]);
valid_session_count = 0;
failed_session_count = 0;

for session_i = 1:numel(meta_array)
    meta = meta_array(session_i);
    session_label = make_session_label(meta);
    fprintf('\n===== Loading %s (%d/%d) =====\n', session_label, session_i, numel(meta_array));

    try
        loaded_states = load_required_states_for_session(root, meta, data_preposts, data_states, data_kernel_indices);

        valid_session_count = valid_session_count + 1;
        loaded_session_data(valid_session_count).meta = meta; %#ok<SAGROW>
        loaded_session_data(valid_session_count).session_label = session_label; %#ok<SAGROW>
        loaded_session_data(valid_session_count).states = loaded_states; %#ok<SAGROW>
    catch ME
        failed_session_count = failed_session_count + 1;
        if skip_failed_sessions
            warning('Skipping %s because loading failed: %s', session_label, ME.message);
            continue;
        else
            rethrow(ME);
        end
    end
end

if valid_session_count == 0
    error('No valid sessions were loaded.');
end

fprintf('\nLoaded sessions: %d. Failed sessions: %d.\n', valid_session_count, failed_session_count);

%% Main plot
% Plot parameters

x_prepost = 'Pre';
x_state = 'RestOpen';
x_kernel_idx = 1;
x_axis_label = sprintf('%s %s K%d', x_prepost, x_state, x_kernel_idx);

y_prepost = 'Pre';
y_state = 'RestClose';
y_kernel_idx = 1;
y_axis_label = sprintf('%s %s K%d', y_prepost, y_state, y_kernel_idx);

% Options: 'all_connection', 'either_significant', 'both_significant', 'same_sign_abs', 'switch_sign_abs'
plot_type = 'same_sign_abs';

x_range = [3, 4];
y_range = [1, 2];

text_out = true;
save_text_out = false;

err_multi = 1;
scatter_marker_size = 8;
scatter_alpha = 0.25;
show_legend = true;
figure_visible = 'on';

colors = struct();
colors.non_sig = [0.65, 0.65, 0.65];
colors.x_only = [0.10, 0.75, 0.10];
colors.y_only = [0.95, 0.55, 0.10];
colors.both_sig = [0.45, 0.10, 0.75];

colors.pos = [1.00, 0.20, 0.20];
colors.neg = [0.10, 0.45, 1.00];
colors.switch_xy = [1.00, 0.50, 0.00];
colors.switch_yx = [0.00, 0.65, 0.80];
colors.switch = [0.45, 0.10, 0.75];

colors.identity_line = [1, 0, 0];
colors.zero_line = [0, 0, 0];
colors.region_box = [1, 0, 0];

% Build selected x/y pair from loaded data

key_x = state_key(x_prepost, x_state, x_kernel_idx);
key_y = state_key(y_prepost, y_state, y_kernel_idx);

pooled_data = empty_connection_data();

for session_i = 1:numel(loaded_session_data)
    meta = loaded_session_data(session_i).meta;
    session_label = loaded_session_data(session_i).session_label;
    loaded_states = loaded_session_data(session_i).states;

    if ~isfield(loaded_states, key_x)
        error('x condition was not loaded: %s. Add it to data_preposts/data_states/data_kernel_indices and rerun Load data.', key_x);
    end
    if ~isfield(loaded_states, key_y)
        error('y condition was not loaded: %s. Add it to data_preposts/data_states/data_kernel_indices and rerun Load data.', key_y);
    end

    state_x = loaded_states.(key_x);
    state_y = loaded_states.(key_y);

    pair_data = make_pair_vectors(state_x, state_y, err_multi, meta, session_label);
    pooled_data = append_connection_data(pooled_data, pair_data);
end

plot_data = prepare_plot_data(pooled_data, plot_type, x_axis_label, y_axis_label);
selection_mask = points_in_region(plot_data.x_plot, plot_data.y_plot, x_range, y_range);
selected_data = subset_plot_data(plot_data, selection_mask);

% Plot full scatter and zoomed selection

f = figure('Color', 'w', 'Visible', figure_visible);
tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

ax = nexttile(1);
plot_scatter_groups(ax, plot_data, scatter_marker_size, scatter_alpha, colors, show_legend);
draw_region_box(ax, x_range, y_range, colors.region_box);
axis_limit = get_axis_limit_for_plot(plot_data, plot_type);
xlim(ax, axis_limit.x);
ylim(ax, axis_limit.y);
axis(ax, 'square');
xlabel(ax, plot_data.x_label);
ylabel(ax, plot_data.y_label);
title(ax, sprintf('%s\nall plotted points, n=%d', plot_data.plot_title, numel(plot_data.x_plot)), 'Interpreter', 'none');

ax = nexttile(2);
plot_scatter_groups(ax, selected_data, scatter_marker_size + 8, min(1, scatter_alpha + 0.35), colors, show_legend);
xlim(ax, sort(x_range));
ylim(ax, sort(y_range));
axis(ax, 'square');
xlabel(ax, plot_data.x_label);
ylabel(ax, plot_data.y_label);
title(ax, sprintf('Zoomed selected region\nselected n=%d', numel(selected_data.x_plot)), 'Interpreter', 'none');

% Optional text output

if text_out
    print_selected_connections(selected_data);
end

if save_text_out
    save_folder = fullfile(root, 'Figures', 'Paper');
    check_path(save_folder);
    text_path = fullfile(save_folder, sprintf('Figure3_region_selection_%s.txt', sanitize_filename(plot_type)));
    write_selected_connections_to_file(selected_data, text_path);
    fprintf('Saved selected connection text to:\n%s\n', text_path);
end

%% Functions
function data = empty_connection_data()
    data = struct();

    data.x_orig = [];
    data.y_orig = [];
    data.x_err_orig = [];
    data.y_err_orig = [];
    data.x_cat_orig = [];
    data.y_cat_orig = [];

    data.session_label = strings(0, 1);
    data.animal_name = strings(0, 1);
    data.injection = strings(0, 1);
    data.session_idx = [];
    data.connection_direction = strings(0, 1);
    data.from_area = strings(0, 1);
    data.to_area = strings(0, 1);

    data.from_cell_area_idx = [];
    data.to_cell_area_idx = [];
    data.from_cell_global_idx = [];
    data.to_cell_global_idx = [];

    % Backward-compatible aliases: these are global cell indices.
    data.from_cell_idx = [];
    data.to_cell_idx = [];
end

function out = append_connection_data(out, incoming)
    fields = fieldnames(out);
    for k = 1:numel(fields)
        field_name = fields{k};
        if isfield(incoming, field_name)
            out.(field_name) = [out.(field_name); incoming.(field_name)];
        end
    end
end

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
                break;
            end
        end

        if has_all_required
            selected_rows(idx) = true;
        end
    end

    fprintf('Default metadata filter selected %d/%d rows.\n', sum(selected_rows), sum(base_filter));

    if ~any(selected_rows)
        warning('Default metadata filter selected no complete Pre/Post x Open/Close sets. Falling back to Pre RestOpen anchors.');
        selected_rows = anchor_filter;
    end
end

function label = make_session_label(meta)
    animal_str = char(string(meta.animal_name));
    injection_str = char(string(meta.injection));
    label = sprintf('%s %s session %s', animal_str, injection_str, char(string(meta.session_idx)));
end


function loaded_states = load_required_states_for_session(root, meta, preposts, states, kernel_indices)
    loaded_states = struct();
    for k_i = 1:numel(kernel_indices)
        kernel_idx = kernel_indices(k_i);
        for pp_i = 1:numel(preposts)
            for st_i = 1:numel(states)
                prepost = preposts{pp_i};
                state = states{st_i};
                key = state_key(prepost, state, kernel_idx);
                loaded_states.(key) = load_state_connectivity(root, meta, prepost, state, kernel_idx);
            end
        end
    end
end

function key = state_key(prepost, state, kernel_idx)
    key = sprintf('k%d_%s_%s', kernel_idx, prepost, state);
end

function state_struct = load_state_connectivity(root, meta, prepost, state, kernel_idx)
    state_struct = struct();
    state_struct.prepost = prepost;
    state_struct.state = state;
    state_struct.kernel_idx = kernel_idx;

    meta.prepost = prepost;
    meta.state = state;

    meta.filename = generate_filename('raster', meta);
    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);

    cell_area = raster_data.data.cell_area;
    filter1 = ismember(cell_area, {'ACC'});
    filter2 = ismember(cell_area, {'VLPFC'});

    meta.filename = generate_filename('GLM', meta);
    GLM_data = load(fullfile(root, 'Data', 'Working', 'GLM', meta.filename));
    fprintf('Loaded GLM data for %s %s %s\n', meta.prepost, meta.state, meta.area);

    N = GLM_data.meta.N;
    J = GLM_data.data.model_par(:, ((2 + N * (kernel_idx - 1)) : (1 + N * kernel_idx)));
    err = GLM_data.data.model_err.total(:, ((2 + N * (kernel_idx - 1)) : (1 + N * kernel_idx)));

    idx1 = find(filter1);
    idx2 = find(filter2);

    state_struct.cell_area = cell_area;
    state_struct.filter1 = filter1;
    state_struct.filter2 = filter2;
    state_struct.idx1 = idx1;
    state_struct.idx2 = idx2;
    state_struct.J = J;
    state_struct.err = err;
    state_struct.J12 = J(filter1, filter2);
    state_struct.J21 = J(filter2, filter1);
    state_struct.err12 = err(filter1, filter2);
    state_struct.err21 = err(filter2, filter1);
end

function data = make_pair_vectors(state_x, state_y, err_multi, meta, session_label)
    validate_matching_filters(state_x, state_y);

    x12 = state_x.J12(:);
    y12 = state_y.J12(:);
    x21 = state_x.J21(:);
    y21 = state_y.J21(:);

    xerr12 = state_x.err12(:);
    yerr12 = state_y.err12(:);
    xerr21 = state_x.err21(:);
    yerr21 = state_y.err21(:);

    [from12_global, to12_global, from12_area, to12_area] = make_connection_index_vectors(state_x.idx1, state_x.idx2);
    [from21_global, to21_global, from21_area, to21_area] = make_connection_index_vectors(state_x.idx2, state_x.idx1);

    x = [x12; x21];
    y = [y12; y21];
    xerr = [xerr12; xerr21];
    yerr = [yerr12; yerr21];

    from_cell_global_idx = [from12_global; from21_global];
    to_cell_global_idx = [to12_global; to21_global];
    from_cell_area_idx = [from12_area; from21_area];
    to_cell_area_idx = [to12_area; to21_area];

    n12 = numel(x12);
    n21 = numel(x21);
    from_area = [repmat("ACC", n12, 1); repmat("VLPFC", n21, 1)];
    to_area = [repmat("VLPFC", n12, 1); repmat("ACC", n21, 1)];
    connection_direction = [repmat("ACC_to_VLPFC", n12, 1); repmat("VLPFC_to_ACC", n21, 1)];

    valid = isfinite(x) & isfinite(y) & isfinite(xerr) & isfinite(yerr);
    x = x(valid);
    y = y(valid);
    xerr = xerr(valid);
    yerr = yerr(valid);
    from_cell_global_idx = from_cell_global_idx(valid);
    to_cell_global_idx = to_cell_global_idx(valid);
    from_cell_area_idx = from_cell_area_idx(valid);
    to_cell_area_idx = to_cell_area_idx(valid);
    from_area = from_area(valid);
    to_area = to_area(valid);
    connection_direction = connection_direction(valid);

    x_pos = x >  err_multi * xerr;
    x_neg = x < -err_multi * xerr;
    y_pos = y >  err_multi * yerr;
    y_neg = y < -err_multi * yerr;

    x_cat = zeros(size(x));
    y_cat = zeros(size(y));
    x_cat(x_pos) = 1;
    x_cat(x_neg) = -1;
    y_cat(y_pos) = 1;
    y_cat(y_neg) = -1;

    n = numel(x);

    data = struct();
    data.x_orig = x;
    data.y_orig = y;
    data.x_err_orig = xerr;
    data.y_err_orig = yerr;
    data.x_cat_orig = x_cat;
    data.y_cat_orig = y_cat;

    data.session_label = repmat(string(session_label), n, 1);
    data.animal_name = repmat(string(meta.animal_name), n, 1);
    data.injection = repmat(string(meta.injection), n, 1);
    data.session_idx = repmat(double(meta.session_idx), n, 1);
    data.connection_direction = connection_direction;
    data.from_area = from_area;
    data.to_area = to_area;

    data.from_cell_area_idx = from_cell_area_idx;
    data.to_cell_area_idx = to_cell_area_idx;
    data.from_cell_global_idx = from_cell_global_idx;
    data.to_cell_global_idx = to_cell_global_idx;

    % Backward-compatible aliases: these are global cell indices.
    data.from_cell_idx = from_cell_global_idx;
    data.to_cell_idx = to_cell_global_idx;
end

function [from_global_vec, to_global_vec, from_area_vec, to_area_vec] = make_connection_index_vectors(from_idx, to_idx)
    [from_area_grid, to_area_grid] = ndgrid((1:numel(from_idx)).', (1:numel(to_idx)).');

    from_area_vec = from_area_grid(:);
    to_area_vec = to_area_grid(:);

    from_global_vec = from_idx(from_area_vec);
    to_global_vec = to_idx(to_area_vec);

    from_global_vec = from_global_vec(:);
    to_global_vec = to_global_vec(:);
end

function plot_data = prepare_plot_data(data, plot_type, x_label_base, y_label_base)
    x_sig = data.x_cat_orig ~= 0;
    y_sig = data.y_cat_orig ~= 0;

    x_pos = data.x_cat_orig == 1;
    x_neg = data.x_cat_orig == -1;
    y_pos = data.y_cat_orig == 1;
    y_neg = data.y_cat_orig == -1;

    switch lower(plot_type)
        case {'all_connection', 'all_connections', 'all'}
            mask = true(size(data.x_orig));
            x_plot = data.x_orig;
            y_plot = data.y_orig;
            group = strings(size(x_plot));
            group(~x_sig & ~y_sig) = "both non-sig";
            group(x_sig & ~y_sig) = sprintf('%s only', x_label_base);
            group(~x_sig & y_sig) = sprintf('%s only', y_label_base);
            group(x_sig & y_sig) = "both sig";
            group_order = ["both non-sig", sprintf('%s only', x_label_base), sprintf('%s only', y_label_base), "both sig"];
            x_label = sprintf('%s J_{ij}', x_label_base);
            y_label = sprintf('%s J_{ij}', y_label_base);
            plot_title = 'All connections';

        case {'either_significant', 'either_sig'}
            mask = x_sig | y_sig;
            x_plot = data.x_orig(mask);
            y_plot = data.y_orig(mask);
            group = strings(size(x_plot));
            x_sig_m = x_sig(mask);
            y_sig_m = y_sig(mask);
            group(x_sig_m & ~y_sig_m) = sprintf('%s only', x_label_base);
            group(~x_sig_m & y_sig_m) = sprintf('%s only', y_label_base);
            group(x_sig_m & y_sig_m) = "both sig";
            group_order = [sprintf('%s only', x_label_base), sprintf('%s only', y_label_base), "both sig"];
            x_label = sprintf('%s J_{ij}', x_label_base);
            y_label = sprintf('%s J_{ij}', y_label_base);
            plot_title = 'Either significant';

        case {'both_significant', 'both_sig'}
            mask = x_sig & y_sig;
            x_plot = data.x_orig(mask);
            y_plot = data.y_orig(mask);
            x_cat = data.x_cat_orig(mask);
            y_cat = data.y_cat_orig(mask);
            group = strings(size(x_plot));
            group(x_cat == 1 & y_cat == 1) = "pos";
            group(x_cat == -1 & y_cat == -1) = "neg";
            group(x_cat ~= y_cat) = "switch";
            group_order = ["pos", "neg", "switch"];
            x_label = sprintf('%s J_{ij}', x_label_base);
            y_label = sprintf('%s J_{ij}', y_label_base);
            plot_title = 'Both significant';

        case {'same_sign_abs', 'same_abs'}
            mask = (x_pos & y_pos) | (x_neg & y_neg);
            x_plot = abs(data.x_orig(mask));
            y_plot = abs(data.y_orig(mask));
            x_cat = data.x_cat_orig(mask);
            group = strings(size(x_plot));
            group(x_cat == 1) = "pos";
            group(x_cat == -1) = "neg";
            group_order = ["pos", "neg"];
            x_label = sprintf('%s |J_{ij}|', x_label_base);
            y_label = sprintf('%s |J_{ij}|', y_label_base);
            plot_title = 'Same-sign abs';

        case {'switch_sign_abs', 'switch_abs'}
            mask = (x_pos & y_neg) | (x_neg & y_pos);
            x_plot = abs(data.x_orig(mask));
            y_plot = abs(data.y_orig(mask));
            x_cat = data.x_cat_orig(mask);
            y_cat = data.y_cat_orig(mask);
            group = strings(size(x_plot));
            group(x_cat == 1 & y_cat == -1) = sprintf('%s pos / %s neg', x_label_base, y_label_base);
            group(x_cat == -1 & y_cat == 1) = sprintf('%s neg / %s pos', x_label_base, y_label_base);
            group_order = [sprintf('%s pos / %s neg', x_label_base, y_label_base), ...
                           sprintf('%s neg / %s pos', x_label_base, y_label_base)];
            x_label = sprintf('%s |J_{ij}|', x_label_base);
            y_label = sprintf('%s |J_{ij}|', y_label_base);
            plot_title = 'Switch-sign abs';

        otherwise
            error('Unknown plot_type: %s', plot_type);
    end

    plot_data = subset_raw_data_for_plot(data, mask);
    plot_data.x_plot = x_plot;
    plot_data.y_plot = y_plot;
    plot_data.group = group;
    plot_data.group_order = group_order;
    plot_data.x_label = x_label;
    plot_data.y_label = y_label;
    plot_data.plot_title = plot_title;
    plot_data.plot_type = string(plot_type);
    plot_data.x_raw = data.x_orig(mask);
    plot_data.y_raw = data.y_orig(mask);
end

function out = subset_raw_data_for_plot(data, mask)
    fields = fieldnames(data);
    out = struct();
    for k = 1:numel(fields)
        field = fields{k};
        values = data.(field);
        if numel(values) == numel(mask)
            out.(field) = values(mask);
        else
            out.(field) = values;
        end
    end
end

function out = subset_plot_data(data, mask)
    fields = fieldnames(data);
    out = struct();
    for k = 1:numel(fields)
        field = fields{k};
        values = data.(field);
        if numel(values) == numel(mask)
            out.(field) = values(mask);
        else
            out.(field) = values;
        end
    end
end

function mask = points_in_region(x, y, x_range, y_range)
    xr = sort(x_range);
    yr = sort(y_range);
    mask = isfinite(x) & isfinite(y) & ...
           x >= xr(1) & x <= xr(2) & ...
           y >= yr(1) & y <= yr(2);
end

function plot_scatter_groups(ax, plot_data, marker_size, marker_alpha, colors, show_legend)
    cla(ax);
    hold(ax, 'on');

    for g = 1:numel(plot_data.group_order)
        group_name = plot_data.group_order(g);
        group_mask = plot_data.group == group_name;
        if ~any(group_mask)
            continue;
        end

        scatter(ax, plot_data.x_plot(group_mask), plot_data.y_plot(group_mask), marker_size, 'filled', ...
            'MarkerFaceColor', group_color(group_name, colors), ...
            'MarkerFaceAlpha', marker_alpha, ...
            'MarkerEdgeAlpha', marker_alpha, ...
            'DisplayName', char(group_name));
    end

    axis_limit = get_axis_limit_for_plot(plot_data, plot_data.plot_type);
    plot(ax, axis_limit.x, axis_limit.y, '--', 'Color', colors.identity_line, 'LineWidth', 1, 'HandleVisibility', 'off');

    if ~contains(lower(plot_data.plot_type), 'abs')
        xline(ax, 0, ':', 'Color', colors.zero_line, 'HandleVisibility', 'off');
        yline(ax, 0, ':', 'Color', colors.zero_line, 'HandleVisibility', 'off');
    end

    hold(ax, 'off');

    if show_legend
        legend(ax, 'Location', 'northeastoutside');
    else
        legend(ax, 'off');
    end
end

function c = group_color(group_name, colors)
    group_name = string(group_name);
    if group_name == "both non-sig"
        c = colors.non_sig;
    elseif group_name == "both sig"
        c = colors.both_sig;
    elseif group_name == "pos"
        c = colors.pos;
    elseif group_name == "neg"
        c = colors.neg;
    elseif group_name == "switch"
        c = colors.switch;
    elseif contains(group_name, "pos /")
        c = colors.switch_xy;
    elseif contains(group_name, "neg /")
        c = colors.switch_yx;
    elseif contains(group_name, "only")
        if contains(group_name, "Pre") || contains(group_name, "Open") || contains(group_name, "K")
            c = colors.x_only;
        else
            c = colors.y_only;
        end
    else
        c = [0.2, 0.2, 0.2];
    end
end

function axis_limit = get_axis_limit_for_plot(plot_data, plot_type)
    if contains(lower(plot_type), 'abs')
        vals = [plot_data.x_plot(:); plot_data.y_plot(:)];
        vals = vals(isfinite(vals));
        if isempty(vals)
            vmax = 1;
        else
            vmax = max(vals);
            if vmax <= 0
                vmax = 1;
            end
        end
        vmax = max(3.5, 1.05 * vmax);
        axis_limit.x = [0, vmax];
        axis_limit.y = [0, vmax];
    else
        vals = [plot_data.x_plot(:); plot_data.y_plot(:)];
        vals = vals(isfinite(vals));
        if isempty(vals)
            vmax = 1;
        else
            vmax = max(abs(vals));
            if vmax <= 0
                vmax = 1;
            end
        end
        vmax = max(3.5, 1.05 * vmax);
        axis_limit.x = [-vmax, vmax];
        axis_limit.y = [-vmax, vmax];
    end
end

function draw_region_box(ax, x_range, y_range, color)
    xr = sort(x_range);
    yr = sort(y_range);
    rectangle(ax, 'Position', [xr(1), yr(1), diff(xr), diff(yr)], ...
        'EdgeColor', color, 'LineWidth', 2);
end

function print_selected_connections(selected_data)
    fprintf('\n===== Selected connections: n = %d =====\n', numel(selected_data.x_plot));
    for i = 1:numel(selected_data.x_plot)
        fprintf('plot_x=%.6g, plot_y=%.6g, raw_x=%.6g, raw_y=%.6g, session="%s", %s(area_idx=%d, global_idx=%d) -> %s(area_idx=%d, global_idx=%d), direction=%s, group=%s\n', ...
            selected_data.x_plot(i), selected_data.y_plot(i), ...
            selected_data.x_raw(i), selected_data.y_raw(i), ...
            selected_data.session_label(i), ...
            selected_data.from_area(i), selected_data.from_cell_area_idx(i), selected_data.from_cell_global_idx(i), ...
            selected_data.to_area(i), selected_data.to_cell_area_idx(i), selected_data.to_cell_global_idx(i), ...
            selected_data.connection_direction(i), selected_data.group(i));
    end
end

function write_selected_connections_to_file(selected_data, text_path)
    fid = fopen(text_path, 'w');
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'plot_x\tplot_y\traw_x\traw_y\tsession\tfrom_area\tfrom_cell_area_idx\tfrom_cell_global_idx\tto_area\tto_cell_area_idx\tto_cell_global_idx\tdirection\tgroup\n');
    for i = 1:numel(selected_data.x_plot)
        fprintf(fid, '%.12g\t%.12g\t%.12g\t%.12g\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n', ...
            selected_data.x_plot(i), selected_data.y_plot(i), ...
            selected_data.x_raw(i), selected_data.y_raw(i), ...
            selected_data.session_label(i), ...
            selected_data.from_area(i), selected_data.from_cell_area_idx(i), selected_data.from_cell_global_idx(i), ...
            selected_data.to_area(i), selected_data.to_cell_area_idx(i), selected_data.to_cell_global_idx(i), ...
            selected_data.connection_direction(i), selected_data.group(i));
    end
end

function validate_matching_filters(state_a, state_b)
    if numel(state_a.filter1) ~= numel(state_b.filter1) || any(state_a.filter1(:) ~= state_b.filter1(:)) || ...
       numel(state_a.filter2) ~= numel(state_b.filter2) || any(state_a.filter2(:) ~= state_b.filter2(:))
        error('State filters do not match between the two comparison members.');
    end

    if ~isequal(size(state_a.J12), size(state_b.J12)) || ~isequal(size(state_a.J21), size(state_b.J21))
        error('State connectivity matrices do not have matching sizes.');
    end
end

function safe_name = sanitize_filename(name)
    safe_name = regexprep(char(string(name)), '[^A-Za-z0-9_-]', '_');
end
