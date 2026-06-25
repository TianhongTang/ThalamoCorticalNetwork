%% J_bidirectional_scatter_hybrid_abs.m - Reciprocal directional J scatter plots.
% Plots paired bidirectional connections: x = area1 -> area2, y = area2 -> area1.
% Edit the configuration block below to change areas, conditions, kernels, or filters.

clear;
set(0, 'DefaultFigureVisible', 'off');

%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
if ~isfolder(fullfile(root, 'Data')) && isfolder(fullfile(pwd, 'Data'))
    root = pwd;
end
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Configuration
% If J(row, col) means row -> col, use 'row_to_col'.
% If J(row, col) means col -> row, use 'col_to_row'.
matrix_direction = 'col_to_row';

area_pairs = struct([]);
% area_pairs(1).area1 = 'ACC';
% area_pairs(1).area2 = 'VLPFC';
% area_pairs(1).label1 = 'ACC';
% area_pairs(1).label2 = 'VLPFC';
area_pairs(1).area1 = 'Thalamus';
area_pairs(1).area2 = 'ACC';
area_pairs(1).label1 = 'MD';
area_pairs(1).label2 = 'ACC';
area_pairs(2).area1 = 'Thalamus';
area_pairs(2).area2 = 'VLPFC';
area_pairs(2).label1 = 'MD';
area_pairs(2).label2 = 'VLPFC';

kernel_indices = 1:3;
% conditions = [ ...
%     make_condition('Pre',  'RestOpen',  'Pre, Eyes Open'), ...
%     make_condition('Pre',  'RestClose', 'Pre, Eyes Closed'), ...
%     make_condition('Post', 'RestOpen',  'Post, Eyes Open'), ...
%     make_condition('Post', 'RestClose', 'Post, Eyes Closed') ...
% ];
% State groups are analysis labels, not metadata states.
% Here, 'All States' pools RestOpen and RestClose into one condition.
state_groups = [ ...
    make_state_group('All States', {'RestOpen', 'RestClose'}) ...
];

conditions = [ ...
    make_condition('Pre',  'RestOpen',  'Pre, Eyes Open'), ...
    make_condition('Pre',  'RestClose', 'Pre, Eyes Closed'), ...
    make_condition('Pre',  'All States', 'Pre, All States'), ...
];

filter_opts = struct();
filter_opts.kernel_name = "DeltaPure";
filter_opts.align = 'Last';
filter_opts.area = "Full";
filter_opts.injection = 'Muscimol';
filter_opts.epoch = 3000;
filter_opts.fold_idx = 0;
filter_opts.shuffle_idx = 0;
filter_opts.resting_dur_threshold = 15;

err_multi = 1;
scatter_marker_size = 8;
scatter_alpha = 0.25;
figure_visible = 'off';
export_pdf = false;
show_legend = false;
show_identity_line = true;
show_fit_line = true;
fit_line_method = 'tls'; % 'ols' or 'tls'.
show_histogram = true;
combine_kernels = false;
by_session = false;
J_lim = 2; % Use NaN to choose the full J matrix color limit automatically.
lim_ratio = 0.99; % When J_lim is NaN, use this central ratio of finite J values to set the color limit.
panel_label_pos = [-0.18, 1.18];
% Scatter modes are composable. Keywords are applied left to right:
% either_sig, both_sig, same_sign, diff_sign, abs, x_flip, y_flip, xy_flip.
scatter_panel_modes = { ...
    'all', ...
    'either_sig', ...
    'both_sig', ...
    'abs', ...
    'either_sig abs', ...
    'both_sig abs', ...
    'xy_flip' ...
    'either_sig xy_flip', ...
    'both_sig xy_flip', ...
};
max_sessions_to_include = inf; % set smaller for debugging.

colors = struct();
colors.non_sig = [0.65, 0.65, 0.65];
colors.x_only = [0.10, 0.75, 0.10];
colors.y_only = [0.95, 0.55, 0.10];
colors.both_sig = [0.45, 0.10, 0.75];
colors.pos = [1.00, 0.20, 0.20];
colors.neg = [0.10, 0.45, 1.00];
colors.switch = [0.45, 0.10, 0.75];
colors.identity_line = [1, 0, 0];
colors.zero_line = [0, 0, 0];
colors.fit_line = [0.20, 0.20, 0.20];

%% Load metadata and select anchor sessions
mt_all = load_meta(root, 'table');
mt = mt_all.GLM;
selected_rows = select_complete_sessions(mt, filter_opts, conditions);
selected_mt = mt(selected_rows, :);
meta_array = table2struct(selected_mt);
if isfinite(max_sessions_to_include)
    meta_array = meta_array(1:min(numel(meta_array), max_sessions_to_include));
end
if isempty(meta_array)
    error('No metadata rows selected.');
end

%% Render figures
if by_session
    render_by_session_figures(root, meta_array, area_pairs, conditions, kernel_indices, state_groups, ...
        matrix_direction, err_multi, colors, scatter_marker_size, scatter_alpha, figure_visible, ...
        show_legend, show_identity_line, show_fit_line, fit_line_method, show_histogram, ...
        combine_kernels, J_lim, lim_ratio, export_pdf, panel_label_pos, scatter_panel_modes);
else
    render_pooled_figures(root, meta_array, area_pairs, conditions, kernel_indices, state_groups, ...
        matrix_direction, err_multi, colors, scatter_marker_size, scatter_alpha, figure_visible, ...
        show_legend, show_identity_line, show_fit_line, fit_line_method, show_histogram, ...
        combine_kernels, export_pdf, panel_label_pos, scatter_panel_modes);
end


function cond = make_condition(prepost, state, title_text)
    cond = struct();
    cond.prepost = prepost;
    cond.state = state;
    cond.title = title_text;
end

function group = make_state_group(name, states)
    group = struct();
    group.name = char(string(name));
    group.states = cellstr(string(states));
end

function render_pooled_figures(root, meta_array, area_pairs, conditions, kernel_indices, state_groups, ...
    matrix_direction, err_multi, colors, scatter_marker_size, scatter_alpha, figure_visible, ...
    show_legend, show_identity_line, show_fit_line, fit_line_method, show_histogram, combine_kernels, export_pdf, panel_label_pos, scatter_panel_modes)

    render_total = count_planned_figures(numel(area_pairs), numel(kernel_indices), combine_kernels);
    render_count = 0;

    for pair_idx = 1:numel(area_pairs)
        area_pair = area_pairs(pair_idx);
        pooled = initialize_pooled_data(conditions, kernel_indices);
        valid_session_counts = zeros(numel(conditions), numel(kernel_indices));

        for session_i = 1:numel(meta_array)
            anchor_meta = meta_array(session_i);
            session_label = make_session_label(anchor_meta);
            fprintf('\n===== Loading %s (%d/%d) =====\n', session_label, session_i, numel(meta_array));

            for condition_idx = 1:numel(conditions)
                cond = conditions(condition_idx);
                for kernel_i = 1:numel(kernel_indices)
                    kernel_idx = kernel_indices(kernel_i);
                    try
                        data = load_condition_reciprocal_data(root, anchor_meta, cond, kernel_idx, area_pair, matrix_direction, err_multi, state_groups);
                        data.J_full = []; % Pooled figures intentionally leave the J-matrix panel blank.
                        data.J_sig_full = [];
                        data.J_borders = [];
                        pooled{condition_idx, kernel_i} = append_reciprocal_data(pooled{condition_idx, kernel_i}, data);
                        valid_session_counts(condition_idx, kernel_i) = valid_session_counts(condition_idx, kernel_i) + 1;
                    catch ME
                        warning('Skipping %s, %s, K%d, %s-%s: %s', ...
                            session_label, cond.title, kernel_idx, ...
                            join_area_names(area_pair.area1), join_area_names(area_pair.area2), ME.message);
                    end
                end
            end
        end

        if combine_kernels
            render_count = render_count + 1;
            if ~has_any_reciprocal_data(pooled)
                print_skipped_figure(render_count, render_total, sprintf('pooled %s-%s', area_pair.label1, area_pair.label2));
                continue;
            end
            render_bidirectional_figure(root, pooled, valid_session_counts, area_pair, matrix_direction, ...
                conditions, kernel_indices, colors, scatter_marker_size, scatter_alpha, figure_visible, ...
                show_legend, show_identity_line, show_fit_line, fit_line_method, show_histogram, ...
                NaN, NaN, '', false, 'J_bidirectional_scatter', export_pdf, panel_label_pos, scatter_panel_modes, render_count, render_total);
        else
            for kernel_i = 1:numel(kernel_indices)
                render_count = render_count + 1;
                if ~has_any_reciprocal_data(pooled(:, kernel_i))
                    print_skipped_figure(render_count, render_total, sprintf('pooled %s-%s k%d', ...
                        area_pair.label1, area_pair.label2, kernel_indices(kernel_i)));
                    continue;
                end
                render_bidirectional_figure(root, pooled(:, kernel_i), valid_session_counts(:, kernel_i), ...
                    area_pair, matrix_direction, conditions, kernel_indices(kernel_i), colors, ...
                    scatter_marker_size, scatter_alpha, figure_visible, show_legend, show_identity_line, ...
                    show_fit_line, fit_line_method, show_histogram, NaN, NaN, sprintf('_k%d', kernel_indices(kernel_i)), ...
                    false, 'J_bidirectional_scatter', export_pdf, panel_label_pos, scatter_panel_modes, render_count, render_total);
            end
        end
    end
end

function render_by_session_figures(root, meta_array, area_pairs, conditions, kernel_indices, state_groups, ...
    matrix_direction, err_multi, colors, scatter_marker_size, scatter_alpha, figure_visible, ...
    show_legend, show_identity_line, show_fit_line, fit_line_method, show_histogram, combine_kernels, J_lim, lim_ratio, export_pdf, panel_label_pos, scatter_panel_modes)

    render_total = count_planned_figures(numel(meta_array) * numel(area_pairs), numel(kernel_indices), combine_kernels);
    render_count = 0;

    for session_i = 1:numel(meta_array)
        anchor_meta = meta_array(session_i);
        session_label = make_session_label(anchor_meta);
        fprintf('\n===== Loading %s (%d/%d) =====\n', session_label, session_i, numel(meta_array));

        for pair_idx = 1:numel(area_pairs)
            area_pair = area_pairs(pair_idx);
            pooled = initialize_pooled_data(conditions, kernel_indices);
            valid_session_counts = zeros(numel(conditions), numel(kernel_indices));

            for condition_idx = 1:numel(conditions)
                cond = conditions(condition_idx);
                for kernel_i = 1:numel(kernel_indices)
                    kernel_idx = kernel_indices(kernel_i);
                    try
                        data = load_condition_reciprocal_data(root, anchor_meta, cond, kernel_idx, area_pair, matrix_direction, err_multi, state_groups);
                        pooled{condition_idx, kernel_i} = append_reciprocal_data(pooled{condition_idx, kernel_i}, data);
                        valid_session_counts(condition_idx, kernel_i) = valid_session_counts(condition_idx, kernel_i) + 1;
                    catch ME
                        warning('Skipping %s, %s, K%d, %s-%s: %s', ...
                            session_label, cond.title, kernel_idx, ...
                            join_area_names(area_pair.area1), join_area_names(area_pair.area2), ME.message);
                    end
                end
            end

            session_suffix = ['_', sanitize_filename(session_label)];
            if combine_kernels
                render_count = render_count + 1;
                if ~has_any_reciprocal_data(pooled)
                    print_skipped_figure(render_count, render_total, sprintf('%s %s-%s', ...
                        session_label, area_pair.label1, area_pair.label2));
                    continue;
                end
                render_bidirectional_figure(root, pooled, valid_session_counts, area_pair, matrix_direction, ...
                    conditions, kernel_indices, colors, scatter_marker_size, scatter_alpha, figure_visible, ...
                    show_legend, show_identity_line, show_fit_line, fit_line_method, show_histogram, ...
                    J_lim, lim_ratio, session_suffix, true, 'J_bidirectional_scatter_by_session', export_pdf, panel_label_pos, scatter_panel_modes, render_count, render_total);
            else
                for kernel_i = 1:numel(kernel_indices)
                    render_count = render_count + 1;
                    if ~has_any_reciprocal_data(pooled(:, kernel_i))
                        print_skipped_figure(render_count, render_total, sprintf('%s %s-%s k%d', ...
                            session_label, area_pair.label1, area_pair.label2, kernel_indices(kernel_i)));
                        continue;
                    end
                    render_bidirectional_figure(root, pooled(:, kernel_i), valid_session_counts(:, kernel_i), ...
                        area_pair, matrix_direction, conditions, kernel_indices(kernel_i), colors, ...
                        scatter_marker_size, scatter_alpha, figure_visible, show_legend, show_identity_line, ...
                        show_fit_line, fit_line_method, show_histogram, J_lim, lim_ratio, ...
                        sprintf('%s_k%d', session_suffix, kernel_indices(kernel_i)), ...
                        true, 'J_bidirectional_scatter_by_session', export_pdf, panel_label_pos, scatter_panel_modes, render_count, render_total);
                end
            end
        end
    end
end

function pooled = initialize_pooled_data(conditions, kernel_indices)
    pooled = cell(numel(conditions), numel(kernel_indices));
    for condition_idx = 1:numel(conditions)
        for kernel_i = 1:numel(kernel_indices)
            pooled{condition_idx, kernel_i} = empty_reciprocal_data();
        end
    end
end

function n_fig = count_planned_figures(n_groups, n_kernel, combine_kernels)
    if combine_kernels
        n_fig = n_groups;
    else
        n_fig = n_groups * n_kernel;
    end
end

function print_skipped_figure(render_count, render_total, label)
    fprintf('Skipping figure %d/%d: %s has no valid reciprocal data.\n', render_count, render_total, label);
end

function data = empty_reciprocal_data()
    data = struct();
    data.x = [];
    data.y = [];
    data.x_err = [];
    data.y_err = [];
    data.x_cat = [];
    data.y_cat = [];
    data.J_full = [];
    data.J_sig_full = [];
    data.J_borders = [];
end

function out = append_reciprocal_data(out, incoming)
    fields = fieldnames(out);
    for k = 1:numel(fields)
        field_name = fields{k};
        if strcmp(field_name, 'J_full') || strcmp(field_name, 'J_sig_full') || strcmp(field_name, 'J_borders')
            if isempty(out.(field_name))
                out.(field_name) = incoming.(field_name);
            end
        else
            out.(field_name) = [out.(field_name); incoming.(field_name)];
        end
    end
end

function has_data = has_any_reciprocal_data(pooled)
    has_data = false;
    for k = 1:numel(pooled)
        if ~isempty(pooled{k}.x)
            has_data = true;
            return;
        end
    end
end

function selected_rows = select_complete_sessions(mt, filter_opts, conditions)
    %#ok<INUSD> Select one anchor row per session; missing condition panels are skipped during loading.
    base_filter = strcmp(mt.kernel_name, filter_opts.kernel_name) & ...
                  strcmp(mt.align, filter_opts.align) & ...
                  strcmp(mt.area, filter_opts.area) & ...
                  strcmp(mt.injection, filter_opts.injection) & ...
                  (mt.epoch == filter_opts.epoch) & ...
                  (mt.fold_idx == filter_opts.fold_idx) & ...
                  (mt.shuffle_idx == filter_opts.shuffle_idx) & ...
                  cellfun(@(x) ~isempty(x) && x == filter_opts.resting_dur_threshold, mt.resting_dur_threshold);

    selected_rows = false(height(mt), 1);
    base_idx = find(base_filter);

    match_fields = {'animal_name', 'injection', 'align', 'session_idx', ...
                    'resting_dur_threshold', 'area', 'kernel_name', ...
                    'reg_name', 'epoch', 'fold_idx', 'shuffle_idx'};

    for k = 1:numel(base_idx)
        idx = base_idx(k);
        already_selected = false;
        selected_idx = find(selected_rows).';

        for prev_idx = selected_idx
            same_as_previous = true;
            for f = 1:numel(match_fields)
                field = match_fields{f};
                if ~ismember(field, mt.Properties.VariableNames)
                    continue;
                end

                if iscell(mt.(field))
                    same_as_previous = same_as_previous && isequal(mt.(field){idx}, mt.(field){prev_idx});
                else
                    same_as_previous = same_as_previous && isequal(mt.(field)(idx), mt.(field)(prev_idx));
                end
            end
            if same_as_previous
                already_selected = true;
                break;
            end
        end

        if ~already_selected
            selected_rows(idx) = true;
        end
    end

    fprintf('Selected %d session anchors.\n', sum(selected_rows));
end

function label = make_session_label(meta)
    label = sprintf('%s %s session %s', char(string(meta.animal_name)), ...
        char(string(meta.injection)), char(string(meta.session_idx)));
end

function data = load_condition_reciprocal_data(root, anchor_meta, cond, kernel_idx, area_pair, matrix_direction, err_multi, state_groups)
    state_list = resolve_condition_states(cond, state_groups);
    data = empty_reciprocal_data();
    loaded_count = 0;
    error_messages = {};

    for state_i = 1:numel(state_list)
        this_cond = cond;
        this_cond.state = state_list{state_i};
        try
            incoming = load_reciprocal_direction_data(root, anchor_meta, this_cond, kernel_idx, area_pair, matrix_direction, err_multi);
            data = append_reciprocal_data(data, incoming);
            loaded_count = loaded_count + 1;
        catch ME
            error_messages{end + 1} = sprintf('%s: %s', state_list{state_i}, ME.message); %#ok<AGROW>
        end
    end

    if loaded_count == 0
        if isempty(error_messages)
            error_detail = 'no states requested';
        else
            error_detail = strjoin(error_messages, ' | ');
        end
        error('No data loaded for condition "%s" from states {%s}: %s', ...
            cond.title, strjoin(state_list, ', '), error_detail);
    elseif ~isempty(error_messages)
        warning('Partially loaded condition "%s" from states {%s}. Missing/failed: %s', ...
            cond.title, strjoin(state_list, ', '), strjoin(error_messages, ' | '));
    end

    % A grouped state pools multiple J matrices into one analysis condition.
    % There is no single representative full N x N matrix, so leave matrix panels blank.
    if numel(state_list) > 1
        data.J_full = [];
        data.J_sig_full = [];
        data.J_borders = [];
    end
end

function state_list = resolve_condition_states(cond, state_groups)
    if iscell(cond.state)
        state_list = cellstr(string(cond.state));
        return;
    end

    state_value = string(cond.state);
    if numel(state_value) > 1
        state_list = cellstr(state_value);
        return;
    end

    requested_state = char(state_value);
    state_list = {requested_state};
    requested_key = normalize_state_group_key(requested_state);
    for group_i = 1:numel(state_groups)
        if strcmp(requested_key, normalize_state_group_key(state_groups(group_i).name))
            state_list = cellstr(string(state_groups(group_i).states));
            return;
        end
    end
end

function key = normalize_state_group_key(name)
    key = lower(regexprep(char(string(name)), '[^a-zA-Z0-9]+', ''));
end

function data = load_reciprocal_direction_data(root, anchor_meta, cond, kernel_idx, area_pair, matrix_direction, err_multi)
    meta = anchor_meta;
    meta.prepost = cond.prepost;
    meta.state = cond.state;

    meta.filename = generate_filename('raster', meta);
    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    cell_area = raster_data.data.cell_area;

    filter1 = ismember(cellstr(string(cell_area)), cellstr(string(area_pair.area1)));
    filter2 = ismember(cellstr(string(cell_area)), cellstr(string(area_pair.area2)));
    if ~any(filter1)
        error('No cells found for area1: %s.', join_area_names(area_pair.area1));
    end
    if ~any(filter2)
        error('No cells found for area2: %s.', join_area_names(area_pair.area2));
    end

    meta.filename = generate_filename('GLM', meta);
    GLM_data = load(fullfile(root, 'Data', 'Working', 'GLM', meta.filename));
    N = GLM_data.meta.N;
    J = GLM_data.data.model_par(:, (2 + N * (kernel_idx - 1)):(1 + N * kernel_idx));
    err = GLM_data.data.model_err.total(:, (2 + N * (kernel_idx - 1)):(1 + N * kernel_idx));

    J12 = J(filter1, filter2);
    J21 = J(filter2, filter1);
    err12 = err(filter1, filter2);
    err21 = err(filter2, filter1);

    switch lower(matrix_direction)
        case 'row_to_col'
            x = J12(:);
            y = J21.'; y = y(:);
            x_err = err12(:);
            y_err = err21.'; y_err = y_err(:);
        case 'col_to_row'
            x = J21.'; x = x(:);
            y = J12(:);
            x_err = err21.'; x_err = x_err(:);
            y_err = err12(:);
        otherwise
            error('Unknown matrix_direction: %s. Use row_to_col or col_to_row.', matrix_direction);
    end

    valid = isfinite(x) & isfinite(y) & isfinite(x_err) & isfinite(y_err);
    x = x(valid);
    y = y(valid);
    x_err = x_err(valid);
    y_err = y_err(valid);

    data = empty_reciprocal_data();
    data.x = x;
    data.y = y;
    data.x_err = x_err;
    data.y_err = y_err;
    data.x_cat = classify_connections(x, x_err, err_multi);
    data.y_cat = classify_connections(y, y_err, err_multi);
    data.J_full = J;
    J_sig_full = J;
    sig_mask_full = isfinite(J) & isfinite(err) & ((J > err_multi * err) | (J < -err_multi * err));
    J_sig_full(~sig_mask_full) = 0;
    data.J_sig_full = J_sig_full;
    data.J_borders = load_area_borders(root, meta, N);
end

function borders = load_area_borders(root, meta, N)
    borders = [];
    border_folder = fullfile(root, 'Data', 'Working', 'border');
    border_file_name = generate_filename('border', meta);
    border_file_path = fullfile(border_folder, border_file_name);
    if ~isfile(border_file_path)
        warning('Border file not found for J matrix: %s', border_file_path);
        return;
    end

    border_data = load(border_file_path);
    if isfield(border_data, 'data') && isfield(border_data.data, 'borders')
        borders = border_data.data.borders;
    elseif isfield(border_data, 'borders')
        borders = border_data.borders;
    else
        warning('No borders variable found in: %s', border_file_path);
        borders = [];
        return;
    end

    borders = double(borders(:).');
    borders = borders(isfinite(borders));
    borders = unique(borders, 'stable');
    if isempty(borders)
        return;
    end
    if borders(1) ~= 1
        borders = [1, borders];
    end
    if borders(end) ~= N + 1
        borders = [borders, N + 1];
    end
end

function cat = classify_connections(J, err, err_multi)
    cat = nan(size(J));
    valid = isfinite(J) & isfinite(err);
    cat(valid) = 0;
    cat(valid & (J >  err_multi * err)) = 1;
    cat(valid & (J < -err_multi * err)) = -1;
end

function render_bidirectional_figure(root, pooled, valid_session_counts, area_pair, matrix_direction, ...
    conditions, kernel_indices, colors, marker_size, marker_alpha, figure_visible, show_legend, ...
    show_identity_line, show_fit_line, fit_line_method, show_histogram, J_lim, lim_ratio, output_suffix, ...
    plot_J_matrix_panel, output_prefix, export_pdf, panel_label_pos, scatter_panel_modes, render_count, render_total)

    [scatter_modes, scatter_titles] = build_scatter_panel_modes(scatter_panel_modes);
    panel_modes = [scatter_modes, {'category_3x3', 'sig_J_matrix'}];
    panel_titles = [scatter_titles, {'3x3 category', 'Significant-only J matrix'}];
    row_multiplier = 1 + double(show_histogram);
    n_row = numel(kernel_indices) * numel(conditions) * row_multiplier;
    n_col = numel(panel_modes);

    f = figure('Color', 'w', 'Visible', figure_visible);
    set(f, 'Visible', figure_visible);
    tiledlayout(n_row, n_col, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    [x_label, y_label] = reciprocal_axis_labels(area_pair, matrix_direction);
    shared_axis_limit = get_shared_symmetric_axis_limit(pooled);

    for kernel_i = 1:numel(kernel_indices)
        for condition_idx = 1:numel(conditions)
            base_row_idx = (kernel_i - 1) * numel(conditions) + condition_idx;
            scatter_row_idx = (base_row_idx - 1) * row_multiplier + 1;
            hist_row_idx = scatter_row_idx + 1;
            row_title = sprintf('K%d, %s, sessions=%d', ...
                kernel_indices(kernel_i), conditions(condition_idx).title, valid_session_counts(condition_idx, kernel_i));

            for panel_i = 1:n_col
                panel_index = (scatter_row_idx - 1) * n_col + panel_i;
                ax = nexttile(panel_index);
                title_text = sprintf('%s\n%s', row_title, panel_titles{panel_i});

                if strcmp(panel_modes{panel_i}, 'category_3x3')
                    plot_category_grid(ax, pooled{condition_idx, kernel_i}, x_label, y_label, title_text);
                elseif strcmp(panel_modes{panel_i}, 'sig_J_matrix')
                    axis(ax, 'off');
                else
                    plot_bidirectional_scatter(ax, pooled{condition_idx, kernel_i}, panel_modes{panel_i}, ...
                        colors, marker_size, marker_alpha, title_text, x_label, y_label, shared_axis_limit, ...
                        show_legend, show_identity_line, show_fit_line, fit_line_method);
                end
                if ~strcmp(panel_modes{panel_i}, 'sig_J_matrix')
                    add_panel_label(ax, panel_index, panel_label_pos);
                end

                if show_histogram
                    hist_panel_index = (hist_row_idx - 1) * n_col + panel_i;
                    hist_ax = nexttile(hist_panel_index);
                    if strcmp(panel_modes{panel_i}, 'category_3x3')
                        if plot_J_matrix_panel
                            full_J_lim = plot_full_J_matrix(hist_ax, pooled{condition_idx, kernel_i}.J_full, ...
                                pooled{condition_idx, kernel_i}.J_borders, J_lim, lim_ratio, ...
                                sprintf('%s\nFull N x N J matrix', row_title));
                        else
                            axis(hist_ax, 'off');
                        end
                    elseif strcmp(panel_modes{panel_i}, 'sig_J_matrix')
                        if plot_J_matrix_panel
                            full_J_lim = compute_J_color_limit(pooled{condition_idx, kernel_i}.J_full, J_lim, lim_ratio);
                            plot_full_J_matrix(hist_ax, pooled{condition_idx, kernel_i}.J_sig_full, ...
                                pooled{condition_idx, kernel_i}.J_borders, full_J_lim, 1, ...
                                sprintf('%s\nSignificant-only N x N J matrix', row_title));
                        else
                            axis(hist_ax, 'off');
                        end
                    else
                        plot_delta_histogram(hist_ax, pooled{condition_idx, kernel_i}, panel_modes{panel_i}, ...
                            colors, sprintf('%s\n%s x-y histogram', row_title, panel_titles{panel_i}), x_label, y_label);
                    end
                    add_panel_label(hist_ax, hist_panel_index, panel_label_pos);
                end
            end
        end
    end

    sgtitle(sprintf('%s <-> %s reciprocal J scatter', area_pair.label1, area_pair.label2), 'Interpreter', 'none');

    save_folder = fullfile(root, 'Figures', 'Paper', 'Analysis', 'bidirectional_scatter');
    check_path(save_folder);

    figWidth = 4.2 * n_col;
    figHeight = 3.8 * n_row;
    resolution = 300;

    set(f, 'Units', 'inches');
    f.Position(3:4) = [figWidth, figHeight];
    set(f, 'PaperUnits', 'inches');
    set(f, 'PaperSize', [figWidth, figHeight]);
    set(f, 'PaperPosition', [0, 0, figWidth, figHeight]);
    set(f, 'Color', 'w');

    output_stub = sprintf('%s_%s_%s%s', output_prefix, ...
        sanitize_filename(area_pair.label1), sanitize_filename(area_pair.label2), output_suffix);
    fprintf('Rendering figure %d/%d: %s\n', render_count, render_total, output_stub);
    preview_filename = fullfile(save_folder, [output_stub, '_preview.jpg']);
    exportgraphics(f, preview_filename, 'ContentType', 'image', 'BackgroundColor', 'white', 'Resolution', resolution);

    if export_pdf
        pdf_filename = fullfile(save_folder, [output_stub, '.pdf']);
        exportgraphics(f, pdf_filename, 'ContentType', 'vector', 'BackgroundColor', 'white', 'Resolution', resolution);
    end

    close(f);
end

function [x_label, y_label] = reciprocal_axis_labels(area_pair, matrix_direction)
    switch lower(matrix_direction)
        case {'row_to_col', 'col_to_row'}
            x_label = sprintf('%s to %s', area_pair.label1, area_pair.label2);
            y_label = sprintf('%s to %s', area_pair.label2, area_pair.label1);
        otherwise
            error('Unknown matrix_direction: %s.', matrix_direction);
    end
end

function [panel_modes, panel_titles] = build_scatter_panel_modes(scatter_panel_modes)
    panel_modes = cell(size(scatter_panel_modes));
    panel_titles = cell(size(scatter_panel_modes));
    for mode_i = 1:numel(scatter_panel_modes)
        panel_modes{mode_i} = normalize_scatter_mode(scatter_panel_modes{mode_i});
        panel_titles{mode_i} = scatter_mode_title(panel_modes{mode_i});
    end
end

function mode = normalize_scatter_mode(mode)
    key = lower(strtrim(char(string(mode))));
    key = regexprep(key, 'same[-_\s]+sign', 'same_sign');
    key = regexprep(key, '(diff|different)[-_\s]+sign', 'diff_sign');
    key = regexprep(key, 'either[-_\s]+(significant|sig)', 'either_sig');
    key = regexprep(key, 'both[-_\s]+(significant|sig)', 'both_sig');
    key = regexprep(key, 'all[-_\s]+connections?', 'all');
    key = regexprep(key, 'x[-_\s]+y[-_\s]+flip', 'xy_flip');
    key = regexprep(key, 'xy[-_\s]+flip', 'xy_flip');
    key = regexprep(key, 'x[-_\s]+flip', 'x_flip');
    key = regexprep(key, 'y[-_\s]+flip', 'y_flip');
    key = regexprep(key, '\s+', ' ');

    legacy_map = containers.Map( ...
        {'either', 'both', 'either abs', 'both abs', 'same abs', 'same_sign abs', 'all absolute'}, ...
        {'either_sig', 'both_sig', 'either_sig abs', 'both_sig abs', ...
         'both_sig same_sign abs', 'both_sig same_sign abs', 'abs'});
    if isKey(legacy_map, key)
        key = legacy_map(key);
    end

    tokens = strsplit(strtrim(key));
    valid_tokens = {'all', 'either_sig', 'both_sig', 'same_sign', 'diff_sign', ...
        'abs', 'x_flip', 'y_flip', 'xy_flip'};
    keep = true(size(tokens));
    for token_i = 1:numel(tokens)
        token = tokens{token_i};
        if isempty(token)
            keep(token_i) = false;
            continue;
        end
        if ~ismember(token, valid_tokens)
            error('Unknown scatter panel mode token "%s" in mode "%s".', token, char(string(mode)));
        end
        if strcmp(token, 'all') && numel(tokens) > 1
            keep(token_i) = false;
        end
    end
    tokens = tokens(keep);
    if isempty(tokens)
        mode = 'all';
    else
        mode = strjoin(tokens, ' ');
    end
end

function title_text = scatter_mode_title(mode)
    if strcmp(mode, 'all')
        title_text = 'All connections';
        return;
    end

    title_text = strrep(mode, '_', ' ');
    title_text = strrep(title_text, 'x flip', 'x-flip');
    title_text = strrep(title_text, 'y flip', 'y-flip');
    title_text = strrep(title_text, 'xy flip', 'xy-flip');
end

function plot_bidirectional_scatter(ax, data, plot_mode, colors, marker_size, marker_alpha, title_text, x_label, y_label, axis_limit, show_legend, show_identity_line, show_fit_line, fit_line_method)
    cla(ax);
    hold(ax, 'on');

    [plot_x, plot_y, plot_x_cat, plot_y_cat, axis_mode] = get_plot_values_for_mode(data, plot_mode);
    draw_mode_scatter(ax, plot_x, plot_y, plot_x_cat, plot_y_cat, plot_mode, colors, ...
        marker_size, marker_alpha, x_label, y_label);

    if strcmp(axis_mode, 'positive')
        axis_limit = [0, max(3.5, get_positive_axis_max(plot_x, plot_y))];
    elseif strcmp(axis_mode, 'auto_signed')
        axis_limit = get_symmetric_axis_limit(plot_x, plot_y);
        axis_mode = 'signed';
    end

    if show_identity_line
        plot(ax, axis_limit, axis_limit, '--', 'Color', colors.identity_line, 'LineWidth', 1, 'DisplayName', 'x=y');
    end
    if show_fit_line
        plot_linear_fit_line(ax, plot_x, plot_y, axis_limit, colors.fit_line, fit_line_method);
    end
    if strcmp(axis_mode, 'signed')
        xline(ax, 0, ':', 'Color', colors.zero_line, 'HandleVisibility', 'off');
        yline(ax, 0, ':', 'Color', colors.zero_line, 'HandleVisibility', 'off');
    end

    xlim(ax, axis_limit);
    ylim(ax, axis_limit);
    axis(ax, 'square');
    if strcmp(axis_mode, 'positive')
        xlabel(ax, sprintf('%s |J|', x_label), 'Interpreter', 'none');
        ylabel(ax, sprintf('%s |J|', y_label), 'Interpreter', 'none');
    else
        xlabel(ax, sprintf('%s J', x_label), 'Interpreter', 'none');
        ylabel(ax, sprintf('%s J', y_label), 'Interpreter', 'none');
    end

    [pearson_r, pearson_p, spearman_rho, spearman_p, n_valid] = correlation_stats(plot_x, plot_y);
    cos_sim = cosine_similarity_omitnan(plot_x, plot_y);
    title(ax, sprintf('%s\nPearson r = %.6f (p = %s, n = %d)\nSpearman rho = %.6f (p = %s)\ncos sim = %.6f', ...
        title_text, pearson_r, format_p_value(pearson_p), n_valid, spearman_rho, format_p_value(spearman_p), cos_sim), ...
        'Interpreter', 'none');

    hold(ax, 'off');
    if show_legend
        legend(ax, 'show', 'Location', 'best');
    else
        legend(ax, 'off');
    end
end

function draw_mode_scatter(ax, x, y, x_cat, y_cat, plot_mode, colors, marker_size, marker_alpha, x_label, y_label)
    x_sig = x_cat ~= 0;
    y_sig = y_cat ~= 0;

    if is_both_style_mode(plot_mode)
        pos_mask = x_cat == 1 & y_cat == 1;
        neg_mask = x_cat == -1 & y_cat == -1;
        switch_mask = x_sig & y_sig & ~(pos_mask | neg_mask);
        scatter(ax, x(pos_mask), y(pos_mask), marker_size, 'filled', ...
            'MarkerFaceColor', colors.pos, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
            'DisplayName', 'pos');
        scatter(ax, x(neg_mask), y(neg_mask), marker_size, 'filled', ...
            'MarkerFaceColor', colors.neg, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
            'DisplayName', 'neg');
        scatter(ax, x(switch_mask), y(switch_mask), marker_size, 'filled', ...
            'MarkerFaceColor', colors.switch, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
            'DisplayName', 'switch');
    else
        none_sig = ~x_sig & ~y_sig;
        x_only = x_sig & ~y_sig;
        y_only = ~x_sig & y_sig;
        both_sig = x_sig & y_sig;
        scatter(ax, x(none_sig), y(none_sig), marker_size, 'filled', ...
            'MarkerFaceColor', colors.non_sig, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
            'DisplayName', 'none sig');
        scatter(ax, x(x_only), y(x_only), marker_size, 'filled', ...
            'MarkerFaceColor', colors.x_only, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
            'DisplayName', sprintf('%s only', x_label));
        scatter(ax, x(y_only), y(y_only), marker_size, 'filled', ...
            'MarkerFaceColor', colors.y_only, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
            'DisplayName', sprintf('%s only', y_label));
        scatter(ax, x(both_sig), y(both_sig), marker_size, 'filled', ...
            'MarkerFaceColor', colors.both_sig, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
            'DisplayName', 'both sig');
    end
end

function tf = is_both_style_mode(plot_mode)
    tokens = mode_tokens(plot_mode);
    tf = any(ismember(tokens, {'both_sig', 'same_sign', 'diff_sign'}));
end

function plot_delta_histogram(ax, data, plot_mode, colors, title_text, x_label, y_label)
    cla(ax);
    hold(ax, 'on');

    [plot_x, plot_y, ~, ~, axis_mode] = get_plot_values_for_mode(data, plot_mode);
    delta = plot_x - plot_y;
    delta = delta(isfinite(delta));

    if isempty(delta)
        text(ax, 0.5, 0.5, 'No valid data', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        axis(ax, 'off');
        return;
    end

    histogram(ax, delta, 'BinMethod', 'fd', 'FaceColor', [0.45, 0.45, 0.45], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.80);
    xline(ax, 0, '-', 'Color', colors.zero_line, 'LineWidth', 1.0, 'DisplayName', 'x-y=0');

    delta_mean = mean(delta);
    delta_median = median(delta);
    p_mean = paired_ttest_pvalue(delta);
    [p_median, median_test_name] = signed_median_pvalue(delta);
    xline(ax, delta_mean, '--', 'Color', colors.fit_line, 'LineWidth', 1.2, 'DisplayName', 'mean');
    xline(ax, delta_median, ':', 'Color', colors.identity_line, 'LineWidth', 1.2, 'DisplayName', 'median');

    hold(ax, 'off');
    if strcmp(axis_mode, 'positive')
        xlabel(ax, sprintf('%s |J| - %s |J|', x_label, y_label), 'Interpreter', 'none');
    else
        xlabel(ax, sprintf('%s J - %s J', x_label, y_label), 'Interpreter', 'none');
    end
    ylabel(ax, 'Count');
    title(ax, sprintf('%s\nmean = %.4f, paired t p = %s %s\nmedian = %.4f, %s p = %s %s\nn = %d', ...
        title_text, delta_mean, format_p_value(p_mean), p_to_stars(p_mean), ...
        delta_median, median_test_name, format_p_value(p_median), p_to_stars(p_median), numel(delta)), ...
        'Interpreter', 'none');
end

function this_lim = plot_full_J_matrix(ax, J_full, J_borders, J_lim, lim_ratio, title_text)
    cla(ax);
    if isempty(J_full)
        text(ax, 0.5, 0.5, 'No J matrix', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        axis(ax, 'off');
        this_lim = NaN;
        return;
    end

    this_lim = compute_J_color_limit(J_full, J_lim, lim_ratio);

    imagesc(ax, J_full, [-this_lim, this_lim]);
    axis(ax, 'square');
    xlabel(ax, 'Column neuron index');
    ylabel(ax, 'Row neuron index');
    cb = colorbar(ax);
    ylabel(cb, 'J');

    try
        cmap = brewermap(256,'*RdBu');
        colormap(ax, cmap);
    catch
        colormap(ax, redblue_fallback_colormap(256));
    end

    draw_area_borders(ax, J_borders, size(J_full, 1));

    title(ax, sprintf('%s\ncolor range = [%.3g, %.3g]', title_text, -this_lim, this_lim), ...
        'Interpreter', 'none');
end

function this_lim = compute_J_color_limit(J_full, J_lim, lim_ratio)
    if ~isnan(J_lim)
        this_lim = J_lim;
        return;
    end

    vals = abs(J_full(isfinite(J_full)));
    vals = vals(vals > 0);
    if isempty(vals)
        this_lim = 1;
        return;
    end

    if nargin < 3 || isnan(lim_ratio)
        lim_ratio = 1;
    end
    lim_ratio = max(0, min(1, lim_ratio));
    this_lim = local_quantile(vals, lim_ratio);
    if ~isfinite(this_lim) || this_lim <= 0
        this_lim = max(vals);
    end
    if ~isfinite(this_lim) || this_lim <= 0
        this_lim = 1;
    end
end

function q = local_quantile(vals, ratio)
    vals = sort(vals(:));
    n = numel(vals);
    if n == 1
        q = vals(1);
        return;
    end
    idx = 1 + (n - 1) * ratio;
    lo = floor(idx);
    hi = ceil(idx);
    if lo == hi
        q = vals(lo);
    else
        q = vals(lo) + (idx - lo) * (vals(hi) - vals(lo));
    end
end

function draw_area_borders(ax, borders, N)
    if isempty(borders)
        return;
    end

    borders = double(borders(:).');
    borders = unique(borders(isfinite(borders)), 'stable');
    boundary_positions = borders(2:end) - 0.5;
    boundary_positions = boundary_positions(boundary_positions > 0.5 & boundary_positions < N + 0.5);
    if isempty(boundary_positions)
        return;
    end

    hold(ax, 'on');
    for boundary_idx = 1:numel(boundary_positions)
        b = boundary_positions(boundary_idx);
        line(ax, [0.5, N + 0.5], [b, b], 'Color', 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
        line(ax, [b, b], [0.5, N + 0.5], 'Color', 'k', 'LineWidth', 1.0, 'HandleVisibility', 'off');
    end
    hold(ax, 'off');
end

function cmap = redblue_fallback_colormap(n)
    if nargin < 1
        n = 256;
    end
    x = linspace(0, 1, n).';
    r = min(1, 2 * x);
    g = 1 - abs(2 * x - 1);
    b = min(1, 2 * (1 - x));
    cmap = [r, g, b];
end

function [plot_x, plot_y, plot_x_cat, plot_y_cat, axis_mode] = get_plot_values_for_mode(data, plot_mode)
    plot_x = data.x;
    plot_y = data.y;
    plot_x_cat = data.x_cat;
    plot_y_cat = data.y_cat;
    axis_mode = 'signed';

    tokens = mode_tokens(plot_mode);
    for token_i = 1:numel(tokens)
        token = tokens{token_i};
        switch token
            case 'all'
                % No filtering or transform.

            case 'either_sig'
                keep = plot_x_cat ~= 0 | plot_y_cat ~= 0;
                [plot_x, plot_y, plot_x_cat, plot_y_cat] = apply_connection_mask( ...
                    plot_x, plot_y, plot_x_cat, plot_y_cat, keep);

            case 'both_sig'
                keep = plot_x_cat ~= 0 & plot_y_cat ~= 0;
                [plot_x, plot_y, plot_x_cat, plot_y_cat] = apply_connection_mask( ...
                    plot_x, plot_y, plot_x_cat, plot_y_cat, keep);

            case 'same_sign'
                keep = (plot_x_cat == 1 & plot_y_cat == 1) | ...
                       (plot_x_cat == -1 & plot_y_cat == -1);
                [plot_x, plot_y, plot_x_cat, plot_y_cat] = apply_connection_mask( ...
                    plot_x, plot_y, plot_x_cat, plot_y_cat, keep);

            case 'diff_sign'
                keep = (plot_x_cat == 1 & plot_y_cat == -1) | ...
                       (plot_x_cat == -1 & plot_y_cat == 1);
                [plot_x, plot_y, plot_x_cat, plot_y_cat] = apply_connection_mask( ...
                    plot_x, plot_y, plot_x_cat, plot_y_cat, keep);

            case 'abs'
                plot_x = abs(plot_x);
                plot_y = abs(plot_y);
                axis_mode = 'positive';

            case 'x_flip'
                flip_mask = plot_x < 0;
                [plot_x, plot_y, plot_x_cat, plot_y_cat] = flip_connections( ...
                    plot_x, plot_y, plot_x_cat, plot_y_cat, flip_mask);
                if ~strcmp(axis_mode, 'positive')
                    axis_mode = 'auto_signed';
                end

            case 'y_flip'
                flip_mask = plot_y < 0;
                [plot_x, plot_y, plot_x_cat, plot_y_cat] = flip_connections( ...
                    plot_x, plot_y, plot_x_cat, plot_y_cat, flip_mask);
                if ~strcmp(axis_mode, 'positive')
                    axis_mode = 'auto_signed';
                end

            case 'xy_flip'
                flip_mask = (plot_x + plot_y) < 0;
                [plot_x, plot_y, plot_x_cat, plot_y_cat] = flip_connections( ...
                    plot_x, plot_y, plot_x_cat, plot_y_cat, flip_mask);
                if ~strcmp(axis_mode, 'positive')
                    axis_mode = 'auto_signed';
                end

            otherwise
                error('Unknown plot mode token "%s" in mode "%s".', token, plot_mode);
        end
    end
end

function tokens = mode_tokens(plot_mode)
    plot_mode = normalize_scatter_mode(plot_mode);
    tokens = strsplit(plot_mode);
end

function [x, y, x_cat, y_cat] = apply_connection_mask(x, y, x_cat, y_cat, keep)
    x = x(keep);
    y = y(keep);
    x_cat = x_cat(keep);
    y_cat = y_cat(keep);
end

function [x, y, x_cat, y_cat] = flip_connections(x, y, x_cat, y_cat, flip_mask)
    x(flip_mask) = -x(flip_mask);
    y(flip_mask) = -y(flip_mask);
    x_cat(flip_mask) = -x_cat(flip_mask);
    y_cat(flip_mask) = -y_cat(flip_mask);
end

function plot_category_grid(ax, data, x_label, y_label, title_text)
    cla(ax);

    class_values = [-1, 0, 1];
    category_labels = {'Negative', 'Non-sig', 'Positive'};
    valid = isfinite(data.x_cat) & isfinite(data.y_cat);
    x_cat = data.x_cat(valid);
    y_cat = data.y_cat(valid);

    counts = zeros(3, 3);
    for x_idx = 1:3
        for y_idx = 1:3
            counts(x_idx, y_idx) = sum(x_cat == class_values(x_idx) & y_cat == class_values(y_idx));
        end
    end

    plot_counts = counts.';
    imagesc(ax, plot_counts);
    set(ax, 'YDir', 'normal');
    axis(ax, 'square');

    xticks(ax, 1:3);
    yticks(ax, 1:3);
    xticklabels(ax, category_labels);
    yticklabels(ax, category_labels);
    xtickangle(ax, 30);
    xlabel(ax, sprintf('%s category', x_label), 'Interpreter', 'none');
    ylabel(ax, sprintf('%s category', y_label), 'Interpreter', 'none');

    cb = colorbar(ax);
    ylabel(cb, 'count');

    max_count = max(plot_counts(:));
    min_count = min(plot_counts(:));
    for y_idx = 1:3
        for x_idx = 1:3
            this_count = plot_counts(y_idx, x_idx);
            this_ratio = (this_count - min_count) / (max_count - min_count + eps);
            if this_ratio < 0.2
                text_color = 'w';
            else
                text_color = 'k';
            end
            text(ax, x_idx, y_idx, sprintf('%d', this_count), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', 'Color', text_color);
        end
    end

    agreement = compute_raw_agreement(counts);
    kappa = compute_cohen_kappa(counts);
    n_valid = sum(counts(:));
    title(ax, sprintf('%s\nAgreement = %s, kappa = %s, n = %d', ...
        title_text, format_stat_value(agreement, '%.3f'), format_stat_value(kappa, '%.3f'), n_valid), ...
        'Interpreter', 'none');
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

function vmax = get_positive_axis_max(x, y)
    vals = [x(:); y(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        vmax = 1;
    else
        vmax = 1.05 * max(vals);
        if vmax <= 0
            vmax = 1;
        end
    end
end

function axis_limit = get_symmetric_axis_limit(x, y)
    vals = [x(:); y(:)];
    vals = vals(isfinite(vals));
    if isempty(vals)
        vmax = 1;
    else
        vmax = 1.05 * max(abs(vals));
        if vmax <= 0
            vmax = 1;
        end
    end
    axis_limit = [-max(3.5, vmax), max(3.5, vmax)];
end

function axis_limit = get_shared_symmetric_axis_limit(pooled)
    vals = [];
    for k = 1:numel(pooled)
        vals = [vals; pooled{k}.x(:); pooled{k}.y(:)]; %#ok<AGROW>
    end
    vals = vals(isfinite(vals));
    if isempty(vals)
        vmax = 1;
    else
        vmax = 1.05 * max(abs(vals));
        if vmax <= 0
            vmax = 1;
        end
    end
    axis_limit = [-max(3.5, vmax), max(3.5, vmax)];
end

function [pearson_r, pearson_p, spearman_rho, spearman_p, n_valid] = correlation_stats(x, y)
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    n_valid = numel(x);

    pearson_r = NaN;
    pearson_p = NaN;
    spearman_rho = NaN;
    spearman_p = NaN;

    if n_valid < 2 || all(x == x(1)) || all(y == y(1))
        return;
    end

    [pearson_r, pearson_p] = corr(x, y, 'Type', 'Pearson', 'Rows', 'complete');
    [spearman_rho, spearman_p] = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');
end

function cos_sim = cosine_similarity_omitnan(x, y)
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    if isempty(x) || norm(x) == 0 || norm(y) == 0
        cos_sim = NaN;
    else
        cos_sim = dot(x, y) / (norm(x) * norm(y));
    end
end

function plot_linear_fit_line(ax, x, y, axis_limit, line_color, fit_line_method)
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if numel(x) < 2
        return;
    end

    switch lower(fit_line_method)
        case {'ols', 'linear'}
            if all(x == x(1))
                return;
            end
            fit_par = polyfit(x, y, 1);
            x_fit = linspace(axis_limit(1), axis_limit(2), 100);
            y_fit = polyval(fit_par, x_fit);

        case {'tls', 'orthogonal', 'pca'}
            [x_fit, y_fit] = total_least_squares_line(x, y, axis_limit);
            if isempty(x_fit)
                return;
            end

        otherwise
            error('Unknown fit_line_method: %s. Use ols or tls.', fit_line_method);
    end

    plot(ax, x_fit, y_fit, '-', 'Color', line_color, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('fit (%s)', fit_line_method));
end

function [x_fit, y_fit] = total_least_squares_line(x, y, axis_limit)
    x_fit = [];
    y_fit = [];

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if numel(x) < 2
        return;
    end

    center = [mean(x), mean(y)];
    XY = [x - center(1), y - center(2)];
    if all(abs(XY(:)) < eps)
        return;
    end

    [~, ~, V] = svd(XY, 0);
    direction = V(:, 1);
    if any(~isfinite(direction)) || norm(direction) == 0
        return;
    end
    direction = direction / norm(direction);

    candidates = [];
    bounds = axis_limit;
    if abs(direction(1)) > eps
        t = (bounds(1) - center(1)) / direction(1);
        y_hit = center(2) + t * direction(2);
        if y_hit >= bounds(1) - eps && y_hit <= bounds(2) + eps
            candidates(end+1, :) = [bounds(1), y_hit]; %#ok<AGROW>
        end
        t = (bounds(2) - center(1)) / direction(1);
        y_hit = center(2) + t * direction(2);
        if y_hit >= bounds(1) - eps && y_hit <= bounds(2) + eps
            candidates(end+1, :) = [bounds(2), y_hit]; %#ok<AGROW>
        end
    end
    if abs(direction(2)) > eps
        t = (bounds(1) - center(2)) / direction(2);
        x_hit = center(1) + t * direction(1);
        if x_hit >= bounds(1) - eps && x_hit <= bounds(2) + eps
            candidates(end+1, :) = [x_hit, bounds(1)]; %#ok<AGROW>
        end
        t = (bounds(2) - center(2)) / direction(2);
        x_hit = center(1) + t * direction(1);
        if x_hit >= bounds(1) - eps && x_hit <= bounds(2) + eps
            candidates(end+1, :) = [x_hit, bounds(2)]; %#ok<AGROW>
        end
    end

    if size(candidates, 1) < 2
        return;
    end
    [~, unique_idx] = unique(round(candidates, 12), 'rows', 'stable');
    candidates = candidates(unique_idx, :);
    if size(candidates, 1) < 2
        return;
    end

    dx = candidates(:, 1) - candidates(:, 1).';
    dy = candidates(:, 2) - candidates(:, 2).';
    D = dx.^2 + dy.^2;
    [~, max_idx] = max(D(:));
    [i1, i2] = ind2sub(size(D), max_idx);
    x_fit = [candidates(i1, 1), candidates(i2, 1)];
    y_fit = [candidates(i1, 2), candidates(i2, 2)];
end

function p_str = format_p_value(p_value)
    if isnan(p_value)
        p_str = 'NaN';
    elseif p_value > 0.001
        p_str = sprintf('%.3f', p_value);
    else
        p_str = sprintf('%.3e', p_value);
    end
end

function stars = p_to_stars(p_value)
    if isnan(p_value)
        stars = '';
    elseif p_value < 0.001
        stars = '***';
    elseif p_value < 0.01
        stars = '**';
    elseif p_value < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end
end

function p_value = paired_ttest_pvalue(delta)
    delta = delta(isfinite(delta));
    if numel(delta) < 2 || all(delta == delta(1))
        p_value = NaN;
        return;
    end

    try
        [~, p_value] = ttest(delta, 0);
    catch
        n = numel(delta);
        t_stat = mean(delta) / (std(delta) / sqrt(n));
        p_value = 2 * tcdf(-abs(t_stat), n - 1);
    end
end

function [p_value, test_name] = signed_median_pvalue(delta)
    test_name = 'signrank';
    delta = delta(isfinite(delta));
    delta = delta(delta ~= 0);
    if isempty(delta)
        p_value = NaN;
        return;
    end

    try
        p_value = signrank(delta, 0);
    catch
        test_name = 'sign test';
        n_pos = sum(delta > 0);
        n = numel(delta);
        p_value = two_sided_sign_test_pvalue(n_pos, n);
    end
end

function p_value = two_sided_sign_test_pvalue(n_pos, n)
    if n == 0
        p_value = NaN;
        return;
    end

    k = min(n_pos, n - n_pos);
    probs = zeros(k + 1, 1);
    for i = 0:k
        probs(i + 1) = exp(gammaln(n + 1) - gammaln(i + 1) - gammaln(n - i + 1) - n * log(2));
    end
    p_value = min(1, 2 * sum(probs));
end

function stat_str = format_stat_value(value, fmt)
    if isnan(value)
        stat_str = 'NaN';
    else
        stat_str = sprintf(fmt, value);
    end
end

function add_panel_label(ax, panel_index, panel_label_pos)
    panel_label = panel_index_to_letters(panel_index);
    text(ax, panel_label_pos(1), panel_label_pos(2), panel_label, 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'Interpreter', 'none', 'Clipping', 'off');
end

function label = panel_index_to_letters(panel_index)
    label = '';
    n = panel_index;
    while n > 0
        rem0 = mod(n - 1, 26);
        label = [char(double('A') + rem0), label]; %#ok<AGROW>
        n = floor((n - 1) / 26);
    end
end

function name = sanitize_filename(name)
    name = regexprep(char(string(name)), '[^A-Za-z0-9_-]+', '_');
    name = regexprep(name, '_+', '_');
    name = regexprep(name, '^_|_$', '');
end

function name = join_area_names(area_names)
    name = strjoin(cellstr(string(area_names)), '+');
end
