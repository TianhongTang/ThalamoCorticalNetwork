%% Supplement Figure: state similarity matrices
% For each kernel, generate one figure with 4 panels:
%   1) 2x2 categorical kappa
%   2) 3x3 categorical kappa
%   3) Pearson r on all J values
%   4) Spearman rho on all J values
%
% Each panel is a 4x4 similarity matrix across:
%   Pre-Open, Pre-Close, Post-Open, Post-Close

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

%% Parameters
kernel_indices = [1, 2, 3];
err_multi = 1;   % same threshold rule as your current script

figure_visible = 'off';
fig_width = 14;   % inches
fig_height = 12;  % inches
resolution = 300;

state_defs = [ ...
    struct('prepost', 'Pre',  'state', 'RestOpen',  'label', 'Pre-Open'),  ...
    struct('prepost', 'Pre',  'state', 'RestClose', 'label', 'Pre-Close'), ...
    struct('prepost', 'Post', 'state', 'RestOpen',  'label', 'Post-Open'), ...
    struct('prepost', 'Post', 'state', 'RestClose', 'label', 'Post-Close') ...
];
state_labels = {state_defs.label};

skip_failed_sessions = false;
max_sessions_to_include = inf;  % reduce if debugging

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

fprintf('Selected %d sessions.\n', numel(meta_array));

%% Run one figure per kernel
for kernel_idx = kernel_indices
    fprintf('\n==============================\n');
    fprintf('Processing kernel %d\n', kernel_idx);
    fprintf('==============================\n');

    [sim, valid_session_count, failed_session_count] = ...
        compute_similarity_for_kernel(root, meta_array, state_defs, kernel_idx, err_multi, skip_failed_sessions);

    if valid_session_count == 0
        warning('No valid sessions for kernel %d. Skipping figure.', kernel_idx);
        continue;
    end

    fprintf('Kernel %d: valid sessions = %d, failed sessions = %d\n', ...
        kernel_idx, valid_session_count, failed_session_count);

    render_similarity_figure(root, kernel_idx, sim, state_labels, valid_session_count, ...
        figure_visible, fig_width, fig_height, resolution);
end

%% ========================= LOCAL FUNCTIONS ==============================

function [sim, valid_session_count, failed_session_count] = ...
    compute_similarity_for_kernel(root, meta_array, state_defs, kernel_idx, err_multi, skip_failed_sessions)

    n_state = numel(state_defs);

    pair_pool = repmat(struct( ...
        'counts3', zeros(3, 3), ...
        'counts2', zeros(2, 2), ...
        'x', [], ...
        'y', []), n_state, n_state);

    valid_session_count = 0;
    failed_session_count = 0;

    for session_i = 1:numel(meta_array)
        meta = meta_array(session_i);
        session_label = make_session_label(meta);

        fprintf('\nLoading %s (%d/%d), kernel %d\n', ...
            session_label, session_i, numel(meta_array), kernel_idx);

        try
            loaded_states = struct();

            % Load all 4 states for this session and this kernel
            for s = 1:n_state
                key = state_key(state_defs(s).prepost, state_defs(s).state, kernel_idx);
                loaded_states.(key) = load_state_connectivity( ...
                    root, meta, state_defs(s).prepost, state_defs(s).state, kernel_idx);
            end

            % Pool pairwise similarities across states
            for i = 1:n_state
                for j = i:n_state
                    key_i = state_key(state_defs(i).prepost, state_defs(i).state, kernel_idx);
                    key_j = state_key(state_defs(j).prepost, state_defs(j).state, kernel_idx);

                    state_i = loaded_states.(key_i);
                    state_j = loaded_states.(key_j);

                    % 3x3 categorical counts
                    [counts3, ~, ~, ~] = make_pair_category_counts(state_i, state_j, err_multi);

                    % 2x2 categorical counts: drop Non-sig, keep only Neg/Pos
                    counts2 = counts3([1, 3], [1, 3]);

                    % Raw J values
                    [x_vals, y_vals] = extract_numeric_pair(state_i, state_j);

                    pair_pool(i, j).counts3 = pair_pool(i, j).counts3 + counts3;
                    pair_pool(i, j).counts2 = pair_pool(i, j).counts2 + counts2;
                    pair_pool(i, j).x = [pair_pool(i, j).x; x_vals];
                    pair_pool(i, j).y = [pair_pool(i, j).y; y_vals];
                end
            end

            valid_session_count = valid_session_count + 1;

        catch ME
            failed_session_count = failed_session_count + 1;
            if skip_failed_sessions
                warning('Skipping %s because processing failed: %s', session_label, ME.message);
                continue;
            else
                rethrow(ME);
            end
        end
    end

    % Convert pooled pair data into symmetric similarity matrices
    sim.kappa2 = nan(n_state, n_state);
    sim.kappa3 = nan(n_state, n_state);
    sim.pearson_r = nan(n_state, n_state);
    sim.spearman_rho = nan(n_state, n_state);

    sim.n2 = zeros(n_state, n_state);
    sim.n3 = zeros(n_state, n_state);
    sim.nr = zeros(n_state, n_state);
    sim.nrho = zeros(n_state, n_state);

    for i = 1:n_state
        for j = i:n_state
            counts2 = pair_pool(i, j).counts2;
            counts3 = pair_pool(i, j).counts3;
            x = pair_pool(i, j).x;
            y = pair_pool(i, j).y;

            kappa2 = compute_cohen_kappa(counts2);
            kappa3 = compute_cohen_kappa(counts3);

            n2 = sum(counts2(:));
            n3 = sum(counts3(:));

            [r, ~, nr] = correlation_stat(x, y, 'Pearson');
            [rho, ~, nrho] = correlation_stat(x, y, 'Spearman');

            sim.kappa2(i, j) = kappa2;
            sim.kappa2(j, i) = kappa2;

            sim.kappa3(i, j) = kappa3;
            sim.kappa3(j, i) = kappa3;

            sim.pearson_r(i, j) = r;
            sim.pearson_r(j, i) = r;

            sim.spearman_rho(i, j) = rho;
            sim.spearman_rho(j, i) = rho;

            sim.n2(i, j) = n2;   sim.n2(j, i) = n2;
            sim.n3(i, j) = n3;   sim.n3(j, i) = n3;
            sim.nr(i, j) = nr;   sim.nr(j, i) = nr;
            sim.nrho(i, j) = nrho; sim.nrho(j, i) = nrho;
        end
    end
end

function render_similarity_figure(root, kernel_idx, sim, state_labels, valid_session_count, ...
    figure_visible, fig_width, fig_height, resolution)

    f = figure('Color', 'w', 'Visible', figure_visible);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(1);
    plot_similarity_heatmap(ax, sim.kappa2, sim.n2, state_labels, ...
        '2x2 categorical kappa (Neg / Pos only)');
    
    ax = nexttile(2);
    plot_similarity_heatmap(ax, sim.kappa3, sim.n3, state_labels, ...
        '3x3 categorical kappa (Neg / Non-sig / Pos)');
    
    ax = nexttile(3);
    plot_similarity_heatmap(ax, sim.pearson_r, sim.nr, state_labels, ...
        'Pearson correlation r on all J');
    
    ax = nexttile(4);
    plot_similarity_heatmap(ax, sim.spearman_rho, sim.nrho, state_labels, ...
        'Spearman correlation \rho on all J');

    sgtitle(sprintf('State similarity summary, Kernel %d (pooled across %d sessions)', ...
        kernel_idx, valid_session_count), 'Interpreter', 'none');

    save_folder = fullfile(root, 'Figures', 'Paper');
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end

    output_stub = sprintf('Figure_state_similarity_k%d', kernel_idx);

    set(f, 'Units', 'inches');
    f.Position(3:4) = [fig_width, fig_height];
    set(f, 'PaperUnits', 'inches');
    set(f, 'PaperSize', [fig_width, fig_height]);
    set(f, 'PaperPosition', [0, 0, fig_width, fig_height]);
    set(f, 'Color', 'w');

    preview_filename = fullfile(save_folder, [output_stub, '_preview.jpg']);
    exportgraphics(f, preview_filename, 'ContentType', 'image', ...
        'BackgroundColor', 'white', 'Resolution', resolution);

    pdf_filename = fullfile(save_folder, [output_stub, '.pdf']);
    exportgraphics(f, pdf_filename, 'ContentType', 'vector', ...
        'BackgroundColor', 'white', 'Resolution', resolution);

    close(f);
end

function plot_similarity_heatmap(ax, mat, nmat, state_labels, title_text)

    n_state = size(mat, 1);

    % Color limits are determined only from finite off-diagonal values.
    off_diag_mask = ~eye(n_state);
    off_diag_vals = mat(off_diag_mask);
    off_diag_vals = off_diag_vals(isfinite(off_diag_vals));

    if isempty(off_diag_vals)
        color_limits = [-1, 1];
    else
        color_min = min(off_diag_vals);
        color_max = max(off_diag_vals);

        % Avoid invalid caxis when all off-diagonal values are identical.
        if color_min == color_max
            padding = max(0.05, 0.05 * abs(color_min));

            if color_min == 0
                padding = 0.05;
            end

            color_limits = [color_min - padding, color_max + padding];
        else
            color_limits = [color_min, color_max];
        end
    end

    % The diagonal values do not affect the colorbar because caxis is
    % explicitly calculated from off-diagonal entries.
    imagesc(ax, mat);
    caxis(ax, color_limits);

    axis(ax, 'square');
    set(ax, 'YDir', 'normal');

    xticks(ax, 1:n_state);
    yticks(ax, 1:n_state);
    xticklabels(ax, state_labels);
    yticklabels(ax, state_labels);
    xtickangle(ax, 30);

    cb = colorbar(ax);
    ylabel(cb, 'similarity');

    title(ax, title_text, 'Interpreter', 'none');

    hold(ax, 'on');

    for row_idx = 1:n_state
        for col_idx = 1:n_state

            if row_idx == col_idx
                % Cover diagonal cells with black squares.
                rectangle(ax, ...
                    'Position', [col_idx - 0.5, row_idx - 0.5, 1, 1], ...
                    'FaceColor', 'k', ...
                    'EdgeColor', 'none');

                text(ax, col_idx, row_idx, '1', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontWeight', 'bold', ...
                    'FontSize', 10, ...
                    'Color', 'w');

                continue;
            end

            value = mat(row_idx, col_idx);
            n_value = nmat(row_idx, col_idx);

            if ~isfinite(value)
                value_text = sprintf('NaN\nn=%d', n_value);
                text_color = 'k';
            else
                value_text = sprintf('%.3f\nn=%d', value, n_value);

                % Choose text color according to its position in the
                % displayed color range.
                color_ratio = ...
                    (value - color_limits(1)) / ...
                    (color_limits(2) - color_limits(1));

                if color_ratio < 0.30
                    text_color = 'w';
                else
                    text_color = 'k';
                end
            end

            text(ax, col_idx, row_idx, value_text, ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', ...
                'FontSize', 9, ...
                'Color', text_color);
        end
    end

    hold(ax, 'off');
end

function [counts, agreement, kappa, n_valid] = make_pair_category_counts(state_x, state_y, err_multi)
    validate_matching_filters(state_x, state_y);

    x_cat12 = classify_connections(state_x.J12(:), state_x.err12(:), err_multi);
    y_cat12 = classify_connections(state_y.J12(:), state_y.err12(:), err_multi);
    x_cat21 = classify_connections(state_x.J21(:), state_x.err21(:), err_multi);
    y_cat21 = classify_connections(state_y.J21(:), state_y.err21(:), err_multi);

    x_cat = [x_cat12; x_cat21];
    y_cat = [y_cat12; y_cat21];

    valid = isfinite(x_cat) & isfinite(y_cat);
    x_cat = x_cat(valid);
    y_cat = y_cat(valid);
    n_valid = numel(x_cat);

    class_values = [-1, 0, 1];
    counts = zeros(3, 3);
    for i_class = 1:3
        for j_class = 1:3
            counts(i_class, j_class) = sum(x_cat == class_values(i_class) & y_cat == class_values(j_class));
        end
    end

    agreement = compute_raw_agreement(counts);
    kappa = compute_cohen_kappa(counts);
end

function cat = classify_connections(J, err, err_multi)
    cat = nan(size(J));
    valid = isfinite(J) & isfinite(err);
    cat(valid) = 0;
    cat(valid & (J >  err_multi * err)) = 1;
    cat(valid & (J < -err_multi * err)) = -1;
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
    col_marginal = sum(counts, 1)' / total_n;
    p_expected = sum(row_marginal .* col_marginal);

    if abs(1 - p_expected) < eps
        kappa = NaN;
    else
        kappa = (p_observed - p_expected) / (1 - p_expected);
    end
end

function [val, pval, n_valid] = correlation_stat(x, y, corr_type)
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    n_valid = numel(x);

    if n_valid < 2
        val = NaN;
        pval = NaN;
        return;
    end

    if all(x == x(1)) || all(y == y(1))
        val = NaN;
        pval = NaN;
        return;
    end

    [R, P] = corr(x, y, 'Type', corr_type, 'Rows', 'complete');
    val = R;
    pval = P;
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
                anchor_name = mt.file_name{idx};
                warning('Session %s is missing required state: %s %s. Skipping this session.', ...
                    anchor_name, required_preposts{q}, required_states{q});
                break;
            end
        end

        if has_all_required
            selected_rows(idx) = true;
        end
    end

    fprintf('Default metadata filter selected %d/%d rows.\n', sum(selected_rows), sum(base_filter));

    if ~any(selected_rows)
        warning('Default metadata filter selected no complete 4-state sets. Falling back to Pre RestOpen anchors.');
        selected_rows = anchor_filter;
    end
end

function label = make_session_label(meta)
    if ischar(meta.animal_name) || isstring(meta.animal_name)
        animal_str = char(meta.animal_name);
    else
        animal_str = char(string(meta.animal_name));
    end

    if ischar(meta.injection) || isstring(meta.injection)
        injection_str = char(meta.injection);
    else
        injection_str = char(string(meta.injection));
    end

    label = sprintf('%s %s session %s', animal_str, injection_str, char(string(meta.session_idx)));
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