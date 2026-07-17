%% Supplement Figure: State similarity matrices
% For each kernel, generate one 2x2-panel figure.
%
% Panels:
%   1. 2x2 categorical Cohen's kappa
%      Negative/Positive only; connections involving Non-sig are excluded.
%   2. 3x3 categorical Cohen's kappa
%      Negative/Non-sig/Positive.
%   3. Pearson correlation r using all raw J values.
%   4. Spearman correlation rho using all raw J values.
%
% Each panel is a 4x4 similarity matrix comparing:
%   Pre-Open, Pre-Close, Post-Open, Post-Close.
%
% Data are pooled across all selected sessions before each similarity metric
% is calculated. ACC->VLPFC and VLPFC->ACC connections are concatenated.
%
% The matrix diagonal is displayed as a black cell with a white "1".
% Diagonal entries are excluded when determining each panel's color limits.

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
kernel_indices = [1, 2];

% A connection is:
%   Positive when J >  err_multi * error
%   Negative when J < -err_multi * error
%   Non-sig otherwise
err_multi = 1;

figure_visible = 'off';
fig_width = 14;       % inches
fig_height = 12;      % inches
resolution = 300;     % dpi

skip_failed_sessions = false;
max_sessions_to_include = inf;  % Set smaller while debugging.

state_defs = [ ...
    struct('prepost', 'Pre',  'state', 'RestOpen',  'label', 'Pre-Open'), ...
    struct('prepost', 'Pre',  'state', 'RestClose', 'label', 'Pre-Close'), ...
    struct('prepost', 'Post', 'state', 'RestOpen',  'label', 'Post-Open'), ...
    struct('prepost', 'Post', 'state', 'RestClose', 'label', 'Post-Close') ...
];

state_labels = {state_defs.label};

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

%% Calculate and plot one figure per kernel
for kernel_idx = kernel_indices
    fprintf('\n==================================================\n');
    fprintf('Processing kernel %d\n', kernel_idx);
    fprintf('==================================================\n');

    [sim, valid_session_count, failed_session_count] = ...
        compute_similarity_for_kernel( ...
            root, meta_array, state_defs, kernel_idx, err_multi, ...
            skip_failed_sessions);

    if valid_session_count == 0
        warning('No valid sessions for kernel %d. Figure was not generated.', ...
            kernel_idx);
        continue;
    end

    fprintf(['Kernel %d complete: valid sessions = %d, ' ...
             'failed sessions = %d.\n'], ...
        kernel_idx, valid_session_count, failed_session_count);

    render_similarity_figure( ...
        root, kernel_idx, sim, state_labels, valid_session_count, ...
        figure_visible, fig_width, fig_height, resolution);
end


%% ========================================================================
% Local functions
% ========================================================================

function [sim, valid_session_count, failed_session_count] = ...
    compute_similarity_for_kernel( ...
        root, meta_array, state_defs, kernel_idx, err_multi, ...
        skip_failed_sessions)

    n_state = numel(state_defs);

    empty_pair = struct( ...
        'counts3', zeros(3, 3), ...
        'counts2', zeros(2, 2), ...
        'x', [], ...
        'y', []);

    pair_pool = repmat(empty_pair, n_state, n_state);

    valid_session_count = 0;
    failed_session_count = 0;

    for session_i = 1:numel(meta_array)
        meta = meta_array(session_i);
        session_label = make_session_label(meta);

        fprintf('\nLoading %s (%d/%d), kernel %d\n', ...
            session_label, session_i, numel(meta_array), kernel_idx);

        try
            loaded_states = struct();

            % Load the four states for this session and kernel.
            for state_i = 1:n_state
                state_def = state_defs(state_i);
                key = state_key( ...
                    state_def.prepost, state_def.state, kernel_idx);

                loaded_states.(key) = load_state_connectivity( ...
                    root, meta, state_def.prepost, ...
                    state_def.state, kernel_idx);
            end

            % Only calculate unique off-diagonal state pairs.
            for row_idx = 1:(n_state - 1)
                for col_idx = (row_idx + 1):n_state
                    row_def = state_defs(row_idx);
                    col_def = state_defs(col_idx);

                    row_key = state_key( ...
                        row_def.prepost, row_def.state, kernel_idx);
                    col_key = state_key( ...
                        col_def.prepost, col_def.state, kernel_idx);

                    state_x = loaded_states.(row_key);
                    state_y = loaded_states.(col_key);

                    % 3x3 categorical counts:
                    % rows = state_x category
                    % columns = state_y category
                    counts3 = make_pair_category_counts( ...
                        state_x, state_y, err_multi);

                    % 2x2 categorical counts:
                    % remove the Non-sig row and column.
                    counts2 = counts3([1, 3], [1, 3]);

                    % Corresponding raw J values from both directions.
                    [x_values, y_values] = extract_numeric_pair( ...
                        state_x, state_y);

                    pair_pool(row_idx, col_idx).counts3 = ...
                        pair_pool(row_idx, col_idx).counts3 + counts3;

                    pair_pool(row_idx, col_idx).counts2 = ...
                        pair_pool(row_idx, col_idx).counts2 + counts2;

                    pair_pool(row_idx, col_idx).x = [ ...
                        pair_pool(row_idx, col_idx).x; x_values];

                    pair_pool(row_idx, col_idx).y = [ ...
                        pair_pool(row_idx, col_idx).y; y_values];
                end
            end

            valid_session_count = valid_session_count + 1;

        catch ME
            failed_session_count = failed_session_count + 1;

            if skip_failed_sessions
                warning('Skipping %s because processing failed: %s', ...
                    session_label, ME.message);
                continue;
            else
                rethrow(ME);
            end
        end
    end

    % Initialize similarity matrices.
    sim.kappa2 = nan(n_state, n_state);
    sim.kappa3 = nan(n_state, n_state);
    sim.pearson_r = nan(n_state, n_state);
    sim.spearman_rho = nan(n_state, n_state);

    % Sample-size matrices.
    sim.n2 = zeros(n_state, n_state);
    sim.n3 = zeros(n_state, n_state);
    sim.nr = zeros(n_state, n_state);
    sim.nrho = zeros(n_state, n_state);

    % The diagonal is defined as perfect self-similarity.
    sim.kappa2(1:(n_state + 1):end) = 1;
    sim.kappa3(1:(n_state + 1):end) = 1;
    sim.pearson_r(1:(n_state + 1):end) = 1;
    sim.spearman_rho(1:(n_state + 1):end) = 1;

    % Calculate each unique off-diagonal pair and mirror it.
    for row_idx = 1:(n_state - 1)
        for col_idx = (row_idx + 1):n_state
            counts2 = pair_pool(row_idx, col_idx).counts2;
            counts3 = pair_pool(row_idx, col_idx).counts3;
            x_values = pair_pool(row_idx, col_idx).x;
            y_values = pair_pool(row_idx, col_idx).y;

            kappa2 = compute_cohen_kappa(counts2);
            kappa3 = compute_cohen_kappa(counts3);

            [pearson_r, ~, n_pearson] = ...
                correlation_stat(x_values, y_values, 'Pearson');

            [spearman_rho, ~, n_spearman] = ...
                correlation_stat(x_values, y_values, 'Spearman');

            n2 = sum(counts2(:));
            n3 = sum(counts3(:));

            sim.kappa2(row_idx, col_idx) = kappa2;
            sim.kappa2(col_idx, row_idx) = kappa2;

            sim.kappa3(row_idx, col_idx) = kappa3;
            sim.kappa3(col_idx, row_idx) = kappa3;

            sim.pearson_r(row_idx, col_idx) = pearson_r;
            sim.pearson_r(col_idx, row_idx) = pearson_r;

            sim.spearman_rho(row_idx, col_idx) = spearman_rho;
            sim.spearman_rho(col_idx, row_idx) = spearman_rho;

            sim.n2(row_idx, col_idx) = n2;
            sim.n2(col_idx, row_idx) = n2;

            sim.n3(row_idx, col_idx) = n3;
            sim.n3(col_idx, row_idx) = n3;

            sim.nr(row_idx, col_idx) = n_pearson;
            sim.nr(col_idx, row_idx) = n_pearson;

            sim.nrho(row_idx, col_idx) = n_spearman;
            sim.nrho(col_idx, row_idx) = n_spearman;
        end
    end
end


function render_similarity_figure( ...
    root, kernel_idx, sim, state_labels, valid_session_count, ...
    figure_visible, fig_width, fig_height, resolution)

    f = figure('Color', 'w', 'Visible', figure_visible);

    tiledlayout(2, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    ax = nexttile(1);
    plot_similarity_heatmap( ...
        ax, sim.kappa2, sim.n2, state_labels, ...
        '2x2 categorical Cohen''s kappa', ...
        'Cohen''s \kappa');

    ax = nexttile(2);
    plot_similarity_heatmap( ...
        ax, sim.kappa3, sim.n3, state_labels, ...
        '3x3 categorical Cohen''s kappa', ...
        'Cohen''s \kappa');

    ax = nexttile(3);
    plot_similarity_heatmap( ...
        ax, sim.pearson_r, sim.nr, state_labels, ...
        'Pearson correlation using all J', ...
        'Pearson r');

    ax = nexttile(4);
    plot_similarity_heatmap( ...
        ax, sim.spearman_rho, sim.nrho, state_labels, ...
        'Spearman correlation using all J', ...
        'Spearman \rho');

    sgtitle(sprintf( ...
        'State similarity, Kernel %d, pooled across %d sessions', ...
        kernel_idx, valid_session_count), ...
        'Interpreter', 'none');

    save_folder = fullfile(root, 'Figures', 'Paper');

    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end

    output_stub = sprintf( ...
        'Supplement_Figure_state_similarity_k%d', kernel_idx);

    set(f, 'Units', 'inches');
    f.Position(3:4) = [fig_width, fig_height];

    set(f, 'PaperUnits', 'inches');
    set(f, 'PaperSize', [fig_width, fig_height]);
    set(f, 'PaperPosition', [0, 0, fig_width, fig_height]);
    set(f, 'Color', 'w');

    jpg_filename = fullfile( ...
        save_folder, [output_stub, '_preview.jpg']);

    exportgraphics( ...
        f, jpg_filename, ...
        'ContentType', 'image', ...
        'BackgroundColor', 'white', ...
        'Resolution', resolution);

    pdf_filename = fullfile( ...
        save_folder, [output_stub, '.pdf']);

    exportgraphics( ...
        f, pdf_filename, ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white', ...
        'Resolution', resolution);

    fprintf('Saved:\n  %s\n  %s\n', jpg_filename, pdf_filename);

    close(f);
end


function plot_similarity_heatmap( ...
    ax, similarity_matrix, n_matrix, state_labels, ...
    title_text, colorbar_label)

    n_state = size(similarity_matrix, 1);

    if size(similarity_matrix, 2) ~= n_state
        error('Similarity matrix must be square.');
    end

    if ~isequal(size(n_matrix), size(similarity_matrix))
        error('Sample-size matrix must match the similarity matrix.');
    end

    % Determine color limits from finite off-diagonal values only.
    off_diagonal_mask = ~eye(n_state);
    off_diagonal_values = similarity_matrix(off_diagonal_mask);
    off_diagonal_values = ...
        off_diagonal_values(isfinite(off_diagonal_values));

    color_limits = calculate_off_diagonal_color_limits( ...
        off_diagonal_values, [-1, 1]);

    image_handle = imagesc(ax, similarity_matrix);
    set(image_handle, 'AlphaData', isfinite(similarity_matrix));

    caxis(ax, color_limits);
    set(ax, 'Color', [0.85, 0.85, 0.85]);
    set(ax, 'YDir', 'normal');

    axis(ax, 'square');
    xlim(ax, [0.5, n_state + 0.5]);
    ylim(ax, [0.5, n_state + 0.5]);

    xticks(ax, 1:n_state);
    yticks(ax, 1:n_state);
    xticklabels(ax, state_labels);
    yticklabels(ax, state_labels);
    xtickangle(ax, 30);

    title(ax, title_text, 'Interpreter', 'none');

    cb = colorbar(ax);
    ylabel(cb, colorbar_label, 'Interpreter', 'tex');

    hold(ax, 'on');

    for row_idx = 1:n_state
        for col_idx = 1:n_state
            if row_idx == col_idx
                % The diagonal is fixed self-similarity and does not
                % participate in the color scale.
                rectangle(ax, ...
                    'Position', [ ...
                        col_idx - 0.5, row_idx - 0.5, 1, 1], ...
                    'FaceColor', 'k', ...
                    'EdgeColor', 'none');

                text(ax, col_idx, row_idx, '1', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'FontWeight', 'bold', ...
                    'FontSize', 11, ...
                    'Color', 'w');

                continue;
            end

            value = similarity_matrix(row_idx, col_idx);
            n_value = n_matrix(row_idx, col_idx);

            if ~isfinite(value)
                value_text = sprintf('NaN\nn=%d', n_value);
                text_color = 'k';
            else
                value_text = sprintf('%.3f\nn=%d', value, n_value);

                normalized_value = ...
                    (value - color_limits(1)) / ...
                    (color_limits(2) - color_limits(1));

                % Appropriate for MATLAB's default parula colormap:
                % lower values are dark, upper values are light.
                if normalized_value < 0.35
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

    % Draw cell boundaries.
    for boundary = 0.5:1:(n_state + 0.5)
        xline(ax, boundary, '-', ...
            'Color', [1, 1, 1], ...
            'LineWidth', 0.5, ...
            'HandleVisibility', 'off');

        yline(ax, boundary, '-', ...
            'Color', [1, 1, 1], ...
            'LineWidth', 0.5, ...
            'HandleVisibility', 'off');
    end

    hold(ax, 'off');
end


function color_limits = calculate_off_diagonal_color_limits( ...
    values, metric_bounds)

    lower_bound = metric_bounds(1);
    upper_bound = metric_bounds(2);

    if isempty(values)
        color_limits = metric_bounds;
        return;
    end

    value_min = min(values);
    value_max = max(values);

    if value_min < value_max
        color_limits = [value_min, value_max];
        return;
    end

    % All finite off-diagonal entries have the same value.
    center_value = value_min;
    padding = max(0.05, 0.05 * abs(center_value));

    color_min = max(lower_bound, center_value - padding);
    color_max = min(upper_bound, center_value + padding);

    % Handle values exactly at one of the theoretical boundaries.
    if color_min == color_max
        if center_value >= upper_bound
            color_limits = [upper_bound - 0.05, upper_bound];
        elseif center_value <= lower_bound
            color_limits = [lower_bound, lower_bound + 0.05];
        else
            color_limits = [center_value - 0.05, center_value + 0.05];
        end
    else
        color_limits = [color_min, color_max];
    end
end


function [x_values, y_values] = extract_numeric_pair(state_x, state_y)
    validate_matching_filters(state_x, state_y);

    x12 = state_x.J12(:);
    y12 = state_y.J12(:);

    x21 = state_x.J21(:);
    y21 = state_y.J21(:);

    x_values = [x12; x21];
    y_values = [y12; y21];

    valid = isfinite(x_values) & isfinite(y_values);
    x_values = x_values(valid);
    y_values = y_values(valid);
end


function counts = make_pair_category_counts( ...
    state_x, state_y, err_multi)

    validate_matching_filters(state_x, state_y);

    x_cat12 = classify_connections( ...
        state_x.J12(:), state_x.err12(:), err_multi);

    y_cat12 = classify_connections( ...
        state_y.J12(:), state_y.err12(:), err_multi);

    x_cat21 = classify_connections( ...
        state_x.J21(:), state_x.err21(:), err_multi);

    y_cat21 = classify_connections( ...
        state_y.J21(:), state_y.err21(:), err_multi);

    x_category = [x_cat12; x_cat21];
    y_category = [y_cat12; y_cat21];

    valid = isfinite(x_category) & isfinite(y_category);
    x_category = x_category(valid);
    y_category = y_category(valid);

    class_values = [-1, 0, 1];
    counts = zeros(3, 3);

    for row_class = 1:3
        for col_class = 1:3
            counts(row_class, col_class) = sum( ...
                x_category == class_values(row_class) & ...
                y_category == class_values(col_class));
        end
    end
end


function category = classify_connections(J, err, err_multi)
    category = nan(size(J));

    valid = isfinite(J) & isfinite(err);
    category(valid) = 0;

    category(valid & (J > err_multi * err)) = 1;
    category(valid & (J < -err_multi * err)) = -1;
end


function kappa = compute_cohen_kappa(counts)
    total_n = sum(counts(:));

    if total_n == 0
        kappa = NaN;
        return;
    end

    observed_agreement = trace(counts) / total_n;

    row_marginal = sum(counts, 2) / total_n;
    column_marginal = sum(counts, 1)' / total_n;

    expected_agreement = sum(row_marginal .* column_marginal);

    if abs(1 - expected_agreement) < eps
        kappa = NaN;
    else
        kappa = ...
            (observed_agreement - expected_agreement) / ...
            (1 - expected_agreement);
    end
end


function [correlation_value, p_value, n_valid] = ...
    correlation_stat(x_values, y_values, correlation_type)

    valid = isfinite(x_values) & isfinite(y_values);
    x_values = x_values(valid);
    y_values = y_values(valid);

    n_valid = numel(x_values);

    if n_valid < 2
        correlation_value = NaN;
        p_value = NaN;
        return;
    end

    if all(x_values == x_values(1)) || ...
       all(y_values == y_values(1))
        correlation_value = NaN;
        p_value = NaN;
        return;
    end

    [R, P] = corr( ...
        x_values, y_values, ...
        'Type', correlation_type, ...
        'Rows', 'complete');

    correlation_value = R;
    p_value = P;
end


function validate_matching_filters(state_a, state_b)
    filter1_matches = ...
        numel(state_a.filter1) == numel(state_b.filter1) && ...
        all(state_a.filter1(:) == state_b.filter1(:));

    filter2_matches = ...
        numel(state_a.filter2) == numel(state_b.filter2) && ...
        all(state_a.filter2(:) == state_b.filter2(:));

    if ~filter1_matches || ~filter2_matches
        error(['State filters do not match between the two ' ...
               'comparison states.']);
    end

    if ~isequal(size(state_a.J12), size(state_b.J12)) || ...
       ~isequal(size(state_a.J21), size(state_b.J21))
        error(['State connectivity matrices do not have ' ...
               'matching sizes.']);
    end
end


function selected_rows = default_metadata_filter(mt)
    base_filter = ...
        strcmp(mt.kernel_name, "DeltaPure") & ...
        strcmp(mt.align, 'Last') & ...
        strcmp(mt.area, "Cortex") & ...
        strcmp(mt.injection, 'Muscimol') & ...
        (mt.epoch == 3000) & ...
        (mt.fold_idx == 0) & ...
        (mt.shuffle_idx == 0) & ...
        cellfun( ...
            @(x) ~isempty(x) && x == 15, ...
            mt.resting_dur_threshold);

    % Use Pre RestOpen as one anchor row per session.
    anchor_filter = ...
        base_filter & ...
        strcmp(mt.prepost, 'Pre') & ...
        strcmp(mt.state, 'RestOpen');

    selected_rows = false(height(mt), 1);
    anchor_indices = find(anchor_filter);

    required_preposts = {'Pre', 'Pre', 'Post', 'Post'};
    required_states = { ...
        'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};

    match_fields = { ...
        'animal_name', ...
        'injection', ...
        'align', ...
        'session_idx', ...
        'resting_dur_threshold', ...
        'area', ...
        'kernel_name', ...
        'reg_name', ...
        'epoch', ...
        'fold_idx', ...
        'shuffle_idx'};

    for anchor_i = 1:numel(anchor_indices)
        anchor_idx = anchor_indices(anchor_i);
        same_session = base_filter;

        for field_i = 1:numel(match_fields)
            field_name = match_fields{field_i};

            if ~ismember(field_name, mt.Properties.VariableNames)
                continue;
            end

            anchor_value = mt.(field_name)(anchor_idx);

            if iscell(mt.(field_name))
                anchor_value = anchor_value{1};

                same_session = same_session & cellfun( ...
                    @(x) isequal(x, anchor_value), ...
                    mt.(field_name));
            else
                same_session = same_session & arrayfun( ...
                    @(x) isequal(x, anchor_value), ...
                    mt.(field_name));
            end
        end

        has_all_required_states = true;

        for required_i = 1:numel(required_preposts)
            has_this_state = any( ...
                same_session & ...
                strcmp(mt.prepost, ...
                    required_preposts{required_i}) & ...
                strcmp(mt.state, ...
                    required_states{required_i}));

            if ~has_this_state
                has_all_required_states = false;

                anchor_name = mt.file_name{anchor_idx};

                warning([ ...
                    'Session %s is missing required state: ' ...
                    '%s %s. Skipping this session.'], ...
                    anchor_name, ...
                    required_preposts{required_i}, ...
                    required_states{required_i});

                break;
            end
        end

        if has_all_required_states
            selected_rows(anchor_idx) = true;
        end
    end

    fprintf( ...
        'Default metadata filter selected %d/%d rows.\n', ...
        sum(selected_rows), sum(base_filter));

    for row_idx = 1:height(mt)
        if selected_rows(row_idx)
            fprintf('  %s\n', mt.file_name{row_idx});
        end
    end

    if ~any(selected_rows)
        warning([ ...
            'No complete Pre/Post x Open/Close state sets were found. ' ...
            'Falling back to Pre RestOpen anchors.']);

        selected_rows = anchor_filter;
    end
end


function label = make_session_label(meta)
    animal_string = char(string(meta.animal_name));
    injection_string = char(string(meta.injection));
    session_string = char(string(meta.session_idx));

    label = sprintf( ...
        '%s %s session %s', ...
        animal_string, injection_string, session_string);
end


function key = state_key(prepost, state, kernel_idx)
    key = sprintf( ...
        'k%d_%s_%s', kernel_idx, prepost, state);
end


function state_struct = load_state_connectivity( ...
    root, meta, prepost, state, kernel_idx)

    state_struct = struct();
    state_struct.prepost = prepost;
    state_struct.state = state;
    state_struct.kernel_idx = kernel_idx;

    meta.prepost = prepost;
    meta.state = state;

    %% Load raster data to identify cell areas
    meta.filename = generate_filename('raster', meta);

    raster_filename = fullfile( ...
        root, 'Data', 'Working', 'raster', meta.filename);

    raster_data = load(raster_filename);

    fprintf('Loaded raster data for %s %s %s\n', ...
        meta.prepost, meta.state, meta.area);

    cell_area = raster_data.data.cell_area;

    filter1 = ismember(cell_area, {'ACC'});
    filter2 = ismember(cell_area, {'VLPFC'});

    %% Load GLM connectivity data
    meta.filename = generate_filename('GLM', meta);

    glm_filename = fullfile( ...
        root, 'Data', 'Working', 'GLM', meta.filename);

    glm_data = load(glm_filename);

    fprintf('Loaded GLM data for %s %s %s\n', ...
        meta.prepost, meta.state, meta.area);

    N = glm_data.meta.N;

    kernel_columns = ...
        (2 + N * (kernel_idx - 1)):(1 + N * kernel_idx);

    J = glm_data.data.model_par(:, kernel_columns);

    err = glm_data.data.model_err.total(:, kernel_columns);

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
