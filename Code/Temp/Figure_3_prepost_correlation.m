%% Figure 3: Aggregated multi-session Pre-Post correlation analysis
% Data from all selected metadata rows are pooled into one figure.
% This version keeps one kernel_idx and compares Pre vs Post within the same state.
% Rows are organized as two rows per state:
%   row 2*i-1: Pre network, all-connection scatter, either-sig scatter, same-sign abs scatter, 3x3 category table.
%   row 2*i  : Post network, all-connection density, both-sig scatter, sign-switch abs scatter, 2x2 category table.

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
compare_states = {'RestOpen', 'RestClose'};
n_compare = numel(compare_states);

%% Parameters

% Keep this as the single GLM kernel column to analyze.
kernel_idx = 2;
kernel_label = sprintf('K%d', kernel_idx);

err_multi = 1; % threshold for significant J, in multiples of the GLM error estimate.
network_err_multi = 2;
density_nbin = 60;
scatter_marker_size = 8;
scatter_alpha = 0.25;
density_clip_percentile = [0.5, 99.5];
density_use_log_count = true;
category_labels = {'Negative', 'Non-sig', 'Positive'};

% Two rows per state comparison.
n_row = 2 * n_compare;
n_col = 5;

show_legend = true; % true or false, applied to all scatter-plot legends.

colors = struct();
colors.non_sig = [0.65, 0.65, 0.65];
colors.x_only = [0.1, 0.75, 0.1];
colors.y_only = [0.95, 0.55, 0.1];
colors.both_sig = [0.45, 0.1, 0.75];

colors.pos = [1, 0.2, 0.2];
colors.neg = [0.1, 0.45, 1];
colors.switch_posneg = [1.0, 0.5, 0.0];
colors.switch_negpos = [0.0, 0.65, 0.8];
colors.switch = [0.45, 0.1, 0.75];

colors.identity_line = [1, 0, 0];
colors.zero_line = [0, 0, 0];

skip_failed_sessions = false;
max_sessions_to_include = inf; % set smaller while debugging.
example_session_idx = 12; % Used for the network example. Falls back to first valid session if unavailable.

%% Load and filter metadata
mt = load_meta(root, 'table');
mt = mt.GLM;

% -------------------------------------------------------------------------
% EDIT THIS FILTER FOR THE FINAL SESSION SET.
% The default filter selects complete Pre/Post x RestOpen/RestClose sets.
% This script then compares Pre vs Post within each state using kernel_idx.
% -------------------------------------------------------------------------
selected_rows = default_metadata_filter(mt); % edit or replace this line for the final session set.

selected_mt = mt(selected_rows, :);
meta_array = table2struct(selected_mt);
if isfinite(max_sessions_to_include)
    meta_array = meta_array(1:min(numel(meta_array), max_sessions_to_include));
end

if isempty(meta_array)
    error('No metadata rows selected.');
end

%% Pool Pre-Post comparison data across all selected sessions
data_all = repmat(empty_connection_data(), n_compare, 1);
cat_counts_total = zeros(numel(category_labels), numel(category_labels), n_compare);

first_valid_pre_state_data = [];
first_valid_post_state_data = [];
first_valid_label = '';
valid_session_count = 0;
failed_session_count = 0;

for session_i = 1:numel(meta_array)
    meta = meta_array(session_i);
    session_label = make_session_label(meta);
    fprintf('\n===== Loading %s (%d/%d) =====\n', session_label, session_i, numel(meta_array));

    try
        pre_state_data = struct([]);
        post_state_data = struct([]);

        for state_i = 1:n_compare
            state = compare_states{state_i};
            sd_pre = load_state_connectivity(root, meta, 'Pre', state, kernel_idx);
            sd_post = load_state_connectivity(root, meta, 'Post', state, kernel_idx);

            if state_i == 1
                pre_state_data = sd_pre;
                post_state_data = sd_post;
            else
                pre_state_data = [pre_state_data; sd_pre]; %#ok<AGROW>
                post_state_data = [post_state_data; sd_post]; %#ok<AGROW>
            end

            data_pair = make_pair_vectors(sd_pre, sd_post, err_multi);
            data_all(state_i) = append_connection_data(data_all(state_i), data_pair);

            [cat_counts, ~, ~, ~] = make_pair_category_counts(sd_pre, sd_post, err_multi);
            cat_counts_total(:, :, state_i) = cat_counts_total(:, :, state_i) + cat_counts;
        end

        % Keep first valid as fallback, but replace with example_session_idx if it exists.
        if isempty(first_valid_pre_state_data) || session_i == example_session_idx
            first_valid_pre_state_data = pre_state_data;
            first_valid_post_state_data = post_state_data;
            first_valid_label = session_label;
        end

        valid_session_count = valid_session_count + 1;
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

%% Pooled statistics by state
stats = struct([]);
category_labels_2x2 = {'Negative', 'Positive'};

for state_i = 1:n_compare
    data = data_all(state_i);

    [stats(state_i).rho_orig, stats(state_i).p_orig, stats(state_i).n_orig] = pearson_stats(data.x_orig, data.y_orig);
    [stats(state_i).rho_abs, stats(state_i).p_abs, stats(state_i).n_abs] = pearson_stats(data.x_abs, data.y_abs);
    [stats(state_i).rho_pos, stats(state_i).p_pos, stats(state_i).n_pos] = pearson_stats(data.xpos_abs, data.ypos_abs);
    [stats(state_i).rho_neg, stats(state_i).p_neg, stats(state_i).n_neg] = pearson_stats(data.xneg_abs, data.yneg_abs);

    counts3 = cat_counts_total(:, :, state_i);
    stats(state_i).cat_counts_3x3 = counts3;
    stats(state_i).agreement_3x3 = compute_raw_agreement(counts3);
    stats(state_i).kappa_3x3 = compute_cohen_kappa(counts3);
    stats(state_i).cat_n_3x3 = sum(counts3(:));

    counts2 = counts3([1, 3], [1, 3]);
    stats(state_i).cat_counts_2x2 = counts2;
    stats(state_i).agreement_2x2 = compute_raw_agreement(counts2);
    stats(state_i).kappa_2x2 = compute_cohen_kappa(counts2);
    stats(state_i).cat_n_2x2 = sum(counts2(:));

    fprintf('Pooled %s Pre-vs-Post %s all signed: rho = %.6f, p = %.3e, n = %d\n', ...
        compare_states{state_i}, kernel_label, stats(state_i).rho_orig, stats(state_i).p_orig, stats(state_i).n_orig);
    fprintf('Pooled %s Pre-vs-Post categorical agreement: agreement = %.6f, kappa = %.6f, n = %d\n', ...
        compare_states{state_i}, stats(state_i).agreement_3x3, stats(state_i).kappa_3x3, stats(state_i).cat_n_3x3);
end

%% Figure
f = figure('Color', 'w', 'Visible', 'off');
tiles = tiledlayout(n_row, n_col, "TileSpacing", "Compact", "Padding", "Compact");
tile_idx = @(row, col) (row - 1) * n_col + col;

x_name = 'Pre';
y_name = 'Post';

for state_i = 1:n_compare
    row_top = 2 * state_i - 1;
    row_bottom = 2 * state_i;
    state_label = compare_states{state_i};
    title_prefix = sprintf('%s: Pre vs Post, %s', state_label, kernel_label);

    % Network examples: top row = Pre, bottom row = Post.
    ax = nexttile(tile_idx(row_top, 1));
    call_plot_network(ax, ...
        first_valid_pre_state_data(state_i).J12, first_valid_pre_state_data(state_i).J21, ...
        first_valid_pre_state_data(state_i).err12, first_valid_pre_state_data(state_i).err21, ...
        network_err_multi, [], []);
    title(ax, sprintf('Example network: %s\n%s, Pre, %s', ...
        first_valid_label, state_label, kernel_label), 'Interpreter', 'none');

    ax = nexttile(tile_idx(row_bottom, 1));
    call_plot_network(ax, ...
        first_valid_post_state_data(state_i).J12, first_valid_post_state_data(state_i).J21, ...
        first_valid_post_state_data(state_i).err12, first_valid_post_state_data(state_i).err21, ...
        network_err_multi, [], []);
    title(ax, sprintf('Example network: %s\n%s, Post, %s', ...
        first_valid_label, state_label, kernel_label), 'Interpreter', 'none');

    % All-connection signed scatter and density.
    ax = nexttile(tile_idx(row_top, 2));
    plot_pair_scatter(ax, data_all(state_i), 'all_signed', ...
        scatter_marker_size, scatter_alpha, colors, show_legend, ...
        sprintf('%s all connections, sessions=%d', title_prefix, valid_session_count), x_name, y_name);

    ax = nexttile(tile_idx(row_bottom, 2));
    plot_pair_density(ax, data_all(state_i).x_orig, data_all(state_i).y_orig, ...
        stats(state_i).rho_orig, stats(state_i).p_orig, stats(state_i).n_orig, ...
        density_nbin, density_clip_percentile, density_use_log_count, ...
        sprintf('%s all-connection density, sessions=%d', title_prefix, valid_session_count), ...
        'signed', colors, x_name, y_name);

    % Significant signed subsets.
    ax = nexttile(tile_idx(row_top, 3));
    plot_pair_scatter(ax, data_all(state_i), 'either_sig_signed', ...
        scatter_marker_size, scatter_alpha, colors, show_legend, ...
        sprintf('%s either significant signed, sessions=%d', title_prefix, valid_session_count), x_name, y_name);

    ax = nexttile(tile_idx(row_bottom, 3));
    plot_pair_scatter(ax, data_all(state_i), 'both_sig_signed', ...
        scatter_marker_size, scatter_alpha, colors, show_legend, ...
        sprintf('%s both significant signed, sessions=%d', title_prefix, valid_session_count), x_name, y_name);

    % Same-sign and sign-switch abs plots.
    ax = nexttile(tile_idx(row_top, 4));
    plot_pair_scatter(ax, data_all(state_i), 'same_abs', ...
        scatter_marker_size, scatter_alpha, colors, show_legend, ...
        sprintf('%s same-sign abs, sessions=%d', title_prefix, valid_session_count), x_name, y_name);

    ax = nexttile(tile_idx(row_bottom, 4));
    plot_pair_scatter(ax, data_all(state_i), 'switch_abs', ...
        scatter_marker_size, scatter_alpha, colors, show_legend, ...
        sprintf('%s sign-switch abs, sessions=%d', title_prefix, valid_session_count), x_name, y_name);

    % Categorical agreement tables.
    ax = nexttile(tile_idx(row_top, 5));
    plot_category_transition_table(ax, stats(state_i).cat_counts_3x3, category_labels, ...
        stats(state_i).agreement_3x3, stats(state_i).kappa_3x3, stats(state_i).cat_n_3x3, ...
        sprintf('%s categorical 3x3, sessions=%d', title_prefix, valid_session_count), x_name, y_name);

    ax = nexttile(tile_idx(row_bottom, 5));
    plot_category_transition_table(ax, stats(state_i).cat_counts_2x2, category_labels_2x2, ...
        stats(state_i).agreement_2x2, stats(state_i).kappa_2x2, stats(state_i).cat_n_2x2, ...
        sprintf('%s categorical 2x2, sessions=%d', title_prefix, valid_session_count), x_name, y_name);
end

%% Export to pdf
fig = gcf;

save_folder = fullfile(root, 'Figures', 'Paper');
check_path(save_folder);

figWidth  = 4 * n_col;  % inches
figHeight = 4 * n_row;  % inches
resolution = 300;  % dpi; mainly affects rasterized components

set(fig, 'Units', 'inches');
fig.Position(3:4) = [figWidth, figHeight];

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [figWidth, figHeight]);
set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);
set(fig, 'Color', 'w');

filename = fullfile(save_folder, sprintf('Figure3_pre_vs_post_k%d_preview.jpg', kernel_idx));
exportgraphics(fig, filename, ...
    'ContentType', 'image', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);
filename = fullfile(save_folder, sprintf('Figure3_pre_vs_post_k%d.pdf', kernel_idx));
exportgraphics(fig, filename, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);

close(fig);

%% functions



function data = empty_connection_data()
    data = struct();
    data.x_orig = [];
    data.y_orig = [];
    data.x_err_orig = [];
    data.y_err_orig = [];
    data.x_cat_orig = [];
    data.y_cat_orig = [];

    data.xpos = [];
    data.ypos = [];
    data.xneg = [];
    data.yneg = [];
    data.xpos_abs = [];
    data.ypos_abs = [];
    data.xneg_abs = [];
    data.yneg_abs = [];
    data.x_abs = [];
    data.y_abs = [];

    data.xswitch_posneg_abs = [];
    data.yswitch_posneg_abs = [];
    data.xswitch_negpos_abs = [];
    data.yswitch_negpos_abs = [];
    data.xswitch_abs = [];
    data.yswitch_abs = [];

    data.x_either_sig = [];
    data.y_either_sig = [];
    data.x_cat_either_sig = [];
    data.y_cat_either_sig = [];

    data.x_both_sig = [];
    data.y_both_sig = [];
    data.x_cat_both_sig = [];
    data.y_cat_both_sig = [];
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
                anchor_name = mt.file_name{idx};
                warning('Session %s is missing required Pre/Post x RestOpen/RestClose combination: %s %s. Skipping this session.', ...
                    anchor_name, required_preposts{q}, required_states{q});
                break;
            end
        end

        if has_all_required
            selected_rows(idx) = true;
        end
    end

    fprintf('Default metadata filter selected %d/%d rows:\n', sum(selected_rows), sum(base_filter));
    for k = 1:height(mt)
        if selected_rows(k)
            fprintf('  %s\n', mt.file_name{k});
        end
    end

    if ~any(selected_rows)
        warning('Default metadata filter selected no complete Pre/Post x RestOpen/RestClose sets. Falling back to Pre RestOpen anchors.');
        selected_rows = anchor_filter;
    end
end



function mask = table_column_match(mt, candidate_names, target_value)
    mask = true(height(mt), 1);
    if ~istable(mt)
        return;
    end

    names = mt.Properties.VariableNames;
    match_name = '';
    for k = 1:numel(candidate_names)
        idx = find(strcmpi(names, candidate_names{k}), 1);
        if ~isempty(idx)
            match_name = names{idx};
            break;
        end
    end

    if isempty(match_name)
        return;
    end

    col = mt.(match_name);
    if iscell(col)
        col_str = string(col);
        target_str = string(target_value);
        mask = strcmp(col_str, target_str);
    elseif isstring(col) || ischar(col) || iscategorical(col)
        mask = strcmp(string(col), string(target_value));
    else
        if isnumeric(target_value) || islogical(target_value)
            mask = col == target_value;
        else
            target_num = str2double(string(target_value));
            if isfinite(target_num)
                mask = col == target_num;
            end
        end
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



function [data, label_text] = make_pair_vectors(x_state, y_state, err_multi)
    validate_matching_filters(x_state, y_state);

    x12 = x_state.J12(:);
    y12 = y_state.J12(:);
    x21 = x_state.J21(:);
    y21 = y_state.J21(:);
    xerr12 = x_state.err12(:);
    yerr12 = y_state.err12(:);
    xerr21 = x_state.err21(:);
    yerr21 = y_state.err21(:);

    x = [x12; x21];
    y = [y12; y21];
    xerr = [xerr12; xerr21];
    yerr = [yerr12; yerr21];

    valid = isfinite(x) & isfinite(y) & isfinite(xerr) & isfinite(yerr);
    x = x(valid);
    y = y(valid);
    xerr = xerr(valid);
    yerr = yerr(valid);

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

    x_sig = x_cat ~= 0;
    y_sig = y_cat ~= 0;

    same_pos = x_pos & y_pos;
    same_neg = x_neg & y_neg;
    switch_posneg = x_pos & y_neg;
    switch_negpos = x_neg & y_pos;
    switch_any = switch_posneg | switch_negpos;

    either_sig = x_sig | y_sig;
    both_sig = x_sig & y_sig;

    data = struct();
    data.x_orig = x;
    data.y_orig = y;
    data.x_err_orig = xerr;
    data.y_err_orig = yerr;
    data.x_cat_orig = x_cat;
    data.y_cat_orig = y_cat;

    % Same-sign significant connections, plotted as abs values.
    data.xpos = x(same_pos);
    data.ypos = y(same_pos);
    data.xneg = x(same_neg);
    data.yneg = y(same_neg);

    data.xpos_abs = data.xpos;
    data.ypos_abs = data.ypos;
    data.xneg_abs = -data.xneg;
    data.yneg_abs = -data.yneg;

    data.x_abs = [data.xpos_abs; data.xneg_abs];
    data.y_abs = [data.ypos_abs; data.yneg_abs];

    % Sign-switch significant connections, plotted as abs values.
    data.xswitch_posneg_abs = abs(x(switch_posneg));
    data.yswitch_posneg_abs = abs(y(switch_posneg));
    data.xswitch_negpos_abs = abs(x(switch_negpos));
    data.yswitch_negpos_abs = abs(y(switch_negpos));

    data.xswitch_abs = abs(x(switch_any));
    data.yswitch_abs = abs(y(switch_any));

    % Signed original-J subsets.
    data.x_either_sig = x(either_sig);
    data.y_either_sig = y(either_sig);
    data.x_cat_either_sig = x_cat(either_sig);
    data.y_cat_either_sig = y_cat(either_sig);

    data.x_both_sig = x(both_sig);
    data.y_both_sig = y(both_sig);
    data.x_cat_both_sig = x_cat(both_sig);
    data.y_cat_both_sig = y_cat(both_sig);

    label_text = 'ACC↔VLPFC J_{ij}';
end


function [counts, agreement, kappa, n_valid] = make_pair_category_counts(x_state, y_state, err_multi)
    validate_matching_filters(x_state, y_state);

    x_cat12 = classify_connections(x_state.J12(:), x_state.err12(:), err_multi);
    y_cat12 = classify_connections(y_state.J12(:), y_state.err12(:), err_multi);
    x_cat21 = classify_connections(x_state.J21(:), x_state.err21(:), err_multi);
    y_cat21 = classify_connections(y_state.J21(:), y_state.err21(:), err_multi);

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



function plot_pair_scatter(ax, data, plot_mode, marker_size, marker_alpha, colors, show_legend, title_text, x_name, y_name)
    cla(ax);
    hold(ax, 'on');

    switch plot_mode
        case 'all_signed'
            x = data.x_orig;
            y = data.y_orig;
            x_sig = data.x_cat_orig ~= 0;
            y_sig = data.y_cat_orig ~= 0;

            none_sig = ~x_sig & ~y_sig;
            x_only = x_sig & ~y_sig;
            y_only = ~x_sig & y_sig;
            both_sig = x_sig & y_sig;

            scatter(ax, x(none_sig), y(none_sig), marker_size, 'filled', ...
                'MarkerFaceColor', colors.non_sig, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'both non-sig');
            scatter(ax, x(x_only), y(x_only), marker_size, 'filled', ...
                'MarkerFaceColor', colors.x_only, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s only', x_name));
            scatter(ax, x(y_only), y(y_only), marker_size, 'filled', ...
                'MarkerFaceColor', colors.y_only, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s only', y_name));
            scatter(ax, x(both_sig), y(both_sig), marker_size, 'filled', ...
                'MarkerFaceColor', colors.both_sig, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'both sig');

            axis_limit = get_symmetric_axis_limit(x, y);
            xlabel_text = sprintf('%s J_{ij}', x_name);
            ylabel_text = sprintf('%s J_{ij}', y_name);
            axis_mode = 'signed';

        case 'same_abs'
            x = data.x_abs;
            y = data.y_abs;
            scatter(ax, data.xpos_abs, data.ypos_abs, marker_size, 'filled', ...
                'MarkerFaceColor', colors.pos, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'pos');
            scatter(ax, data.xneg_abs, data.yneg_abs, marker_size, 'filled', ...
                'MarkerFaceColor', colors.neg, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'neg');

            axis_limit = [0, max(3.5, get_positive_axis_max(x, y))];
            xlabel_text = sprintf('%s |J_{ij}|', x_name);
            ylabel_text = sprintf('%s |J_{ij}|', y_name);
            axis_mode = 'positive';

        case 'switch_abs'
            x = data.xswitch_abs;
            y = data.yswitch_abs;
            scatter(ax, data.xswitch_posneg_abs, data.yswitch_posneg_abs, marker_size, 'filled', ...
                'MarkerFaceColor', colors.switch_posneg, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s pos / %s neg', x_name, y_name));
            scatter(ax, data.xswitch_negpos_abs, data.yswitch_negpos_abs, marker_size, 'filled', ...
                'MarkerFaceColor', colors.switch_negpos, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s neg / %s pos', x_name, y_name));

            axis_limit = [0, max(3.5, get_positive_axis_max(x, y))];
            xlabel_text = sprintf('%s |J_{ij}|', x_name);
            ylabel_text = sprintf('%s |J_{ij}|', y_name);
            axis_mode = 'positive';

        case 'either_sig_signed'
            x = data.x_either_sig;
            y = data.y_either_sig;
            x_sig = data.x_cat_either_sig ~= 0;
            y_sig = data.y_cat_either_sig ~= 0;
            x_only = x_sig & ~y_sig;
            y_only = ~x_sig & y_sig;
            both = x_sig & y_sig;

            scatter(ax, x(x_only), y(x_only), marker_size, 'filled', ...
                'MarkerFaceColor', colors.x_only, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s only', x_name));
            scatter(ax, x(y_only), y(y_only), marker_size, 'filled', ...
                'MarkerFaceColor', colors.y_only, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s only', y_name));
            scatter(ax, x(both), y(both), marker_size, 'filled', ...
                'MarkerFaceColor', colors.both_sig, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'both sig');

            axis_limit = get_symmetric_axis_limit(x, y);
            xlabel_text = sprintf('%s J_{ij}', x_name);
            ylabel_text = sprintf('%s J_{ij}', y_name);
            axis_mode = 'signed';

        case 'both_sig_signed'
            x = data.x_both_sig;
            y = data.y_both_sig;
            x_cat = data.x_cat_both_sig;
            y_cat = data.y_cat_both_sig;

            pos_mask = x_cat == 1 & y_cat == 1;
            neg_mask = x_cat == -1 & y_cat == -1;
            switch_mask = (x_cat ~= 0) & (y_cat ~= 0) & (x_cat ~= y_cat);

            scatter(ax, x(pos_mask), y(pos_mask), marker_size, 'filled', ...
                'MarkerFaceColor', colors.pos, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'pos');
            scatter(ax, x(neg_mask), y(neg_mask), marker_size, 'filled', ...
                'MarkerFaceColor', colors.neg, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'neg');
            scatter(ax, x(switch_mask), y(switch_mask), marker_size, 'filled', ...
                'MarkerFaceColor', colors.switch, ...
                'MarkerFaceAlpha', marker_alpha, ...
                'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'switch');

            axis_limit = get_symmetric_axis_limit(x, y);
            xlabel_text = sprintf('%s J_{ij}', x_name);
            ylabel_text = sprintf('%s J_{ij}', y_name);
            axis_mode = 'signed';

        otherwise
            error('Unknown plot_mode: %s', plot_mode);
    end

    [rho, pval, n_valid] = pearson_stats(x, y);
    cos_sim = cosine_similarity_omitnan(x, y);

    plot(ax, axis_limit, axis_limit, '--', 'Color', colors.identity_line, 'LineWidth', 1, 'HandleVisibility', 'off');
    if strcmp(axis_mode, 'signed')
        xline(ax, 0, ':', 'Color', colors.zero_line, 'HandleVisibility', 'off');
        yline(ax, 0, ':', 'Color', colors.zero_line, 'HandleVisibility', 'off');
    end

    xlim(ax, axis_limit);
    ylim(ax, axis_limit);
    hold(ax, 'off');

    axis(ax, 'square');
    xlabel(ax, xlabel_text);
    ylabel(ax, ylabel_text);
    title(ax, sprintf('%s\nPearson r = %.6f (p = %.3e, n = %d)\ncos sim = %.6f', ...
        title_text, rho, pval, n_valid, cos_sim), 'Interpreter', 'none');

    if show_legend
        legend(ax, 'Location', 'northeastoutside');
    else
        legend(ax, 'off');
    end
    add_stats_text(ax, rho, pval, n_valid);
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




function plot_pair_density(ax, x, y, rho, pval, n_valid, nbin, clip_percentile, use_log_count, title_text, axis_mode, colors, x_name, y_name)
    if nargin < 10 || isempty(title_text)
        title_text = '';
    end
    if nargin < 11 || isempty(axis_mode)
        axis_mode = 'abs';
    end
    if nargin < 12 || isempty(colors)
        colors = struct();
        colors.identity_line = [1, 0, 0];
        colors.zero_line = [0, 0, 0];
    end
    if nargin < 13 || isempty(x_name)
        x_name = 'X';
    end
    if nargin < 14 || isempty(y_name)
        y_name = 'Y';
    end

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
    plot_identity_line_with_limits(ax, edges_x(1), edges_x(end), edges_y(1), edges_y(end), colors.identity_line);
    if strcmp(axis_mode, 'signed')
        xline(ax, 0, ':', 'Color', colors.zero_line, 'HandleVisibility', 'off');
        yline(ax, 0, ':', 'Color', colors.zero_line, 'HandleVisibility', 'off');
    end
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
    if strcmp(axis_mode, 'signed')
        xlabel(ax, sprintf('%s J_{ij}', x_name));
        ylabel(ax, sprintf('%s J_{ij}', y_name));
    else
        xlabel(ax, sprintf('%s |J_{ij}|', x_name));
        ylabel(ax, sprintf('%s |J_{ij}|', y_name));
    end
    title(ax, title_text, 'Interpreter', 'none');
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

    in_range = x_valid >= edges_x(1) & x_valid <= edges_x(end) & y_valid >= edges_y(1) & y_valid <= edges_y(end);
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



function plot_identity_line_with_limits(ax, xmin, xmax, ymin, ymax, line_color)
    vmin = max(xmin, ymin);
    vmax = min(xmax, ymax);
    if nargin < 6 || isempty(line_color)
        line_color = [1, 0, 0];
    end
    if isfinite(vmin) && isfinite(vmax) && vmin < vmax
        plot(ax, [vmin, vmax], [vmin, vmax], '--', 'Color', line_color, 'LineWidth', 1, 'HandleVisibility', 'off');
    else
        vals_min = min([xmin, ymin]);
        vals_max = max([xmax, ymax]);
        plot(ax, [vals_min, vals_max], [vals_min, vals_max], '--', 'Color', line_color, 'LineWidth', 1, 'HandleVisibility', 'off');
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
        'BackgroundColor', 'w', 'Margin', 4, 'FontSize', 8);
end



function plot_category_transition_table(ax, counts, category_labels, agreement, kappa, n_valid, title_text, x_name, y_name)
    % counts rows = kernel 1 category, columns = kernel 2 category.
    % Plot with x = kernel 1 and y = kernel 2 to match scatter/density axes.
    if nargin < 8 || isempty(x_name)
        x_name = 'X';
    end
    if nargin < 9 || isempty(y_name)
        y_name = 'Y';
    end

    plot_counts = counts.';
    imagesc(ax, plot_counts);
    set(ax, 'YDir', 'normal');
    axis(ax, 'square');

    n_cat = numel(category_labels);
    xticks(ax, 1:n_cat);
    yticks(ax, 1:n_cat);
    xticklabels(ax, category_labels);
    yticklabels(ax, category_labels);
    xtickangle(ax, 30);

    xlabel(ax, sprintf('%s category', x_name));
    ylabel(ax, sprintf('%s category', y_name));

    cb = colorbar(ax);
    ylabel(cb, 'count');

    max_count = max(plot_counts(:));
    if max_count == 0
        max_count = 1;
    end

    for row = 1:n_cat
        for col = 1:n_cat
            value = plot_counts(row, col);
            if value > 0.5 * max_count
                text_color = 'w';
            else
                text_color = 'k';
            end

            text(ax, col, row, sprintf('%d', value), ...
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

    title(ax, sprintf('%s\nAgreement = %s, kappa = %s, n = %d', ...
        title_text, agreement_str, kappa_str, n_valid), 'Interpreter', 'none');
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
