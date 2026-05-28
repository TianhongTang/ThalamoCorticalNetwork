%% Figure 3: Aggregated multi-session Open-Close network/correlation analysis
% Data from all selected metadata rows are pooled into one figure.
% Row 1: example network plots from the first valid selected session.
% Row 2: pooled Open-Close density plots, separated by positive/negative connections.
% Row 3: pooled scatter plots and categorical transition tables.

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

%% Base metadata template
% meta_template = struct();
% meta_template.animal_name = 'Slayer';
% meta_template.injection = 'Muscimol';
% meta_template.align = 'Last';
% meta_template.session_idx = 6;
% meta_template.resting_dur_threshold = 15;
% meta_template.area = 'Cortex';
% meta_template.shuffle_idx = 0;
% meta_template.kernel_name = 'DeltaPure';
% meta_template.reg_name = 'L2=0_2';
% meta_template.epoch = 3000;
% meta_template.fold_idx = 0;

preposts = {'Pre', 'Pre', 'Post', 'Post'};
states   = {'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};
n_state = numel(states);

%% Parameters
kernel_idx = 1;
err_multi = 1; % threshold for significant J, in multiples of the GLM error estimate.
network_err_multi = 2;
density_nbin = 60;
scatter_marker_size = 8;
scatter_alpha = 0.25;
density_clip_percentile = [0.5, 99.5];
density_use_log_count = true;
category_labels = {'Negative', 'Non-sig', 'Positive'};

skip_failed_sessions = false;
max_sessions_to_include = inf; % set smaller while debugging.

%% Load and filter metadata
mt = load_meta(root, 'table');
mt = mt.GLM;

% -------------------------------------------------------------------------
% EDIT THIS FILTER FOR THE FINAL SESSION SET.
% Examples:
% selected_rows = strcmp(string(mt.animal_name), "Slayer") & strcmp(string(mt.injection), "Muscimol");
% selected_rows = mt.session_idx >= 1 & strcmp(string(mt.injection), "Muscimol");
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

%% Pool data across all selected sessions
data_pre_all = empty_connection_data();
data_post_all = empty_connection_data();

cat_counts_pre_total = zeros(numel(category_labels), numel(category_labels));
cat_counts_post_total = zeros(numel(category_labels), numel(category_labels));

first_valid_state_data = [];
first_valid_label = '';
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

        data_pre_all = append_connection_data(data_pre_all, data_pre);
        data_post_all = append_connection_data(data_post_all, data_post);

        [cat_counts_pre, ~, ~, ~] = make_open_close_category_counts(pre_open, pre_close, err_multi);
        [cat_counts_post, ~, ~, ~] = make_open_close_category_counts(post_open, post_close, err_multi);

        cat_counts_pre_total = cat_counts_pre_total + cat_counts_pre;
        cat_counts_post_total = cat_counts_post_total + cat_counts_post;

        % Not use first but 12th session.
        if isempty(first_valid_state_data) && session_i==12
            first_valid_state_data = state_data;
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

%% Pooled statistics
[rho_pre_all, p_pre_all, n_pre_all] = pearson_stats(data_pre_all.x_abs, data_pre_all.y_abs);
[rho_post_all, p_post_all, n_post_all] = pearson_stats(data_post_all.x_abs, data_post_all.y_abs);

[rho_pre_pos, p_pre_pos, n_pre_pos] = pearson_stats(data_pre_all.xpos_abs, data_pre_all.ypos_abs);
[rho_pre_neg, p_pre_neg, n_pre_neg] = pearson_stats(data_pre_all.xneg_abs, data_pre_all.yneg_abs);
[rho_post_pos, p_post_pos, n_post_pos] = pearson_stats(data_post_all.xpos_abs, data_post_all.ypos_abs);
[rho_post_neg, p_post_neg, n_post_neg] = pearson_stats(data_post_all.xneg_abs, data_post_all.yneg_abs);

agreement_pre = compute_raw_agreement(cat_counts_pre_total);
kappa_pre = compute_cohen_kappa(cat_counts_pre_total);
cat_n_pre = sum(cat_counts_pre_total(:));

agreement_post = compute_raw_agreement(cat_counts_post_total);
kappa_post = compute_cohen_kappa(cat_counts_post_total);
cat_n_post = sum(cat_counts_post_total(:));

fprintf('Pooled Pre all same-sign significant: rho = %.6f, p = %.3e, n = %d\n', rho_pre_all, p_pre_all, n_pre_all);
fprintf('Pooled Post all same-sign significant: rho = %.6f, p = %.3e, n = %d\n', rho_post_all, p_post_all, n_post_all);
fprintf('Pooled Pre categorical agreement: agreement = %.6f, kappa = %.6f, n = %d\n', agreement_pre, kappa_pre, cat_n_pre);
fprintf('Pooled Post categorical agreement: agreement = %.6f, kappa = %.6f, n = %d\n', agreement_post, kappa_post, cat_n_post);

%% Figure
f = figure('Color', 'w');
tiles = tiledlayout(3, 4, "TileSpacing", "Compact", "Padding", "Compact");

%% Row 1: example network plot from the first valid session
for state_i = 1:n_state
    ax = nexttile(state_i);
    call_plot_network(ax, ...
        first_valid_state_data(state_i).J12, first_valid_state_data(state_i).J21, ...
        first_valid_state_data(state_i).err12, first_valid_state_data(state_i).err21, ...
        network_err_multi, [], []);
    title(ax, sprintf('Example network: %s, %s\n%s', ...
        first_valid_state_data(state_i).prepost, first_valid_state_data(state_i).state, first_valid_label), ...
        'Interpreter', 'none');
end

%% Row 2: pooled density plots, positive and negative separated
ax = nexttile(5);
plot_open_close_density(ax, data_pre_all.xpos_abs, data_pre_all.ypos_abs, ...
    rho_pre_pos, p_pre_pos, n_pre_pos, density_nbin, density_clip_percentile, density_use_log_count, ...
    sprintf('Pooled Pre positive density, sessions=%d', valid_session_count));

ax = nexttile(6);
plot_open_close_density(ax, data_pre_all.xneg_abs, data_pre_all.yneg_abs, ...
    rho_pre_neg, p_pre_neg, n_pre_neg, density_nbin, density_clip_percentile, density_use_log_count, ...
    sprintf('Pooled Pre negative density, sessions=%d', valid_session_count));

ax = nexttile(7);
plot_open_close_density(ax, data_post_all.xpos_abs, data_post_all.ypos_abs, ...
    rho_post_pos, p_post_pos, n_post_pos, density_nbin, density_clip_percentile, density_use_log_count, ...
    sprintf('Pooled Post positive density, sessions=%d', valid_session_count));

ax = nexttile(8);
plot_open_close_density(ax, data_post_all.xneg_abs, data_post_all.yneg_abs, ...
    rho_post_neg, p_post_neg, n_post_neg, density_nbin, density_clip_percentile, density_use_log_count, ...
    sprintf('Pooled Post negative density, sessions=%d', valid_session_count));

%% Row 3: pooled scatter and categorical transition tables
ax = nexttile(9);
plot_open_close_scatter(ax, data_pre_all, rho_pre_all, p_pre_all, n_pre_all, ...
    scatter_marker_size, scatter_alpha, sprintf('Pooled Pre scatter, sessions=%d', valid_session_count));

ax = nexttile(10);
plot_category_transition_table(ax, cat_counts_pre_total, category_labels, agreement_pre, kappa_pre, cat_n_pre, ...
    sprintf('Pooled Pre categorical, sessions=%d', valid_session_count));

ax = nexttile(11);
plot_open_close_scatter(ax, data_post_all, rho_post_all, p_post_all, n_post_all, ...
    scatter_marker_size, scatter_alpha, sprintf('Pooled Post scatter, sessions=%d', valid_session_count));

ax = nexttile(12);
plot_category_transition_table(ax, cat_counts_post_total, category_labels, agreement_post, kappa_post, cat_n_post, ...
    sprintf('Pooled Post categorical, sessions=%d', valid_session_count));

%% Export to pdf
fig = gcf;

save_folder = fullfile(root, 'Figures', 'Paper');
check_path(save_folder);

figWidth  = 16.0;   % inches
figHeight = 15.0;   % inches
resolution = 300;  % dpi; mainly affects rasterized components

set(fig, 'Units', 'inches');
fig.Position(3:4) = [figWidth, figHeight];

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [figWidth, figHeight]);
set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);
set(fig, 'Color', 'w');

filename = fullfile(save_folder, sprintf('Figure3_k%d_preview.jpg', kernel_idx));
exportgraphics(fig, filename, ...
    'ContentType', 'image', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);
filename = fullfile(save_folder, sprintf('Figure3_k%d.pdf', kernel_idx));
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
                warning('Session %s is missing required Pre/Post x Open/Close combination: %s %s. Skipping this session.', ...
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
        warning('Default metadata filter selected no complete Pre/Post x Open/Close sets. Falling back to Pre RestOpen anchors.');
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

function plot_failed_session_block(row_offset, session_label, message)
    for row = 1:3
        for col = 1:4
            ax = nexttile((row_offset + row - 1) * 4 + col);
            axis(ax, 'off');
            if row == 1 && col == 1
                text(ax, 0, 0.5, sprintf('Skipped: %s\n%s', session_label, message), ...
                    'Interpreter', 'none', 'Color', 'r');
            end
        end
    end
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

    pos = x > err_multi * xerr & y > err_multi * yerr;
    neg = x < -err_multi * xerr & y < -err_multi * yerr;

    data = struct();
    data.x_orig = x;
    data.y_orig = y;
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
            counts(i_class, j_class) = sum(open_cat == class_values(i_class) & close_cat == class_values(j_class));
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

function plot_open_close_scatter(ax, data, rho, pval, n_valid, marker_size, marker_alpha, title_text)
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
    plot_identity_line(ax, data.x_abs, data.y_abs);
    hold(ax, 'off');

    % calculate cos similarity for orig and abs values, and add to title
    cos_sim_orig = dot(data.x_orig, data.y_orig) / (norm(data.x_orig) * norm(data.y_orig));
    cos_sim_abs = dot(data.x_abs, data.y_abs) / (norm(data.x_abs) * norm(data.y_abs));
    title_text = sprintf('%s\nPearson r = %.6f (p = %.3e, n = %d)\ncos sim orig = %.6f, cos sim abs = %.6f', ...
        title_text, rho, pval, n_valid, cos_sim_orig, cos_sim_abs);

    axis(ax, 'square');
    xlim([0, 3.5]);
    ylim([0, 3.5]);
    xlabel(ax, 'Open |J_{ij}|');
    ylabel(ax, 'Close |J_{ij}|');
    title(ax, title_text, 'Interpreter', 'none');
    legend(ax, 'Location', 'northeast');
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
    xlabel(ax, 'Open |J_{ij}|');
    ylabel(ax, 'Close |J_{ij}|');
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
        'BackgroundColor', 'w', 'Margin', 4, 'FontSize', 8);
end

function plot_category_transition_table(ax, counts, category_labels, agreement, kappa, n_valid, title_text)
    % counts rows = Open category, columns = Close category.
    % Plot with x = Open and y = Close to match scatter/density axes.
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

    xlabel(ax, 'Open category');
    ylabel(ax, 'Close category');

    cb = colorbar(ax);
    ylabel(cb, 'count');

    max_count = max(plot_counts(:));
    for y_idx = 1:n_cat
        for x_idx = 1:n_cat
            this_count = plot_counts(y_idx, x_idx);
            if max_count > 0 && this_count > 0.55 * max_count
                text_color = 'w';
            else
                text_color = 'k';
            end
            text(ax, x_idx, y_idx, sprintf('%d', this_count), ...
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
