%% Supplement Figure: Open/Close and Pre/Post compound-pattern transition matrices
% This script generates 6 separate 2x6 figures and saves each as PDF + JPG.
%
% Figure versions:
%   1. Kernel 1, Pre -> Post transition of Open/Close compound patterns.
%   2. Kernel 2, Pre -> Post transition of Open/Close compound patterns.
%   3. Kernel 3, Pre -> Post transition of Open/Close compound patterns.
%   4. Kernel 1, Open vs Close comparison of Pre/Post compound patterns.
%   5. Kernel 2, Open vs Close comparison of Pre/Post compound patterns.
%   6. Kernel 3, Open vs Close comparison of Pre/Post compound patterns.
%
% Layout for every saved figure:
%   Row 1 = 9x9 compound-category table including non-significant category.
%   Row 2 = 4x4 compound-category table excluding any non-significant member.
%
% Columns:
%   1. Raw count.
%   2. Normalized by row.
%   3. Normalized by column.
%   4. Normalized by total.
%   5. O/E ratio under independence baseline.
%   6. Standardized residual under independence baseline.
%
% Compound category labels use +, o, -.
%   For Open/Close patterns: first symbol = RestOpen, second symbol = RestClose.
%   For Pre/Post patterns:   first symbol = Pre,      second symbol = Post.

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

err_multi = 1; % threshold for significant J, in multiples of the GLM error estimate.
category_labels_9 = {'++', '+o', '+-', 'o+', 'oo', 'o-', '-+', '-o', '--'};
category_labels_4 = {'++', '+-', '-+', '--'};

n_row = 2;
n_col = 6;
figure_visible = 'off';

skip_failed_sessions = false;
max_sessions_to_include = inf; % set smaller while debugging.

params = struct();
params.err_multi = err_multi;
params.category_labels_9 = category_labels_9;
params.category_labels_4 = category_labels_4;
params.n_row = n_row;
params.n_col = n_col;
params.figure_visible = figure_visible;

figure_configs = build_pattern_figure_configs(kernel_indices);

%% Load and filter metadata
mt = load_meta(root, 'table');
mt = mt.GLM;

% -------------------------------------------------------------------------
% EDIT THIS FILTER FOR THE FINAL SESSION SET.
% -------------------------------------------------------------------------
selected_rows = default_metadata_filter(mt);

selected_mt = mt(selected_rows, :);
meta_array = table2struct(selected_mt);
if isfinite(max_sessions_to_include)
    meta_array = meta_array(1:min(numel(meta_array), max_sessions_to_include));
end

if isempty(meta_array)
    error('No metadata rows selected.');
end

%% Initialize pooled storage
n_fig = numel(figure_configs);
pooled = struct([]);
for fig_i = 1:n_fig
    pooled(fig_i).counts9 = zeros(numel(category_labels_9), numel(category_labels_9));
    pooled(fig_i).counts4 = zeros(numel(category_labels_4), numel(category_labels_4));
    pooled(fig_i).valid_session_count = 0;
end

valid_session_count = 0;
failed_session_count = 0;

%% Pool data across all selected sessions
for session_i = 1:numel(meta_array)
    meta = meta_array(session_i);
    session_label = make_session_label(meta);
    fprintf('\n===== Loading %s (%d/%d) =====\n', session_label, session_i, numel(meta_array));

    try
        loaded_states = load_all_required_states(root, meta, kernel_indices);

        for fig_i = 1:n_fig
            cfg = figure_configs(fig_i);
            k = cfg.kernel_idx;

            pre_open_state   = loaded_states.(state_key('Pre',  'RestOpen',  k));
            pre_close_state  = loaded_states.(state_key('Pre',  'RestClose', k));
            post_open_state  = loaded_states.(state_key('Post', 'RestOpen',  k));
            post_close_state = loaded_states.(state_key('Post', 'RestClose', k));

            switch cfg.analysis_type
                case 'pre_to_post_oc_pattern'
                    [counts9, counts4] = make_pre_to_post_oc_pattern_transition_counts( ...
                        pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi);

                case 'open_vs_close_prepost_pattern'
                    [counts9, counts4] = make_open_vs_close_prepost_pattern_counts( ...
                        pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi);

                otherwise
                    error('Unknown analysis_type: %s', cfg.analysis_type);
            end

            pooled(fig_i).counts9 = pooled(fig_i).counts9 + counts9;
            pooled(fig_i).counts4 = pooled(fig_i).counts4 + counts4;
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

for fig_i = 1:n_fig
    pooled(fig_i).valid_session_count = valid_session_count;
end

%% Render and save all 6 figures
for fig_i = 1:n_fig
    render_pattern_transition_figure(root, figure_configs(fig_i), pooled(fig_i), params);
end


function figure_configs = build_pattern_figure_configs(kernel_indices)
    figure_configs = struct([]);
    fig_i = 0;

    for k_i = 1:numel(kernel_indices)
        k = kernel_indices(k_i);
        fig_i = fig_i + 1;
        figure_configs(fig_i).kernel_idx = k;
        figure_configs(fig_i).analysis_type = 'pre_to_post_oc_pattern';
        figure_configs(fig_i).output_stub = sprintf('SuppFig_OC_pattern_pre_to_post_k%d', k);
        figure_configs(fig_i).figure_title = sprintf('Kernel %d: Pre to Post transition of Open/Close patterns', k);
        figure_configs(fig_i).x_label = 'Pre Open/Close pattern';
        figure_configs(fig_i).y_label = 'Post Open/Close pattern';
        figure_configs(fig_i).row1_title = '9x9 Open/Close pattern transition';
        figure_configs(fig_i).row2_title = '4x4 Open/Close pattern transition, significant only';
    end

    for k_i = 1:numel(kernel_indices)
        k = kernel_indices(k_i);
        fig_i = fig_i + 1;
        figure_configs(fig_i).kernel_idx = k;
        figure_configs(fig_i).analysis_type = 'open_vs_close_prepost_pattern';
        figure_configs(fig_i).output_stub = sprintf('SuppFig_PrePost_pattern_open_vs_close_k%d', k);
        figure_configs(fig_i).figure_title = sprintf('Kernel %d: Open vs Close comparison of Pre/Post transition patterns', k);
        figure_configs(fig_i).x_label = 'Open Pre/Post pattern';
        figure_configs(fig_i).y_label = 'Close Pre/Post pattern';
        figure_configs(fig_i).row1_title = '9x9 Pre/Post pattern comparison';
        figure_configs(fig_i).row2_title = '4x4 Pre/Post pattern comparison, significant only';
    end
end

function loaded_states = load_all_required_states(root, meta, kernel_indices)
    loaded_states = struct();
    preposts = {'Pre', 'Post'};
    states = {'RestOpen', 'RestClose'};
    for k_i = 1:numel(kernel_indices)
        kernel_idx = kernel_indices(k_i);
        for pp_i = 1:numel(preposts)
            for st_i = 1:numel(states)
                pp = preposts{pp_i};
                st = states{st_i};
                key = state_key(pp, st, kernel_idx);
                loaded_states.(key) = load_state_connectivity(root, meta, pp, st, kernel_idx);
            end
        end
    end
end

function key = state_key(prepost, state, kernel_idx)
    key = sprintf('k%d_%s_%s', kernel_idx, prepost, state);
end

function [counts9, counts4] = make_pre_to_post_oc_pattern_transition_counts( ...
    pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi)

    [pre_open_cat, pre_close_cat, post_open_cat, post_close_cat] = make_all_condition_categories( ...
        pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi);

    % 9x9: Pre Open/Close compound pattern -> Post Open/Close compound pattern.
    pre_code9  = compound_category_code_9(pre_open_cat,  pre_close_cat);
    post_code9 = compound_category_code_9(post_open_cat, post_close_cat);
    counts9 = make_transition_counts(pre_code9, post_code9, 9);

    % 4x4: same transition, but only connections with no non-significant member
    % in either the Pre or Post Open/Close pair are retained.
    pre_code4  = compound_category_code_4(pre_open_cat,  pre_close_cat);
    post_code4 = compound_category_code_4(post_open_cat, post_close_cat);
    counts4 = make_transition_counts(pre_code4, post_code4, 4);
end

function [counts9, counts4] = make_open_vs_close_prepost_pattern_counts( ...
    pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi)

    [pre_open_cat, pre_close_cat, post_open_cat, post_close_cat] = make_all_condition_categories( ...
        pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi);

    % 9x9: Open Pre/Post compound pattern -> Close Pre/Post compound pattern.
    open_code9  = compound_category_code_9(pre_open_cat,  post_open_cat);
    close_code9 = compound_category_code_9(pre_close_cat, post_close_cat);
    counts9 = make_transition_counts(open_code9, close_code9, 9);

    % 4x4: same comparison, but only connections with no non-significant member
    % in either the Open or Close Pre/Post pair are retained.
    open_code4  = compound_category_code_4(pre_open_cat,  post_open_cat);
    close_code4 = compound_category_code_4(pre_close_cat, post_close_cat);
    counts4 = make_transition_counts(open_code4, close_code4, 4);
end

function [pre_open_cat, pre_close_cat, post_open_cat, post_close_cat] = make_all_condition_categories( ...
    pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi)

    validate_matching_filters(pre_open_state, pre_close_state);
    validate_matching_filters(pre_open_state, post_open_state);
    validate_matching_filters(pre_open_state, post_close_state);

    pre_open_cat12   = classify_connections(pre_open_state.J12(:),   pre_open_state.err12(:),   err_multi);
    pre_close_cat12  = classify_connections(pre_close_state.J12(:),  pre_close_state.err12(:),  err_multi);
    post_open_cat12  = classify_connections(post_open_state.J12(:),  post_open_state.err12(:),  err_multi);
    post_close_cat12 = classify_connections(post_close_state.J12(:), post_close_state.err12(:), err_multi);

    pre_open_cat21   = classify_connections(pre_open_state.J21(:),   pre_open_state.err21(:),   err_multi);
    pre_close_cat21  = classify_connections(pre_close_state.J21(:),  pre_close_state.err21(:),  err_multi);
    post_open_cat21  = classify_connections(post_open_state.J21(:),  post_open_state.err21(:),  err_multi);
    post_close_cat21 = classify_connections(post_close_state.J21(:), post_close_state.err21(:), err_multi);

    pre_open_cat   = [pre_open_cat12;   pre_open_cat21];
    pre_close_cat  = [pre_close_cat12;  pre_close_cat21];
    post_open_cat  = [post_open_cat12;  post_open_cat21];
    post_close_cat = [post_close_cat12; post_close_cat21];

    valid = isfinite(pre_open_cat) & isfinite(pre_close_cat) & ...
            isfinite(post_open_cat) & isfinite(post_close_cat);

    pre_open_cat   = pre_open_cat(valid);
    pre_close_cat  = pre_close_cat(valid);
    post_open_cat  = post_open_cat(valid);
    post_close_cat = post_close_cat(valid);
end

function code = compound_category_code_9(first_cat, second_cat)
    % Category order:
    %   ++, +o, +-, o+, oo, o-, -+, -o, --
    % first_cat and second_cat must use: +1 = positive, 0 = non-significant, -1 = negative.
    code = nan(size(first_cat));

    code(first_cat ==  1 & second_cat ==  1) = 1;
    code(first_cat ==  1 & second_cat ==  0) = 2;
    code(first_cat ==  1 & second_cat == -1) = 3;

    code(first_cat ==  0 & second_cat ==  1) = 4;
    code(first_cat ==  0 & second_cat ==  0) = 5;
    code(first_cat ==  0 & second_cat == -1) = 6;

    code(first_cat == -1 & second_cat ==  1) = 7;
    code(first_cat == -1 & second_cat ==  0) = 8;
    code(first_cat == -1 & second_cat == -1) = 9;
end

function code = compound_category_code_4(first_cat, second_cat)
    % Category order:
    %   ++, +-, -+, --
    % Non-significant members are returned as NaN and excluded from counts.
    code = nan(size(first_cat));

    code(first_cat ==  1 & second_cat ==  1) = 1;
    code(first_cat ==  1 & second_cat == -1) = 2;
    code(first_cat == -1 & second_cat ==  1) = 3;
    code(first_cat == -1 & second_cat == -1) = 4;
end

function counts = make_transition_counts(x_code, y_code, n_cat)
    valid = isfinite(x_code) & isfinite(y_code);
    x_code = x_code(valid);
    y_code = y_code(valid);

    counts = zeros(n_cat, n_cat);
    for i_cat = 1:n_cat
        for j_cat = 1:n_cat
            counts(i_cat, j_cat) = sum(x_code == i_cat & y_code == j_cat);
        end
    end
end

function render_pattern_transition_figure(root, cfg, pooled_one, params)
    f = figure('Color', 'w', 'Visible', params.figure_visible);
    tiledlayout(params.n_row, params.n_col, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    tile_idx = @(row, col) rowcol_to_panel_index(row, col, params.n_col);

    plot_modes = {'raw', 'row', 'column', 'total', 'oe_ratio', 'standardized_residual'};
    column_titles = {'raw count', 'normalized by row', 'normalized by column', ...
                     'normalized by total', 'O/E ratio', 'standardized residual'};
    colorbar_labels = {'count', 'row fraction', 'column fraction', ...
                       'total fraction', 'O/E ratio', 'std residual'};
    value_modes = {'integer', 'fraction', 'fraction', 'fraction', 'ratio', 'signed'};

    for col_i = 1:numel(plot_modes)
        mode = plot_modes{col_i};

        counts = pooled_one.counts9;
        [agreement, kappa, n_valid] = summarize_counts(counts);
        mat = transform_transition_matrix(counts, mode);
        ax = nexttile(tile_idx(1, col_i));
        plot_transition_table_generic(ax, mat, params.category_labels_9, agreement, kappa, n_valid, ...
            cfg.row1_title, cfg.x_label, cfg.y_label, column_titles{col_i}, colorbar_labels{col_i}, value_modes{col_i});
        add_panel_label(ax, 1, col_i, params.n_col);

        counts = pooled_one.counts4;
        [agreement, kappa, n_valid] = summarize_counts(counts);
        mat = transform_transition_matrix(counts, mode);
        ax = nexttile(tile_idx(2, col_i));
        plot_transition_table_generic(ax, mat, params.category_labels_4, agreement, kappa, n_valid, ...
            cfg.row2_title, cfg.x_label, cfg.y_label, column_titles{col_i}, colorbar_labels{col_i}, value_modes{col_i});
        add_panel_label(ax, 2, col_i, params.n_col);
    end

    sgtitle(sprintf('Supplement Fig: %s; valid sessions = %d', cfg.figure_title, pooled_one.valid_session_count), ...
        'Interpreter', 'none');

    save_folder = fullfile(root, 'Figures', 'Paper');
    check_path(save_folder);

    figWidth = 24.0;  % inches.
    figHeight = 8.5;  % inches.
    resolution = 300; % dpi.

    set(f, 'Units', 'inches');
    f.Position(3:4) = [figWidth, figHeight];
    set(f, 'PaperUnits', 'inches');
    set(f, 'PaperSize', [figWidth, figHeight]);
    set(f, 'PaperPosition', [0, 0, figWidth, figHeight]);
    set(f, 'Color', 'w');

    preview_filename = fullfile(save_folder, [cfg.output_stub, '_supp_preview.jpg']);
    exportgraphics(f, preview_filename, 'ContentType', 'image', 'BackgroundColor', 'white', 'Resolution', resolution);

    % pdf_filename = fullfile(save_folder, [cfg.output_stub, '_supp.pdf']);
    % exportgraphics(f, pdf_filename, 'ContentType', 'vector', 'BackgroundColor', 'white', 'Resolution', resolution);

    close(f);
end

function [agreement, kappa, n_valid] = summarize_counts(counts)
    agreement = compute_raw_agreement(counts);
    kappa = compute_cohen_kappa(counts);
    n_valid = sum(counts(:));
end

function mat_out = transform_transition_matrix(counts, mode)
    counts = double(counts);
    total_n = sum(counts(:));

    switch lower(mode)
        case 'raw'
            mat_out = counts;

        case 'row'
            denom = sum(counts, 2);
            mat_out = zeros(size(counts));
            valid = denom > 0;
            mat_out(valid, :) = bsxfun(@rdivide, counts(valid, :), denom(valid));

        case 'column'
            denom = sum(counts, 1);
            mat_out = zeros(size(counts));
            valid = denom > 0;
            mat_out(:, valid) = bsxfun(@rdivide, counts(:, valid), denom(valid));

        case 'total'
            if total_n > 0
                mat_out = counts / total_n;
            else
                mat_out = zeros(size(counts));
            end

        case 'oe_ratio'
            expected_prob = expected_probability_under_independence(counts);
            observed_prob = zeros(size(counts));
            if total_n > 0
                observed_prob = counts / total_n;
            end
            mat_out = nan(size(counts));
            valid = expected_prob > 0;
            mat_out(valid) = observed_prob(valid) ./ expected_prob(valid);

        case 'standardized_residual'
            expected_counts = expected_counts_under_independence(counts);
            mat_out = nan(size(counts));
            valid = expected_counts > 0;
            mat_out(valid) = (counts(valid) - expected_counts(valid)) ./ sqrt(expected_counts(valid));

        otherwise
            error('Unknown transform mode: %s', mode);
    end
end

function expected_prob = expected_probability_under_independence(counts)
    counts = double(counts);
    total_n = sum(counts(:));
    if total_n <= 0
        expected_prob = zeros(size(counts));
        return;
    end
    x_marginal_prob = sum(counts, 2) / total_n;
    y_marginal_prob = sum(counts, 1) / total_n;
    expected_prob = x_marginal_prob * y_marginal_prob;
end

function expected_counts = expected_counts_under_independence(counts)
    counts = double(counts);
    total_n = sum(counts(:));
    if total_n <= 0
        expected_counts = zeros(size(counts));
        return;
    end
    expected_counts = total_n * expected_probability_under_independence(counts);
end

function plot_transition_table_generic(ax, mat, category_labels, agreement, kappa, n_valid, title_text, xlabel_text, ylabel_text, normalization_label, colorbar_label, value_mode)
    plot_mat = mat.';
    imagesc(ax, plot_mat);
    set(ax, 'YDir', 'normal');
    axis(ax, 'square');

    finite_vals = plot_mat(isfinite(plot_mat));
    if isempty(finite_vals)
        caxis(ax, [0, 1]);
    elseif strcmp(value_mode, 'signed')
        max_abs = max(abs(finite_vals));
        if max_abs > 0
            caxis(ax, [-max_abs, max_abs]);
        end
    end

    n_cat = numel(category_labels);
    xticks(ax, 1:n_cat);
    yticks(ax, 1:n_cat);
    xticklabels(ax, category_labels);
    yticklabels(ax, category_labels);
    xtickangle(ax, 30);

    xlabel(ax, xlabel_text);
    ylabel(ax, ylabel_text);

    cb = colorbar(ax);
    ylabel(cb, colorbar_label);

    max_val = max(finite_vals);
    min_val = min(finite_vals);
    for y_idx = 1:n_cat
        for x_idx = 1:n_cat
            this_val = plot_mat(y_idx, x_idx);
            if ~isfinite(this_val)
                value_str = 'NaN';
                text_color = 'k';
            else
                this_ratio = (this_val - min_val) / (max_val - min_val + eps);
                if this_ratio < 0.2
                    text_color = 'w';
                else
                    text_color = 'k';
                end

                switch value_mode
                    case 'integer'
                        value_str = sprintf('%d', round(this_val));
                    case 'fraction'
                        value_str = sprintf('%.2f', this_val);
                    case 'ratio'
                        value_str = sprintf('%.2f', this_val);
                    case 'signed'
                        value_str = sprintf('%.2f', this_val);
                    otherwise
                        value_str = sprintf('%.2f', this_val);
                end
            end

            text(ax, x_idx, y_idx, value_str, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', 'Color', text_color, 'FontSize', choose_cell_font_size(n_cat));
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

    title(ax, sprintf('%s\n%s\nAgreement = %s, kappa = %s, n = %d', ...
        title_text, normalization_label, agreement_str, kappa_str, n_valid), 'Interpreter', 'none');
end

function font_size = choose_cell_font_size(n_cat)
    if n_cat <= 4
        font_size = 8;
    else
        font_size = 6;
    end
end

function add_panel_label(ax, row, col, n_col)
    panel_index = rowcol_to_panel_index(row, col, n_col);
    panel_label = panel_index_to_letters(panel_index);
    text(ax, -0.10, 1.08, panel_label, 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'Interpreter', 'none', 'Clipping', 'off');
end

function panel_index = rowcol_to_panel_index(row, col, n_col)
    panel_index = (row - 1) * n_col + col;
end

function label = panel_index_to_letters(panel_index)
    if panel_index < 1 || panel_index ~= floor(panel_index)
        error('panel_index must be a positive integer.');
    end
    label = '';
    n = panel_index;
    while n > 0
        rem0 = mod(n - 1, 26);
        label = [char(double('A') + rem0), label]; %#ok<AGROW>
        n = floor((n - 1) / 26);
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
    state_struct.kernel_idx = kernel_idx;

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
    col_marginal = sum(counts, 1) / total_n;
    p_expected = row_marginal.' * col_marginal.';

    if abs(1 - p_expected) < eps
        kappa = NaN;
    else
        kappa = (p_observed - p_expected) / (1 - p_expected);
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
