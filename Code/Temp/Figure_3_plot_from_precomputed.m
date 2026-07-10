%% Figure 3: pooled connection-pair correlation analyses
% Plot-only version.
% Loads precomputed MAT files and renders figures without reloading metadata
% or model data. Default: export JPG only, PDF export off.

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

%% Plot parameters
figure_visible = 'off';
show_legend = true;
show_inset_stat = false;
show_identity_line = false;
show_fit_line = true;
fit_line_method = 'tls';
export_jpg = true;
export_pdf = false;

base_plot_params = struct();
base_plot_params.n_row = 4;
base_plot_params.n_col = 4;
base_plot_params.figure_visible = figure_visible;
base_plot_params.show_legend = show_legend;
base_plot_params.show_inset_stat = show_inset_stat;
base_plot_params.show_identity_line = show_identity_line;
base_plot_params.show_fit_line = show_fit_line;
base_plot_params.fit_line_method = fit_line_method;
base_plot_params.export_jpg = export_jpg;
base_plot_params.export_pdf = export_pdf;

%% Hyperparameter configurations
hyperparam_configs = build_auto_hyperparam_configs();

cache_root = fullfile(root, 'Data', 'Working', 'FigureCache', 'Figure3_Main');
output_root = fullfile(root, 'Figures/Paper/fig3');
check_path(output_root);
log_filename = fullfile(output_root, 'run_issue_log_plot_from_precomputed.txt');
initialize_run_log(log_filename, numel(hyperparam_configs));

for hp_i = 1:numel(hyperparam_configs)
    hp = hyperparam_configs(hp_i);
    fprintf('\n\n##### Plot hyperparameter configuration %d/%d #####\n%s\n', hp_i, numel(hyperparam_configs), make_hyperparam_title(hp));
    try
        cache_file = fullfile(cache_root, make_hyperparam_folder_name(hp), 'precomputed_figure3.mat');
        if ~exist(cache_file, 'file')
            error('Missing precomputed MAT: %s', cache_file);
        end
        loaded = load(cache_file, 'precomputed');
        precomputed = loaded.precomputed;

        params = precomputed.params;
        params.n_row = base_plot_params.n_row;
        params.n_col = base_plot_params.n_col;
        params.figure_visible = base_plot_params.figure_visible;
        params.show_legend = base_plot_params.show_legend;
        params.show_inset_stat = base_plot_params.show_inset_stat;
        params.show_identity_line = base_plot_params.show_identity_line;
        params.show_fit_line = base_plot_params.show_fit_line;
        params.fit_line_method = base_plot_params.fit_line_method;
        params.export_jpg = base_plot_params.export_jpg;
        params.export_pdf = base_plot_params.export_pdf;
        params.output_folder = fullfile(output_root, make_hyperparam_folder_name(precomputed.hp));
        params.hyperparam_title = make_hyperparam_title(precomputed.hp);
        check_path(params.output_folder);

        figure_configs = precomputed.figure_configs;
        pooled = precomputed.pooled;
        for fig_i = 1:numel(figure_configs)
            render_comparison_figure(root, figure_configs(fig_i), pooled(fig_i), params);
        end
    catch ME
        log_hyperparam_issue(log_filename, hp_i, numel(hyperparam_configs), hp, ME);
        warning('Skipping plot hyperparameter configuration %d/%d after error: %s', hp_i, numel(hyperparam_configs), ME.message);
        continue;
    end
end

function run_one_hyperparam_config(root, mt, hp, params, base_output_folder)
    kernel_indices = 1:hp.kernel_num;
    figure_configs = build_figure_configs(kernel_indices);

    params.output_folder = fullfile(base_output_folder, make_hyperparam_folder_name(hp));
    params.hyperparam_title = make_hyperparam_title(hp);
    check_path(params.output_folder);

    %% Load and filter metadata
    selected_rows = metadata_filter(mt, hp);

    selected_mt = mt(selected_rows, :);
    meta_array = table2struct(selected_mt);
    if isfinite(params.max_sessions_to_include)
        meta_array = meta_array(1:min(numel(meta_array), params.max_sessions_to_include));
    end

    if isempty(meta_array)
        error('No metadata rows selected for this hyperparameter configuration.');
    end

    %% Initialize pooled outputs
    n_fig = numel(figure_configs);
    pooled = struct([]);
    for fig_i = 1:n_fig
        pooled(fig_i).contexts = struct([]);
        for ctx_i = 1:numel(figure_configs(fig_i).contexts)
            pooled(fig_i).contexts(ctx_i).data = empty_connection_data();
            pooled(fig_i).contexts(ctx_i).cat_counts = zeros(numel(params.category_labels), numel(params.category_labels));
        end
        pooled(fig_i).example_states = cell(4, 1);
        pooled(fig_i).example_label = '';
    end

    valid_session_count = 0;
    failed_session_count = 0;
    first_valid_loaded = [];
    first_valid_label = '';
    preferred_loaded = [];
    preferred_label = '';

    %% Pool data across selected sessions
    for session_i = 1:numel(meta_array)
        meta = meta_array(session_i);
        session_label = make_session_label(meta);
        fprintf('\n===== Loading %s (%d/%d) =====\n', session_label, session_i, numel(meta_array));

        try
            loaded_states = load_all_required_states(root, meta, kernel_indices);

            if isempty(first_valid_loaded)
                first_valid_loaded = loaded_states;
                first_valid_label = session_label;
            end
            if session_i == params.preferred_example_session_index
                preferred_loaded = loaded_states;
                preferred_label = session_label;
            end

            for fig_i = 1:n_fig
                cfg = figure_configs(fig_i);
                for ctx_i = 1:numel(cfg.contexts)
                    cond_x = cfg.contexts(ctx_i).x;
                    cond_y = cfg.contexts(ctx_i).y;

                    state_x = loaded_states.(state_key(cond_x.prepost, cond_x.state, cond_x.kernel_idx));
                    state_y = loaded_states.(state_key(cond_y.prepost, cond_y.state, cond_y.kernel_idx));

                    [pair_data, ~] = make_pair_vectors(state_x, state_y, params.err_multi);
                    pooled(fig_i).contexts(ctx_i).data = append_connection_data(pooled(fig_i).contexts(ctx_i).data, pair_data);

                    [pair_counts, ~, ~, ~] = make_pair_category_counts(state_x, state_y, params.err_multi);
                    pooled(fig_i).contexts(ctx_i).cat_counts = pooled(fig_i).contexts(ctx_i).cat_counts + pair_counts;
                end
            end

            valid_session_count = valid_session_count + 1;
        catch ME
            failed_session_count = failed_session_count + 1;
            if params.skip_failed_sessions
                warning('Skipping %s because loading/processing failed: %s', session_label, ME.message);
                continue;
            else
                rethrow(ME);
            end
        end
    end

    if valid_session_count == 0
        error('No valid sessions were loaded for this hyperparameter configuration.');
    end

    fprintf('\nValid sessions: %d. Failed sessions: %d.\n', valid_session_count, failed_session_count);

    if isempty(preferred_loaded)
        preferred_loaded = first_valid_loaded;
        preferred_label = first_valid_label;
    end

    for fig_i = 1:n_fig
        cfg = figure_configs(fig_i);
        pooled(fig_i).example_label = preferred_label;
        for ctx_i = 1:numel(cfg.contexts)
            row_top = (ctx_i - 1) * 2 + 1;
            row_bottom = row_top + 1;
            pooled(fig_i).example_states{row_top} = preferred_loaded.(state_key(cfg.contexts(ctx_i).x.prepost, cfg.contexts(ctx_i).x.state, cfg.contexts(ctx_i).x.kernel_idx));
            pooled(fig_i).example_states{row_bottom} = preferred_loaded.(state_key(cfg.contexts(ctx_i).y.prepost, cfg.contexts(ctx_i).y.state, cfg.contexts(ctx_i).y.kernel_idx));
        end
        pooled(fig_i).valid_session_count = valid_session_count;
    end

    %% Render and save figures
    for fig_i = 1:n_fig
        render_comparison_figure(root, figure_configs(fig_i), pooled(fig_i), params);
    end
end


function figure_configs = build_figure_configs(kernel_indices)
    figure_configs = struct([]);
    fig_i = 0;

    % Open vs Close for each kernel, split by Pre/Post.
    for k_i = 1:numel(kernel_indices)
        k = kernel_indices(k_i);

        x_pre_open   = make_condition('Pre',  'RestOpen',  k, 'Open');
        y_pre_close  = make_condition('Pre',  'RestClose', k, 'Close');
        x_post_open  = make_condition('Post', 'RestOpen',  k, 'Open');
        y_post_close = make_condition('Post', 'RestClose', k, 'Close');

        fig_i = fig_i + 1;
        figure_configs(fig_i).output_stub = sprintf('Figure3_k%d_open_vs_close', k);
        figure_configs(fig_i).figure_title = sprintf('Kernel %d: Open vs Close', k);
        figure_configs(fig_i).contexts = [ ...
            make_context('Pre',  x_pre_open,  y_pre_close), ...
            make_context('Post', x_post_open, y_post_close) ...
        ];
    end

    % Pre vs Post for each kernel, split by Open/Close.
    for k_i = 1:numel(kernel_indices)
        k = kernel_indices(k_i);

        x_open_pre   = make_condition('Pre',  'RestOpen',  k, 'Pre');
        y_open_post  = make_condition('Post', 'RestOpen',  k, 'Post');
        x_close_pre  = make_condition('Pre',  'RestClose', k, 'Pre');
        y_close_post = make_condition('Post', 'RestClose', k, 'Post');

        fig_i = fig_i + 1;
        figure_configs(fig_i).output_stub = sprintf('Figure3_k%d_pre_vs_post', k);
        figure_configs(fig_i).figure_title = sprintf('Kernel %d: Pre vs Post', k);
        figure_configs(fig_i).contexts = [ ...
            make_context('RestOpen',  x_open_pre,  y_open_post), ...
            make_context('RestClose', x_close_pre, y_close_post) ...
        ];
    end

    % Kernel-vs-kernel comparisons for every kernel pair, split by Open/Close.
    for a_i = 1:numel(kernel_indices)
        for b_i = (a_i + 1):numel(kernel_indices)
            k_a = kernel_indices(a_i);
            k_b = kernel_indices(b_i);

            x_pre_open_a  = make_condition('Pre', 'RestOpen',  k_a, sprintf('K%d', k_a));
            y_pre_open_b  = make_condition('Pre', 'RestOpen',  k_b, sprintf('K%d', k_b));
            x_pre_close_a = make_condition('Pre', 'RestClose', k_a, sprintf('K%d', k_a));
            y_pre_close_b = make_condition('Pre', 'RestClose', k_b, sprintf('K%d', k_b));

            fig_i = fig_i + 1;
            figure_configs(fig_i).output_stub = sprintf('Figure3_pre_k%d_vs_k%d', k_a, k_b);
            figure_configs(fig_i).figure_title = sprintf('Pre: Kernel %d vs Kernel %d', k_a, k_b);
            figure_configs(fig_i).contexts = [ ...
                make_context('RestOpen',  x_pre_open_a,  y_pre_open_b), ...
                make_context('RestClose', x_pre_close_a, y_pre_close_b) ...
            ];

            x_post_open_a  = make_condition('Post', 'RestOpen',  k_a, sprintf('K%d', k_a));
            y_post_open_b  = make_condition('Post', 'RestOpen',  k_b, sprintf('K%d', k_b));
            x_post_close_a = make_condition('Post', 'RestClose', k_a, sprintf('K%d', k_a));
            y_post_close_b = make_condition('Post', 'RestClose', k_b, sprintf('K%d', k_b));

            fig_i = fig_i + 1;
            figure_configs(fig_i).output_stub = sprintf('Figure3_post_k%d_vs_k%d', k_a, k_b);
            figure_configs(fig_i).figure_title = sprintf('Post: Kernel %d vs Kernel %d', k_a, k_b);
            figure_configs(fig_i).contexts = [ ...
                make_context('RestOpen',  x_post_open_a,  y_post_open_b), ...
                make_context('RestClose', x_post_close_a, y_post_close_b) ...
            ];
        end
    end
end

function cond = make_condition(prepost, state, kernel_idx, axis_label)
    cond = struct();
    cond.prepost = prepost;
    cond.state = state;
    cond.kernel_idx = kernel_idx;
    cond.axis_label = axis_label;
end

function ctx = make_context(name, x_cond, y_cond)
    ctx = struct();
    ctx.name = name;
    ctx.x = x_cond;
    ctx.y = y_cond;
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

function render_comparison_figure(root, cfg, pooled_one, params)
    f = figure('Color', 'w', 'Visible', params.figure_visible);
    tiledlayout(params.n_row, params.n_col, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    tile_idx = @(row, col) rowcol_to_panel_index(row, col, params.n_col);

    for ctx_i = 1:numel(cfg.contexts)
        ctx = cfg.contexts(ctx_i);
        row_top = (ctx_i - 1) * 2 + 1;
        row_bottom = row_top + 1;
        data = pooled_one.contexts(ctx_i).data;
        counts = pooled_one.contexts(ctx_i).cat_counts;

        [pearson_r_orig, pearson_p_orig, spearman_rho_orig, spearman_p_orig, n_orig] = correlation_stats(data.x_orig, data.y_orig);
        orig_axis_limit = get_symmetric_axis_limit(data.x_orig, data.y_orig);
        agreement = compute_raw_agreement(counts);
        kappa = compute_cohen_kappa(counts);
        n_cat = sum(counts(:));

        counts_2x2 = counts([1, 3], [1, 3]);
        agreement_2x2 = compute_raw_agreement(counts_2x2);
        kappa_2x2 = compute_cohen_kappa(counts_2x2);
        n_cat_2x2 = sum(counts_2x2(:));
        category_labels_2x2 = {'Negative', 'Positive'};

        ax = nexttile(tile_idx(row_top, 1));
        plot_pair_scatter(ax, data, 'all_signed', params.scatter_marker_size, params.scatter_alpha, ...
            params.colors, params.show_legend, sprintf('%s all connections, sessions=%d', ctx.name, pooled_one.valid_session_count), ...
            ctx.x.axis_label, ctx.y.axis_label, orig_axis_limit, ...
            params.show_inset_stat, params.show_identity_line, params.show_fit_line, params.fit_line_method);
        add_panel_label(ax, row_top, 1, params.n_col);

        ax = nexttile(tile_idx(row_bottom, 1));
        plot_pair_density(ax, data.x_orig, data.y_orig, pearson_r_orig, pearson_p_orig, spearman_rho_orig, spearman_p_orig, n_orig, params.density_nbin, params.density_clip_percentile, ...
            params.density_use_log_count, sprintf('%s all-connection density, sessions=%d', ctx.name, pooled_one.valid_session_count), ...
            'signed', params.colors, ctx.x.axis_label, ctx.y.axis_label, orig_axis_limit);
        add_panel_label(ax, row_bottom, 1, params.n_col);

        ax = nexttile(tile_idx(row_top, 2));
        plot_pair_scatter(ax, data, 'either_sig_signed', params.scatter_marker_size, params.scatter_alpha, ...
            params.colors, params.show_legend, sprintf('%s either significant signed, sessions=%d', ctx.name, pooled_one.valid_session_count), ...
            ctx.x.axis_label, ctx.y.axis_label, [], ...
            params.show_inset_stat, params.show_identity_line, params.show_fit_line, params.fit_line_method);
        add_panel_label(ax, row_top, 2, params.n_col);

        ax = nexttile(tile_idx(row_bottom, 2));
        plot_pair_scatter(ax, data, 'both_sig_signed', params.scatter_marker_size, params.scatter_alpha, ...
            params.colors, params.show_legend, sprintf('%s both significant signed, sessions=%d', ctx.name, pooled_one.valid_session_count), ...
            ctx.x.axis_label, ctx.y.axis_label, [], ...
            params.show_inset_stat, params.show_identity_line, params.show_fit_line, params.fit_line_method);
        add_panel_label(ax, row_bottom, 2, params.n_col);

        ax = nexttile(tile_idx(row_top, 3));
        plot_pair_scatter(ax, data, 'same_abs', params.scatter_marker_size, params.scatter_alpha, ...
            params.colors, params.show_legend, sprintf('%s same-sign abs, sessions=%d', ctx.name, pooled_one.valid_session_count), ...
            ctx.x.axis_label, ctx.y.axis_label, [], ...
            params.show_inset_stat, params.show_identity_line, params.show_fit_line, params.fit_line_method);
        add_panel_label(ax, row_top, 3, params.n_col);

        ax = nexttile(tile_idx(row_bottom, 3));
        plot_pair_scatter(ax, data, 'switch_abs', params.scatter_marker_size, params.scatter_alpha, ...
            params.colors, params.show_legend, sprintf('%s sign-switch abs, sessions=%d', ctx.name, pooled_one.valid_session_count), ...
            ctx.x.axis_label, ctx.y.axis_label, [], ...
            params.show_inset_stat, params.show_identity_line, params.show_fit_line, params.fit_line_method);
        add_panel_label(ax, row_bottom, 3, params.n_col);

        ax = nexttile(tile_idx(row_top, 4));
        plot_category_transition_table(ax, counts, params.category_labels, agreement, kappa, n_cat, ...
            sprintf('%s categorical 3x3, sessions=%d', ctx.name, pooled_one.valid_session_count), ...
            sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label));
        add_panel_label(ax, row_top, 4, params.n_col);

        ax = nexttile(tile_idx(row_bottom, 4));
        plot_category_transition_table(ax, counts_2x2, category_labels_2x2, agreement_2x2, kappa_2x2, n_cat_2x2, ...
            sprintf('%s categorical 2x2, sessions=%d', ctx.name, pooled_one.valid_session_count), ...
            sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label));
        add_panel_label(ax, row_bottom, 4, params.n_col);
    end

    if isfield(params, 'hyperparam_title') && ~isempty(params.hyperparam_title)
        sgtitle(sprintf('%s\n%s', cfg.figure_title, params.hyperparam_title), 'Interpreter', 'none');
    else
        sgtitle(cfg.figure_title, 'Interpreter', 'none');
    end

    if isfield(params, 'output_folder') && ~isempty(params.output_folder)
        save_folder = params.output_folder;
    else
        save_folder = fullfile(root, 'Figures', 'Paper', 'fig3');
    end
    check_path(save_folder);

    figWidth = 16.0;  % inches.
    figHeight = 16.0; % inches.
    resolution = 300;

    set(f, 'Units', 'inches');
    f.Position(3:4) = [figWidth, figHeight];
    set(f, 'PaperUnits', 'inches');
    set(f, 'PaperSize', [figWidth, figHeight]);
    set(f, 'PaperPosition', [0, 0, figWidth, figHeight]);
    set(f, 'Color', 'w');

    if ~isfield(params, 'export_jpg') || params.export_jpg
        preview_filename = fullfile(save_folder, [cfg.output_stub, '_preview.jpg']);
        exportgraphics(f, preview_filename, 'ContentType', 'image', 'BackgroundColor', 'white', 'Resolution', resolution);
    end

    if isfield(params, 'export_pdf') && params.export_pdf
        pdf_filename = fullfile(save_folder, [cfg.output_stub, '.pdf']);
        exportgraphics(f, pdf_filename, 'ContentType', 'vector', 'BackgroundColor', 'white', 'Resolution', resolution);
    end

    close(f);
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

function title_str = make_condition_title(cond)
    title_str = sprintf('%s, %s, K%d', cond.prepost, cond.state, cond.kernel_idx);
end

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

    data.xswitch_xy_abs = [];
    data.yswitch_xy_abs = [];
    data.xswitch_yx_abs = [];
    data.yswitch_yx_abs = [];
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

function selected_rows = metadata_filter(mt, hp)
    base_filter = metadata_field_equals(mt, 'kernel_name', hp.kernel_name) & ...
                  metadata_field_equals(mt, 'align', hp.align) & ...
                  metadata_field_equals(mt, 'area', "Cortex") & ...
                  metadata_field_equals(mt, 'injection', hp.injection) & ...
                  metadata_field_equals(mt, 'epoch', 3000) & ...
                  metadata_field_equals(mt, 'fold_idx', 0) & ...
                  metadata_field_equals(mt, 'shuffle_idx', 0) & ...
                  metadata_field_equals(mt, 'resting_dur_threshold', hp.resting_dur_threshold);

    if isfield(hp, 'reg_name') && ~isempty(hp.reg_name)
        base_filter = base_filter & metadata_field_equals(mt, 'reg_name', hp.reg_name);
    end

    anchor_filter = base_filter & ...
                    metadata_field_equals(mt, 'prepost', 'Pre') & ...
                    metadata_field_equals(mt, 'state', 'RestOpen');

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
            same_session = same_session & metadata_matches_anchor(mt, field, idx);
        end

        has_all_required = true;
        for q = 1:numel(required_preposts)
            has_this = any(same_session & ...
                           metadata_field_equals(mt, 'prepost', required_preposts{q}) & ...
                           metadata_field_equals(mt, 'state', required_states{q}));

            if ~has_this
                has_all_required = false;
                anchor_name = get_metadata_value_as_char(mt, 'file_name', idx);
                warning('Session %s is missing required Pre/Post x Open/Close combination: %s %s. Skipping this session.', ...
                    anchor_name, required_preposts{q}, required_states{q});
                break;
            end
        end

        if has_all_required
            selected_rows(idx) = true;
        end
    end

    fprintf('Metadata filter selected %d/%d rows for %s:
', ...
        sum(selected_rows), sum(base_filter), make_hyperparam_title(hp));
    for k = 1:height(mt)
        if selected_rows(k)
            fprintf('  %s
', get_metadata_value_as_char(mt, 'file_name', k));
        end
    end

    if ~any(selected_rows)
        warning('Metadata filter selected no complete Pre/Post x Open/Close sets for %s.', make_hyperparam_title(hp));
    end
end

function mask = metadata_field_equals(mt, field, target)
    if isempty(target)
        mask = true(height(mt), 1);
        return;
    end

    if ~ismember(field, mt.Properties.VariableNames)
        error('Metadata table is missing required field: %s', field);
    end

    col = mt.(field);
    mask = false(height(mt), 1);

    if iscell(col)
        for i = 1:numel(col)
            mask(i) = metadata_value_equals(col{i}, target);
        end
    else
        for i = 1:numel(col)
            mask(i) = metadata_value_equals(col(i), target);
        end
    end

    mask = mask(:);
end

function tf = metadata_value_equals(value, target)
    if iscell(value)
        if isempty(value)
            tf = isempty(target);
            return;
        end
        value = value{1};
    end

    if isempty(value)
        tf = isempty(target);
        return;
    end

    if isstring(target) || ischar(target)
        tf = strcmp(char(string(value)), char(string(target)));
    elseif isnumeric(target)
        if isnumeric(value) || islogical(value)
            tf = isequal(value, target) || (isscalar(value) && isscalar(target) && value == target);
        else
            tf = false;
        end
    else
        tf = isequal(value, target);
    end
end

function mask = metadata_matches_anchor(mt, field, anchor_idx)
    col = mt.(field);
    if iscell(col)
        anchor_value = col{anchor_idx};
    else
        anchor_value = col(anchor_idx);
    end
    mask = metadata_field_equals(mt, field, anchor_value);
end

function value_str = get_metadata_value_as_char(mt, field, idx)
    if ~ismember(field, mt.Properties.VariableNames)
        value_str = '';
        return;
    end

    value = mt.(field)(idx);
    if iscell(mt.(field))
        value = value{1};
    end
    value_str = value_to_char(value);
end


function hyperparam_configs = build_auto_hyperparam_configs()
    kernel_specs = build_kernel_specs();
    injections = ["Muscimol", "Saline"];
    aligns = ["Last", "Longest"];

    config_i = 0;
    hyperparam_configs = struct([]);

    for k_i = 1:numel(kernel_specs)
        for inj_i = 1:numel(injections)
            for align_i = 1:numel(aligns)
                config_i = config_i + 1;
                hyperparam_configs(config_i).kernel_name = kernel_specs(k_i).kernel_name;
                hyperparam_configs(config_i).kernel_num = kernel_specs(k_i).kernel_num;
                hyperparam_configs(config_i).reg_name = "L2=0_2";
                hyperparam_configs(config_i).resting_dur_threshold = 15;
                hyperparam_configs(config_i).injection = injections(inj_i);
                hyperparam_configs(config_i).align = aligns(align_i);
            end
        end
    end
end

function kernel_specs = build_kernel_specs()
    kernel_specs = struct([]);
    spec_i = 0;

    first_names = ["Exp5", "Exp40", "Exp100", "Exp40Mix", "Exp100Mix"];
    first_kernel_nums = [1, 1, 1, 2, 3];
    for i = 1:numel(first_names)
        spec_i = spec_i + 1;
        kernel_specs(spec_i).kernel_name = first_names(i);
        kernel_specs(spec_i).kernel_num = first_kernel_nums(i);
    end

    suffixes = ["5", "10", "20", "40", "80", "160", "320", "640", "1000"];
    prefixes = ["LongExp", "LongGaussC", "LongGDeriv"];
    for p_i = 1:numel(prefixes)
        for s_i = 1:numel(suffixes)
            spec_i = spec_i + 1;
            kernel_specs(spec_i).kernel_name = prefixes(p_i) + suffixes(s_i);
            kernel_specs(spec_i).kernel_num = 1;
        end
    end

    step_suffixes = ["0_5", "5_10", "10_20", "20_40", "40_80", ...
                     "80_160", "160_320", "320_640", "640_1000", "1000_3000"];
    for s_i = 1:numel(step_suffixes)
        spec_i = spec_i + 1;
        kernel_specs(spec_i).kernel_name = "LongStepB" + step_suffixes(s_i);
        kernel_specs(spec_i).kernel_num = 1;
    end
end

function title_str = make_hyperparam_title(hp)
    title_str = sprintf('kernel_name=%s; kernel_num=%d; reg_name=%s; resting_dur_threshold=%s; injection=%s; align=%s', ...
        value_to_char(hp.kernel_name), hp.kernel_num, value_to_char(hp.reg_name), ...
        value_to_char(hp.resting_dur_threshold), value_to_char(hp.injection), value_to_char(hp.align));
end

function folder_name = make_hyperparam_folder_name(hp)
    folder_name = sprintf('kernel_%s__kernelNum_%d__reg_%s__restDur_%s__inj_%s__align_%s', ...
        value_to_char(hp.kernel_name), hp.kernel_num, value_to_char(hp.reg_name), ...
        value_to_char(hp.resting_dur_threshold), value_to_char(hp.injection), value_to_char(hp.align));
    folder_name = regexprep(folder_name, '[^\w=\.-]', '_');
end

function out = value_to_char(value)
    if iscell(value)
        if isempty(value)
            out = '[]';
        else
            out = value_to_char(value{1});
        end
    elseif isempty(value)
        out = '[]';
    elseif isstring(value) || ischar(value)
        out = char(string(value));
    elseif isnumeric(value) || islogical(value)
        out = strrep(mat2str(value), ' ', '_');
    else
        out = char(string(value));
    end
end

function initialize_run_log(log_filename, n_configs)
    fid = fopen(log_filename, 'w');
    if fid < 0
        warning('Could not initialize log file: %s', log_filename);
        return;
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, 'Figure 3 batch run issue log\n');
    fprintf(fid, 'Total hyperparameter configurations: %d\n', n_configs);
    fprintf(fid, 'Started: %s\n\n', datestr(now));
end

function log_hyperparam_issue(log_filename, hp_i, n_configs, hp, ME)
    fid = fopen(log_filename, 'a');
    if fid < 0
        warning('Could not write to log file: %s', log_filename);
        return;
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '----- Hyperparameter configuration %d/%d failed -----\n', hp_i, n_configs);
    fprintf(fid, '%s\n', make_hyperparam_title(hp));
    fprintf(fid, 'Error: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf(fid, 'Stack:\n');
        for s_i = 1:numel(ME.stack)
            fprintf(fid, '  %s, line %d, function %s\n', ME.stack(s_i).file, ME.stack(s_i).line, ME.stack(s_i).name);
        end
    end
    fprintf(fid, '\n');
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

function [data, label_text] = make_pair_vectors(state_x, state_y, err_multi)
    validate_matching_filters(state_x, state_y);

    x12 = state_x.J12(:);
    y12 = state_y.J12(:);
    x21 = state_x.J21(:);
    y21 = state_y.J21(:);
    xerr12 = state_x.err12(:);
    yerr12 = state_y.err12(:);
    xerr21 = state_x.err21(:);
    yerr21 = state_y.err21(:);

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
    switch_xy = x_pos & y_neg;
    switch_yx = x_neg & y_pos;
    switch_any = switch_xy | switch_yx;

    either_sig = x_sig | y_sig;
    both_sig = x_sig & y_sig;

    data = struct();
    data.x_orig = x;
    data.y_orig = y;
    data.x_err_orig = xerr;
    data.y_err_orig = yerr;
    data.x_cat_orig = x_cat;
    data.y_cat_orig = y_cat;

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

    data.xswitch_xy_abs = abs(x(switch_xy));
    data.yswitch_xy_abs = abs(y(switch_xy));
    data.xswitch_yx_abs = abs(x(switch_yx));
    data.yswitch_yx_abs = abs(y(switch_yx));
    data.xswitch_abs = abs(x(switch_any));
    data.yswitch_abs = abs(y(switch_any));

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

function [pearson_r, pearson_p, spearman_rho, spearman_p, n_valid] = correlation_stats(x, y)
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    n_valid = numel(x);

    pearson_r = NaN;
    pearson_p = NaN;
    spearman_rho = NaN;
    spearman_p = NaN;

    if n_valid < 2
        return;
    end

    if all(x == x(1)) || all(y == y(1))
        return;
    end

    [pearson_r, pearson_p] = corr(x, y, 'Type', 'Pearson', 'Rows', 'complete');
    [spearman_rho, spearman_p] = corr(x, y, 'Type', 'Spearman', 'Rows', 'complete');
end

function [rho, pval, n_valid] = pearson_stats(x, y)
    % Compatibility wrapper.
    [rho, pval, ~, ~, n_valid] = correlation_stats(x, y);
end


function plot_pair_scatter(ax, data, plot_mode, marker_size, marker_alpha, colors, show_legend, title_text, x_label, y_label, axis_limit_override, show_inset_stat, show_identity_line, show_fit_line, fit_line_method)
    cla(ax);
    hold(ax, 'on');

    if nargin < 11
        axis_limit_override = [];
    end
    if nargin < 12
        show_inset_stat = true;
    end
    if nargin < 13
        show_identity_line = true;
    end
    if nargin < 14
        show_fit_line = false;
    end
    if nargin < 15
        fit_line_method = 'ols';
    end

    switch plot_mode
        case 'all_signed'
            x = data.x_orig;
            y = data.y_orig;
            x_sig = data.x_cat_orig ~= 0;
            y_sig = data.y_cat_orig ~= 0;

            none_sig = ~x_sig & ~y_sig;
            x_only = x_sig & ~y_sig;
            y_only = ~x_sig & y_sig;
            both = x_sig & y_sig;

            scatter(ax, x(none_sig), y(none_sig), marker_size, 'filled', ...
                'MarkerFaceColor', colors.non_sig, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'none sig');
            scatter(ax, x(x_only), y(x_only), marker_size, 'filled', ...
                'MarkerFaceColor', colors.x_only, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s only', x_label));
            scatter(ax, x(y_only), y(y_only), marker_size, 'filled', ...
                'MarkerFaceColor', colors.y_only, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s only', y_label));
            scatter(ax, x(both), y(both), marker_size, 'filled', ...
                'MarkerFaceColor', colors.both_sig, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'both sig');

            axis_limit = get_symmetric_axis_limit(x, y);
            xlabel_text = sprintf('%s J_{ij}', x_label);
            ylabel_text = sprintf('%s J_{ij}', y_label);
            axis_mode = 'signed';

        case 'either_sig_signed'
            x = data.x_either_sig;
            y = data.y_either_sig;
            x_sig = data.x_cat_either_sig ~= 0;
            y_sig = data.y_cat_either_sig ~= 0;

            x_only = x_sig & ~y_sig;
            y_only = ~x_sig & y_sig;
            both = x_sig & y_sig;

            scatter(ax, x(x_only), y(x_only), marker_size, 'filled', ...
                'MarkerFaceColor', colors.x_only, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s only', x_label));
            scatter(ax, x(y_only), y(y_only), marker_size, 'filled', ...
                'MarkerFaceColor', colors.y_only, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s only', y_label));
            scatter(ax, x(both), y(both), marker_size, 'filled', ...
                'MarkerFaceColor', colors.both_sig, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'both sig');

            axis_limit = get_symmetric_axis_limit(x, y);
            xlabel_text = sprintf('%s J_{ij}', x_label);
            ylabel_text = sprintf('%s J_{ij}', y_label);
            axis_mode = 'signed';

        case 'both_sig_signed'
            x = data.x_both_sig;
            y = data.y_both_sig;
            x_cat = data.x_cat_both_sig;
            y_cat = data.y_cat_both_sig;

            pos_mask = x_cat == 1 & y_cat == 1;
            neg_mask = x_cat == -1 & y_cat == -1;
            switch_mask = (x_cat == 1 & y_cat == -1) | (x_cat == -1 & y_cat == 1);

            scatter(ax, x(pos_mask), y(pos_mask), marker_size, 'filled', ...
                'MarkerFaceColor', colors.pos, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'pos');
            scatter(ax, x(neg_mask), y(neg_mask), marker_size, 'filled', ...
                'MarkerFaceColor', colors.neg, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'neg');
            scatter(ax, x(switch_mask), y(switch_mask), marker_size, 'filled', ...
                'MarkerFaceColor', colors.switch, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'switch');

            axis_limit = get_symmetric_axis_limit(x, y);
            xlabel_text = sprintf('%s J_{ij}', x_label);
            ylabel_text = sprintf('%s J_{ij}', y_label);
            axis_mode = 'signed';

        case 'same_abs'
            x = data.x_abs;
            y = data.y_abs;

            scatter(ax, data.xpos_abs, data.ypos_abs, marker_size, 'filled', ...
                'MarkerFaceColor', colors.pos, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'pos');
            scatter(ax, data.xneg_abs, data.yneg_abs, marker_size, 'filled', ...
                'MarkerFaceColor', colors.neg, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'neg');

            axis_limit = [0, max(3.5, get_positive_axis_max(x, y))];
            xlabel_text = sprintf('%s |J_{ij}|', x_label);
            ylabel_text = sprintf('%s |J_{ij}|', y_label);
            axis_mode = 'positive';

        case 'switch_abs'
            x = data.xswitch_abs;
            y = data.yswitch_abs;

            scatter(ax, data.xswitch_xy_abs, data.yswitch_xy_abs, marker_size, 'filled', ...
                'MarkerFaceColor', colors.switch_xy, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s pos / %s neg', x_label, y_label));
            scatter(ax, data.xswitch_yx_abs, data.yswitch_yx_abs, marker_size, 'filled', ...
                'MarkerFaceColor', colors.switch_yx, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s neg / %s pos', x_label, y_label));

            axis_limit = [0, max(3.5, get_positive_axis_max(x, y))];
            xlabel_text = sprintf('%s |J_{ij}|', x_label);
            ylabel_text = sprintf('%s |J_{ij}|', y_label);
            axis_mode = 'positive';

        otherwise
            error('Unknown plot_mode: %s', plot_mode);
    end

    if ~isempty(axis_limit_override)
        axis_limit = axis_limit_override;
    end

    [pearson_r, pearson_p, spearman_rho, spearman_p, n_valid] = correlation_stats(x, y);
    cos_sim = normalized_cosine_similarity_omitnan(x, y);

    if show_identity_line
        plot(ax, axis_limit, axis_limit, '--', 'Color', colors.identity_line, 'LineWidth', 1, 'DisplayName', 'x=y');
    end
    if show_fit_line
        plot_linear_fit_line(ax, x, y, axis_limit, colors.fit_line, fit_line_method);
    end
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

    if pearson_p > 0.001
        pearson_p_str = sprintf('%.3f', pearson_p);
    else
        pearson_p_str = sprintf('%.3e', pearson_p);
    end
    if spearman_p > 0.001
        spearman_p_str = sprintf('%.3f', spearman_p);
    else
        spearman_p_str = sprintf('%.3e', spearman_p);
    end

    title(ax, sprintf('%s\nPearson r = %.6f (p = %s, n = %d)\nSpearman rho = %.6f (p = %s)\ncos sim = %.6f', ...
        title_text, pearson_r, pearson_p_str, n_valid, spearman_rho, spearman_p_str, cos_sim), 'Interpreter', 'none');

    if show_legend
        legend(ax, 'Location', 'northeastoutside');
    else
        legend(ax, 'off');
    end
    if show_inset_stat
        add_stats_text(ax, pearson_r, pearson_p, spearman_rho, spearman_p, n_valid);
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
            error('Unknown fit_line_method: %s. Use ''ols'' or ''tls''.', fit_line_method);
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
    x_min = axis_limit(1);
    x_max = axis_limit(2);
    y_min = axis_limit(1);
    y_max = axis_limit(2);

    candidates = [];
    if abs(direction(1)) > eps
        t = (x_min - center(1)) / direction(1);
        y_hit = center(2) + t * direction(2);
        if y_hit >= y_min - eps && y_hit <= y_max + eps
            candidates(end+1, :) = [x_min, y_hit]; %#ok<AGROW>
        end
        t = (x_max - center(1)) / direction(1);
        y_hit = center(2) + t * direction(2);
        if y_hit >= y_min - eps && y_hit <= y_max + eps
            candidates(end+1, :) = [x_max, y_hit]; %#ok<AGROW>
        end
    end

    if abs(direction(2)) > eps
        t = (y_min - center(2)) / direction(2);
        x_hit = center(1) + t * direction(1);
        if x_hit >= x_min - eps && x_hit <= x_max + eps
            candidates(end+1, :) = [x_hit, y_min]; %#ok<AGROW>
        end
        t = (y_max - center(2)) / direction(2);
        x_hit = center(1) + t * direction(1);
        if x_hit >= x_min - eps && x_hit <= x_max + eps
            candidates(end+1, :) = [x_hit, y_max]; %#ok<AGROW>
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

function cos_sim = normalized_cosine_similarity_omitnan(x, y)
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    x = x - mean(x);
    y = y - mean(y);
    if isempty(x) || norm(x) == 0 || norm(y) == 0
        cos_sim = NaN;
    else
        cos_sim = dot(x, y) / (norm(x) * norm(y));
    end
end

function plot_pair_density(ax, x, y, pearson_r, pearson_p, spearman_rho, spearman_p, n_valid, nbin, clip_percentile, use_log_count, title_text, axis_mode, colors, x_label, y_label, axis_limit_override)
    if nargin < 12 || isempty(title_text)
        title_text = '';
    end
    if nargin < 13 || isempty(axis_mode)
        axis_mode = 'abs';
    end
    if nargin < 14 || isempty(colors)
        colors = struct();
        colors.identity_line = [1, 0, 0];
        colors.zero_line = [0, 0, 0];
    end
    if nargin < 17
        axis_limit_override = [];
    end

    [edges_x, edges_y, n_in_range] = make_density_edges(x, y, nbin, clip_percentile, axis_limit_override);
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
        xlabel(ax, sprintf('%s J_{ij}', x_label));
        ylabel(ax, sprintf('%s J_{ij}', y_label));
    else
        xlabel(ax, sprintf('%s |J_{ij}|', x_label));
        ylabel(ax, sprintf('%s |J_{ij}|', y_label));
    end
    title(ax, title_text, 'Interpreter', 'none');
    add_stats_text(ax, pearson_r, pearson_p, spearman_rho, spearman_p, n_valid);

    fprintf('%s density: total n = %d, in density range = %d, max bin count = %d\n', ...
        title_text, n_valid, n_in_range, max(count_mat(:)));
end

function [edges_x, edges_y, n_in_range] = make_density_edges(x, y, nbin, clip_percentile, axis_limit_override)
    valid = isfinite(x) & isfinite(y);
    xv = x(valid);
    yv = y(valid);

    if isempty(xv) || isempty(yv)
        edges_x = linspace(-1, 1, nbin + 1);
        edges_y = linspace(-1, 1, nbin + 1);
        n_in_range = 0;
        return;
    end

    if ~isempty(axis_limit_override)
        x_limits = axis_limit_override;
        y_limits = axis_limit_override;
    else
        x_limits = local_percentile(xv, clip_percentile);
        y_limits = local_percentile(yv, clip_percentile);
    end;

    if ~all(isfinite(x_limits)) || x_limits(1) == x_limits(2)
        delta = max(abs(xv));
        if delta == 0, delta = 1; end
        x_limits = [-delta, delta];
    end
    if ~all(isfinite(y_limits)) || y_limits(1) == y_limits(2)
        delta = max(abs(yv));
        if delta == 0, delta = 1; end
        y_limits = [-delta, delta];
    end

    edges_x = linspace(x_limits(1), x_limits(2), nbin + 1);
    edges_y = linspace(y_limits(1), y_limits(2), nbin + 1);

    in_range = xv >= edges_x(1) & xv <= edges_x(end) & yv >= edges_y(1) & yv <= edges_y(end);
    n_in_range = sum(in_range);
end

function q = local_percentile(vals, pct)
    vals = vals(isfinite(vals));
    if isempty(vals)
        q = [NaN, NaN];
        return;
    end
    vals = sort(vals(:));
    pct = pct(:).';
    q = nan(size(pct));
    n = numel(vals);
    for i = 1:numel(pct)
        p = pct(i);
        if p <= 0
            q(i) = vals(1);
        elseif p >= 100
            q(i) = vals(end);
        else
            idx = 1 + (n - 1) * p / 100;
            lo = floor(idx);
            hi = ceil(idx);
            if lo == hi
                q(i) = vals(lo);
            else
                q(i) = vals(lo) + (idx - lo) * (vals(hi) - vals(lo));
            end
        end
    end
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

function add_stats_text(ax, pearson_r, pearson_p, spearman_rho, spearman_p, n_valid)
    pearson_r_str = format_stat_value(pearson_r, '%.4f');
    pearson_p_str = format_stat_value(pearson_p, '%.3e');
    spearman_rho_str = format_stat_value(spearman_rho, '%.4f');
    spearman_p_str = format_stat_value(spearman_p, '%.3e');

    text(ax, 0.02, 0.98, sprintf('Pearson r = %s\np = %s\nSpearman \\rho = %s\np = %s\nn = %d', ...
        pearson_r_str, pearson_p_str, spearman_rho_str, spearman_p_str, n_valid), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'w', 'Margin', 2, 'FontSize', 8, 'Interpreter', 'tex');
end

function stat_str = format_stat_value(value, fmt)
    if isnan(value)
        stat_str = 'NaN';
    else
        stat_str = sprintf(fmt, value);
    end
end

function plot_category_transition_table(ax, counts, category_labels, agreement, kappa, n_valid, title_text, xlabel_text, ylabel_text)
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

    xlabel(ax, xlabel_text);
    ylabel(ax, ylabel_text);

    cb = colorbar(ax);
    ylabel(cb, 'count');

    max_count = max(plot_counts(:));
    min_count = min(plot_counts(:));
    for y_idx = 1:n_cat
        for x_idx = 1:n_cat
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
