%% Supplement Figure 3: Combined pooled comparison script (auto hyperparameter version)
% This script generates pooled categorical comparison figures for each
% hyperparameter configuration and saves each figure as PDF + JPG.
%
% Figure versions for each hyperparameter configuration:
%   For each kernel index:
%     1. Kernel k, Open vs Close, split by Pre/Post.
%     2. Kernel k, Pre vs Post, split by Open/Close.
%   For each kernel-index pair:
%     3. Pre, Kernel ki vs Kernel kj, split by Open/Close.
%     4. Post, Kernel ki vs Kernel kj, split by Open/Close.
%
% Output organization:
%   Figures/Paper/Fig3_Supp/<hyperparameter-folder>/
%
% Layout for every saved figure:
%   Row 1 = context A, 3x3 tables.
%   Row 2 = context A, 2x2 tables.
%   Row 3 = context B, 3x3 tables.
%   Row 4 = context B, 2x2 tables.
%   Columns = raw count, row-normalized, column-normalized, total-normalized, O/E ratio, standardized residual.

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

%% Hyperparameter configurations
% Automatically generate all requested kernel_name x injection x align combinations.
% Fixed filters for this batch:
%   reg_name = "L2=0_2"
%   resting_dur_threshold = 15
hyperparam_configs = build_auto_hyperparam_configs();

%% Plot and analysis parameters
err_multi = 1; % threshold for significant J, in multiples of the GLM error estimate.
network_err_multi = 2;
density_nbin = 60;
scatter_marker_size = 8;
scatter_alpha = 0.25;
density_clip_percentile = [0.5, 99.5];
density_use_log_count = true;
category_labels = {'Negative', 'Non-sig', 'Positive'};

n_row = 4;
n_col = 6;
figure_visible = 'off';
show_legend = true; % true or false, applied to all scatter-plot legends.

skip_failed_sessions = false;
max_sessions_to_include = inf; % set smaller while debugging.

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

base_params = struct();
base_params.err_multi = err_multi;
base_params.network_err_multi = network_err_multi;
base_params.density_nbin = density_nbin;
base_params.scatter_marker_size = scatter_marker_size;
base_params.scatter_alpha = scatter_alpha;
base_params.density_clip_percentile = density_clip_percentile;
base_params.density_use_log_count = density_use_log_count;
base_params.category_labels = category_labels;
base_params.n_row = n_row;
base_params.n_col = n_col;
base_params.figure_visible = figure_visible;
base_params.show_legend = show_legend;
base_params.colors = colors;
base_params.skip_failed_sessions = skip_failed_sessions;
base_params.max_sessions_to_include = max_sessions_to_include;

%% Load metadata
mt = load_meta(root, 'table');
mt = mt.GLM;

fig3_root_folder = fullfile(root, 'Figures', 'Paper', 'Fig3_Supp');
check_path(fig3_root_folder);
run_issue_log = empty_run_issue_log();
run_issue_log_filename = fullfile(fig3_root_folder, 'run_issue_log_category_pairwise.txt');

total_rendered_figure_count = 0;

%% Run each hyperparameter configuration
for hp_i = 1:numel(hyperparam_configs)
    hp_title_for_log = sprintf('hyperparameter set %d/%d', hp_i, numel(hyperparam_configs));

    try
        hp = fill_default_hyperparams(hyperparam_configs(hp_i));
        hp_title_for_log = make_hyperparam_title(hp);
        kernel_indices = 1:hp.kernel_num;

        params = base_params;
        params.hyperparams = hp;
        params.hyperparam_title = hp_title_for_log;
        params.output_folder = fullfile(fig3_root_folder, make_hyperparam_folder_name(hp));

        fprintf('\n===== Hyperparameter set %d/%d: %s =====\n', hp_i, numel(hyperparam_configs), params.hyperparam_title);
        fprintf('Output folder: %s\n', params.output_folder);

        figure_configs = build_figure_configs(kernel_indices);

        %% Load and filter metadata for this hyperparameter set
        selected_rows = metadata_filter(mt, hp);

        selected_mt = mt(selected_rows, :);
        meta_array = table2struct(selected_mt);
        if isfinite(params.max_sessions_to_include)
            meta_array = meta_array(1:min(numel(meta_array), params.max_sessions_to_include));
        end

        if isempty(meta_array)
            message = 'No metadata rows selected. Skipping this hyperparameter set.';
            warning('%s Hyperparameter set %d: %s', message, hp_i, params.hyperparam_title);
            run_issue_log = append_run_issue(run_issue_log, hp_i, params.hyperparam_title, 'NoMetadata', message);
            write_run_issue_log(run_issue_log_filename, run_issue_log);
            continue;
        end

        %% Initialize pooled storage
        n_fig = numel(figure_configs);
        pooled = struct([]);
        for fig_i = 1:n_fig
            pooled(fig_i).contexts = struct([]);
            for ctx_i = 1:numel(figure_configs(fig_i).contexts)
                pooled(fig_i).contexts(ctx_i).data = empty_connection_data();
                pooled(fig_i).contexts(ctx_i).cat_counts = zeros(numel(category_labels), numel(category_labels));
            end
            pooled(fig_i).valid_session_count = 0;
        end

        valid_session_count = 0;
        failed_session_count = 0;

        %% Pool data across all selected sessions for this hyperparameter set
        for session_i = 1:numel(meta_array)
            meta = meta_array(session_i);
            session_label = make_session_label(meta);
            fprintf('\n===== Loading %s (%d/%d) =====\n', session_label, session_i, numel(meta_array));

            try
                loaded_states = load_all_required_states(root, meta, kernel_indices);

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
            catch ME_session
                failed_session_count = failed_session_count + 1;
                if params.skip_failed_sessions
                    warning('Skipping %s because loading/processing failed: %s', session_label, ME_session.message);
                    continue;
                else
                    rethrow(ME_session);
                end
            end
        end

        if valid_session_count == 0
            message = sprintf('No valid sessions were loaded. Failed sessions: %d. Skipping rendering for this hyperparameter set.', failed_session_count);
            warning('%s Hyperparameter set %d: %s', message, hp_i, params.hyperparam_title);
            run_issue_log = append_run_issue(run_issue_log, hp_i, params.hyperparam_title, 'NoValidSessions', message);
            write_run_issue_log(run_issue_log_filename, run_issue_log);
            continue;
        end

        fprintf('\nValid sessions: %d. Failed sessions: %d.\n', valid_session_count, failed_session_count);

        for fig_i = 1:n_fig
            pooled(fig_i).valid_session_count = valid_session_count;
        end

        %% Render and save figures for this hyperparameter set
        for fig_i = 1:n_fig
            render_comparison_figure(root, figure_configs(fig_i), pooled(fig_i), params);
            total_rendered_figure_count = total_rendered_figure_count + 1;
        end

    catch ME_hp
        message = sprintf('%s: %s', ME_hp.identifier, ME_hp.message);
        warning('Skipping hyperparameter set %d/%d after error: %s\n%s', ...
            hp_i, numel(hyperparam_configs), hp_title_for_log, message);
        run_issue_log = append_run_issue(run_issue_log, hp_i, hp_title_for_log, 'Error', message);
        write_run_issue_log(run_issue_log_filename, run_issue_log);
        continue;
    end
end

if total_rendered_figure_count == 0
    warning('No figures were rendered. Check hyperparameter values, metadata coverage, and run_issue_log_category_pairwise.txt.');
else
    fprintf('\nRendered %d figures across all hyperparameter sets.\n', total_rendered_figure_count);
end

if ~isempty(run_issue_log)
    write_run_issue_log(run_issue_log_filename, run_issue_log);
    fprintf('Run issue log written to: %s\n', run_issue_log_filename);
end

function figure_configs = build_figure_configs(kernel_indices)
    figure_configs = struct([]);
    fig_i = 0;

    % Per-kernel Open vs Close comparison, split by Pre/Post.
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

    % Per-kernel Pre vs Post comparison, split by Open/Close.
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

    % Cross-kernel comparisons. For kernel_num = 2, this reproduces the
    % original K1 vs K2 Pre and Post figures. For kernel_num > 2, it uses all
    % unique kernel-index pairs.
    for i = 1:numel(kernel_indices)
        for j = (i + 1):numel(kernel_indices)
            kx = kernel_indices(i);
            ky = kernel_indices(j);

            x_pre_open  = make_condition('Pre',  'RestOpen',  kx, sprintf('K%d', kx));
            y_pre_open  = make_condition('Pre',  'RestOpen',  ky, sprintf('K%d', ky));
            x_pre_close = make_condition('Pre',  'RestClose', kx, sprintf('K%d', kx));
            y_pre_close = make_condition('Pre',  'RestClose', ky, sprintf('K%d', ky));

            fig_i = fig_i + 1;
            figure_configs(fig_i).output_stub = sprintf('Figure3_pre_k%d_vs_k%d', kx, ky);
            figure_configs(fig_i).figure_title = sprintf('Pre: Kernel %d vs Kernel %d', kx, ky);
            figure_configs(fig_i).contexts = [ ...
                make_context('RestOpen',  x_pre_open,  y_pre_open), ...
                make_context('RestClose', x_pre_close, y_pre_close) ...
            ];

            x_post_open  = make_condition('Post', 'RestOpen',  kx, sprintf('K%d', kx));
            y_post_open  = make_condition('Post', 'RestOpen',  ky, sprintf('K%d', ky));
            x_post_close = make_condition('Post', 'RestClose', kx, sprintf('K%d', kx));
            y_post_close = make_condition('Post', 'RestClose', ky, sprintf('K%d', ky));

            fig_i = fig_i + 1;
            figure_configs(fig_i).output_stub = sprintf('Figure3_post_k%d_vs_k%d', kx, ky);
            figure_configs(fig_i).figure_title = sprintf('Post: Kernel %d vs Kernel %d', kx, ky);
            figure_configs(fig_i).contexts = [ ...
                make_context('RestOpen',  x_post_open,  y_post_open), ...
                make_context('RestClose', x_post_close, y_post_close) ...
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
        counts = pooled_one.contexts(ctx_i).cat_counts;

        agreement = compute_raw_agreement(counts);
        kappa = compute_cohen_kappa(counts);
        n_cat = sum(counts(:));

        counts_2x2 = counts([1, 3], [1, 3]);
        agreement_2x2 = compute_raw_agreement(counts_2x2);
        kappa_2x2 = compute_cohen_kappa(counts_2x2);
        n_cat_2x2 = sum(counts_2x2(:));
        category_labels_2x2 = {'Negative', 'Positive'};

        counts_row_norm = normalize_transition_matrix(counts, 'row');
        counts_col_norm = normalize_transition_matrix(counts, 'column');
        counts_total_norm = normalize_transition_matrix(counts, 'total');
        counts_oe_ratio = normalize_transition_matrix(counts, 'oe_ratio');
        counts_standardized_residual = normalize_transition_matrix(counts, 'standardized_residual');

        counts_2x2_row_norm = normalize_transition_matrix(counts_2x2, 'row');
        counts_2x2_col_norm = normalize_transition_matrix(counts_2x2, 'column');
        counts_2x2_total_norm = normalize_transition_matrix(counts_2x2, 'total');
        counts_2x2_oe_ratio = normalize_transition_matrix(counts_2x2, 'oe_ratio');
        counts_2x2_standardized_residual = normalize_transition_matrix(counts_2x2, 'standardized_residual');

        % Column 1: raw counts
        ax = nexttile(tile_idx(row_top, 1));
        plot_transition_table_generic(ax, counts, params.category_labels, agreement, kappa, n_cat, ...
            sprintf('%s categorical 3x3', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'raw count', 'count', 'integer');
        add_panel_label(ax, row_top, 1, params.n_col);

        ax = nexttile(tile_idx(row_bottom, 1));
        plot_transition_table_generic(ax, counts_2x2, category_labels_2x2, agreement_2x2, kappa_2x2, n_cat_2x2, ...
            sprintf('%s categorical 2x2', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'raw count', 'count', 'integer');
        add_panel_label(ax, row_bottom, 1, params.n_col);

        % Column 2: normalized by row
        ax = nexttile(tile_idx(row_top, 2));
        plot_transition_table_generic(ax, counts_row_norm, params.category_labels, agreement, kappa, n_cat, ...
            sprintf('%s categorical 3x3', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'normalized by row', 'row fraction', 'fraction');
        add_panel_label(ax, row_top, 2, params.n_col);

        ax = nexttile(tile_idx(row_bottom, 2));
        plot_transition_table_generic(ax, counts_2x2_row_norm, category_labels_2x2, agreement_2x2, kappa_2x2, n_cat_2x2, ...
            sprintf('%s categorical 2x2', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'normalized by row', 'row fraction', 'fraction');
        add_panel_label(ax, row_bottom, 2, params.n_col);

        % Column 3: normalized by column
        ax = nexttile(tile_idx(row_top, 3));
        plot_transition_table_generic(ax, counts_col_norm, params.category_labels, agreement, kappa, n_cat, ...
            sprintf('%s categorical 3x3', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'normalized by column', 'column fraction', 'fraction');
        add_panel_label(ax, row_top, 3, params.n_col);

        ax = nexttile(tile_idx(row_bottom, 3));
        plot_transition_table_generic(ax, counts_2x2_col_norm, category_labels_2x2, agreement_2x2, kappa_2x2, n_cat_2x2, ...
            sprintf('%s categorical 2x2', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'normalized by column', 'column fraction', 'fraction');
        add_panel_label(ax, row_bottom, 3, params.n_col);

        % Column 4: normalized by total
        ax = nexttile(tile_idx(row_top, 4));
        plot_transition_table_generic(ax, counts_total_norm, params.category_labels, agreement, kappa, n_cat, ...
            sprintf('%s categorical 3x3', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'normalized by total', 'total fraction', 'fraction');
        add_panel_label(ax, row_top, 4, params.n_col);

        ax = nexttile(tile_idx(row_bottom, 4));
        plot_transition_table_generic(ax, counts_2x2_total_norm, category_labels_2x2, agreement_2x2, kappa_2x2, n_cat_2x2, ...
            sprintf('%s categorical 2x2', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'normalized by total', 'total fraction', 'fraction');
        add_panel_label(ax, row_bottom, 4, params.n_col);

        % Column 5: O/E ratio under the independent-marginal chance baseline
        ax = nexttile(tile_idx(row_top, 5));
        plot_transition_table_generic(ax, counts_oe_ratio, params.category_labels, agreement, kappa, n_cat, ...
            sprintf('%s categorical 3x3', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'O/E ratio', 'observed / expected', 'fraction');
        add_panel_label(ax, row_top, 5, params.n_col);

        ax = nexttile(tile_idx(row_bottom, 5));
        plot_transition_table_generic(ax, counts_2x2_oe_ratio, category_labels_2x2, agreement_2x2, kappa_2x2, n_cat_2x2, ...
            sprintf('%s categorical 2x2', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'O/E ratio', 'observed / expected', 'fraction');
        add_panel_label(ax, row_bottom, 5, params.n_col);

        % Column 6: standardized residual relative to the independent-marginal expected counts
        ax = nexttile(tile_idx(row_top, 6));
        plot_transition_table_generic(ax, counts_standardized_residual, params.category_labels, agreement, kappa, n_cat, ...
            sprintf('%s categorical 3x3', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'standardized residual', '(O - E) / sqrt(E)', 'fraction');
        add_panel_label(ax, row_top, 6, params.n_col);

        ax = nexttile(tile_idx(row_bottom, 6));
        plot_transition_table_generic(ax, counts_2x2_standardized_residual, category_labels_2x2, agreement_2x2, kappa_2x2, n_cat_2x2, ...
            sprintf('%s categorical 2x2', ctx.name), sprintf('%s category', ctx.x.axis_label), sprintf('%s category', ctx.y.axis_label), ...
            'standardized residual', '(O - E) / sqrt(E)', 'fraction');
        add_panel_label(ax, row_bottom, 6, params.n_col);
    end

    sgtitle(sprintf('Supplement Fig 3: %s\n%s; valid sessions = %d', ...
        cfg.figure_title, params.hyperparam_title, pooled_one.valid_session_count), 'Interpreter', 'none');

    if isfield(params, 'output_folder') && ~isempty(params.output_folder)
        save_folder = params.output_folder;
    else
        save_folder = fullfile(root, 'Figures', 'Paper');
    end
    check_path(save_folder);

    figWidth = 24.0;  % inches. Wider layout for 6 columns.
    figHeight = 16.0; % inches.
    resolution = 300; % dpi.

    set(f, 'Units', 'inches');
    f.Position(3:4) = [figWidth, figHeight];
    set(f, 'PaperUnits', 'inches');
    set(f, 'PaperSize', [figWidth, figHeight]);
    set(f, 'PaperPosition', [0, 0, figWidth, figHeight]);
    set(f, 'Color', 'w');

    preview_filename = fullfile(save_folder, [cfg.output_stub, '_supp_preview.jpg']);
    exportgraphics(f, preview_filename, 'ContentType', 'image', 'BackgroundColor', 'white', 'Resolution', resolution);

    pdf_filename = fullfile(save_folder, [cfg.output_stub, '_supp.pdf']);
    exportgraphics(f, pdf_filename, 'ContentType', 'vector', 'BackgroundColor', 'white', 'Resolution', resolution);

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

function hyperparam_configs = build_auto_hyperparam_configs()
    reg_name = "L2=0_2";
    resting_dur_threshold = 15;
    injections = {"Muscimol", "Saline"};
    aligns = {"Last", "Longest"};

    kernel_specs = struct('kernel_name', {}, 'kernel_num', {});

    first_kernel_names = {"Exp5", "Exp40", "Exp100", "Exp40Mix", "Exp100Mix"};
    first_kernel_nums = [1, 1, 1, 2, 3];
    for i = 1:numel(first_kernel_names)
        kernel_specs(end + 1).kernel_name = first_kernel_names{i}; %#ok<AGROW>
        kernel_specs(end).kernel_num = first_kernel_nums(i);
    end

    suffixes = {'5', '10', '20', '40', '80', '160', '320', '640', '1000'};
    prefixes = {'LongExp', 'LongGaussC', 'LongGDeriv'};
    for p = 1:numel(prefixes)
        for s = 1:numel(suffixes)
            kernel_specs(end + 1).kernel_name = string(sprintf('%s%s', prefixes{p}, suffixes{s})); %#ok<AGROW>
            kernel_specs(end).kernel_num = 1;
        end
    end

    step_suffixes = {'0_5', '5_10', '10_20', '20_40', '40_80', ...
                     '80_160', '160_320', '320_640', '640_1000', '1000_3000'};
    for s = 1:numel(step_suffixes)
        kernel_specs(end + 1).kernel_name = string(sprintf('LongStepB%s', step_suffixes{s})); %#ok<AGROW>
        kernel_specs(end).kernel_num = 1;
    end

    hyperparam_configs = struct([]);
    cfg_i = 0;
    for k = 1:numel(kernel_specs)
        for inj_i = 1:numel(injections)
            for align_i = 1:numel(aligns)
                cfg_i = cfg_i + 1;
                hyperparam_configs(cfg_i).kernel_name = kernel_specs(k).kernel_name; %#ok<AGROW>
                hyperparam_configs(cfg_i).kernel_num = kernel_specs(k).kernel_num;
                hyperparam_configs(cfg_i).reg_name = reg_name;
                hyperparam_configs(cfg_i).resting_dur_threshold = resting_dur_threshold;
                hyperparam_configs(cfg_i).injection = injections{inj_i};
                hyperparam_configs(cfg_i).align = aligns{align_i};
            end
        end
    end

    fprintf('Generated %d hyperparameter configurations from %d kernel specifications.\n', ...
        numel(hyperparam_configs), numel(kernel_specs));
end

function issue_log = empty_run_issue_log()
    issue_log = struct('hp_index', {}, 'hyperparam_title', {}, 'issue_type', {}, 'message', {});
end

function issue_log = append_run_issue(issue_log, hp_index, hyperparam_title, issue_type, message)
    entry_i = numel(issue_log) + 1;
    issue_log(entry_i).hp_index = hp_index;
    issue_log(entry_i).hyperparam_title = char(string(hyperparam_title));
    issue_log(entry_i).issue_type = char(string(issue_type));
    issue_log(entry_i).message = char(string(message));
end

function write_run_issue_log(filename, issue_log)
    fid = fopen(filename, 'w');
    if fid < 0
        warning('Could not open run issue log for writing: %s', filename);
        return;
    end
    cleanup_obj = onCleanup(@() fclose(fid));

    fprintf(fid, 'hp_index\tissue_type\thyperparam_title\tmessage\n');
    for i = 1:numel(issue_log)
        fprintf(fid, '%d\t%s\t%s\t%s\n', ...
            issue_log(i).hp_index, ...
            escape_log_field(issue_log(i).issue_type), ...
            escape_log_field(issue_log(i).hyperparam_title), ...
            escape_log_field(issue_log(i).message));
    end
end

function out = escape_log_field(value)
    out = char(string(value));
    out = strrep(out, sprintf('\t'), ' ');
    out = strrep(out, sprintf('\n'), ' ');
    out = strrep(out, sprintf('\r'), ' ');
end

function hp = fill_default_hyperparams(hp)
    default_hp = struct();
    default_hp.kernel_name = "DeltaPure";
    default_hp.kernel_num = 3;
    default_hp.reg_name = [];
    default_hp.resting_dur_threshold = 15;
    default_hp.injection = "Muscimol";
    default_hp.align = "Last";

    fields = fieldnames(default_hp);
    for f = 1:numel(fields)
        field = fields{f};
        if ~isfield(hp, field)
            hp.(field) = default_hp.(field);
        end
    end

    if isempty(hp.kernel_num) || hp.kernel_num < 1 || hp.kernel_num ~= floor(hp.kernel_num)
        error('kernel_num must be a positive integer.');
    end
end

function selected_rows = metadata_filter(mt, hp)
    base_filter = true(height(mt), 1);

    base_filter = base_filter & metadata_matches(mt, 'kernel_name', hp.kernel_name);
    base_filter = base_filter & metadata_matches(mt, 'align', hp.align);
    base_filter = base_filter & metadata_matches(mt, 'area', "Cortex");
    base_filter = base_filter & metadata_matches(mt, 'injection', hp.injection);
    base_filter = base_filter & metadata_matches(mt, 'epoch', 3000);
    base_filter = base_filter & metadata_matches(mt, 'fold_idx', 0);
    base_filter = base_filter & metadata_matches(mt, 'shuffle_idx', 0);
    base_filter = base_filter & metadata_matches(mt, 'resting_dur_threshold', hp.resting_dur_threshold);

    if ~is_empty_filter_value(hp.reg_name)
        base_filter = base_filter & metadata_matches(mt, 'reg_name', hp.reg_name);
    end

    anchor_filter = base_filter & ...
                    metadata_matches(mt, 'prepost', 'Pre') & ...
                    metadata_matches(mt, 'state', 'RestOpen');

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

            anchor_value = get_table_scalar(mt, field, idx);
            same_session = same_session & metadata_matches(mt, field, anchor_value);
        end

        has_all_required = true;
        for q = 1:numel(required_preposts)
            has_this = any(same_session & ...
                           metadata_matches(mt, 'prepost', required_preposts{q}) & ...
                           metadata_matches(mt, 'state', required_states{q}));

            if ~has_this
                has_all_required = false;
                anchor_name = get_optional_table_string(mt, 'file_name', idx, sprintf('row_%d', idx));
                warning('Session %s is missing required Pre/Post x Open/Close combination: %s %s. Skipping this session.', ...
                    anchor_name, required_preposts{q}, required_states{q});
                break;
            end
        end

        if has_all_required
            selected_rows(idx) = true;
        end
    end

    fprintf('Metadata filter selected %d/%d anchor rows for %s:\n', ...
        sum(selected_rows), sum(base_filter), make_hyperparam_title(hp));
    for k = 1:height(mt)
        if selected_rows(k)
            fprintf('  %s\n', get_optional_table_string(mt, 'file_name', k, sprintf('row_%d', k)));
        end
    end

    if ~any(selected_rows)
        warning('Metadata filter selected no complete Pre/Post x Open/Close sets for %s.', make_hyperparam_title(hp));
    end
end

function mask = metadata_matches(mt, field, target)
    if is_empty_filter_value(target)
        mask = true(height(mt), 1);
        return;
    end

    if ~ismember(field, mt.Properties.VariableNames)
        error('Metadata table is missing required field: %s', field);
    end

    mask = false(height(mt), 1);
    for idx = 1:height(mt)
        value = get_table_scalar(mt, field, idx);
        mask(idx) = metadata_scalar_matches(value, target);
    end
end

function value = get_table_scalar(mt, field, idx)
    col = mt.(field);
    if iscell(col)
        value = col{idx};
    else
        value = col(idx);
    end
end

function out = get_optional_table_string(mt, field, idx, fallback)
    if ~ismember(field, mt.Properties.VariableNames)
        out = fallback;
        return;
    end
    out = value_to_display_string(get_table_scalar(mt, field, idx), fallback);
end

function tf = metadata_scalar_matches(value, target)
    if is_empty_filter_value(target)
        tf = true;
        return;
    end

    if iscell(target)
        tf = false;
        for i = 1:numel(target)
            if metadata_scalar_matches(value, target{i})
                tf = true;
                return;
            end
        end
        return;
    end

    if isstring(target) && numel(target) > 1
        tf = false;
        for i = 1:numel(target)
            if metadata_scalar_matches(value, target(i))
                tf = true;
                return;
            end
        end
        return;
    end

    if isnumeric(target) && numel(target) > 1
        tf = false;
        for i = 1:numel(target)
            if metadata_scalar_matches(value, target(i))
                tf = true;
                return;
            end
        end
        return;
    end

    value = unwrap_cell_scalar(value);
    target = unwrap_cell_scalar(target);

    if ischar(value) || isstring(value) || ischar(target) || isstring(target) || iscategorical_safe(value) || iscategorical_safe(target)
        tf = strcmp(string(value), string(target));
    elseif isnumeric(value) && isnumeric(target)
        tf = isequal(value, target);
    elseif islogical(value) && islogical(target)
        tf = isequal(value, target);
    else
        tf = isequal(value, target);
    end

    if numel(tf) > 1
        tf = any(tf(:));
    end
end

function value = unwrap_cell_scalar(value)
    while iscell(value) && isscalar(value)
        value = value{1};
    end
end

function tf = iscategorical_safe(value)
    tf = false;
    try
        tf = iscategorical(value);
    catch
        tf = false;
    end
end

function tf = is_empty_filter_value(value)
    if isempty(value)
        tf = true;
        return;
    end
    if isstring(value) && all(strlength(value) == 0)
        tf = true;
        return;
    end
    if iscell(value) && isempty(value)
        tf = true;
        return;
    end
    tf = false;
end

function title_str = make_hyperparam_title(hp)
    title_str = sprintf('kernel_name=%s, kernel_num=%d, reg_name=%s, resting_dur_threshold=%s, injection=%s, align=%s', ...
        value_to_display_string(hp.kernel_name, 'all'), ...
        hp.kernel_num, ...
        value_to_display_string(hp.reg_name, 'all'), ...
        value_to_display_string(hp.resting_dur_threshold, 'all'), ...
        value_to_display_string(hp.injection, 'all'), ...
        value_to_display_string(hp.align, 'all'));
end

function folder_name = make_hyperparam_folder_name(hp)
    parts = { ...
        ['kernel_', sanitize_for_path(value_to_display_string(hp.kernel_name, 'all'))], ...
        ['kernelNum_', sanitize_for_path(value_to_display_string(hp.kernel_num, 'all'))], ...
        ['reg_', sanitize_for_path(value_to_display_string(hp.reg_name, 'all'))], ...
        ['restDur_', sanitize_for_path(value_to_display_string(hp.resting_dur_threshold, 'all'))], ...
        ['inj_', sanitize_for_path(value_to_display_string(hp.injection, 'all'))], ...
        ['align_', sanitize_for_path(value_to_display_string(hp.align, 'all'))] ...
    };
    folder_name = strjoin(parts, '__');
end

function out = value_to_display_string(value, empty_label)
    if nargin < 2
        empty_label = 'all';
    end

    if is_empty_filter_value(value)
        out = empty_label;
        return;
    end

    value = unwrap_cell_scalar(value);

    if iscell(value)
        parts = cell(size(value));
        for i = 1:numel(value)
            parts{i} = value_to_display_string(value{i}, empty_label);
        end
        out = strjoin(parts(:).', '+');
    elseif isstring(value)
        parts = cell(1, numel(value));
        for i = 1:numel(value)
            parts{i} = char(value(i));
        end
        out = strjoin(parts, '+');
    elseif ischar(value)
        out = value;
    elseif isnumeric(value) || islogical(value)
        if isscalar(value)
            out = num2str(value);
        else
            parts = cell(1, numel(value));
            for i = 1:numel(value)
                parts{i} = num2str(value(i));
            end
            out = strjoin(parts, '+');
        end
    else
        out = char(string(value));
    end
end

function out = sanitize_for_path(str_in)
    out = char(str_in);
    out = regexprep(out, '[^A-Za-z0-9._=-]+', '_');
    out = regexprep(out, '^_+|_+$', '');
    if isempty(out)
        out = 'empty';
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

function [rho, pval, n_valid] = pearson_stats(x, y)
    if nargin < 5
        axis_limit_override = [];
    end

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    n_valid = numel(x);

    if n_valid < 2
        rho = NaN;
        pval = NaN;
        return;
    end

    if all(x == x(1)) || all(y == y(1))
        rho = NaN;
        pval = NaN;
        return;
    end

    [R, P] = corr(x, y, 'Type', 'Pearson', 'Rows', 'complete');
    rho = R;
    pval = P;
end

function plot_pair_scatter(ax, data, plot_mode, marker_size, marker_alpha, colors, show_legend, title_text, x_label, y_label, axis_limit_override)
    cla(ax);
    hold(ax, 'on');

    if nargin < 11
        axis_limit_override = [];
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
            both_sig = x_sig & y_sig;

            scatter(ax, x(none_sig), y(none_sig), marker_size, 'filled', ...
                'MarkerFaceColor', colors.non_sig, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'both non-sig');
            scatter(ax, x(x_only), y(x_only), marker_size, 'filled', ...
                'MarkerFaceColor', colors.x_only, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s only', x_label));
            scatter(ax, x(y_only), y(y_only), marker_size, 'filled', ...
                'MarkerFaceColor', colors.y_only, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', sprintf('%s only', y_label));
            scatter(ax, x(both_sig), y(both_sig), marker_size, 'filled', ...
                'MarkerFaceColor', colors.both_sig, 'MarkerFaceAlpha', marker_alpha, 'MarkerEdgeAlpha', marker_alpha, ...
                'DisplayName', 'both sig');

            axis_limit = get_symmetric_axis_limit(x, y);
            if ~isempty(axis_limit_override)
                axis_limit = axis_limit_override;
            end;
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
            switch_mask = (x_cat ~= 0) & (y_cat ~= 0) & (x_cat ~= y_cat);

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

function plot_pair_density(ax, x, y, rho, pval, n_valid, nbin, clip_percentile, use_log_count, title_text, axis_mode, colors, x_label, y_label, axis_limit_override)
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
    if nargin < 15
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
    add_stats_text(ax, rho, pval, n_valid);

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

function add_stats_text(ax, rho, pval, n_valid)
    if isnan(rho)
        rho_str = 'NaN';
    else
        rho_str = sprintf('%.4f', rho);
    end
    if isnan(pval)
        p_str = 'NaN';
    else
        p_str = sprintf('%.3e', pval);
    end
    text(ax, 0.02, 0.98, sprintf('r = %s\np = %s\nn = %d', rho_str, p_str, n_valid), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'w', 'Margin', 2, 'FontSize', 8);
end

function plot_transition_table_generic(ax, mat, category_labels, agreement, kappa, n_valid, title_text, xlabel_text, ylabel_text, normalization_label, colorbar_label, value_mode)
    plot_mat = mat.';
    imagesc(ax, plot_mat);
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
    ylabel(cb, colorbar_label);

    finite_vals = plot_mat(isfinite(plot_mat));
    if isempty(finite_vals)
        max_val = 1;
        min_val = 0;
    else
        max_val = max(finite_vals);
        min_val = min(finite_vals);
    end
    for y_idx = 1:n_cat
        for x_idx = 1:n_cat
            this_val = plot_mat(y_idx, x_idx);
            if isfinite(this_val)
                this_ratio = (this_val - min_val) / (max_val - min_val + eps);
            else
                this_ratio = 1;
            end
            if this_ratio < 0.2
                text_color = 'w';
            else
                text_color = 'k';
            end

            if ~isfinite(this_val)
                value_str = 'NaN';
            elseif strcmp(value_mode, 'integer')
                value_str = sprintf('%d', round(this_val));
            else
                value_str = sprintf('%.2f', this_val);
            end

            text(ax, x_idx, y_idx, value_str, ...
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

    title(ax, sprintf('%s\n%s\nAgreement = %s, kappa = %s, n = %d', ...
        title_text, normalization_label, agreement_str, kappa_str, n_valid), 'Interpreter', 'none');
end

function mat_norm = normalize_transition_matrix(counts, mode)
    counts = double(counts);
    switch lower(mode)
        case 'row'
            denom = sum(counts, 2);
            mat_norm = zeros(size(counts));
            valid = denom > 0;
            mat_norm(valid, :) = bsxfun(@rdivide, counts(valid, :), denom(valid));
        case 'column'
            denom = sum(counts, 1);
            mat_norm = zeros(size(counts));
            valid = denom > 0;
            mat_norm(:, valid) = bsxfun(@rdivide, counts(:, valid), denom(valid));
        case 'total'
            denom = sum(counts(:));
            if denom > 0
                mat_norm = counts / denom;
            else
                mat_norm = zeros(size(counts));
            end
        case 'oe_ratio'
            total_n = sum(counts(:));
            mat_norm = nan(size(counts));
            if total_n > 0
                observed_prob = counts / total_n;
                x_marginal_prob = sum(counts, 2) / total_n;
                y_marginal_prob = sum(counts, 1) / total_n;
                expected_prob = x_marginal_prob * y_marginal_prob;
                valid = expected_prob > 0;
                mat_norm(valid) = observed_prob(valid) ./ expected_prob(valid);
            end
        case 'standardized_residual'
            total_n = sum(counts(:));
            mat_norm = nan(size(counts));
            if total_n > 0
                row_counts = sum(counts, 2);
                col_counts = sum(counts, 1);
                expected_counts = row_counts * col_counts / total_n;
                valid = expected_counts > 0;
                mat_norm(valid) = (counts(valid) - expected_counts(valid)) ./ sqrt(expected_counts(valid));
            end
        otherwise
            error('Unknown normalization mode: %s', mode);
    end
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
