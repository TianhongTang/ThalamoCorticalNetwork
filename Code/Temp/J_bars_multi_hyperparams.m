%% J_bars_multi_hyperparams.m - Significant-J bar plots for multiple GLM configurations.
%
% Each configuration is analyzed and plotted separately.
%
% A configuration is:
%   kernel_name x reg_name x align x resting_dur_threshold
%
% Injection, state, and Pre/Post remain comparison dimensions within each
% configuration, preserving the figures and pooled-effect analyses from the
% original script.
%
% Output organization:
%   Figures/Paper/J_bars_multi_config/
%       <configuration-folder>/
%           primary_comparisons/
%           pooled_effects/
%
% Metadata behavior:
%   1. The first run loads the full GLM metadata as a struct, applies the
%      fixed J-bars filters, and saves a compact metadata index.
%   2. Later runs load only the compact index.
%
% Edit build_default_kernel_specs(), reg_names, alignments, and
% resting_dur_thresholds to control the configurations.

clear;
run_tic = tic;
progress_log('SCRIPT', 'Started.');

%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end

addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));
progress_log('SCRIPT', 'Project root: %s', root);

%% Configuration sweep
% kernel_specs supports both one-block kernels and composite models with
% several internal kernel blocks.
kernel_specs = build_default_kernel_specs();

% Current analysis is fixed at L2=0.2. Add entries here for future sweeps,
% for example:
%   reg_names = {'L2=0_2', 'L2=0_4', 'L1=0_2', 'L12=0_2'};
reg_names = {'L2=0_2'};

% Add 'Longest' here when needed.
alignments = {'Last'};

% Kept as a vector so future threshold sweeps need no structural changes.
resting_dur_thresholds = 15;

analysis_configs = build_analysis_configs( ...
    kernel_specs, reg_names, alignments, resting_dur_thresholds);
progress_log('CONFIG', 'Generated %d configurations from %d kernel specifications, %d reg values, %d alignments, and %d duration thresholds.', ...
    numel(analysis_configs), numel(kernel_specs), numel(reg_names), ...
    numel(alignments), numel(resting_dur_thresholds));

%% Dataset and comparison parameters
% Empty animal_names means all animals in the compact metadata index.
animal_names = {};

injection_types = {'Muscimol', 'Saline'};
states = {'RestOpen', 'RestClose'};
state_labels = {'Eyes Open', 'Eyes Closed'};
posneg_labels = {'Positive J', 'Negative J'};
prepost_labels = {'Pre', 'Post'};

bar_width = 0.30;
bar_offset = 0.18;
y_limit = [0, 20];

% Pooled connection-level effect figures.
effect_injection_types = {'Saline', 'Muscimol'};
effect_injection_labels = {'Saline', 'Muscimol'};
effect_y_limit = [];
effect_diff_y_limit = [];
effect_alpha = 0.05;

% 'difference' = Post - Pre in percentage points.
% 'ratio'      = Post / Pre significant-J ratio.
effect_metric_mode = 'difference';
make_combined_effect_figure = true;
make_single_injection_effect_figures = true;
make_prepost_percent_figures = true;
make_difference_only_figure = false;

%% Execution and output controls
skip_failed_configs = false;
max_sessions_per_condition = inf; % Set smaller while debugging.

figure_visible = 'off';           % 'on' or 'off'.
export_jpg = true;
export_pdf = false;
export_resolution = 300;
close_figures_after_export = true;

output_root = fullfile(root, 'Figures', 'Paper', 'J_bars_multi_config');
check_path(output_root);

run_issue_log = empty_run_issue_log();
run_issue_log_filename = fullfile(output_root, 'run_issue_log.txt');

%% Compact metadata index
metadata_index_version = 1;
force_rebuild_metadata_index = false;
metadata_index_filename = fullfile(root, 'Data', 'Working', ...
    'metadata_index', 'GLM_epoch3000_Jbars_index.mat');

progress_log('INDEX', 'Index file: %s', metadata_index_filename);
[mt_glm, metadata_index_info] = load_or_build_jbars_metadata_index( ...
    root, metadata_index_filename, metadata_index_version, ...
    force_rebuild_metadata_index);
progress_log('INDEX', 'Index ready: %d rows; created %s.', ...
    height(mt_glm), metadata_index_info.created_at);

animal_filter = build_animal_filter(mt_glm, animal_names);
progress_log('CONFIG', 'Animal selection retained %d/%d metadata rows (%s).', ...
    sum(animal_filter), height(mt_glm), ...
    value_to_display_string(animal_names, 'all animals'));

%% Runtime output context shared by plotting helpers
setappdata(0, 'JBarsFigureVisible', figure_visible);
setappdata(0, 'JBarsExportJPG', export_jpg);
setappdata(0, 'JBarsExportPDF', export_pdf);
setappdata(0, 'JBarsExportResolution', export_resolution);
setappdata(0, 'JBarsCloseFigures', close_figures_after_export);
set(0, 'DefaultFigureVisible', figure_visible);

progress_log('CONFIG', 'Figure output: visible=%s, JPG=%d, PDF=%d, resolution=%d dpi.', ...
    figure_visible, export_jpg, export_pdf, export_resolution);

total_rendered_figure_count = 0;
successful_config_count = 0;
failed_config_count = 0;
skipped_config_count = 0;

%% Run each configuration separately
for config_i = 1:numel(analysis_configs)
    config_tic = tic;
    cfg = analysis_configs(config_i);
    config_title = make_config_title(cfg);
    config_folder = fullfile(output_root, make_config_folder_name(cfg));
    primary_output_folder = fullfile(config_folder, 'primary_comparisons');
    effect_output_folder = fullfile(config_folder, 'pooled_effects');

    progress_log('CONFIG', '[%d/%d] START: %s', ...
        config_i, numel(analysis_configs), config_title);

    try
        config_meta = select_config_metadata(mt_glm, animal_filter, cfg);
        progress_log('CONFIG', '[%d/%d][1/4] Selected %d metadata rows.', ...
            config_i, numel(analysis_configs), height(config_meta));

        if isempty(config_meta)
            message = 'No metadata rows matched this configuration.';
            warning('%s %s', message, config_title);
            run_issue_log = append_run_issue(run_issue_log, config_i, ...
                config_title, 'NoMetadata', message);
            write_run_issue_log(run_issue_log_filename, run_issue_log);
            skipped_config_count = skipped_config_count + 1;
            progress_log('CONFIG', '[%d/%d] SKIPPED: no matching metadata.', ...
                config_i, numel(analysis_configs));
            continue;
        end

        kernel_indices = 1:cfg.kernel_num;
        kernel_labels = cfg.kernel_labels;
        if numel(kernel_labels) ~= numel(kernel_indices)
            kernel_labels = arrayfun(@(k) sprintf('Kernel %d', k), ...
                kernel_indices, 'UniformOutput', false);
        end

        setappdata(0, 'JBarsConfigTitle', config_title);

        progress_log('CONFIG', '[%d/%d][2/4] Pooling all condition counts. Internal kernels=%s.', ...
            config_i, numel(analysis_configs), mat2str(kernel_indices));
        condition_counts = collect_config_condition_counts( ...
            root, config_meta, states, injection_types, prepost_labels, ...
            kernel_indices, max_sessions_per_condition);

        if ~any(condition_counts.total(:) > 0)
            message = 'No valid GLM connections were counted for this configuration.';
            warning('%s %s', message, config_title);
            run_issue_log = append_run_issue(run_issue_log, config_i, ...
                config_title, 'NoValidCounts', message);
            write_run_issue_log(run_issue_log_filename, run_issue_log);
            skipped_config_count = skipped_config_count + 1;
            progress_log('CONFIG', '[%d/%d] SKIPPED: no valid counts.', ...
                config_i, numel(analysis_configs));
            continue;
        end

        %% Original primary figures: Pre/Post and Open/Closed
        progress_log('CONFIG', '[%d/%d][3/4] Rendering primary comparison figures.', ...
            config_i, numel(analysis_configs));
        setappdata(0, 'JBarsOutputFolder', primary_output_folder);

        primary_figure_count = render_primary_comparison_figures( ...
            root, condition_counts, injection_types, states, state_labels, ...
            posneg_labels, prepost_labels, kernel_labels, ...
            bar_offset, bar_width, y_limit);
        total_rendered_figure_count = total_rendered_figure_count + ...
            primary_figure_count;

        %% Original extended pooled-effect figures
        progress_log('CONFIG', '[%d/%d][4/4] Rendering pooled-effect figures.', ...
            config_i, numel(analysis_configs));
        setappdata(0, 'JBarsOutputFolder', effect_output_folder);

        effect_stats_grid = build_effect_stats_grid_from_counts( ...
            condition_counts, injection_types, effect_injection_types, ...
            prepost_labels);

        effect_figure_count = 0;
        if make_combined_effect_figure
            plot_effect_metric_combined_pooled_ext( ...
                root, effect_stats_grid, states, state_labels, posneg_labels, ...
                kernel_indices, kernel_labels, effect_injection_types, ...
                effect_injection_labels, bar_offset, bar_width, ...
                effect_metric_mode, effect_y_limit, effect_alpha);
            effect_figure_count = effect_figure_count + 1;
        end

        if make_single_injection_effect_figures
            for effect_injection_idx = 1:numel(effect_injection_types)
                plot_effect_metric_single_injection_pooled_ext( ...
                    root, effect_stats_grid, states, state_labels, ...
                    posneg_labels, kernel_indices, kernel_labels, ...
                    effect_injection_types, effect_injection_labels, ...
                    effect_injection_idx, effect_metric_mode, ...
                    effect_y_limit, effect_alpha);
                effect_figure_count = effect_figure_count + 1;
            end
        end

        if make_prepost_percent_figures
            for effect_injection_idx = 1:numel(effect_injection_types)
                plot_prepost_percent_injection_pooled_ext( ...
                    root, effect_stats_grid, states, state_labels, ...
                    posneg_labels, kernel_indices, kernel_labels, ...
                    effect_injection_types, effect_injection_labels, ...
                    effect_injection_idx, bar_offset, bar_width, ...
                    effect_alpha);
                effect_figure_count = effect_figure_count + 1;
            end
        end

        if make_difference_only_figure
            plot_muscimol_minus_saline_effect_pooled_ext( ...
                root, effect_stats_grid, states, state_labels, ...
                posneg_labels, kernel_indices, kernel_labels, ...
                effect_injection_types, effect_metric_mode, ...
                effect_diff_y_limit, effect_alpha);
            effect_figure_count = effect_figure_count + 1;
        end

        total_rendered_figure_count = total_rendered_figure_count + ...
            effect_figure_count;
        successful_config_count = successful_config_count + 1;

        progress_log('CONFIG', ...
            '[%d/%d] DONE in %.1f s. Metadata rows=%d, figures=%d.', ...
            config_i, numel(analysis_configs), toc(config_tic), ...
            height(config_meta), primary_figure_count + effect_figure_count);

    catch ME_config
        failed_config_count = failed_config_count + 1;
        message = sprintf('%s: %s', ME_config.identifier, ME_config.message);
        warning('Configuration %d/%d failed: %s\n%s', ...
            config_i, numel(analysis_configs), config_title, message);
        run_issue_log = append_run_issue(run_issue_log, config_i, ...
            config_title, 'Error', message);
        write_run_issue_log(run_issue_log_filename, run_issue_log);
        progress_log('CONFIG', '[%d/%d] FAILED after %.1f s: %s', ...
            config_i, numel(analysis_configs), toc(config_tic), message);

        close_created_figures();
        if ~skip_failed_configs
            rethrow(ME_config);
        end
    end
end

if ~isempty(run_issue_log)
    write_run_issue_log(run_issue_log_filename, run_issue_log);
    progress_log('SUMMARY', 'Issue log written to: %s', run_issue_log_filename);
end

progress_log('SUMMARY', ...
    ['Finished in %.1f s. Configurations: successful=%d, skipped=%d, ' ...
     'failed=%d. Figures=%d.'], ...
    toc(run_tic), successful_config_count, skipped_config_count, ...
    failed_config_count, total_rendered_figure_count);


function figure_count = render_primary_comparison_figures( ...
    root, condition_counts, injection_types, states, state_labels, ...
    posneg_labels, prepost_labels, kernel_labels, ...
    bar_offset, bar_width, y_limit)

    n_kernel = numel(kernel_labels);
    figure_count = 0;

    for injection_idx = 1:numel(injection_types)
        injection = injection_types{injection_idx};

        %% Figure 1: Pre vs Post
        progress_log('FIGURE', 'Rendering Pre vs Post for injection=%s.', injection);
        f1 = figure('Position', [100, 100, 1200, 800], ...
            'Visible', get_jbars_figure_visible());
        t1 = tiledlayout(f1, 2, 2, ...
            'TileSpacing', 'Compact', 'Padding', 'Compact');

        for state_idx = 1:numel(states)
            ratios_all = nan(n_kernel, 2, numel(prepost_labels));
            low_errs_all = nan(n_kernel, 2, numel(prepost_labels));
            high_errs_all = nan(n_kernel, 2, numel(prepost_labels));
            sig_labels_all = cell(n_kernel, 2);

            for kernel_i = 1:n_kernel
                pos_counts = reshape( ...
                    condition_counts.pos(state_idx, injection_idx, :, kernel_i), ...
                    1, []);
                neg_counts = reshape( ...
                    condition_counts.neg(state_idx, injection_idx, :, kernel_i), ...
                    1, []);
                total_vec = reshape( ...
                    condition_counts.total(state_idx, injection_idx, :, kernel_i), ...
                    1, []);

                J_counts = [pos_counts; neg_counts];
                total_counts = [total_vec; total_vec];

                [ratios, low_errs, high_errs, sig_labels] = ...
                    compute_bar_stats(J_counts, total_counts);
                ratios_all(kernel_i, :, :) = ratios;
                low_errs_all(kernel_i, :, :) = low_errs;
                high_errs_all(kernel_i, :, :) = high_errs;
                sig_labels_all(kernel_i, :) = sig_labels;
            end

            for posneg_idx = 1:2
                tile_idx = (state_idx - 1) * 2 + posneg_idx;
                tile = nexttile(t1, tile_idx);
                plot_kernel_grouped_bars( ...
                    tile, squeeze(ratios_all(:, posneg_idx, :)), ...
                    squeeze(low_errs_all(:, posneg_idx, :)), ...
                    squeeze(high_errs_all(:, posneg_idx, :)), ...
                    sig_labels_all(:, posneg_idx), prepost_labels, ...
                    kernel_labels, bar_offset, bar_width);
                ylabel(tile, 'Significant J %');
                title(tile, sprintf('%s, %s', ...
                    state_labels{state_idx}, posneg_labels{posneg_idx}));
                if ~isempty(y_limit)
                    ylim(tile, y_limit);
                end
            end
        end

        jbars_sgtitle(f1, sprintf('%s, Pre vs Post', injection));
        export_figure(root, f1, sprintf( ...
            'J_bars_prepost_kernel_groups_by_sign_%s', injection));
        figure_count = figure_count + 1;

        %% Figure 2: Open vs Closed
        progress_log('FIGURE', 'Rendering Open vs Closed for injection=%s.', injection);
        f2 = figure('Position', [100, 100, 1200, 800], ...
            'Visible', get_jbars_figure_visible());
        t2 = tiledlayout(f2, 2, 2, ...
            'TileSpacing', 'Compact', 'Padding', 'Compact');

        for prepost_idx = 1:numel(prepost_labels)
            ratios_all = nan(n_kernel, 2, numel(states));
            low_errs_all = nan(n_kernel, 2, numel(states));
            high_errs_all = nan(n_kernel, 2, numel(states));
            sig_labels_all = cell(n_kernel, 2);

            for kernel_i = 1:n_kernel
                pos_counts = reshape( ...
                    condition_counts.pos(:, injection_idx, prepost_idx, kernel_i), ...
                    1, []);
                neg_counts = reshape( ...
                    condition_counts.neg(:, injection_idx, prepost_idx, kernel_i), ...
                    1, []);
                total_vec = reshape( ...
                    condition_counts.total(:, injection_idx, prepost_idx, kernel_i), ...
                    1, []);

                J_counts = [pos_counts; neg_counts];
                total_counts = [total_vec; total_vec];

                [ratios, low_errs, high_errs, sig_labels] = ...
                    compute_bar_stats(J_counts, total_counts);
                ratios_all(kernel_i, :, :) = ratios;
                low_errs_all(kernel_i, :, :) = low_errs;
                high_errs_all(kernel_i, :, :) = high_errs;
                sig_labels_all(kernel_i, :) = sig_labels;
            end

            for posneg_idx = 1:2
                tile_idx = (prepost_idx - 1) * 2 + posneg_idx;
                tile = nexttile(t2, tile_idx);
                plot_kernel_grouped_bars( ...
                    tile, squeeze(ratios_all(:, posneg_idx, :)), ...
                    squeeze(low_errs_all(:, posneg_idx, :)), ...
                    squeeze(high_errs_all(:, posneg_idx, :)), ...
                    sig_labels_all(:, posneg_idx), state_labels, ...
                    kernel_labels, bar_offset, bar_width);
                ylabel(tile, 'Significant J %');
                title(tile, sprintf('%s, %s', ...
                    prepost_labels{prepost_idx}, ...
                    posneg_labels{posneg_idx}));
                if ~isempty(y_limit)
                    ylim(tile, y_limit);
                end
            end
        end

        jbars_sgtitle(f2, sprintf('%s, Open vs Closed', injection));
        export_figure(root, f2, sprintf( ...
            'J_bars_openclosed_kernel_groups_by_sign_%s', injection));
        figure_count = figure_count + 1;
    end
end

function counts = collect_config_condition_counts( ...
    root, config_meta, states, injection_types, prepost_labels, ...
    kernel_indices, max_sessions_per_condition)

    n_state = numel(states);
    n_injection = numel(injection_types);
    n_prepost = numel(prepost_labels);
    n_kernel = numel(kernel_indices);

    counts = struct();
    counts.pos = zeros(n_state, n_injection, n_prepost, n_kernel);
    counts.neg = zeros(n_state, n_injection, n_prepost, n_kernel);
    counts.total = zeros(n_state, n_injection, n_prepost, n_kernel);
    counts.session_count = zeros(n_state, n_injection, n_prepost);

    condition_total = n_state * n_injection * n_prepost;
    condition_i = 0;

    for state_idx = 1:n_state
        state = states{state_idx};
        for injection_idx = 1:n_injection
            injection = injection_types{injection_idx};
            for prepost_idx = 1:n_prepost
                prepost = prepost_labels{prepost_idx};
                condition_i = condition_i + 1;

                condition_meta = filter_condition_meta( ...
                    config_meta, state, injection, prepost);
                if isfinite(max_sessions_per_condition)
                    condition_meta = condition_meta( ...
                        1:min(height(condition_meta), max_sessions_per_condition), :);
                end

                counts.session_count(state_idx, injection_idx, prepost_idx) = ...
                    height(condition_meta);

                progress_log('COUNT', ...
                    '[%d/%d] state=%s, injection=%s, prepost=%s: %d metadata rows.', ...
                    condition_i, condition_total, state, injection, prepost, ...
                    height(condition_meta));

                if isempty(condition_meta)
                    warning('No metadata for state=%s, injection=%s, prepost=%s.', ...
                        state, injection, prepost);
                    continue;
                end

                [pos_count, neg_count, total_count] = ...
                    count_condition_J_all_kernels( ...
                        root, condition_meta, kernel_indices);

                counts.pos(state_idx, injection_idx, prepost_idx, :) = ...
                    reshape(pos_count, 1, 1, 1, []);
                counts.neg(state_idx, injection_idx, prepost_idx, :) = ...
                    reshape(neg_count, 1, 1, 1, []);
                counts.total(state_idx, injection_idx, prepost_idx, :) = ...
                    reshape(total_count, 1, 1, 1, []);

                progress_log('COUNT', ...
                    ['Completed state=%s, injection=%s, prepost=%s. ' ...
                     'Totals by kernel=%s.'], ...
                    state, injection, prepost, mat2str(total_count(:).'));
            end
        end
    end
end

function stats_grid = build_effect_stats_grid_from_counts( ...
    counts, injection_types, effect_injection_types, prepost_labels)

    pre_idx = find(strcmp(prepost_labels, 'Pre'), 1);
    post_idx = find(strcmp(prepost_labels, 'Post'), 1);
    if isempty(pre_idx) || isempty(post_idx)
        error('prepost_labels must contain Pre and Post.');
    end

    n_state = size(counts.pos, 1);
    n_kernel = size(counts.pos, 4);
    stats_grid = cell(n_state, n_kernel, numel(effect_injection_types));

    for effect_injection_idx = 1:numel(effect_injection_types)
        injection = effect_injection_types{effect_injection_idx};
        source_injection_idx = find(strcmp(injection_types, injection), 1);
        if isempty(source_injection_idx)
            error('Effect injection %s is not present in injection_types.', injection);
        end

        for state_idx = 1:n_state
            for kernel_i = 1:n_kernel
                stats = struct();
                stats.pre_counts = [ ...
                    counts.pos(state_idx, source_injection_idx, pre_idx, kernel_i); ...
                    counts.neg(state_idx, source_injection_idx, pre_idx, kernel_i)];
                stats.post_counts = [ ...
                    counts.pos(state_idx, source_injection_idx, post_idx, kernel_i); ...
                    counts.neg(state_idx, source_injection_idx, post_idx, kernel_i)];
                stats.pre_total = ...
                    counts.total(state_idx, source_injection_idx, pre_idx, kernel_i);
                stats.post_total = ...
                    counts.total(state_idx, source_injection_idx, post_idx, kernel_i);
                stats_grid{state_idx, kernel_i, effect_injection_idx} = stats;
            end
        end
    end
end

function condition_meta = filter_condition_meta( ...
    config_meta, state, injection, prepost)

    condition_filter = ...
        metadata_table_matches(config_meta, 'state', state) & ...
        metadata_table_matches(config_meta, 'injection', injection) & ...
        metadata_table_matches(config_meta, 'prepost', prepost);
    condition_meta = config_meta(condition_filter, :);
end

function [pos_count, neg_count, total_count] = count_condition_J( ...
    root, condition_meta, kernel_idx)
    [pos_count, neg_count, total_count] = ...
        count_condition_J_all_kernels(root, condition_meta, kernel_idx);
    pos_count = pos_count(1);
    neg_count = neg_count(1);
    total_count = total_count(1);
end

function [pos_counts, neg_counts, total_counts] = ...
    count_condition_J_all_kernels(root, condition_meta, kernel_indices)

    kernel_indices = kernel_indices(:).';
    n_kernel = numel(kernel_indices);
    pos_counts = zeros(n_kernel, 1);
    neg_counts = zeros(n_kernel, 1);
    total_counts = zeros(n_kernel, 1);

    progress_log('LOAD', 'Counting %d sessions for kernel blocks %s.', ...
        height(condition_meta), mat2str(kernel_indices));

    for session_i = 1:height(condition_meta)
        meta = table2struct(condition_meta(session_i, :));

        if ~isfield(meta, 'file_name') || isempty(meta.file_name)
            warning('Metadata row %d has no file_name. Skipping.', session_i);
            continue;
        end

        data_folder = fullfile(root, 'Data', 'Working', 'GLM');
        file_path = fullfile(data_folder, char(string(meta.file_name)));
        progress_log('LOAD', '[%d/%d] GLM: %s', ...
            session_i, height(condition_meta), file_path);

        if ~isfile(file_path)
            warning('GLM file not found: %s. Skipping.', file_path);
            continue;
        end

        session_data = load(file_path, 'meta', 'data');
        if ~isfield(session_data, 'data') || ...
                ~isfield(session_data.data, 'model_par') || ...
                ~isfield(session_data.data, 'model_err') || ...
                ~isfield(session_data.data.model_err, 'total')
            warning('GLM file is missing model_par/model_err.total: %s. Skipping.', ...
                file_path);
            continue;
        end

        model_par = session_data.data.model_par;
        model_err = session_data.data.model_err.total;
        if ~isequal(size(model_par), size(model_err))
            warning('model_par and model_err.total sizes differ: %s. Skipping.', ...
                file_path);
            continue;
        end

        border_folder = fullfile(root, 'Data', 'Working', 'border');
        border_file_name = generate_filename('border', meta);
        border_file_path = fullfile(border_folder, border_file_name);
        if ~isfile(border_file_path)
            warning('Border file not found: %s. Skipping.', border_file_path);
            continue;
        end

        border_data = load(border_file_path, 'data', 'meta');
        if ~isfield(border_data, 'data') || ...
                ~isfield(border_data.data, 'borders') || ...
                ~isfield(border_data, 'meta') || ...
                ~isfield(border_data.meta, 'N')
            warning('Border file is missing borders/N: %s. Skipping.', ...
                border_file_path);
            continue;
        end

        borders = border_data.data.borders;
        N = border_data.meta.N;
        if numel(borders) ~= 2
            warning('Expected two area borders but found %d: %s. Skipping.', ...
                numel(borders), border_file_path);
            continue;
        end
        borders = [borders, N + 1];

        required_last_col = 1 + N * max(kernel_indices);
        if required_last_col > size(model_par, 2)
            warning(['Kernel block %d requires columns through %d, but GLM has ' ...
                     '%d columns: %s. Skipping.'], ...
                max(kernel_indices), required_last_col, ...
                size(model_par, 2), file_path);
            continue;
        end

        selected_areas = {[1, 2], [2, 1]};
        for kernel_i = 1:n_kernel
            kernel_idx = kernel_indices(kernel_i);
            [pos, neg, total] = J_count( ...
                model_par, model_err, kernel_idx, borders, selected_areas);
            pos_counts(kernel_i) = pos_counts(kernel_i) + pos;
            neg_counts(kernel_i) = neg_counts(kernel_i) + neg;
            total_counts(kernel_i) = total_counts(kernel_i) + total;
        end
    end
end

function kernel_specs = build_default_kernel_specs()
    % Edit this function to add, remove, or reorder kernel configurations.
    %
    % A one-block model uses kernel_num=1. Composite models can use more
    % blocks and optionally provide custom x-axis labels.

    kernel_specs = struct('kernel_name', {}, 'kernel_num', {}, ...
        'kernel_labels', {});

    kernel_specs(end + 1) = make_kernel_spec( ...
        'DeltaPure', 3, {'Kernel 1', 'Kernel 2', 'Kernel 3'});

    first_kernel_names = {'Exp5', 'Exp40', 'Exp100', ...
        'Exp40Mix', 'Exp100Mix'};
    first_kernel_nums = [1, 1, 1, 2, 3];
    for i = 1:numel(first_kernel_names)
        kernel_specs(end + 1) = make_kernel_spec( ...
            first_kernel_names{i}, first_kernel_nums(i), {});
    end

    suffixes = {'5', '10', '20', '40', '80', '160', ...
        '320', '640', '1000'};
    prefixes = {'LongExp', 'LongGaussC', 'LongGDeriv'};
    for prefix_i = 1:numel(prefixes)
        for suffix_i = 1:numel(suffixes)
            kernel_specs(end + 1) = make_kernel_spec( ...
                sprintf('%s%s', prefixes{prefix_i}, suffixes{suffix_i}), ...
                1, {});
        end
    end

    step_suffixes = {'0_5', '5_10', '10_20', '20_40', ...
        '40_80', '80_160', '160_320', '320_640', ...
        '640_1000', '1000_3000'};
    for suffix_i = 1:numel(step_suffixes)
        kernel_specs(end + 1) = make_kernel_spec( ...
            sprintf('LongStepB%s', step_suffixes{suffix_i}), 1, {});
    end
end

function spec = make_kernel_spec(kernel_name, kernel_num, kernel_labels)
    if nargin < 3 || isempty(kernel_labels)
        kernel_labels = arrayfun(@(k) sprintf('Kernel %d', k), ...
            1:kernel_num, 'UniformOutput', false);
    end
    if kernel_num < 1 || kernel_num ~= floor(kernel_num)
        error('kernel_num must be a positive integer for %s.', kernel_name);
    end
    if numel(kernel_labels) ~= kernel_num
        error('kernel_labels count must equal kernel_num for %s.', kernel_name);
    end

    spec = struct();
    spec.kernel_name = char(string(kernel_name));
    spec.kernel_num = kernel_num;
    spec.kernel_labels = kernel_labels;
end

function configs = build_analysis_configs( ...
    kernel_specs, reg_names, alignments, resting_dur_thresholds)

    configs = struct('kernel_name', {}, 'kernel_num', {}, ...
        'kernel_labels', {}, 'reg_name', {}, 'align', {}, ...
        'resting_dur_threshold', {});

    config_i = 0;
    for kernel_i = 1:numel(kernel_specs)
        for reg_i = 1:numel(reg_names)
            for align_i = 1:numel(alignments)
                for duration_i = 1:numel(resting_dur_thresholds)
                    config_i = config_i + 1;
                    configs(config_i).kernel_name = ...
                        kernel_specs(kernel_i).kernel_name;
                    configs(config_i).kernel_num = ...
                        kernel_specs(kernel_i).kernel_num;
                    configs(config_i).kernel_labels = ...
                        kernel_specs(kernel_i).kernel_labels;
                    configs(config_i).reg_name = reg_names{reg_i};
                    configs(config_i).align = alignments{align_i};
                    configs(config_i).resting_dur_threshold = ...
                        resting_dur_thresholds(duration_i);
                end
            end
        end
    end
end

function selected = select_config_metadata(mt_glm, animal_filter, cfg)
    mask = animal_filter & ...
        metadata_table_matches(mt_glm, 'kernel_name', cfg.kernel_name) & ...
        metadata_table_matches(mt_glm, 'align', cfg.align) & ...
        metadata_table_matches(mt_glm, 'resting_dur_threshold', ...
            cfg.resting_dur_threshold);

    if ~is_empty_filter_value(cfg.reg_name)
        mask = mask & metadata_table_matches( ...
            mt_glm, 'reg_name', cfg.reg_name);
    end

    selected = mt_glm(mask, :);
end

function mask = build_animal_filter(mt_glm, animal_names)
    if is_empty_filter_value(animal_names)
        mask = true(height(mt_glm), 1);
    else
        mask = metadata_table_matches(mt_glm, 'animal_name', animal_names);
    end
end

function title_str = make_config_title(cfg)
    title_str = sprintf( ...
        'kernel=%s, kernel_num=%d, reg=%s, align=%s, restDur=%s', ...
        value_to_display_string(cfg.kernel_name, 'all'), ...
        cfg.kernel_num, value_to_display_string(cfg.reg_name, 'all'), ...
        value_to_display_string(cfg.align, 'all'), ...
        value_to_display_string(cfg.resting_dur_threshold, 'all'));
end

function folder_name = make_config_folder_name(cfg)
    parts = { ...
        ['kernel_', sanitize_for_path(value_to_display_string( ...
            cfg.kernel_name, 'all'))], ...
        ['kernelNum_', sanitize_for_path(value_to_display_string( ...
            cfg.kernel_num, 'all'))], ...
        ['reg_', sanitize_for_path(value_to_display_string( ...
            cfg.reg_name, 'all'))], ...
        ['align_', sanitize_for_path(value_to_display_string( ...
            cfg.align, 'all'))], ...
        ['restDur_', sanitize_for_path(value_to_display_string( ...
            cfg.resting_dur_threshold, 'all'))] ...
    };
    folder_name = strjoin(parts, '__');
end

function [mt_glm, index_info] = load_or_build_jbars_metadata_index( ...
    root, index_filename, expected_version, force_rebuild)

    check_path(fileparts(index_filename));

    if isfile(index_filename) && ~force_rebuild
        load_tic = tic;
        progress_log('INDEX', 'Loading existing compact J-bars index.');
        loaded = load(index_filename, 'metadata_index', 'index_info');

        if ~isfield(loaded, 'metadata_index') || ...
                ~isfield(loaded, 'index_info')
            error(['J-bars metadata index is missing metadata_index or ' ...
                   'index_info. Rebuild the index.']);
        end

        if ~isfield(loaded.index_info, 'version') || ...
                loaded.index_info.version ~= expected_version
            error(['J-bars metadata index version mismatch. Expected %d. ' ...
                   'Set force_rebuild_metadata_index=true.'], ...
                expected_version);
        end

        metadata_index = loaded.metadata_index;
        index_info = loaded.index_info;
        if istable(metadata_index)
            mt_glm = metadata_index;
        elseif isstruct(metadata_index)
            if isempty(metadata_index)
                mt_glm = table();
            else
                mt_glm = struct2table(metadata_index);
            end
        else
            error('metadata_index must be a table or struct array.');
        end

        if ~isfield(index_info, 'created_at')
            index_info.created_at = 'unknown';
        end
        progress_log('INDEX', ...
            'Loaded %d compact rows in %.1f s without full metadata.', ...
            height(mt_glm), toc(load_tic));
        return;
    end

    progress_log('INDEX', ...
        'Building compact J-bars index from full struct metadata.');
    [metadata_index, index_info] = ...
        build_jbars_metadata_index(root, expected_version);

    progress_log('INDEX', 'Saving compact index: %s', index_filename);
    save(index_filename, 'metadata_index', 'index_info', '-v7.3');

    if isempty(metadata_index)
        mt_glm = table();
    else
        mt_glm = struct2table(metadata_index);
    end
end

function [metadata_index, index_info] = ...
    build_jbars_metadata_index(root, index_version)

    target_epoch = 3000;
    chunk_size = 100000;

    load_tic = tic;
    loaded_meta = load_meta(root, 'struct');
    if ~isfield(loaded_meta, 'GLM') || ~isstruct(loaded_meta.GLM)
        error('load_meta(root, ''struct'') did not return GLM structs.');
    end
    records = loaded_meta.GLM;
    clear loaded_meta;
    source_row_count = numel(records);
    progress_log('INDEX-BUILD', ...
        'Loaded %d GLM metadata records in %.1f s.', ...
        source_row_count, toc(load_tic));

    required_fields = {'epoch', 'area', 'fold_idx', 'shuffle_idx', ...
        'state', 'prepost', 'resting_dur_threshold', 'animal_name', ...
        'injection', 'align', 'kernel_name', 'reg_name', 'file_name'};
    missing_fields = required_fields(~isfield(records, required_fields));
    if ~isempty(missing_fields)
        error('Metadata is missing required fields: %s', ...
            strjoin(missing_fields, ', '));
    end

    all_idx = (1:source_row_count).';
    epoch = struct_field_numeric(records, all_idx, 'epoch', chunk_size);
    epoch_idx = find(epoch == target_epoch);
    clear epoch;
    if isempty(epoch_idx)
        error('No metadata records found for epoch=%d.', target_epoch);
    end

    progress_log('INDEX-BUILD', ...
        'Extracting fixed-filter fields for %d epoch candidates.', ...
        numel(epoch_idx));
    area = struct_field_string(records, epoch_idx, 'area', chunk_size);
    fold_idx = struct_field_numeric(records, epoch_idx, 'fold_idx', chunk_size);
    shuffle_idx = struct_field_numeric(records, epoch_idx, 'shuffle_idx', chunk_size);
    state = struct_field_string(records, epoch_idx, 'state', chunk_size);
    prepost = struct_field_string(records, epoch_idx, 'prepost', chunk_size);

    keep = area == "Cortex" & ...
        fold_idx == 0 & ...
        shuffle_idx == 0 & ...
        ismember(state, ["RestOpen", "RestClose"]) & ...
        ismember(prepost, ["Pre", "Post"]);

    selected_source_idx = epoch_idx(keep);
    metadata_index = records(selected_source_idx);
    clear records;

    index_info = struct();
    index_info.version = index_version;
    index_info.format = 'struct';
    index_info.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    index_info.target_epoch = target_epoch;
    index_info.fixed_area = 'Cortex';
    index_info.fixed_fold_idx = 0;
    index_info.fixed_shuffle_idx = 0;
    index_info.source_row_count = source_row_count;
    index_info.epoch_candidate_count = numel(epoch_idx);
    index_info.index_row_count = numel(metadata_index);

    progress_log('INDEX-BUILD', ...
        'Compact index retained %d/%d source records.', ...
        numel(metadata_index), source_row_count);
end

function values = struct_field_numeric(records, indices, field, chunk_size)
    if nargin < 4 || isempty(chunk_size)
        chunk_size = 100000;
    end
    indices = indices(:);
    values = nan(numel(indices), 1);

    for first_i = 1:chunk_size:numel(indices)
        last_i = min(first_i + chunk_size - 1, numel(indices));
        source_idx = indices(first_i:last_i);
        target_idx = first_i:last_i;

        used_fast_path = false;
        try
            raw = [records(source_idx).(field)];
            if numel(raw) == numel(source_idx) && ...
                    (isnumeric(raw) || islogical(raw))
                values(target_idx) = double(raw(:));
                used_fast_path = true;
            end
        catch
            used_fast_path = false;
        end

        if ~used_fast_path
            for j = 1:numel(source_idx)
                value = unwrap_cell_scalar(records(source_idx(j)).(field));
                if isempty(value)
                    continue;
                elseif (isnumeric(value) || islogical(value)) && isscalar(value)
                    values(target_idx(j)) = double(value);
                elseif (ischar(value) || isstring(value))
                    parsed = str2double(string(value));
                    if isscalar(parsed) && ~isnan(parsed)
                        values(target_idx(j)) = parsed;
                    end
                end
            end
        end
    end
end

function values = struct_field_string(records, indices, field, chunk_size)
    if nargin < 4 || isempty(chunk_size)
        chunk_size = 100000;
    end
    indices = indices(:);
    values = strings(numel(indices), 1);

    for first_i = 1:chunk_size:numel(indices)
        last_i = min(first_i + chunk_size - 1, numel(indices));
        source_idx = indices(first_i:last_i);
        target_idx = first_i:last_i;

        for j = 1:numel(source_idx)
            value = unwrap_cell_scalar(records(source_idx(j)).(field));
            if isempty(value)
                values(target_idx(j)) = missing;
            else
                converted = string(value);
                if ~isscalar(converted)
                    error('Field %s is nonscalar at metadata row %d.', ...
                        field, source_idx(j));
                end
                values(target_idx(j)) = converted;
            end
        end
    end
end

function mask = metadata_table_matches(tbl, field, target)
    if ~ismember(field, tbl.Properties.VariableNames)
        error('Metadata table is missing required field: %s', field);
    end

    if is_empty_filter_value(target)
        mask = true(height(tbl), 1);
        return;
    end

    if iscell(target)
        mask = false(height(tbl), 1);
        for i = 1:numel(target)
            mask = mask | metadata_table_matches(tbl, field, target{i});
        end
        return;
    end

    if isstring(target) && numel(target) > 1
        mask = false(height(tbl), 1);
        for i = 1:numel(target)
            mask = mask | metadata_table_matches(tbl, field, target(i));
        end
        return;
    end

    mask = false(height(tbl), 1);
    for row_i = 1:height(tbl)
        value = get_table_scalar(tbl, field, row_i);
        mask(row_i) = metadata_scalar_matches(value, target);
    end
end

function value = get_table_scalar(tbl, field, row_i)
    column = tbl.(field);
    if iscell(column)
        value = column{row_i};
    elseif ischar(column) && size(column, 1) == height(tbl)
        value = column(row_i, :);
    else
        value = column(row_i);
    end
end

function tf = metadata_scalar_matches(value, target)
    value = unwrap_cell_scalar(value);
    target = unwrap_cell_scalar(target);

    if ischar(value) || isstring(value) || ...
            ischar(target) || isstring(target) || ...
            iscategorical_safe(value) || iscategorical_safe(target)
        tf = strcmp(string(value), string(target));
    elseif isnumeric(value) && isnumeric(target)
        tf = isequal(double(value), double(target));
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
    elseif isstring(value)
        tf = all(strlength(value) == 0);
    elseif iscell(value)
        tf = isempty(value);
    else
        tf = false;
    end
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
        out = strjoin(cellstr(value(:).'), '+');
    elseif ischar(value)
        out = value;
    elseif isnumeric(value) || islogical(value)
        parts = arrayfun(@num2str, value(:).', 'UniformOutput', false);
        out = strjoin(parts, '+');
    else
        out = char(string(value));
    end
end

function out = sanitize_for_path(value)
    out = char(string(value));
    out = regexprep(out, '[^A-Za-z0-9._=-]+', '_');
    out = regexprep(out, '^_+|_+$', '');
    if isempty(out)
        out = 'empty';
    end
end

function progress_log(stage, format_string, varargin)
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    fprintf('[%s][%s] %s\n', timestamp, stage, ...
        sprintf(format_string, varargin{:}));
end

function issue_log = empty_run_issue_log()
    issue_log = struct('config_index', {}, 'config_title', {}, ...
        'issue_type', {}, 'message', {});
end

function issue_log = append_run_issue( ...
    issue_log, config_index, config_title, issue_type, message)
    entry_i = numel(issue_log) + 1;
    issue_log(entry_i).config_index = config_index;
    issue_log(entry_i).config_title = char(string(config_title));
    issue_log(entry_i).issue_type = char(string(issue_type));
    issue_log(entry_i).message = char(string(message));
end

function write_run_issue_log(filename, issue_log)
    fid = fopen(filename, 'w');
    if fid < 0
        warning('Could not open issue log: %s', filename);
        return;
    end
    cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fprintf(fid, 'config_index\tissue_type\tconfig_title\tmessage\n');
    for i = 1:numel(issue_log)
        fprintf(fid, '%d\t%s\t%s\t%s\n', ...
            issue_log(i).config_index, ...
            escape_log_field(issue_log(i).issue_type), ...
            escape_log_field(issue_log(i).config_title), ...
            escape_log_field(issue_log(i).message));
    end
end

function out = escape_log_field(value)
    out = char(string(value));
    out = strrep(out, sprintf('\t'), ' ');
    out = strrep(out, sprintf('\n'), ' ');
    out = strrep(out, sprintf('\r'), ' ');
end

function close_created_figures()
    figs = findall(0, 'Type', 'figure');
    for i = 1:numel(figs)
        if isgraphics(figs(i))
            close(figs(i));
        end
    end
end

function jbars_sgtitle(fig, base_title)
    config_title = '';
    if isappdata(0, 'JBarsConfigTitle')
        config_title = getappdata(0, 'JBarsConfigTitle');
    end

    if isempty(config_title)
        full_title = char(string(base_title));
    else
        full_title = sprintf('%s\n%s', char(string(base_title)), ...
            char(string(config_title)));
    end

    sgtitle(fig, full_title, 'FontSize', 14, 'Interpreter', 'none');
end

function plot_kernel_grouped_bars(tile, ratios, low_errs, high_errs, sig_labels, condition_labels, kernel_labels, bar_offset, bar_width)
    % ratios: kernel x condition
    n_kernel = size(ratios, 1);
    x = 1:n_kernel;

    hold(tile, "on");

    bar(tile, x - bar_offset, ratios(:, 1), 'FaceColor', 'b', 'BarWidth', bar_width);
    bar(tile, x + bar_offset, ratios(:, 2), 'FaceColor', 'r', 'BarWidth', bar_width);

    errorbar(tile, x - bar_offset, ratios(:, 1), low_errs(:, 1), high_errs(:, 1), ...
        'k', 'LineStyle', 'none');
    errorbar(tile, x + bar_offset, ratios(:, 2), low_errs(:, 2), high_errs(:, 2), ...
        'k', 'LineStyle', 'none');

    for kernel_i = 1:n_kernel
        y_text = max(ratios(kernel_i, :) + high_errs(kernel_i, :), [], 'omitnan') + 2;
        if ~isfinite(y_text)
            y_text = 2;
        end
        text(tile, x(kernel_i), y_text, sig_labels{kernel_i}, ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
    end

    hold(tile, "off");

    xticks(tile, x);
    xticklabels(tile, kernel_labels);
    legend(tile, condition_labels, 'Location', 'Best');
end

function [ratios, low_errs, high_errs, sig_labels] = compute_bar_stats(J_counts, total_counts)
    ratios_raw = J_counts ./ total_counts;
    CIs = zeros(2, 2, 2); % pos/neg x condition x low/high
    sig_labels = cell(1, 2);

    for posneg_idx = 1:2
        count = J_counts(posneg_idx, :);
        total = total_counts(posneg_idx, :);

        for condition_idx = 1:2
            [CI_low, CI_high] = wilsonCI(count(condition_idx), total(condition_idx), 0.05);
            CIs(posneg_idx, condition_idx, :) = [CI_low, CI_high];
        end

        p_val = twoProportionPValue(count(1), total(1), count(2), total(2));
        if p_val < 0.001
            sig_labels{posneg_idx} = '***';
        elseif p_val < 0.01
            sig_labels{posneg_idx} = '**';
        elseif p_val < 0.05
            sig_labels{posneg_idx} = '*';
        else
            sig_labels{posneg_idx} = 'N.S.';
        end
    end

    ratios = ratios_raw * 100;
    low_errs = ratios - squeeze(CIs(:, :, 1)) * 100;
    high_errs = squeeze(CIs(:, :, 2)) * 100 - ratios;
end

function [effect_pct, ci_low_pct, ci_high_pct] = prepost_effect_ci(stats, posneg_idx, alpha)
    if stats.pre_total == 0 || stats.post_total == 0
        effect_pct = NaN;
        ci_low_pct = NaN;
        ci_high_pct = NaN;
        return;
    end

    p_pre = stats.pre_counts(posneg_idx) / stats.pre_total;
    p_post = stats.post_counts(posneg_idx) / stats.post_total;
    effect = p_post - p_pre;
    se = sqrt(p_post * (1 - p_post) / stats.post_total + ...
              p_pre * (1 - p_pre) / stats.pre_total);
    zcrit = norminv(1 - alpha / 2);

    effect_pct = 100 * effect;
    ci_low_pct = 100 * (effect - zcrit * se);
    ci_high_pct = 100 * (effect + zcrit * se);
end

function [diff_pct, ci_low_pct, ci_high_pct, p_val] = prepost_effect_difference_ci(stats_a, stats_b, posneg_idx, alpha)
    % stats_a - stats_b, where each effect is Post - Pre.
    if stats_a.pre_total == 0 || stats_a.post_total == 0 || stats_b.pre_total == 0 || stats_b.post_total == 0
        diff_pct = NaN;
        ci_low_pct = NaN;
        ci_high_pct = NaN;
        p_val = NaN;
        return;
    end

    a_pre = stats_a.pre_counts(posneg_idx) / stats_a.pre_total;
    a_post = stats_a.post_counts(posneg_idx) / stats_a.post_total;
    b_pre = stats_b.pre_counts(posneg_idx) / stats_b.pre_total;
    b_post = stats_b.post_counts(posneg_idx) / stats_b.post_total;

    diff = (a_post - a_pre) - (b_post - b_pre);
    se = sqrt(a_post * (1 - a_post) / stats_a.post_total + ...
              a_pre * (1 - a_pre) / stats_a.pre_total + ...
              b_post * (1 - b_post) / stats_b.post_total + ...
              b_pre * (1 - b_pre) / stats_b.pre_total);

    if se == 0
        if diff == 0
            p_val = 1;
        else
            p_val = 0;
        end
        ci_low = diff;
        ci_high = diff;
    else
        z = diff / se;
        p_val = twoSidedNormalPValue(z);
        zcrit = norminv(1 - alpha / 2);
        ci_low = diff - zcrit * se;
        ci_high = diff + zcrit * se;
    end

    diff_pct = 100 * diff;
    ci_low_pct = 100 * ci_low;
    ci_high_pct = 100 * ci_high;
end

function p_val = twoSidedNormalPValue(z)
    if ~isfinite(z)
        p_val = NaN;
        return;
    end
    abs_z = abs(z);
    Phi = 0.5 * erfc(-abs_z / sqrt(2));
    p_val = 2 * (1 - Phi);
end

function label = p_to_sig_label(p_val)
    if ~isfinite(p_val)
        label = 'n/a';
    elseif p_val < 0.001
        label = '***';
    elseif p_val < 0.01
        label = '**';
    elseif p_val < 0.05
        label = '*';
    else
        label = 'N.S.';
    end
end

function plot_effect_grouped_bars(tile, effects, low_errs, high_errs, sig_labels, condition_labels, kernel_labels, bar_offset, bar_width, y_limit)
    n_kernel = size(effects, 1);
    x = 1:n_kernel;

    hold(tile, 'on');
    bar(tile, x - bar_offset, effects(:, 1), 'FaceColor', 'b', 'BarWidth', bar_width);
    bar(tile, x + bar_offset, effects(:, 2), 'FaceColor', 'r', 'BarWidth', bar_width);

    errorbar(tile, x - bar_offset, effects(:, 1), low_errs(:, 1), high_errs(:, 1), ...
        'k', 'LineStyle', 'none');
    errorbar(tile, x + bar_offset, effects(:, 2), low_errs(:, 2), high_errs(:, 2), ...
        'k', 'LineStyle', 'none');

    set_effect_ylim(tile, effects, low_errs, high_errs, y_limit);
    yl = ylim(tile);
    y_offset = 0.04 * diff(yl);
    for kernel_i = 1:n_kernel
        y_text = max([effects(kernel_i, :) + high_errs(kernel_i, :), 0], [], 'omitnan') + y_offset;
        if ~isfinite(y_text)
            y_text = yl(2) - y_offset;
        end
        text(tile, x(kernel_i), y_text, sig_labels{kernel_i}, ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
    end

    hold(tile, 'off');
    xticks(tile, x);
    xticklabels(tile, kernel_labels);
    legend(tile, condition_labels, 'Location', 'Best');
end

function plot_effect_difference_bars(tile, diff_effects, low_errs, high_errs, sig_labels, kernel_labels, y_limit)
    n_kernel = numel(diff_effects);
    x = 1:n_kernel;

    hold(tile, 'on');
    bar(tile, x, diff_effects, 'FaceColor', [0.4, 0.4, 0.4], 'BarWidth', 0.45);
    errorbar(tile, x, diff_effects, low_errs, high_errs, 'k', 'LineStyle', 'none');

    set_effect_ylim(tile, diff_effects(:), low_errs(:), high_errs(:), y_limit);
    yl = ylim(tile);
    y_offset = 0.04 * diff(yl);
    for kernel_i = 1:n_kernel
        y_text = max([diff_effects(kernel_i) + high_errs(kernel_i), 0], [], 'omitnan') + y_offset;
        if ~isfinite(y_text)
            y_text = yl(2) - y_offset;
        end
        text(tile, x(kernel_i), y_text, sig_labels{kernel_i}, ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
    end

    hold(tile, 'off');
    xticks(tile, x);
    xticklabels(tile, kernel_labels);
end

function set_effect_ylim(tile, values, low_errs, high_errs, y_limit)
    if ~isempty(y_limit)
        ylim(tile, y_limit);
        return;
    end

    y_lower = values(:) - low_errs(:);
    upper = values(:) + high_errs(:);
    vals = [y_lower; upper; 0];
    vals = vals(isfinite(vals));
    if isempty(vals)
        ylim(tile, [-1, 1]);
        return;
    end
    yr = max(vals) - min(vals);
    if yr == 0
        yr = max(abs(vals));
    end
    if yr == 0
        yr = 1;
    end
    ylim(tile, [min(vals) - 0.15 * yr, max(vals) + 0.25 * yr]);
end

function draw_zero_line(tile)
    xl = xlim(tile);
    hold(tile, 'on');
    plot(tile, xl, [0, 0], 'k-', 'LineWidth', 0.75);
    xlim(tile, xl);
    hold(tile, 'off');
end

function [p_low, p_high] = wilsonCI(M, N, alpha)
    if nargin < 3
        alpha = 0.05;
    end

    if N == 0
        p_low = NaN;
        p_high = NaN;
        return;
    end

    p = M / N;
    z = norminv(1 - alpha/2);

    denominator = 1 + (z^2)/N;
    center = p + (z^2)/(2*N);
    radius = z * sqrt((p*(1-p)/N) + (z^2)/(4*N^2));

    p_low = (center - radius) / denominator;
    p_high = (center + radius) / denominator;
end

function p_val = twoProportionPValue(success1, trials1, success2, trials2)
    if trials1 == 0 || trials2 == 0
        p_val = NaN;
        return;
    end

    p1 = success1 / trials1;
    p2 = success2 / trials2;
    pooled = (success1 + success2) / (trials1 + trials2);
    denom = sqrt(pooled * (1 - pooled) * (1 / trials1 + 1 / trials2));
    if denom == 0
        p_val = 1;
        return;
    end

    z = (p1 - p2) / denom;
    abs_z = abs(z);
    Phi = 0.5 * erfc(-abs_z / sqrt(2));
    p_val = 2 * (1 - Phi);
end


%% Extended pooled effect plotting functions.

function plot_effect_metric_combined_pooled_ext(root, stats_grid, states, state_labels, posneg_labels, ...
    kernel_indices, kernel_labels, injection_types, injection_labels, bar_offset, bar_width, metric_mode, y_limit, alpha)

    if numel(injection_types) ~= 2 || ~strcmp(injection_types{1}, 'Saline') || ~strcmp(injection_types{2}, 'Muscimol')
        error('Combined effect figure expects injection_types = {''Saline'', ''Muscimol''}.');
    end

    [values, low_errs, high_errs, prepost_sig_labels, did_values, did_low_errs, did_high_errs, did_sig_labels] = ...
        compute_metric_arrays_from_stats_grid_ext(stats_grid, states, posneg_labels, kernel_indices, injection_types, metric_mode, alpha);

    if isempty(y_limit)
        y_limit = common_ylim_for_metric_ext(values, low_errs, high_errs, metric_mode);
    end

    f = figure('Position', [100, 100, 1200, 800], 'Visible', get_jbars_figure_visible());
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for state_idx = 1:numel(states)
        for posneg_idx = 1:2
            tile = nexttile(t, (state_idx - 1) * 2 + posneg_idx);
            plot_metric_grouped_bars_ext(tile, squeeze(values(state_idx, posneg_idx, :, :)), ...
                squeeze(low_errs(state_idx, posneg_idx, :, :)), squeeze(high_errs(state_idx, posneg_idx, :, :)), ...
                squeeze(prepost_sig_labels(state_idx, posneg_idx, :, :)), squeeze(did_sig_labels(state_idx, posneg_idx, :)), ...
                injection_labels, kernel_labels, bar_offset, bar_width, y_limit, metric_mode);
            ylabel(tile, metric_axis_label_ext(metric_mode));
            title(tile, sprintf('%s, %s', state_labels{state_idx}, posneg_labels{posneg_idx}));
        end
    end

    jbars_sgtitle(f, sprintf('Pooled %s: Saline vs Muscimol', metric_axis_label_ext(metric_mode)));
    export_figure(root, f, sprintf('J_bars_%s_effect_pooled_saline_muscimol_combined', metric_output_stub_ext(metric_mode)));
end

function plot_effect_metric_single_injection_pooled_ext(root, stats_grid, states, state_labels, posneg_labels, ...
    kernel_indices, kernel_labels, injection_types, injection_labels, injection_idx, metric_mode, y_limit, alpha)

    [values, low_errs, high_errs, prepost_sig_labels] = ...
        compute_metric_arrays_from_stats_grid_ext(stats_grid, states, posneg_labels, kernel_indices, injection_types, metric_mode, alpha);

    values_one = values(:, :, :, injection_idx);
    low_one = low_errs(:, :, :, injection_idx);
    high_one = high_errs(:, :, :, injection_idx);
    labels_one = prepost_sig_labels(:, :, :, injection_idx);

    if isempty(y_limit)
        y_limit = common_ylim_for_metric_ext(values_one, low_one, high_one, metric_mode);
    end

    f = figure('Position', [100, 100, 1200, 800], 'Visible', get_jbars_figure_visible());
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    injection_label = injection_labels{injection_idx};

    for state_idx = 1:numel(states)
        for posneg_idx = 1:2
            tile = nexttile(t, (state_idx - 1) * 2 + posneg_idx);
            plot_metric_single_bars_ext(tile, squeeze(values_one(state_idx, posneg_idx, :)), ...
                squeeze(low_one(state_idx, posneg_idx, :)), squeeze(high_one(state_idx, posneg_idx, :)), ...
                squeeze(labels_one(state_idx, posneg_idx, :)), kernel_labels, y_limit, metric_mode, injection_label);
            ylabel(tile, metric_axis_label_ext(metric_mode));
            title(tile, sprintf('%s, %s', state_labels{state_idx}, posneg_labels{posneg_idx}));
        end
    end

    jbars_sgtitle(f, sprintf('Pooled %s: %s only', metric_axis_label_ext(metric_mode), injection_label));
    export_figure(root, f, sprintf('J_bars_%s_effect_pooled_%s_only', metric_output_stub_ext(metric_mode), lower(injection_types{injection_idx})));
end

function plot_prepost_percent_injection_pooled_ext(root, stats_grid, states, state_labels, posneg_labels, ...
    kernel_indices, kernel_labels, injection_types, injection_labels, injection_idx, bar_offset, bar_width, alpha)

    n_state = numel(states);
    n_sign = 2;
    n_kernel = numel(kernel_indices);
    ratios = nan(n_state, n_sign, n_kernel, 2);
    low_errs = nan(n_state, n_sign, n_kernel, 2);
    high_errs = nan(n_state, n_sign, n_kernel, 2);
    sig_labels = cell(n_state, n_sign, n_kernel);

    for state_idx = 1:n_state
        for kernel_i = 1:n_kernel
            stats = stats_grid{state_idx, kernel_i, injection_idx};
            for posneg_idx = 1:n_sign
                counts = [stats.pre_counts(posneg_idx), stats.post_counts(posneg_idx)];
                totals = [stats.pre_total, stats.post_total];
                for pp_idx = 1:2
                    [ci_low, ci_high] = wilsonCI(counts(pp_idx), totals(pp_idx), alpha);
                    p = counts(pp_idx) / totals(pp_idx);
                    ratios(state_idx, posneg_idx, kernel_i, pp_idx) = 100 * p;
                    low_errs(state_idx, posneg_idx, kernel_i, pp_idx) = 100 * (p - ci_low);
                    high_errs(state_idx, posneg_idx, kernel_i, pp_idx) = 100 * (ci_high - p);
                end
                p_val = twoProportionPValue(stats.pre_counts(posneg_idx), stats.pre_total, stats.post_counts(posneg_idx), stats.post_total);
                sig_labels{state_idx, posneg_idx, kernel_i} = p_to_sig_label(p_val);
            end
        end
    end

    y_limit = common_ylim_percent_ext(ratios, low_errs, high_errs);
    f = figure('Position', [100, 100, 1200, 800], 'Visible', get_jbars_figure_visible());
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    injection_label = injection_labels{injection_idx};

    for state_idx = 1:n_state
        for posneg_idx = 1:n_sign
            tile = nexttile(t, (state_idx - 1) * 2 + posneg_idx);
            plot_prepost_grouped_bars_ext(tile, squeeze(ratios(state_idx, posneg_idx, :, :)), ...
                squeeze(low_errs(state_idx, posneg_idx, :, :)), squeeze(high_errs(state_idx, posneg_idx, :, :)), ...
                squeeze(sig_labels(state_idx, posneg_idx, :)), {'Pre', 'Post'}, kernel_labels, bar_offset, bar_width, y_limit);
            ylabel(tile, 'Significant J %');
            title(tile, sprintf('%s, %s', state_labels{state_idx}, posneg_labels{posneg_idx}));
        end
    end

    jbars_sgtitle(f, sprintf('Pooled Pre vs Post: %s', injection_label));
    export_figure(root, f, sprintf('J_bars_prepost_percent_pooled_%s', lower(injection_types{injection_idx})));
end

function plot_muscimol_minus_saline_effect_pooled_ext(root, stats_grid, states, state_labels, posneg_labels, ...
    kernel_indices, kernel_labels, injection_types, metric_mode, y_limit, alpha)

    if numel(injection_types) ~= 2 || ~strcmp(injection_types{1}, 'Saline') || ~strcmp(injection_types{2}, 'Muscimol')
        error('DID effect figure expects injection_types = {''Saline'', ''Muscimol''}.');
    end

    [~, ~, ~, ~, did_values, did_low_errs, did_high_errs, did_sig_labels] = ...
        compute_metric_arrays_from_stats_grid_ext(stats_grid, states, posneg_labels, kernel_indices, injection_types, metric_mode, alpha);

    if isempty(y_limit)
        y_limit = common_ylim_for_metric_ext(did_values, did_low_errs, did_high_errs, metric_mode);
    end

    f = figure('Position', [100, 100, 1200, 800], 'Visible', get_jbars_figure_visible());
    t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    for state_idx = 1:numel(states)
        for posneg_idx = 1:2
            tile = nexttile(t, (state_idx - 1) * 2 + posneg_idx);
            plot_metric_single_bars_ext(tile, squeeze(did_values(state_idx, posneg_idx, :)), ...
                squeeze(did_low_errs(state_idx, posneg_idx, :)), squeeze(did_high_errs(state_idx, posneg_idx, :)), ...
                squeeze(did_sig_labels(state_idx, posneg_idx, :)), kernel_labels, y_limit, metric_mode, 'DID');
            ylabel(tile, did_axis_label_ext(metric_mode));
            title(tile, sprintf('%s, %s', state_labels{state_idx}, posneg_labels{posneg_idx}));
        end
    end
    jbars_sgtitle(f, sprintf('Pooled DID: Muscimol vs Saline (%s)', metric_axis_label_ext(metric_mode)));
    export_figure(root, f, sprintf('J_bars_did_%s_pooled_muscimol_vs_saline', metric_output_stub_ext(metric_mode)));
end

function [values, low_errs, high_errs, prepost_sig_labels, did_values, did_low_errs, did_high_errs, did_sig_labels] = ...
    compute_metric_arrays_from_stats_grid_ext(stats_grid, states, posneg_labels, kernel_indices, injection_types, metric_mode, alpha)

    n_state = numel(states);
    n_sign = numel(posneg_labels);
    n_kernel = numel(kernel_indices);
    n_injection = numel(injection_types);

    values = nan(n_state, n_sign, n_kernel, n_injection);
    low_errs = nan(n_state, n_sign, n_kernel, n_injection);
    high_errs = nan(n_state, n_sign, n_kernel, n_injection);
    prepost_sig_labels = cell(n_state, n_sign, n_kernel, n_injection);
    did_values = nan(n_state, n_sign, n_kernel);
    did_low_errs = nan(n_state, n_sign, n_kernel);
    did_high_errs = nan(n_state, n_sign, n_kernel);
    did_sig_labels = cell(n_state, n_sign, n_kernel);

    for state_idx = 1:n_state
        for posneg_idx = 1:n_sign
            for kernel_i = 1:n_kernel
                for injection_idx = 1:n_injection
                    stats = stats_grid{state_idx, kernel_i, injection_idx};
                    [v, ci_low, ci_high, p_val] = prepost_metric_ci_ext(stats, posneg_idx, metric_mode, alpha);
                    values(state_idx, posneg_idx, kernel_i, injection_idx) = v;
                    low_errs(state_idx, posneg_idx, kernel_i, injection_idx) = v - ci_low;
                    high_errs(state_idx, posneg_idx, kernel_i, injection_idx) = ci_high - v;
                    prepost_sig_labels{state_idx, posneg_idx, kernel_i, injection_idx} = p_to_sig_label(p_val);
                end
                [v, ci_low, ci_high, p_val] = did_metric_ci_ext(stats_grid{state_idx, kernel_i, 2}, stats_grid{state_idx, kernel_i, 1}, posneg_idx, metric_mode, alpha);
                did_values(state_idx, posneg_idx, kernel_i) = v;
                did_low_errs(state_idx, posneg_idx, kernel_i) = v - ci_low;
                did_high_errs(state_idx, posneg_idx, kernel_i) = ci_high - v;
                did_sig_labels{state_idx, posneg_idx, kernel_i} = p_to_sig_label(p_val);
            end
        end
    end
end

function [value, ci_low, ci_high, p_val] = prepost_metric_ci_ext(stats, posneg_idx, metric_mode, alpha)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            [value, ci_low, ci_high] = prepost_effect_ci(stats, posneg_idx, alpha);
            p_val = twoProportionPValue(stats.pre_counts(posneg_idx), stats.pre_total, stats.post_counts(posneg_idx), stats.post_total);
        case {'ratio', 'post/pre', 'post_over_pre'}
            [value, ci_low, ci_high, p_val] = prepost_ratio_ci_ext(stats, posneg_idx, alpha);
        otherwise
            error('Unknown metric_mode: %s', metric_mode);
    end
end

function [rr, ci_low, ci_high, p_val] = prepost_ratio_ci_ext(stats, posneg_idx, alpha)
    [post_count, post_total, pre_count, pre_total] = safe_counts_for_log_ratio_ext( ...
        stats.post_counts(posneg_idx), stats.post_total, stats.pre_counts(posneg_idx), stats.pre_total);
    if pre_total == 0 || post_total == 0 || pre_count <= 0 || post_count <= 0
        rr = NaN; ci_low = NaN; ci_high = NaN; p_val = NaN; return;
    end
    p_post = post_count / post_total;
    p_pre = pre_count / pre_total;
    rr = p_post / p_pre;
    log_rr = log(rr);
    se = sqrt(1 / post_count - 1 / post_total + 1 / pre_count - 1 / pre_total);
    if se == 0
        p_val = double(log_rr == 0);
        ci_low = rr; ci_high = rr; return;
    end
    zcrit = norminv(1 - alpha / 2);
    ci_low = exp(log_rr - zcrit * se);
    ci_high = exp(log_rr + zcrit * se);
    p_val = twoSidedNormalPValue(log_rr / se);
end

function [value, ci_low, ci_high, p_val] = did_metric_ci_ext(stats_a, stats_b, posneg_idx, metric_mode, alpha)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            [value, ci_low, ci_high, p_val] = prepost_effect_difference_ci(stats_a, stats_b, posneg_idx, alpha);
        case {'ratio', 'post/pre', 'post_over_pre'}
            [value, ci_low, ci_high, p_val] = did_ratio_ci_ext(stats_a, stats_b, posneg_idx, alpha);
        otherwise
            error('Unknown metric_mode: %s', metric_mode);
    end
end

function [ror, ci_low, ci_high, p_val] = did_ratio_ci_ext(stats_a, stats_b, posneg_idx, alpha)
    [a_post_count, a_post_total, a_pre_count, a_pre_total] = safe_counts_for_log_ratio_ext(stats_a.post_counts(posneg_idx), stats_a.post_total, stats_a.pre_counts(posneg_idx), stats_a.pre_total);
    [b_post_count, b_post_total, b_pre_count, b_pre_total] = safe_counts_for_log_ratio_ext(stats_b.post_counts(posneg_idx), stats_b.post_total, stats_b.pre_counts(posneg_idx), stats_b.pre_total);
    if any([a_post_count, a_pre_count, b_post_count, b_pre_count] <= 0) || any([a_post_total, a_pre_total, b_post_total, b_pre_total] == 0)
        ror = NaN; ci_low = NaN; ci_high = NaN; p_val = NaN; return;
    end
    log_rr_a = log((a_post_count / a_post_total) / (a_pre_count / a_pre_total));
    log_rr_b = log((b_post_count / b_post_total) / (b_pre_count / b_pre_total));
    log_ror = log_rr_a - log_rr_b;
    se = sqrt(1 / a_post_count - 1 / a_post_total + 1 / a_pre_count - 1 / a_pre_total + ...
              1 / b_post_count - 1 / b_post_total + 1 / b_pre_count - 1 / b_pre_total);
    ror = exp(log_ror);
    if se == 0
        p_val = double(log_ror == 0);
        ci_low = ror; ci_high = ror; return;
    end
    zcrit = norminv(1 - alpha / 2);
    ci_low = exp(log_ror - zcrit * se);
    ci_high = exp(log_ror + zcrit * se);
    p_val = twoSidedNormalPValue(log_ror / se);
end

function [num1, den1, num0, den0] = safe_counts_for_log_ratio_ext(num1, den1, num0, den0)
    if num1 == 0 || num0 == 0 || num1 == den1 || num0 == den0
        num1 = num1 + 0.5;
        num0 = num0 + 0.5;
        den1 = den1 + 1;
        den0 = den0 + 1;
    end
end

function plot_metric_grouped_bars_ext(tile, values, low_errs, high_errs, bar_sig_labels, did_sig_labels, condition_labels, kernel_labels, bar_offset, bar_width, y_limit, metric_mode)
    n_kernel = size(values, 1);
    x = 1:n_kernel;
    hold(tile, 'on');
    bar(tile, x - bar_offset, values(:, 1), 'FaceColor', 'b', 'BarWidth', bar_width, 'DisplayName', condition_labels{1});
    bar(tile, x + bar_offset, values(:, 2), 'FaceColor', 'r', 'BarWidth', bar_width, 'DisplayName', condition_labels{2});
    errorbar(tile, x - bar_offset, values(:, 1), low_errs(:, 1), high_errs(:, 1), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    errorbar(tile, x + bar_offset, values(:, 2), low_errs(:, 2), high_errs(:, 2), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    ylim(tile, y_limit);
    draw_metric_reference_line_ext(tile, metric_mode);
    yr = diff(y_limit);
    label_offset = 0.025 * yr;
    bracket_offset = 0.11 * yr;
    bracket_height = 0.025 * yr;
    for kernel_i = 1:n_kernel
        xpos = [x(kernel_i) - bar_offset, x(kernel_i) + bar_offset];
        for injection_idx = 1:2
            y_bar = values(kernel_i, injection_idx) + high_errs(kernel_i, injection_idx) + label_offset;
            text(tile, xpos(injection_idx), y_bar, bar_sig_labels{kernel_i, injection_idx}, ...
                'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', 'HandleVisibility', 'off');
        end
        y_bracket = max(values(kernel_i, :) + high_errs(kernel_i, :), [], 'omitnan') + bracket_offset;
        draw_sig_bracket_ext(tile, xpos(1), xpos(2), y_bracket, bracket_height, did_sig_labels{kernel_i});
    end
    hold(tile, 'off');
    xticks(tile, x);
    xticklabels(tile, kernel_labels);
    legend(tile, condition_labels, 'Location', 'Best');
end

function plot_metric_single_bars_ext(tile, values, low_errs, high_errs, sig_labels, kernel_labels, y_limit, metric_mode, bar_label)
    n_kernel = numel(values);
    x = 1:n_kernel;
    hold(tile, 'on');
    bar(tile, x, values, 'BarWidth', 0.45, 'DisplayName', char(bar_label));
    errorbar(tile, x, values, low_errs, high_errs, 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    ylim(tile, y_limit);
    draw_metric_reference_line_ext(tile, metric_mode);
    yr = diff(y_limit);
    label_offset = 0.035 * yr;
    for kernel_i = 1:n_kernel
        y_text = values(kernel_i) + high_errs(kernel_i) + label_offset;
        text(tile, x(kernel_i), y_text, sig_labels{kernel_i}, 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'HandleVisibility', 'off');
    end
    hold(tile, 'off');
    xticks(tile, x);
    xticklabels(tile, kernel_labels);
end

function plot_prepost_grouped_bars_ext(tile, ratios, low_errs, high_errs, sig_labels, condition_labels, kernel_labels, bar_offset, bar_width, y_limit)
    n_kernel = size(ratios, 1);
    x = 1:n_kernel;
    hold(tile, 'on');
    bar(tile, x - bar_offset, ratios(:, 1), 'FaceColor', 'b', 'BarWidth', bar_width, 'DisplayName', condition_labels{1});
    bar(tile, x + bar_offset, ratios(:, 2), 'FaceColor', 'r', 'BarWidth', bar_width, 'DisplayName', condition_labels{2});
    errorbar(tile, x - bar_offset, ratios(:, 1), low_errs(:, 1), high_errs(:, 1), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    errorbar(tile, x + bar_offset, ratios(:, 2), low_errs(:, 2), high_errs(:, 2), 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    ylim(tile, y_limit);
    yr = diff(y_limit);
    bracket_height = 0.025 * yr;
    bracket_offset = 0.06 * yr;
    for kernel_i = 1:n_kernel
        y_bracket = max(ratios(kernel_i, :) + high_errs(kernel_i, :), [], 'omitnan') + bracket_offset;
        draw_sig_bracket_ext(tile, x(kernel_i) - bar_offset, x(kernel_i) + bar_offset, y_bracket, bracket_height, sig_labels{kernel_i});
    end
    hold(tile, 'off');
    xticks(tile, x);
    xticklabels(tile, kernel_labels);
    legend(tile, condition_labels, 'Location', 'Best');
end

function draw_sig_bracket_ext(tile, x1, x2, y, h, label)
    plot(tile, [x1, x1, x2, x2], [y, y + h, y + h, y], 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(tile, mean([x1, x2]), y + h, label, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold', 'HandleVisibility', 'off');
end

function y_limit = common_ylim_for_metric_ext(values, low_errs, high_errs, metric_mode)
    y_lower = values(:) - low_errs(:);
    upper = values(:) + high_errs(:);
    switch lower(metric_mode)
        case {'ratio', 'post/pre', 'post_over_pre'}
            ref = 1;
        otherwise
            ref = 0;
    end
    vals = [y_lower; upper; ref];
    vals = vals(isfinite(vals));
    if isempty(vals)
        y_limit = [ref - 1, ref + 1];
        return;
    end
    y_min = min(vals);
    y_max = max(vals);
    yr = y_max - y_min;
    if yr == 0
        yr = max(abs([y_min, y_max, 1]));
    end
    if yr == 0
        yr = 1;
    end
    y_limit = [y_min - 0.18 * yr, y_max + 0.45 * yr];
    if strcmpi(metric_mode, 'ratio') || strcmpi(metric_mode, 'post/pre')
        y_limit(1) = max(0, y_limit(1));
    end
end

function y_limit = common_ylim_percent_ext(values, low_errs, high_errs)
    vals = [values(:) - low_errs(:); values(:) + high_errs(:); 0];
    vals = vals(isfinite(vals));
    if isempty(vals)
        y_limit = [0, 1];
        return;
    end
    y_min = min(vals);
    y_max = max(vals);
    yr = y_max - y_min;
    if yr == 0, yr = max(abs(vals)); end
    if yr == 0, yr = 1; end
    y_limit = [max(0, y_min - 0.10 * yr), y_max + 0.35 * yr];
end

function draw_metric_reference_line_ext(tile, metric_mode)
    xl = xlim(tile);
    hold(tile, 'on');
    switch lower(metric_mode)
        case {'ratio', 'post/pre', 'post_over_pre'}
            plot(tile, xl, [1, 1], 'k--', 'LineWidth', 0.9, 'HandleVisibility', 'off');
        otherwise
            plot(tile, xl, [0, 0], 'k-', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    end
    xlim(tile, xl);
end

function label = metric_axis_label_ext(metric_mode)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            label = 'Post - Pre significant J (percentage points)';
        case {'ratio', 'post/pre', 'post_over_pre'}
            label = 'Post / Pre significant J ratio';
        otherwise
            label = metric_mode;
    end
end

function label = did_axis_label_ext(metric_mode)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            label = '(Muscimol Post-Pre) - (Saline Post-Pre) (percentage points)';
        case {'ratio', 'post/pre', 'post_over_pre'}
            label = '(Muscimol Post/Pre) / (Saline Post/Pre)';
        otherwise
            label = metric_mode;
    end
end

function stub = metric_output_stub_ext(metric_mode)
    switch lower(metric_mode)
        case {'difference', 'post-pre', 'post_minus_pre'}
            stub = 'post_minus_pre';
        case {'ratio', 'post/pre', 'post_over_pre'}
            stub = 'post_over_pre';
        otherwise
            stub = regexprep(metric_mode, '[^A-Za-z0-9_-]', '_');
    end
end


function export_figure(root, fig, output_stub) %#ok<INUSD>
    output_folder = get_jbars_output_folder();
    check_path(output_folder);

    export_jpg = get_appdata_bool('JBarsExportJPG', true);
    export_pdf = get_appdata_bool('JBarsExportPDF', false);
    close_after = get_appdata_bool('JBarsCloseFigures', true);
    resolution = get_appdata_numeric('JBarsExportResolution', 300);

    fig_width = 12.0;
    fig_height = 8.0;

    set(fig, 'Units', 'inches');
    fig.Position(3:4) = [fig_width, fig_height];
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperSize', [fig_width, fig_height]);
    set(fig, 'PaperPosition', [0, 0, fig_width, fig_height]);
    set(fig, 'Color', 'w');

    if export_jpg
        jpg_filename = fullfile(output_folder, ...
            [output_stub, '_preview.jpg']);
        progress_log('EXPORT', 'Writing JPG: %s', jpg_filename);
        exportgraphics(fig, jpg_filename, ...
            'ContentType', 'image', ...
            'BackgroundColor', 'white', ...
            'Resolution', resolution);
    end

    if export_pdf
        pdf_filename = fullfile(output_folder, [output_stub, '.pdf']);
        progress_log('EXPORT', 'Writing PDF: %s', pdf_filename);
        exportgraphics(fig, pdf_filename, ...
            'ContentType', 'vector', ...
            'BackgroundColor', 'white');
    end

    if ~export_jpg && ~export_pdf
        progress_log('EXPORT', ...
            'JPG and PDF exports are disabled for %s.', output_stub);
    end

    if close_after
        if isgraphics(fig)
            close(fig);
        end
    else
        set(fig, 'Visible', 'on');
    end
end

function figure_visible = get_jbars_figure_visible()
    if isappdata(0, 'JBarsFigureVisible')
        figure_visible = getappdata(0, 'JBarsFigureVisible');
    else
        figure_visible = 'off';
    end
end

function output_folder = get_jbars_output_folder()
    if isappdata(0, 'JBarsOutputFolder')
        output_folder = getappdata(0, 'JBarsOutputFolder');
    else
        error('JBarsOutputFolder has not been configured.');
    end
end

function value = get_appdata_bool(name, default_value)
    if isappdata(0, name)
        value = logical(getappdata(0, name));
    else
        value = default_value;
    end
end

function value = get_appdata_numeric(name, default_value)
    if isappdata(0, name)
        value = double(getappdata(0, name));
    else
        value = default_value;
    end
end
