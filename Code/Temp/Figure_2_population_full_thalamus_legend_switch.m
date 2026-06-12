%% Figure 2: Full-area Thalamus-ACC / Thalamus-VLPFC population analyses.
% Similar to the all-trial Figure 2 population script, but using Full-area
% rasters that include ACC, VLPFC, and Thalamus.
%
% Layout:
%   Rows 1-2: Thalamus vs ACC, RestOpen / RestClose
%   Rows 3-4: Thalamus vs VLPFC, RestOpen / RestClose
%   Only Pre is used (Full-area data do not have Post here).

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

%% Session definition
% Change these fields to generate a different session.
meta = struct();
meta.animal_name = 'Slayer';
meta.injection = 'Muscimol';
meta.align = 'Longest';
meta.session_idx = 6;
meta.resting_dur_threshold = 15;

%% Conditions shown in the figure
meta_area = 'Full';
pair_left_areas  = {'Thalamus', 'Thalamus', 'Thalamus', 'Thalamus'};
pair_right_areas = {'ACC',      'ACC',      'VLPFC',    'VLPFC'};
pair_labels      = {'Thalamus-ACC', 'Thalamus-ACC', 'Thalamus-VLPFC', 'Thalamus-VLPFC'};
preposts         = {'Pre',      'Pre',       'Pre',      'Pre'};
states           = {'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};
n_state = numel(states);

%% Figure layout
% Column 1 = pair population firing-rate traces.
% Column 2 = population cross-correlogram.
% Column 3 = population auto-correlograms.
% Column 4 = population cross-PSD.
% Column 5 = pair population PSDs.
% Column 6 = pair spectral coherence.
n_column = 6;
figure_visible = 'off';
show_legend = false;

%% Parameters
shuffle_N = 2;
std_multiplier = 2;
sig_min_run_corr = 50;
sig_min_run_spec = 15;

analysis_trial_mode = 'all'; % 'all' or 'single'.
analysis_trial_idx = 1;      % used only when analysis_trial_mode = 'single'.
analysis_t_range = [];       % [] means full length of each selected trial.
lag_weight_correction = true; % true: pooled Pearson per lag; false: duration-weighted trial-wise correlogram.

display_trial_idx = 1;
display_t_range = 1:2000; % [] means full length of display_trial_idx.

corr_range = 1500; % ms.
smooth_window = 25; % ms.
sample_rate = 1000; % Hz.
freqs = linspace(0, 150, 301); % Hz.
spec_smooth_window = 15;
psd_window_sec = 10;
psd_overlap_frac = 0.5;

smooth_kernel = 1 - abs(-smooth_window:smooth_window) / smooth_window;
smooth_kernel = smooth_kernel / sum(smooth_kernel);

%% Figure
f = figure('Color', 'w', 'Visible', figure_visible);
tiledlayout(n_state, n_column, "TileSpacing", "Compact", "Padding", "Compact");
col_axes = gobjects(n_state, n_column);

%% Column 1-6: population firing rate, correlation, and spectra
for i = 1:n_state
    meta.area = meta_area;
    meta.prepost = preposts{i};
    meta.state = states{i};
    meta.filename = generate_filename('raster', meta);

    left_area = pair_left_areas{i};
    right_area = pair_right_areas{i};
    pair_label = pair_labels{i};

    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s | %s | %s (%s)\n', pair_label, meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    rasters = raster_data.data.rasters;
    cell_area = raster_data.data.cell_area;
    filter_left = get_area_filter(cell_area, left_area);
    filter_right = get_area_filter(cell_area, right_area);
    check_pair_area_filters(filter_left, filter_right, left_area, right_area, meta);

    N_left = sum(filter_left);
    N_right = sum(filter_right);

    [left_trials, right_trials] = get_population_trial_segments( ...
        rasters, filter_left, filter_right, ...
        analysis_trial_mode, analysis_trial_idx, analysis_t_range, corr_range);

    [display_time_ms, left_display_plot, right_display_plot] = get_population_display_trace( ...
        rasters, filter_left, filter_right, display_trial_idx, display_t_range, ...
        sample_rate, smooth_kernel);

    fprintf('Population firing rates for %s | %s %s %s:\n', pair_label, meta.prepost, meta.state, meta.area);
    fprintf('%s population mean: %.3f Hz/neuron, N=%d, trials=%d\n', ...
        left_area, weighted_signal_mean(left_trials, sample_rate), N_left, numel(left_trials));
    fprintf('%s population mean: %.3f Hz/neuron, N=%d, trials=%d\n', ...
        right_area, weighted_signal_mean(right_trials, sample_rate), N_right, numel(right_trials));

    %% Column 1: pair population firing-rate traces.
    column = 1;
    tile = nexttile(column + n_column*(i-1));
    col_axes(i, column) = tile;
    plot(tile, display_time_ms, left_display_plot, 'r-', 'LineWidth', 1, ...
        'DisplayName', sprintf('%s, N=%d', left_area, N_left));
    hold(tile, 'on');
    plot(tile, display_time_ms, right_display_plot, 'b-', 'LineWidth', 1, ...
        'DisplayName', sprintf('%s, N=%d', right_area, N_right));
    hold(tile, 'off');
    title(tile, sprintf('%s firing rate: %s', pair_label, meta.state));
    xlabel(tile, 'Time (ms)');
    ylabel(tile, 'Firing rate (Hz/neuron)');
    set_legend_visibility(tile, show_legend, 'northeast');

    %% Compute all-trial population correlations and spectra.
    [correlogram, lags] = aggregate_norm_xcorr_trials(left_trials, right_trials, corr_range, lag_weight_correction);
    [auto_left, ~] = aggregate_norm_xcorr_trials(left_trials, left_trials, corr_range, lag_weight_correction);
    [auto_right, ~] = aggregate_norm_xcorr_trials(right_trials, right_trials, corr_range, lag_weight_correction);

    smooth_correlogram = same_conv(correlogram, smooth_kernel);
    smooth_auto_left = same_conv(auto_left, smooth_kernel);
    smooth_auto_right = same_conv(auto_right, smooth_kernel);

    [psd_left, psd_right, cpsd, coherence, spec_freqs] = aggregate_pair_spectra_trials( ...
        left_trials, right_trials, sample_rate, freqs, psd_window_sec, psd_overlap_frac);

    coherence_plot = smooth_spectrum_for_plot(coherence, spec_smooth_window);
    cpsd_plot = smooth_spectrum_for_plot(cpsd, spec_smooth_window);
    psd_left_plot = smooth_spectrum_for_plot(psd_left, spec_smooth_window);
    psd_right_plot = smooth_spectrum_for_plot(psd_right, spec_smooth_window);

    %% Compute both shuffled and shifted controls.
    shuffle_correlograms = zeros(shuffle_N, length(correlogram));
    shuffle_auto_left_controls = zeros(shuffle_N, length(auto_left));
    shuffle_auto_right_controls = zeros(shuffle_N, length(auto_right));
    shuffle_cpsd = zeros(shuffle_N, numel(spec_freqs));
    shuffle_coh = zeros(shuffle_N, numel(spec_freqs));
    shuffle_psd_left = zeros(shuffle_N, numel(spec_freqs));
    shuffle_psd_right = zeros(shuffle_N, numel(spec_freqs));

    shifted_correlograms = zeros(shuffle_N, length(correlogram));
    shifted_auto_left_controls = zeros(shuffle_N, length(auto_left));
    shifted_auto_right_controls = zeros(shuffle_N, length(auto_right));
    shifted_cpsd = zeros(shuffle_N, numel(spec_freqs));
    shifted_coh = zeros(shuffle_N, numel(spec_freqs));

    for j = 1:shuffle_N
        fprintf('%d/%d population controls started\n', j, shuffle_N);

        [shuffled_corr, shuffled_auto_left, shuffled_auto_right, ...
            shuffled_psd_left, shuffled_psd_right, shuffled_cpsd, shuffled_coh_j] = ...
            aggregate_control_metrics_trials(left_trials, right_trials, corr_range, lag_weight_correction, ...
                sample_rate, freqs, psd_window_sec, psd_overlap_frac, 'shuffled');

        shuffle_correlograms(j, :) = same_conv(shuffled_corr, smooth_kernel);
        shuffle_auto_left_controls(j, :) = same_conv(shuffled_auto_left, smooth_kernel);
        shuffle_auto_right_controls(j, :) = same_conv(shuffled_auto_right, smooth_kernel);
        shuffle_cpsd(j, :) = smooth_spectrum_for_plot(shuffled_cpsd, spec_smooth_window);
        shuffle_coh(j, :) = smooth_spectrum_for_plot(shuffled_coh_j, spec_smooth_window);
        shuffle_psd_left(j, :) = smooth_spectrum_for_plot(shuffled_psd_left, spec_smooth_window);
        shuffle_psd_right(j, :) = smooth_spectrum_for_plot(shuffled_psd_right, spec_smooth_window);

        [shifted_corr, shifted_auto_left, shifted_auto_right, ...
            ~, ~, shifted_cpsd_j, shifted_coh_j] = ...
            aggregate_control_metrics_trials(left_trials, right_trials, corr_range, lag_weight_correction, ...
                sample_rate, freqs, psd_window_sec, psd_overlap_frac, 'shifted');

        shifted_correlograms(j, :) = same_conv(shifted_corr, smooth_kernel);
        shifted_auto_left_controls(j, :) = same_conv(shifted_auto_left, smooth_kernel);
        shifted_auto_right_controls(j, :) = same_conv(shifted_auto_right, smooth_kernel);
        shifted_cpsd(j, :) = smooth_spectrum_for_plot(shifted_cpsd_j, spec_smooth_window);
        shifted_coh(j, :) = smooth_spectrum_for_plot(shifted_coh_j, spec_smooth_window);

        fprintf('%d/%d population controls finished\n', j, shuffle_N);
    end

    %% Control summary.
    shuffle_auto_left_mean = mean_omitnan(shuffle_auto_left_controls);
    shuffle_auto_left_std = std_omitnan(shuffle_auto_left_controls);

    shuffle_auto_right_mean = mean_omitnan(shuffle_auto_right_controls);
    shuffle_auto_right_std = std_omitnan(shuffle_auto_right_controls);

    cpsd_shuffle_mean = mean_omitnan(shuffle_cpsd);
    cpsd_shuffle_std = std_omitnan(shuffle_cpsd);

    coh_shuffle_mean = mean_omitnan(shuffle_coh);
    coh_shuffle_std = std_omitnan(shuffle_coh);

    psd_left_shuffle_mean = mean_omitnan(shuffle_psd_left);
    psd_left_shuffle_std = std_omitnan(shuffle_psd_left);

    psd_right_shuffle_mean = mean_omitnan(shuffle_psd_right);
    psd_right_shuffle_std = std_omitnan(shuffle_psd_right);

    shifted_ccg_mean = mean_omitnan(shifted_correlograms);
    shifted_ccg_std = std_omitnan(shifted_correlograms);

    shifted_auto_left_mean = mean_omitnan(shifted_auto_left_controls);
    shifted_auto_left_std = std_omitnan(shifted_auto_left_controls);

    shifted_auto_right_mean = mean_omitnan(shifted_auto_right_controls);
    shifted_auto_right_std = std_omitnan(shifted_auto_right_controls);

    cpsd_shifted_mean = mean_omitnan(shifted_cpsd);
    cpsd_shifted_std = std_omitnan(shifted_cpsd);

    shifted_coh_mean = mean_omitnan(shifted_coh);
    shifted_coh_std = std_omitnan(shifted_coh);

    %% Column 2: population cross-correlogram with shifted control.
    column = 2;
    tile = nexttile(column + n_column*(i-1));
    col_axes(i, column) = tile;
    shifted_ccg_upper = shifted_ccg_mean + std_multiplier * shifted_ccg_std;
    shifted_ccg_lower = shifted_ccg_mean - std_multiplier * shifted_ccg_std;

    fill_control_band(tile, lags, shifted_ccg_mean, shifted_ccg_std, std_multiplier, ...
        [0.8, 0.8, 0.8], 'Shifted mean ± 2SD');
    hold(tile, 'on');
    xline(tile, 0, 'k--', 'HandleVisibility', 'off');
    yline(tile, 0, 'k--', 'HandleVisibility', 'off');
    plot(tile, lags, smooth_correlogram, 'm-', 'LineWidth', 1, ...
        'DisplayName', sprintf('%s cross-corr', pair_label));

    ccg_sig_mask = (smooth_correlogram > shifted_ccg_upper) | ...
                   (smooth_correlogram < shifted_ccg_lower);
    plot_significant_segments(tile, lags, smooth_correlogram, ccg_sig_mask, ...
        sig_min_run_corr, 'k', 'Significant');
    hold(tile, 'off');

    title(tile, sprintf('%s correlogram: %s', pair_label, meta.state));
    xlabel(tile, 'Lag (ms)');
    ylabel(tile, 'Normalized correlation');
    set_legend_visibility(tile, show_legend, 'northeast');
    xlim(tile, [-corr_range, corr_range]);

    %% Column 3: pair population auto-correlograms with shifted self-controls.
    column = 3;
    tile = nexttile(column + n_column*(i-1));
    col_axes(i, column) = tile;

    fill_control_band(tile, lags, shifted_auto_left_mean, shifted_auto_left_std, ...
        std_multiplier, [1.0, 0.85, 0.85], '');
    hold(tile, 'on');
    fill_control_band(tile, lags, shifted_auto_right_mean, shifted_auto_right_std, ...
        std_multiplier, [0.85, 0.85, 1.0], '');

    plot(tile, lags, smooth_correlogram, 'm-', 'LineWidth', 2, ...
        'DisplayName', 'Cross-corr', 'Color', [1, 0, 1, 0.2]);
    xline(tile, 0, 'k--', 'HandleVisibility', 'off');
    yline(tile, 0, 'k--', 'HandleVisibility', 'off');
    plot(tile, lags, smooth_auto_left, 'r-', 'LineWidth', 1, ...
        'DisplayName', sprintf('%s auto-corr, N=%d', left_area, N_left));
    plot(tile, lags, smooth_auto_right, 'b-', 'LineWidth', 1, ...
        'DisplayName', sprintf('%s auto-corr, N=%d', right_area, N_right));

    auto_left_upper = shifted_auto_left_mean + std_multiplier * shifted_auto_left_std;
    auto_left_lower = shifted_auto_left_mean - std_multiplier * shifted_auto_left_std;
    auto_right_upper = shifted_auto_right_mean + std_multiplier * shifted_auto_right_std;
    auto_right_lower = shifted_auto_right_mean - std_multiplier * shifted_auto_right_std;

    auto_left_sig_mask = (smooth_auto_left > auto_left_upper) | ...
                         (smooth_auto_left < auto_left_lower);
    auto_right_sig_mask = (smooth_auto_right > auto_right_upper) | ...
                          (smooth_auto_right < auto_right_lower);

    plot_significant_segments(tile, lags, smooth_auto_left, auto_left_sig_mask, ...
        sig_min_run_corr, [0.5, 0, 0], '');
    plot_significant_segments(tile, lags, smooth_auto_right, auto_right_sig_mask, ...
        sig_min_run_corr, [0, 0, 0.5], '');
    hold(tile, 'off');

    title(tile, sprintf('%s auto-correlograms: %s', pair_label, meta.state));
    xlabel(tile, 'Lag (ms)');
    ylabel(tile, 'Normalized correlation');
    set_legend_visibility(tile, show_legend, 'southeast');
    xlim(tile, [-corr_range, corr_range]);

    %% Column 4: pair population cross-PSD with shifted control.
    column = 4;
    selected_control = 'shifted'; % 'shuffled' or 'shifted'
    if strcmp(selected_control, 'shuffled')
        cpsd_control_mean = cpsd_shuffle_mean;
        cpsd_control_std = cpsd_shuffle_std;
        control_label = 'Shuffled';
    elseif strcmp(selected_control, 'shifted')
        cpsd_control_mean = cpsd_shifted_mean;
        cpsd_control_std = cpsd_shifted_std;
        control_label = 'Shifted';
    else
        error('Invalid control type selected. Use "shuffled" or "shifted".');
    end

    tile = nexttile(column + n_column*(i-1));
    col_axes(i, column) = tile;
    fill_control_band(tile, spec_freqs, cpsd_control_mean, cpsd_control_std, ...
        std_multiplier, [0.8, 0.8, 0.8], sprintf('%s mean ± 2SD', control_label));
    hold(tile, 'on');
    plot(tile, spec_freqs, cpsd_plot, 'm-', 'LineWidth', 1, ...
        'DisplayName', 'Cross-PSD');

    cpsd_upper = cpsd_control_mean + std_multiplier * cpsd_control_std;
    cpsd_sig_mask = cpsd_plot > cpsd_upper;
    plot_significant_segments(tile, spec_freqs, cpsd_plot, cpsd_sig_mask, ...
        sig_min_run_spec, 'k', 'Significant');
    hold(tile, 'off');

    title(tile, sprintf('%s cross-PSD: %s', pair_label, meta.state));
    xlabel(tile, 'Frequency (Hz)');
    ylabel(tile, 'Cross-PSD');
    set_legend_visibility(tile, show_legend, 'northeast');

    %% Column 5: pair population PSDs with shuffled controls.
    column = 5;
    tile = nexttile(column + n_column*(i-1));
    col_axes(i, column) = tile;

    fill_control_band(tile, spec_freqs, psd_left_shuffle_mean, psd_left_shuffle_std, ...
        std_multiplier, [1.0, 0.85, 0.85], '');
    hold(tile, 'on');
    fill_control_band(tile, spec_freqs, psd_right_shuffle_mean, psd_right_shuffle_std, ...
        std_multiplier, [0.85, 0.85, 1.0], '');

    plot(tile, spec_freqs, psd_left_plot, 'r-', 'LineWidth', 1, ...
        'DisplayName', sprintf('%s population PSD', left_area));
    plot(tile, spec_freqs, psd_right_plot, 'b-', 'LineWidth', 1, ...
        'DisplayName', sprintf('%s population PSD', right_area));

    psd_left_upper = psd_left_shuffle_mean + std_multiplier * psd_left_shuffle_std;
    psd_right_upper = psd_right_shuffle_mean + std_multiplier * psd_right_shuffle_std;
    psd_left_sig_mask = psd_left_plot > psd_left_upper;
    psd_right_sig_mask = psd_right_plot > psd_right_upper;

    plot_significant_segments(tile, spec_freqs, psd_left_plot, psd_left_sig_mask, ...
        sig_min_run_spec, [0.5, 0, 0], '');
    plot_significant_segments(tile, spec_freqs, psd_right_plot, psd_right_sig_mask, ...
        sig_min_run_spec, [0, 0, 0.5], '');
    hold(tile, 'off');

    title(tile, sprintf('%s signal PSD: %s', pair_label, meta.state));
    xlabel(tile, 'Frequency (Hz)');
    ylabel(tile, 'PSD');
    set_legend_visibility(tile, show_legend, 'northeast');

    %% Column 6: pair spectral coherence.
    column = 6;
    selected_control = 'shifted'; % 'shuffled' or 'shifted'
    if strcmp(selected_control, 'shuffled')
        coh_control_mean = coh_shuffle_mean;
        coh_control_std = coh_shuffle_std;
        control_label = 'Shuffled';
    elseif strcmp(selected_control, 'shifted')
        coh_control_mean = shifted_coh_mean;
        coh_control_std = shifted_coh_std;
        control_label = 'Shifted';
    else
        error('Invalid control type selected. Use "shuffled" or "shifted".');
    end

    tile = nexttile(column + n_column*(i-1));
    col_axes(i, column) = tile;
    fill_control_band(tile, spec_freqs, coh_control_mean, coh_control_std, ...
        std_multiplier, [0.8, 0.8, 0.8], sprintf('%s mean ± 2SD', control_label));
    hold(tile, 'on');
    plot(tile, spec_freqs, coherence_plot, 'm-', 'LineWidth', 1, ...
        'DisplayName', 'Coherence');

    coh_upper = coh_control_mean + std_multiplier * coh_control_std;
    coh_sig_mask = coherence_plot > coh_upper;
    plot_significant_segments(tile, spec_freqs, coherence_plot, coh_sig_mask, ...
        sig_min_run_spec, 'k', 'Significant');
    hold(tile, 'off');

    title(tile, sprintf('%s spectral coherence: %s', pair_label, meta.state));
    xlabel(tile, 'Frequency (Hz)');
    ylabel(tile, 'Coherence');
    set_legend_visibility(tile, show_legend, 'northeast');
end

% Harmonize y-limits within each analysis column across rows.
set_common_ylim_by_column(col_axes, 2:n_column);

%% Export to pdf and preview image
fig = gcf;

save_folder = fullfile(root, 'Figures', 'Paper');
check_path(save_folder);

figWidth  = 4*n_column;
figHeight = 4*n_state;
resolution = 300;

set(fig, 'Units', 'inches');
fig.Position(3:4) = [figWidth, figHeight];

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [figWidth, figHeight]);
set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);
set(fig, 'Color', 'w');

output_stub = sprintf('Figure2_population_full_preonly_%s_%s_session%s_thal_acc_vlpfc', ...
    sanitize_filename(meta.animal_name), sanitize_filename(meta.injection), ...
    sanitize_filename(char(string(meta.session_idx))));

sgtitle(sprintf('Population Full-area | %s %s session %s | Thalamus-ACC / Thalamus-VLPFC | %s trials', ...
    meta.animal_name, meta.injection, char(string(meta.session_idx)), analysis_trial_mode), ...
    'Interpreter', 'none');

preview_filename = fullfile(save_folder, [output_stub, '_preview.jpg']);
exportgraphics(fig, preview_filename, ...
    'ContentType', 'image', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);

pdf_filename = fullfile(save_folder, [output_stub, '.pdf']);
exportgraphics(fig, pdf_filename, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);

close(fig);

%% functions


function [acc_trials, vlpfc_trials] = get_population_trial_segments(rasters, filter_acc, filter_vlpfc, trial_mode, trial_idx, t_range, corr_range)
    if strcmpi(trial_mode, 'single')
        trial_indices = trial_idx;
    elseif strcmpi(trial_mode, 'all')
        trial_indices = 1:numel(rasters);
    else
        error('analysis_trial_mode must be ''single'' or ''all''.');
    end

    acc_trials = {};
    vlpfc_trials = {};

    for k = 1:numel(trial_indices)
        this_trial = trial_indices(k);
        if this_trial < 1 || this_trial > numel(rasters)
            error('analysis_trial_idx contains an out-of-range trial.');
        end

        raster = rasters{this_trial};
        idx = resolve_t_range(size(raster, 2), t_range);

        if numel(idx) <= corr_range
            warning('Skipping trial %d because selected length (%d) <= corr_range (%d).', ...
                this_trial, numel(idx), corr_range);
            continue;
        end

        acc_trials{end+1} = mean(raster(filter_acc, idx), 1); %#ok<AGROW>
        vlpfc_trials{end+1} = mean(raster(filter_vlpfc, idx), 1); %#ok<AGROW>
    end

    if isempty(acc_trials)
        error('No valid analysis trials after applying analysis_t_range and corr_range.');
    end
end

function [display_time_ms, pop_acc_display_plot, pop_vlpfc_display_plot] = get_population_display_trace(rasters, filter_acc, filter_vlpfc, trial_idx, t_range, sample_rate, smooth_kernel)
    if trial_idx < 1 || trial_idx > numel(rasters)
        error('display_trial_idx is out of range.');
    end

    raster = rasters{trial_idx};
    idx = resolve_t_range(size(raster, 2), t_range);

    pop_acc_display = mean(raster(filter_acc, idx), 1) * sample_rate;
    pop_vlpfc_display = mean(raster(filter_vlpfc, idx), 1) * sample_rate;

    pop_acc_display_plot = same_conv(pop_acc_display, smooth_kernel);
    pop_vlpfc_display_plot = same_conv(pop_vlpfc_display, smooth_kernel);
    display_time_ms = idx - idx(1);
end

function idx = resolve_t_range(T, t_range)
    if isempty(t_range)
        idx = 1:T;
    else
        idx = t_range(t_range >= 1 & t_range <= T);
    end

    if isempty(idx)
        error('Selected time range is empty after clipping to trial length.');
    end
end

function fr = weighted_signal_mean(signal_trials, sample_rate)
    total_value = 0;
    total_samples = 0;
    for k = 1:numel(signal_trials)
        x = signal_trials{k};
        total_value = total_value + sum(x);
        total_samples = total_samples + numel(x);
    end
    fr = total_value / total_samples * sample_rate;
end

function [corr_avg, lags] = aggregate_norm_xcorr_trials(r1_trials, r2_trials, max_lag, lag_weight_correction)
    if nargin < 4
        lag_weight_correction = true;
    end

    lags = -max_lag:max_lag;

    if lag_weight_correction
        corr_avg = pooled_norm_xcorr_trials(r1_trials, r2_trials, max_lag);
    else
        corr_avg = duration_weighted_norm_xcorr_trials(r1_trials, r2_trials, max_lag);
    end
end

function corr_avg = pooled_norm_xcorr_trials(r1_trials, r2_trials, max_lag)
    lags = -max_lag:max_lag;
    corr_avg = nan(size(lags));

    for lag_i = 1:numel(lags)
        lag = lags(lag_i);
        x_all = [];
        y_all = [];

        for trial_i = 1:numel(r1_trials)
            r1 = r1_trials{trial_i}(:);
            r2 = r2_trials{trial_i}(:);
            N = min(numel(r1), numel(r2));
            r1 = r1(1:N);
            r2 = r2(1:N);

            if max_lag >= N
                continue;
            end

            if lag >= 0
                x = r1((1+lag):N);
                y = r2(1:(N-lag));
            else
                x = r1(1:(N+lag));
                y = r2((1-lag):N);
            end

            x_all = [x_all; x(:)]; %#ok<AGROW>
            y_all = [y_all; y(:)]; %#ok<AGROW>
        end

        valid = isfinite(x_all) & isfinite(y_all);
        x_all = x_all(valid);
        y_all = y_all(valid);

        if numel(x_all) < 2
            continue;
        end

        x_all = x_all - mean(x_all);
        y_all = y_all - mean(y_all);
        denom = sqrt(sum(x_all.^2) * sum(y_all.^2));

        if denom > 0
            corr_avg(lag_i) = sum(x_all .* y_all) / denom;
        end
    end
end

function corr_avg = duration_weighted_norm_xcorr_trials(r1_trials, r2_trials, max_lag)
    n_trial = numel(r1_trials);
    lags = -max_lag:max_lag;
    corr_mat = nan(n_trial, numel(lags));
    weights = zeros(n_trial, 1);

    for trial_i = 1:n_trial
        r1 = r1_trials{trial_i};
        r2 = r2_trials{trial_i};
        N = min(numel(r1), numel(r2));
        if N <= max_lag
            continue;
        end

        [corr_k, ~] = norm_xcorr(r1, r2, max_lag);
        corr_mat(trial_i, :) = corr_k;
        weights(trial_i) = N;
    end

    corr_avg = weighted_average_rows(corr_mat, weights);
end

function [pxx, pyy, pxy_abs, coh, f] = aggregate_pair_spectra_trials(x_trials, y_trials, sample_rate, freqs, window_sec, overlap_frac)
    n_trial = numel(x_trials);
    pxx_mat = nan(n_trial, numel(freqs));
    pyy_mat = nan(n_trial, numel(freqs));
    sxy_mat = nan(n_trial, numel(freqs));
    weights = zeros(n_trial, 1);

    f = freqs(:).';
    for trial_i = 1:n_trial
        x = x_trials{trial_i};
        y = y_trials{trial_i};
        N = min(numel(x), numel(y));
        if N < 4
            continue;
        end

        [pxx_k, pyy_k, ~, ~, f, sxy_k] = compute_pair_spectra(x, y, sample_rate, freqs, window_sec, overlap_frac);
        pxx_mat(trial_i, :) = pxx_k;
        pyy_mat(trial_i, :) = pyy_k;
        sxy_mat(trial_i, :) = sxy_k;
        weights(trial_i) = N;
    end

    pxx = weighted_average_rows(pxx_mat, weights);
    pyy = weighted_average_rows(pyy_mat, weights);
    sxy = weighted_average_rows(sxy_mat, weights);

    pxy_abs = abs(sxy);
    denom = pxx .* pyy;
    coh = abs(sxy).^2 ./ denom;
    coh(~isfinite(coh) | denom <= 0) = NaN;
    coh = min(max(real(coh), 0), 1);
end

function [corr_avg, auto1_avg, auto2_avg, psd1_avg, psd2_avg, cpsd_avg, coh_avg] = aggregate_control_metrics_trials( ...
    r1_trials, r2_trials, corr_range, lag_weight_correction, sample_rate, freqs, window_sec, overlap_frac, control_type)

    n_trial = numel(r1_trials);
    c1_trials = cell(1, n_trial);
    c2_trials = cell(1, n_trial);
    a1_ref_trials = cell(1, n_trial);
    a2_ref_trials = cell(1, n_trial);

    for k = 1:n_trial
        r1 = r1_trials{k};
        r2 = r2_trials{k};

        if strcmpi(control_type, 'shuffled')
            c1 = r1(randperm(numel(r1)));
            c2 = r2(randperm(numel(r2)));
        elseif strcmpi(control_type, 'shifted')
            [shift1, shift2] = random_distinct_circular_shifts(min(numel(r1), numel(r2)), corr_range);
            c1 = circshift(r1, shift1);
            c2 = circshift(r2, shift2);
        else
            error('Unknown control_type: %s', control_type);
        end

        c1_trials{k} = c1;
        c2_trials{k} = c2;
        a1_ref_trials{k} = r1;
        a2_ref_trials{k} = r2;
    end

    [corr_avg, ~] = aggregate_norm_xcorr_trials(c1_trials, c2_trials, corr_range, lag_weight_correction);
    [auto1_avg, ~] = aggregate_norm_xcorr_trials(a1_ref_trials, c1_trials, corr_range, lag_weight_correction);
    [auto2_avg, ~] = aggregate_norm_xcorr_trials(a2_ref_trials, c2_trials, corr_range, lag_weight_correction);
    [psd1_avg, psd2_avg, cpsd_avg, coh_avg, ~] = aggregate_pair_spectra_trials(c1_trials, c2_trials, sample_rate, freqs, window_sec, overlap_frac);
end

function avg = weighted_average_rows(values, weights)
    if isvector(weights)
        weights = weights(:);
        weights = repmat(weights, 1, size(values, 2));
    end

    valid = isfinite(values) & isfinite(weights) & weights > 0;
    weighted_values = values;
    weighted_values(~valid) = 0;
    weight_values = weights;
    weight_values(~valid) = 0;

    denom = sum(weight_values, 1);
    avg = sum(weighted_values .* weight_values, 1) ./ denom;
    avg(denom == 0) = NaN;
end

function set_common_ylim_by_column(ax_mat, columns)
    for c = columns
        axes_this_col = ax_mat(:, c);
        y_min = inf;
        y_max = -inf;

        for r = 1:numel(axes_this_col)
            ax = axes_this_col(r);
            if ~isgraphics(ax)
                continue;
            end
            lim = ylim(ax);
            if all(isfinite(lim))
                y_min = min(y_min, lim(1));
                y_max = max(y_max, lim(2));
            end
        end

        if isfinite(y_min) && isfinite(y_max) && y_max > y_min
            pad = 0.05 * (y_max - y_min);
            if y_min >= 0
                y_min_new = 0;
            else
                y_min_new = y_min - pad;
            end
            y_max_new = y_max + pad;

            for r = 1:numel(axes_this_col)
                ax = axes_this_col(r);
                if isgraphics(ax)
                    ylim(ax, [y_min_new, y_max_new]);
                end
            end
        end
    end
end



function safe_name = sanitize_filename(name)
    safe_name = regexprep(char(string(name)), '[^A-Za-z0-9_-]', '_');
end


function set_legend_visibility(ax, show_legend, location)
    if nargin < 3 || isempty(location)
        location = 'best';
    end

    if show_legend
        legend(ax, 'show', 'Location', location);
    else
        legend(ax, 'off');
    end
end

function filter_area = get_area_filter(cell_area, area_name)
    filter_area = ismember(cell_area, {char(string(area_name))});
end

function check_pair_area_filters(filter_1, filter_2, area_name_1, area_name_2, meta)
    if ~any(filter_1)
        error('No %s neurons found for %s %s %s.', char(string(area_name_1)), meta.prepost, meta.state, meta.area);
    end
    if ~any(filter_2)
        error('No %s neurons found for %s %s %s.', char(string(area_name_2)), meta.prepost, meta.state, meta.area);
    end
end

function t_idx = safe_time_range(t_range, raster_len, meta, range_name)
    if nargin < 4
        range_name = 't_range';
    end
    t_idx = t_range(t_range >= 1 & t_range <= raster_len);
    if isempty(t_idx)
        error('Requested %s is outside raster length for %s %s %s.', ...
            range_name, meta.prepost, meta.state, meta.area);
    end
    if numel(t_idx) < numel(t_range)
        warning('Truncated %s for %s %s %s to fit raster length %d.', ...
            range_name, meta.prepost, meta.state, meta.area, raster_len);
    end
end

function [correlogram, lags] = norm_xcorr(r1, r2, max_lag)
%NORM_XCORR Normalized cross-correlation between two signals.
%
%   [correlogram, lags] = norm_xcorr(r1, r2, max_lag)
%
%   Computes Pearson-normalized cross-correlation for lags:
%       -max_lag : max_lag
%
%   Positive lag means r1 is delayed relative to r2.
%
%   Inputs:
%       r1      - first signal, vector
%       r2      - second signal, vector
%       max_lag - maximum lag in samples
%
%   Outputs:
%       correlogram - normalized cross-correlation values
%       lags        - lag values in samples

    r1 = r1(:);
    r2 = r2(:);

    N = min(length(r1), length(r2));
    r1 = r1(1:N);
    r2 = r2(1:N);

    if max_lag >= N
        error('max_lag must be smaller than the signal length.');
    end

    lags = -max_lag:max_lag;
    correlogram = nan(size(lags));

    for i = 1:length(lags)
        lag = lags(i);

        if lag >= 0
            x = r1((1+lag):N);
            y = r2(1:(N-lag));
        else
            x = r1(1:(N+lag));
            y = r2((1-lag):N);
        end

        % Remove mean
        x = x - mean(x);
        y = y - mean(y);

        % Normalize by energy
        denom = sqrt(sum(x.^2) * sum(y.^2));

        if denom > 0
            correlogram(i) = sum(x .* y) / denom;
        else
            correlogram(i) = NaN;
        end
    end
end

function [correlogram, lags] = norm_xcorr2(r1, r2, max_lag)
    % Compute normalized cross-correlogram, not only by the number of overlapping points, but also by the overall firing rates of the two neurons.
    % Compute raw cross-correlogram
    [raw_correlogram, lags] = xcorr(r1, r2, max_lag, 'unbiased');
    ave_rate = zeros(size(lags));
    for i = 1:length(lags)
        lag = lags(i);
        if lag < 0
            ave_rate(i) = mean(r1(1:end+lag)) * mean(r2(1-lag:end));
        else
            ave_rate(i) = mean(r1(1+lag:end)) * mean(r2(1:end-lag));
        end
    end
    valid_filter = ave_rate > 0;
    correlogram = zeros(size(raw_correlogram));
    correlogram(valid_filter) = raw_correlogram(valid_filter) ./ ave_rate(valid_filter);

    if any(~valid_filter)
        correlogram(~valid_filter) = 0;
        warning('Some lags have zero average firing rate product, resulting in 0 in the normalized correlogram.');
    end
end

function [p_low, p_high] = wilsonCI(M, N, alpha)
    % Wilson score interval
    % M = successes, N = trials
    % alpha = significance level (e.g., 0.05 for 95% CI)

    if nargin < 3
        alpha = 0.05;
    end

    p = M / N;
    z = norminv(1 - alpha/2); % Z for CI

    denominator = 1 + (z^2)/N;
    center = p + (z^2)/(2*N);
    radius = z * sqrt( (p*(1-p)/N) + (z^2)/(4*N^2) );

    p_low = (center - radius) / denominator;
    p_high = (center + radius) / denominator;
end

function p_val = twoProportionPValue(success1, trials1, success2, trials2)
    % Two-proportion z-test (two-tailed) between pre and post proportions

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

function convolved = same_conv(signal, kernel)
    % Convolve signal and kernel, return the central part. Normalize the border part by only including the overlapping part of the kernel.
    convolved = conv(signal, kernel, 'same');
    weight = conv(ones(size(signal)), kernel, 'same');
    convolved = convolved ./ weight;
end


function signal = remove_dc(signal)
    % Remove the finite-sample mean before FFT and replace non-finite values with 0.
    signal = signal(:);
    finite_filter = isfinite(signal);
    if any(finite_filter)
        signal(finite_filter) = signal(finite_filter) - mean(signal(finite_filter));
    end
    signal(~finite_filter) = 0;
end

function spec_plot = smooth_spectrum_for_plot(spec, smooth_window)
    % Smooth only the frequency-domain amplitude for visualization.
    if smooth_window > 1
        spec_plot = movmean(spec, smooth_window);
    else
        spec_plot = spec;
    end
end

function [pxx, pyy, pxy_abs, coh, f, sxy] = compute_pair_spectra(x, y, sample_rate, freqs, window_sec, overlap_frac)
    % Estimate Sxx, Syy, Sxy, and coherence using one consistent Welch
    % cross-spectral framework. This requires Signal Processing Toolbox for cpsd.
    % pxx: power spectral density of x (Sxx)
    % pyy: power spectral density of y (Syy)
    % pxy_abs: absolute value of cross-spectrum (|Sxy|)
    % coh: coherence = |Sxy|^2 / (Sxx*Syy)

    x = remove_dc(x(:));
    y = remove_dc(y(:));
    freqs = freqs(:);

    N = min(length(x), length(y));
    x = x(1:N);
    y = y(1:N);

    if N < 4
        f = freqs.';
        pxx = nan(size(f));
        pyy = nan(size(f));
        pxy_abs = nan(size(f));
        coh = nan(size(f));
        sxy = nan(size(f));
        return;
    end

    [win, noverlap] = make_welch_window(N, sample_rate, window_sec, overlap_frac);

    try
        [sxx_col, f_col] = cpsd(x, x, win, noverlap, freqs, sample_rate);
        [syy_col, ~] = cpsd(y, y, win, noverlap, freqs, sample_rate);
        [sxy_col, ~] = cpsd(x, y, win, noverlap, freqs, sample_rate);

        pxx = real(sxx_col(:).');
        pyy = real(syy_col(:).');
        sxy = sxy_col(:).';

        pxy_abs = abs(sxy);
        denom = pxx .* pyy;
        coh = abs(sxy).^2 ./ denom;
        coh(~isfinite(coh) | denom <= 0) = NaN;
        coh = min(max(real(coh), 0), 1);

        f = f_col(:).';
    catch ME
        warning('cpsd failed (%s). Falling back to direct power-like spectra; coherence is set to NaN.', ME.message);
        amp_x = myFT(x, sample_rate, freqs);
        amp_y = myFT(y, sample_rate, freqs);
        pxx = amp_x(:).'.^2;
        pyy = amp_y(:).'.^2;
        pxy_abs = nan(size(pxx));
        coh = nan(size(pxx));
        sxy = nan(size(pxx));
        f = freqs(:).';
    end
end

function [pxx, f] = compute_signal_psd(signal, sample_rate, freqs, window_sec, overlap_frac)
    % Estimate a single-signal PSD using the same cpsd(x,x) framework used
    % for the pair spectra.
    signal = remove_dc(signal(:));
    freqs = freqs(:);
    N = length(signal);

    if N < 4
        f = freqs.';
        pxx = nan(size(f));
        return;
    end

    [win, noverlap] = make_welch_window(N, sample_rate, window_sec, overlap_frac);

    try
        [sxx_col, f_col] = cpsd(signal, signal, win, noverlap, freqs, sample_rate);
        pxx = real(sxx_col(:).');
        f = f_col(:).';
    catch ME
        warning('cpsd failed for signal PSD (%s). Falling back to direct power-like spectrum.', ME.message);
        amp = myFT(signal, sample_rate, freqs);
        pxx = amp(:).'.^2;
        f = freqs(:).';
    end
end

function [win, noverlap] = make_welch_window(N, sample_rate, window_sec, overlap_frac)
    % Hann window used by all Welch spectral estimates.
    win_len = min(max(round(window_sec * sample_rate), 4), N);
    if win_len <= 1
        win = 1;
    else
        win = 0.5 - 0.5*cos(2*pi*(0:win_len-1).'/(win_len-1));
    end
    noverlap = min(floor(overlap_frac * win_len), win_len - 1);
end

function fill_control_band(tile, x, control_mean, control_std, multiplier, color, display_name)
    % Plot mean ± multiplier*SD as a shaded control band.
    x = x(:).';
    control_mean = control_mean(:).';
    control_std = control_std(:).';

    upper = control_mean + multiplier * control_std;
    lower = control_mean - multiplier * control_std;

    valid = isfinite(x) & isfinite(upper) & isfinite(lower);
    if ~any(valid)
        return;
    end

    fill(tile, [x(valid), fliplr(x(valid))], ...
        [upper(valid), fliplr(lower(valid))], ...
        color, 'EdgeColor', 'none', 'FaceAlpha', 0.35, ...
        'DisplayName', display_name);
end

function plot_significant_segments(tile, x, y, sig_mask, min_run_bins, color, display_name)
    % Overlay only significant runs that last for at least min_run_bins.
    x = x(:).';
    y = y(:).';
    sig_mask = sig_mask(:).' & isfinite(x) & isfinite(y);

    runs = find_logical_runs(sig_mask);
    has_label = false;

    for r = 1:size(runs, 1)
        idx = runs(r, 1):runs(r, 2);
        if numel(idx) < min_run_bins
            continue;
        end

        if ~has_label
            plot(tile, x(idx), y(idx), '-', 'Color', color, ...
                'LineWidth', 3, 'DisplayName', display_name);
            has_label = true;
        else
            plot(tile, x(idx), y(idx), '-', 'Color', color, ...
                'LineWidth', 3, 'HandleVisibility', 'off');
        end
    end
end

function runs = find_logical_runs(mask)
    % Return [start_idx, end_idx] for each contiguous true run in a logical vector.
    mask = logical(mask(:).');
    edge = diff([false, mask, false]);
    starts = find(edge == 1);
    ends = find(edge == -1) - 1;
    runs = [starts(:), ends(:)];
end

function m = mean_omitnan(A)
    % Column-wise mean, ignoring NaNs/Infs, without requiring toolbox support.
    m = nan(1, size(A, 2));
    for k = 1:size(A, 2)
        vals = A(:, k);
        vals = vals(isfinite(vals));
        if ~isempty(vals)
            m(k) = mean(vals);
        end
    end
end

function s = std_omitnan(A)
    % Column-wise sample standard deviation, ignoring NaNs/Infs.
    s = nan(1, size(A, 2));
    for k = 1:size(A, 2)
        vals = A(:, k);
        vals = vals(isfinite(vals));
        if numel(vals) >= 2
            s(k) = std(vals);
        elseif numel(vals) == 1
            s(k) = 0;
        end
    end
end


function shift_val = random_circular_shift(signal_len, min_shift)
    % Draw a circular shift that avoids near-zero shifts when possible.
    if signal_len <= 1
        shift_val = 0;
        return;
    end

    min_shift = min(max(round(min_shift), 1), max(signal_len - 1, 1));
    max_shift = max(signal_len - min_shift, 1);

    if max_shift >= min_shift
        shift_val = randi([min_shift, max_shift]);
    else
        shift_val = randi(signal_len - 1);
    end
end

function [shift1, shift2] = random_distinct_circular_shifts(signal_len, min_relative_shift)
    % Draw two non-identical circular shifts with a relative shift that is
    % not close to zero when possible.
    if signal_len <= 1
        shift1 = 0;
        shift2 = 0;
        return;
    end

    min_relative_shift = min(max(round(min_relative_shift), 1), floor(signal_len / 2));
    max_tries = 1000;

    for attempt = 1:max_tries %#ok<NASGU>
        shift1 = random_circular_shift(signal_len, min_relative_shift);
        shift2 = random_circular_shift(signal_len, min_relative_shift);

        rel_shift = mod(shift2 - shift1, signal_len);
        circular_distance = min(rel_shift, signal_len - rel_shift);

        if shift1 ~= shift2 && circular_distance >= min_relative_shift
            return;
        end
    end

    shift1 = random_circular_shift(signal_len, min_relative_shift);
    shift2 = mod(shift1 + min_relative_shift - 1, signal_len) + 1;
end


function amp = myFT(signal, sample_rate, freq)
%MYFT Real amplitude spectrum at specified frequencies.
%
%   amp = myFT(signal, sample_rate, freq)
%
%   Inputs:
%       signal      - signal vector
%       sample_rate - sampling rate in Hz
%       freq        - frequency or vector of frequencies in Hz
%
%   Output:
%       amp         - real amplitude spectrum at requested frequencies

    signal = signal(:);
    freq_shape = size(freq);
    freq = freq(:);

    N = length(signal);
    n = (0:N-1).';

    if any(abs(freq) > sample_rate / 2)
        warning('Some frequencies exceed the Nyquist frequency.');
    end

    basis = exp(-1i * 2*pi/sample_rate * (n * freq.'));
    spec = signal.' * basis;

    amp = 2 * abs(spec) / N;

    % Do not double DC component
    dc_idx = (freq == 0);
    amp(dc_idx) = abs(spec(dc_idx)) / N;

    amp = reshape(amp, freq_shape);
end