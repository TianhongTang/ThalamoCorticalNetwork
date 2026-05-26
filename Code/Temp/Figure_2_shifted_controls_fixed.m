%% Figure 2: Example session population-level ACC-VLPFC analyses.
% Positive lag in correlogram: signal1(ACC population) is later than signal2(VLPFC population).
% Same layout and spectral-control logic as Figure 1, replacing the example
% neuron pair with population-average firing-rate traces.

clear;
%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end

% include code folder and utils
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));
% addpath(fullfile(root, 'Code', 'Utils', 'HELPER_GENERAL'));

%% Main

% Example session: Slayer Mus 6.
meta = struct();
meta.animal_name = 'Slayer';
meta.injection = 'Muscimol';
meta.align = 'Longest';
meta.session_idx = 6;
meta.resting_dur_threshold = 15;

areas    = {'Cortex',   'Cortex',    'Cortex',   'Cortex'};
preposts = {'Pre',      'Pre',       'Post',     'Post'};
states   = {'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};
n_state = numel(areas);

% figure: same row structure as Figure 1.
% Row 1 = ACC/VLPFC population firing-rate traces.
% Row 2 = population cross-correlogram.
% Row 3 = population auto-correlograms.
% Row 4 = population cross-PSD.
% Row 5 = ACC/VLPFC population PSDs.
% Row 6 = ACC-VLPFC population spectral coherence.
f = figure('Color', 'w', 'Visible', 'off');
tiles = tiledlayout(6, n_state, "TileSpacing", "Compact", "Padding", "Compact");

%% parameters
shuffle_N = 10;
std_multiplier = 2; % threshold for shuffled controls, in multiples of shuffled SD.
sig_min_run_bins = 5; % only mark significance if it lasts for at least this many consecutive bins.
t_range = 1:60000;
display_t_range = 1:2000; % time range for row 1 firing-rate traces.
corr_range = 5000; % ms.
smooth_window = 25; % ms.
sample_rate = 1000; % Hz.
freqs = linspace(0, 150, 301); % Hz.
spec_smooth_window = 15; % frequency-bin smoothing for spectral visualization only.
psd_window_sec = 10; % Welch PSD window length for row 5 signal spectra.
psd_overlap_frac = 0.5; % Welch PSD fractional overlap.

% Same triangular smoothing kernel as Figure 1.
smooth_kernel = 1 - abs(-smooth_window:smooth_window) / smooth_window;
smooth_kernel = smooth_kernel / sum(smooth_kernel); % normalize kernel

%% Rows 1-6: population firing rate, correlation, and spectra
for i = 1:n_state
    % Load raster data.
    meta.area = areas{i};
    meta.prepost = preposts{i};
    meta.state = states{i};
    meta.filename = generate_filename('raster', meta);

    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    raster = raster_data.data.rasters{1};
    t_idx = safe_time_range(t_range, size(raster, 2), meta, 't_range');
    display_idx = safe_time_range(display_t_range, size(raster, 2), meta, 'display_t_range');

    cell_area = raster_data.data.cell_area;
    [filter_acc, filter_vlpfc] = get_area_filters(cell_area);
    check_area_filters(filter_acc, filter_vlpfc, meta);

    r_acc = raster(filter_acc, t_idx);
    r_vlpfc = raster(filter_vlpfc, t_idx);

    N_acc = sum(filter_acc);
    N_vlpfc = sum(filter_vlpfc);

    % Population-average activity: mean spikes/bin across neurons.
    % Multiplying by sample_rate converts it to average firing rate in Hz/neuron.
    pop_acc = mean(r_acc, 1);
    pop_vlpfc = mean(r_vlpfc, 1);

    % Display-only population firing-rate traces for row 1.
    pop_acc_display = mean(raster(filter_acc, display_idx), 1) * sample_rate;
    pop_vlpfc_display = mean(raster(filter_vlpfc, display_idx), 1) * sample_rate;
    pop_acc_display_plot = same_conv(pop_acc_display, smooth_kernel);
    pop_vlpfc_display_plot = same_conv(pop_vlpfc_display, smooth_kernel);
    display_time_ms = display_idx - display_idx(1);

    fprintf('Population firing rates for %s %s %s:\n', meta.prepost, meta.state, meta.area);
    fprintf('ACC population mean: %.3f Hz/neuron, N=%d\n', mean(pop_acc) * sample_rate, N_acc);
    fprintf('VLPFC population mean: %.3f Hz/neuron, N=%d\n', mean(pop_vlpfc) * sample_rate, N_vlpfc);

    %% Row 1: population firing-rate traces.
    row = 1;
    tile = nexttile(i + n_state*(row-1));
    plot(tile, display_time_ms, pop_acc_display_plot, 'r-', 'LineWidth', 1, ...
        'DisplayName', sprintf('ACC, N=%d', N_acc));
    hold(tile, 'on');
    plot(tile, display_time_ms, pop_vlpfc_display_plot, 'b-', 'LineWidth', 1, ...
        'DisplayName', sprintf('VLPFC, N=%d', N_vlpfc));
    hold(tile, 'off');
    title(tile, sprintf('Population firing rate: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Time (ms)');
    ylabel(tile, 'Firing rate (Hz/neuron)');
    legend(tile, 'show', 'Location', 'northeast');

    %% Compute population correlations.
    [correlogram, lags] = norm_xcorr(pop_acc, pop_vlpfc, corr_range);
    [auto_acc, ~] = norm_xcorr(pop_acc, pop_acc, corr_range);
    [auto_vlpfc, ~] = norm_xcorr(pop_vlpfc, pop_vlpfc, corr_range);

    smooth_correlogram = same_conv(correlogram, smooth_kernel);
    smooth_auto_acc = same_conv(auto_acc, smooth_kernel);
    smooth_auto_vlpfc = same_conv(auto_vlpfc, smooth_kernel);

    % Row 4/5 spectra are computed from the original population signals with
    % the same Welch/cross-spectral framework as Figure 1.
    [psd_acc, psd_vlpfc, cpsd, coherence, spec_freqs] = compute_pair_spectra( ...
        pop_acc, pop_vlpfc, sample_rate, freqs, psd_window_sec, psd_overlap_frac);

    coherence_plot = smooth_spectrum_for_plot(coherence, spec_smooth_window);
    cpsd_plot = smooth_spectrum_for_plot(cpsd, spec_smooth_window);
    psd_acc_plot = smooth_spectrum_for_plot(psd_acc, spec_smooth_window);
    psd_vlpfc_plot = smooth_spectrum_for_plot(psd_vlpfc, spec_smooth_window);

    %% Compute both shuffled and shifted controls.
    % - Shuffled controls use random time permutation and are used for PSD.
    % - Shifted controls use independent circular shifts, preserving each signal's
    %   temporal structure, and are used for Xcorr, autocorr, cross-PSD, and coherence.
    shuffle_correlograms = zeros(shuffle_N, length(correlogram));
    shuffle_auto_acc_controls = zeros(shuffle_N, length(auto_acc));
    shuffle_auto_vlpfc_controls = zeros(shuffle_N, length(auto_vlpfc));
    shuffle_cpsd = zeros(shuffle_N, numel(spec_freqs));
    shuffle_coh = zeros(shuffle_N, numel(spec_freqs));
    shuffle_psd_acc = zeros(shuffle_N, numel(spec_freqs));
    shuffle_psd_vlpfc = zeros(shuffle_N, numel(spec_freqs));

    shifted_correlograms = zeros(shuffle_N, length(correlogram));
    shifted_auto_acc_controls = zeros(shuffle_N, length(auto_acc));
    shifted_auto_vlpfc_controls = zeros(shuffle_N, length(auto_vlpfc));
    shifted_cpsd = zeros(shuffle_N, numel(spec_freqs));
    shifted_coh = zeros(shuffle_N, numel(spec_freqs));

    for j = 1:shuffle_N
        fprintf('%d/%d population controls started\n', j, shuffle_N);

        % Shuffled controls: destroy temporal structure.
        shuffled_acc = pop_acc(randperm(length(pop_acc)));
        shuffled_vlpfc = pop_vlpfc(randperm(length(pop_vlpfc)));

        [shuffle_corr, ~] = norm_xcorr(shuffled_acc, shuffled_vlpfc, corr_range);
        [shuffle_auto_acc, ~] = norm_xcorr(pop_acc, shuffled_acc, corr_range);
        [shuffle_auto_vlpfc, ~] = norm_xcorr(pop_vlpfc, shuffled_vlpfc, corr_range);

        shuffle_correlograms(j, :) = same_conv(shuffle_corr, smooth_kernel);
        shuffle_auto_acc_controls(j, :) = same_conv(shuffle_auto_acc, smooth_kernel);
        shuffle_auto_vlpfc_controls(j, :) = same_conv(shuffle_auto_vlpfc, smooth_kernel);

        [shuffled_psd_acc, shuffled_psd_vlpfc, shuffled_cpsd, shuffled_coh_j, ~] = compute_pair_spectra( ...
            shuffled_acc, shuffled_vlpfc, sample_rate, freqs, psd_window_sec, psd_overlap_frac);

        shuffle_cpsd(j, :) = smooth_spectrum_for_plot(shuffled_cpsd, spec_smooth_window);
        shuffle_coh(j, :) = smooth_spectrum_for_plot(shuffled_coh_j, spec_smooth_window);
        shuffle_psd_acc(j, :) = smooth_spectrum_for_plot(shuffled_psd_acc, spec_smooth_window);
        shuffle_psd_vlpfc(j, :) = smooth_spectrum_for_plot(shuffled_psd_vlpfc, spec_smooth_window);

        % Shifted controls: preserve each population trace's temporal structure
        % while disrupting the original alignment.
        [shift_acc, shift_vlpfc] = random_distinct_circular_shifts(length(pop_acc), corr_range);
        shifted_acc = circshift(pop_acc, shift_acc);
        shifted_vlpfc = circshift(pop_vlpfc, shift_vlpfc);

        [shifted_corr, ~] = norm_xcorr(shifted_acc, shifted_vlpfc, corr_range);
        [shifted_auto_acc, ~] = norm_xcorr(pop_acc, shifted_acc, corr_range);
        [shifted_auto_vlpfc, ~] = norm_xcorr(pop_vlpfc, shifted_vlpfc, corr_range);

        shifted_correlograms(j, :) = same_conv(shifted_corr, smooth_kernel);
        shifted_auto_acc_controls(j, :) = same_conv(shifted_auto_acc, smooth_kernel);
        shifted_auto_vlpfc_controls(j, :) = same_conv(shifted_auto_vlpfc, smooth_kernel);

        [~, ~, shifted_cpsd_j, shifted_coh_j, ~] = compute_pair_spectra( ...
            shifted_acc, shifted_vlpfc, sample_rate, freqs, psd_window_sec, psd_overlap_frac);
        shifted_cpsd(j, :) = smooth_spectrum_for_plot(shifted_cpsd_j, spec_smooth_window);
        shifted_coh(j, :) = smooth_spectrum_for_plot(shifted_coh_j, spec_smooth_window);

        fprintf('%d/%d population controls finished\n', j, shuffle_N);
    end

    % Shuffled controls.
    shuffle_mean = mean_omitnan(shuffle_correlograms);
    shuffle_std = std_omitnan(shuffle_correlograms);

    shuffle_auto_acc_mean = mean_omitnan(shuffle_auto_acc_controls);
    shuffle_auto_acc_std = std_omitnan(shuffle_auto_acc_controls);

    shuffle_auto_vlpfc_mean = mean_omitnan(shuffle_auto_vlpfc_controls);
    shuffle_auto_vlpfc_std = std_omitnan(shuffle_auto_vlpfc_controls);

    cpsd_shuffle_mean = mean_omitnan(shuffle_cpsd);
    cpsd_shuffle_std = std_omitnan(shuffle_cpsd);

    coh_shuffle_mean = mean_omitnan(shuffle_coh);
    coh_shuffle_std = std_omitnan(shuffle_coh);

    psd_acc_shuffle_mean = mean_omitnan(shuffle_psd_acc);
    psd_acc_shuffle_std = std_omitnan(shuffle_psd_acc);

    psd_vlpfc_shuffle_mean = mean_omitnan(shuffle_psd_vlpfc);
    psd_vlpfc_shuffle_std = std_omitnan(shuffle_psd_vlpfc);

    % Shifted controls used for Xcorr, autocorr, and cross-PSD plotting.
    shifted_ccg_mean = mean_omitnan(shifted_correlograms);
    shifted_ccg_std = std_omitnan(shifted_correlograms);

    shifted_auto_acc_mean = mean_omitnan(shifted_auto_acc_controls);
    shifted_auto_acc_std = std_omitnan(shifted_auto_acc_controls);

    shifted_auto_vlpfc_mean = mean_omitnan(shifted_auto_vlpfc_controls);
    shifted_auto_vlpfc_std = std_omitnan(shifted_auto_vlpfc_controls);

    cpsd_shifted_mean = mean_omitnan(shifted_cpsd);
    cpsd_shifted_std = std_omitnan(shifted_cpsd);

    coh_shifted_mean = mean_omitnan(shifted_coh);
    coh_shifted_std = std_omitnan(shifted_coh);

    %% Row 2: population cross-correlogram with shifted control and significant segments.
    row = 2;
    tile = nexttile(i + n_state*(row-1));
    shifted_ccg_upper = shifted_ccg_mean + std_multiplier * shifted_ccg_std;
    shifted_ccg_lower = shifted_ccg_mean - std_multiplier * shifted_ccg_std;

    fill_control_band(tile, lags, shifted_ccg_mean, shifted_ccg_std, std_multiplier, ...
        [0.8, 0.8, 0.8], 'Shifted mean ± 2SD');
    hold(tile, 'on');
    xline(tile, 0, 'k--', 'HandleVisibility', 'off');
    yline(tile, 0, 'k--', 'HandleVisibility', 'off');
    plot(tile, lags, smooth_correlogram, 'm-', 'LineWidth', 1, ...
        'DisplayName', 'ACC-VLPFC cross-corr');

    ccg_sig_mask = (smooth_correlogram > shifted_ccg_upper) | ...
                   (smooth_correlogram < shifted_ccg_lower);
    plot_significant_segments(tile, lags, smooth_correlogram, ccg_sig_mask, ...
        sig_min_run_bins, 'k', 'Significant');
    hold(tile, 'off');

    title(tile, sprintf('Population correlogram: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Lag (ms)');
    ylabel(tile, 'Normalized correlation');
    legend(tile, 'show', 'Location', 'southeast');
    ylim(tile, [-0.05, 0.05]);

    %% Row 3: population auto-correlograms with shifted self-controls.
    row = 3;
    tile = nexttile(i + n_state*(row-1));

    fill_control_band(tile, lags, shifted_auto_acc_mean, shifted_auto_acc_std, ...
        std_multiplier, [1.0, 0.85, 0.85], 'ACC shifted ± 2SD');
    hold(tile, 'on');
    fill_control_band(tile, lags, shifted_auto_vlpfc_mean, shifted_auto_vlpfc_std, ...
        std_multiplier, [0.85, 0.85, 1.0], 'VLPFC shifted ± 2SD');

    plot(tile, lags, smooth_correlogram, 'm-', 'LineWidth', 2, ...
        'DisplayName', 'Cross-corr', 'Color', [1, 0, 1, 0.2]);
    xline(tile, 0, 'k--', 'HandleVisibility', 'off');
    yline(tile, 0, 'k--', 'HandleVisibility', 'off');
    plot(tile, lags, smooth_auto_acc, 'r-', 'LineWidth', 1, ...
        'DisplayName', sprintf('ACC auto-corr, N=%d', N_acc));
    plot(tile, lags, smooth_auto_vlpfc, 'b-', 'LineWidth', 1, ...
        'DisplayName', sprintf('VLPFC auto-corr, N=%d', N_vlpfc));

    auto_acc_upper = shifted_auto_acc_mean + std_multiplier * shifted_auto_acc_std;
    auto_acc_lower = shifted_auto_acc_mean - std_multiplier * shifted_auto_acc_std;
    auto_vlpfc_upper = shifted_auto_vlpfc_mean + std_multiplier * shifted_auto_vlpfc_std;
    auto_vlpfc_lower = shifted_auto_vlpfc_mean - std_multiplier * shifted_auto_vlpfc_std;

    auto_acc_sig_mask = (smooth_auto_acc > auto_acc_upper) | ...
                        (smooth_auto_acc < auto_acc_lower);
    auto_vlpfc_sig_mask = (smooth_auto_vlpfc > auto_vlpfc_upper) | ...
                          (smooth_auto_vlpfc < auto_vlpfc_lower);

    plot_significant_segments(tile, lags, smooth_auto_acc, auto_acc_sig_mask, ...
        sig_min_run_bins, [0.5, 0, 0], 'ACC significant');
    plot_significant_segments(tile, lags, smooth_auto_vlpfc, auto_vlpfc_sig_mask, ...
        sig_min_run_bins, [0, 0, 0.5], 'VLPFC significant');

    hold(tile, 'off');

    title(tile, sprintf('Population auto-correlograms: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Lag (ms)');
    ylabel(tile, 'Normalized correlation');
    legend(tile, 'show', 'Location', 'southeast');
    ylim(tile, [-0.05, 0.05]);

    %% Row 4: population cross-PSD with shifted control.
    row = 4;

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

    tile = nexttile(i + n_state*(row-1));
    fill_control_band(tile, spec_freqs, cpsd_control_mean, cpsd_control_std, ...
        std_multiplier, [0.8, 0.8, 0.8], sprintf('%s mean ± 2SD', control_label));
    hold(tile, 'on');
    plot(tile, spec_freqs, cpsd_plot, 'm-', 'LineWidth', 1, ...
        'DisplayName', 'Cross-PSD');

    cpsd_upper = cpsd_control_mean + std_multiplier * cpsd_control_std;
    cpsd_sig_mask = cpsd_plot > cpsd_upper;
    plot_significant_segments(tile, spec_freqs, cpsd_plot, cpsd_sig_mask, ...
        sig_min_run_bins, 'k', 'Significant');
    hold(tile, 'off');

    title(tile, sprintf('Population cross-PSD: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Frequency (Hz)');
    ylabel(tile, 'Cross-PSD');
    ylim([0, 6e-7]);
    legend(tile, 'show', 'Location', 'northeast');

    %% Row 5: population PSDs with shuffled controls.
    row = 5;
    tile = nexttile(i + n_state*(row-1));

    fill_control_band(tile, spec_freqs, psd_acc_shuffle_mean, psd_acc_shuffle_std, ...
        std_multiplier, [1.0, 0.85, 0.85], 'ACC shuffled ± 2SD');
    hold(tile, 'on');
    fill_control_band(tile, spec_freqs, psd_vlpfc_shuffle_mean, psd_vlpfc_shuffle_std, ...
        std_multiplier, [0.85, 0.85, 1.0], 'VLPFC shuffled ± 2SD');

    plot(tile, spec_freqs, psd_acc_plot, 'r-', 'LineWidth', 1, ...
        'DisplayName', 'ACC population PSD');
    plot(tile, spec_freqs, psd_vlpfc_plot, 'b-', 'LineWidth', 1, ...
        'DisplayName', 'VLPFC population PSD');

    psd_acc_upper = psd_acc_shuffle_mean + std_multiplier * psd_acc_shuffle_std;
    psd_vlpfc_upper = psd_vlpfc_shuffle_mean + std_multiplier * psd_vlpfc_shuffle_std;
    psd_acc_sig_mask = psd_acc_plot > psd_acc_upper;
    psd_vlpfc_sig_mask = psd_vlpfc_plot > psd_vlpfc_upper;

    plot_significant_segments(tile, spec_freqs, psd_acc_plot, psd_acc_sig_mask, ...
        sig_min_run_bins, [0.5, 0, 0], 'ACC significant');
    plot_significant_segments(tile, spec_freqs, psd_vlpfc_plot, psd_vlpfc_sig_mask, ...
        sig_min_run_bins, [0, 0, 0.5], 'VLPFC significant');
    hold(tile, 'off');

    title(tile, sprintf('Population signal PSD: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Frequency (Hz)');
    ylabel(tile, 'PSD');
    ylim(tile, [0, 2e-6]);
    legend(tile, 'show', 'Location', 'northeast');

    %% Row 6: population spectral coherence with shuffled control.
    % Coherence is the normalized cross-spectrum: |Sxy|^2 / (Sxx*Syy).
    % Global circular shifts only rotate cross-spectral phase and do not
    % reliably change coherence magnitude, so use shuffled control here.
    row = 6;
    tile = nexttile(i + n_state*(row-1));
    fill_control_band(tile, spec_freqs, coh_shuffle_mean, coh_shuffle_std, ...
        std_multiplier, [0.8, 0.8, 0.8], 'Shuffled mean ± 2SD');
    hold(tile, 'on');
    plot(tile, spec_freqs, coherence_plot, 'm-', 'LineWidth', 1, ...
        'DisplayName', 'Coherence');

    coh_upper = coh_shuffle_mean + std_multiplier * coh_shuffle_std;
    coh_sig_mask = coherence_plot > coh_upper;
    plot_significant_segments(tile, spec_freqs, coherence_plot, coh_sig_mask, ...
        sig_min_run_bins, 'k', 'Significant');
    hold(tile, 'off');

    title(tile, sprintf('Population spectral coherence: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Frequency (Hz)');
    ylabel(tile, 'Coherence');
    ylim(tile, [0, 1]);
    legend(tile, 'show', 'Location', 'northeast');
end

%% Export to pdf and preview image
fig = gcf;

save_folder = fullfile(root, 'Figures', 'Paper');
check_path(save_folder);

figWidth  = 16.0;   % inches.
figHeight = 24.0;   % inches.
resolution = 300;   % dpi; mainly affects rasterized components.

set(fig, 'Units', 'inches');
fig.Position(3:4) = [figWidth, figHeight];

set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [figWidth, figHeight]);
set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);
set(fig, 'Color', 'w');

pdf_filename = fullfile(save_folder, 'Figure2.pdf');
exportgraphics(fig, pdf_filename, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);

preview_filename = fullfile(save_folder, 'Figure2_preview.jpg');
exportgraphics(fig, preview_filename, ...
    'ContentType', 'image', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);

close(fig);

%% functions
function [filter_acc, filter_vlpfc] = get_area_filters(cell_area)
    filter_acc = ismember(cell_area, {'ACC'});
    filter_vlpfc = ismember(cell_area, {'VLPFC'});
end

function check_area_filters(filter_acc, filter_vlpfc, meta)
    if ~any(filter_acc)
        error('No ACC neurons found for %s %s %s.', meta.prepost, meta.state, meta.area);
    end
    if ~any(filter_vlpfc)
        error('No VLPFC neurons found for %s %s %s.', meta.prepost, meta.state, meta.area);
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

function [pxx, pyy, pxy_abs, coh, f] = compute_pair_spectra(x, y, sample_rate, freqs, window_sec, overlap_frac)
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