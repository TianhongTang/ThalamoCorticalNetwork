%% Figure 1: An example neuron pair showing state-dependent changes in firing and correlation.
% Positive lag in correlogram: signal1(ACC) is later than signal2(VLPFC). 

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

% areas    = {'Full',     'Full',      'Cortex',   'Cortex'};
areas    = {'Cortex',   'Cortex',    'Cortex',   'Cortex'};
preposts = {'Pre',      'Pre',       'Post',     'Post'};
states   = {'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};
% areas    = {'Cortex',   'Cortex',  };
% preposts = {'Pre',      'Pre',     };
% states   = {'RestOpen', 'RestClose'};
n_state = numel(areas);

% figure
f = figure();
tiles = tiledlayout(5, n_state, "TileSpacing", "Compact", "Padding", "Compact");

%% parameters
shuffle_N = 10;
std_multiplier = 2; % threshold for shuffled controls, in multiples of shuffled SD.
sig_min_run_bins = 5; % only mark significance if it lasts for at least this many consecutive bins.
err_multi = 2; % threshold for significant J, in multiples of the error estimate from GLM.
t_range = 1:60000;
corr_range = 500; % ms
smooth_window = 15; % ms
sample_rate = 1000; % Hz
freqs = linspace(0, 150, 301); % Hz
spec_smooth_window = 15; % frequency-bin smoothing for spectral visualization only
psd_window_sec = 10; % Welch PSD window length for row 5 signal spectra
psd_overlap_frac = 0.5; % Welch PSD fractional overlap
mygauss = @(size, sigma) exp(-(-floor(size/2):floor(size/2)).^2/(2*sigma^2));
% smooth_kernel = mygauss(10*smooth_window+1, smooth_window); % gaussian kernel
% smooth_kernel = ones(1, smooth_window); % moving average kernel
smooth_kernel = 1 - abs(-smooth_window:smooth_window) / smooth_window;
smooth_kernel = smooth_kernel / sum(smooth_kernel); % normalize kernel

%% Row 1: Raster of selected neurons

% Example session: Slayer Mus 6. Selected t_range. Justification: Best sleep period for Pre-eyeclose.
% Remove this range in population analysis, or add proper filters to select sleep periods.

% selected_neurons = 1:95;
% selected_neurons = [16:21, 39, 44, 61:67, 83:93];
% selected_neurons = [14:19, 33:40, 42, 44, 46, 49, 58:73, 83:95];
% selected_neurons = [16, 61];
selected_neurons = [2, 42];
% selected_neurons = [29, 81];
display_t_range = 00001:2000;

for i = 1:length(areas)
% for i = 2:2
    % Load data
    meta.area = areas{i};
    meta.prepost = preposts{i};
    meta.state = states{i};
    meta.filename = generate_filename('raster', meta);

    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    % Plot
    % N = raster_data.meta.N;

    cell_area = raster_data.data.cell_area;
    % tile = nexttile(i*2-1);
    tile = nexttile(i);
    raster = raster_data.data.rasters{1}(:, display_t_range);

    cell_area = cell_area(selected_neurons);
    raster = raster(selected_neurons, :);

    N = numel(selected_neurons);
    colors = zeros(N, 3);
    for j = 1:N
        switch cell_area{j}
            case 'Thalamus'
                colors(j, :) = [1, 0, 1];
            case 'ACC'
                colors(j, :) = [0, 0, 1]; % blue
            case 'VLPFC'
                colors(j, :) = [1, 0, 0]; % red
            otherwise
                colors(j, :) = [0, 0, 0]; % black
        end
    end
    cla(tile);
    raster_visualization_plot(tile, raster, colors)
    title(sprintf('Example pair, %s, %s', meta.prepost, meta.state));
    xlabel('Time (ms)');
    % ylabel('Neuron No.');
    % ylabel('');
    % yticks([1, 2]);
    % yticklabels({'ACC #2', 'VLPFC #1'});
    % ytickangle(0);
    % ylim([-1.5, 4.5]);
    % ylim([-1.5, 4.5]);
end

%% Row 2: Pairwise correlogram between selected neurons. Row 3: Auto-corrogram for the same neurons.
% f = figure();
% tiles = tiledlayout(2, 2, "TileSpacing", "Compact", "Padding", "Compact");

% selected_neurons = [85, 86];
% selected_neurons = [33, 61];
% selected_neurons = [2, 42];
% selected_neurons = [35, 60];

for i = 1:length(areas)
    % Load data
    meta.area = areas{i};
    meta.prepost = preposts{i};
    meta.state = states{i};
    meta.filename = generate_filename('raster', meta);

    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    raster = raster_data.data.rasters{1};

    r1 = raster(selected_neurons(1), t_range);
    r2 = raster(selected_neurons(2), t_range);
    % smooth_r1 = movmean(r1, smooth_window);
    % smooth_r2 = movmean(r2, smooth_window);

    % Compute correlogram
    [correlogram, lags] = norm_xcorr(r1, r2, corr_range);
    [auto1, ~] = norm_xcorr(r1, r1, corr_range);
    [auto2, ~] = norm_xcorr(r2, r2, corr_range);
    % [correlogram, lags] = norm_xcorr(smooth_r1, smooth_r2, corr_range);
    % [auto1, ~] = norm_xcorr(smooth_r1, smooth_r1, corr_range);
    % [auto2, ~] = norm_xcorr(smooth_r2, smooth_r2, corr_range);

    % smooth_correlogram = movmean(correlogram, smooth_window);
    % smooth_auto1 = movmean(auto1, smooth_window);
    % smooth_auto2 = movmean(auto2, smooth_window);
    smooth_correlogram = same_conv(correlogram, smooth_kernel);
    smooth_auto1 = same_conv(auto1, smooth_kernel);
    smooth_auto2 = same_conv(auto2, smooth_kernel);

    % Row 4/5 spectra are computed from the original signals with the same
    % Welch/cross-spectral framework. This avoids comparing FFT(normalized
    % correlogram) with PSD(original signal), which have different scales.
    [psd_r1, psd_r2, cross_spectrum_abs, coherence, spec_freqs] = compute_pair_spectra( ...
        r1, r2, sample_rate, freqs, psd_window_sec, psd_overlap_frac);

    coherence_plot = smooth_spectrum_for_plot(coherence, spec_smooth_window);
    psd_r1_plot = smooth_spectrum_for_plot(psd_r1, spec_smooth_window);
    psd_r2_plot = smooth_spectrum_for_plot(psd_r2, spec_smooth_window);

    % Compute shuffled controls.
    % - CCG and coherence: shuffle signal 2 relative to signal 1.
    % - Signal PSDs: independently time-permute each signal to preserve the
    %   value distribution / spike count but destroy temporal structure.
    shuffle_correlograms = zeros(shuffle_N, length(correlogram));
    shuffle_coherences = zeros(shuffle_N, numel(spec_freqs));
    shuffle_psd_r1 = zeros(shuffle_N, numel(spec_freqs));
    shuffle_psd_r2 = zeros(shuffle_N, numel(spec_freqs));

    for j = 1:shuffle_N
        fprintf('%d/%d shuffles started\n', j, shuffle_N);

        shuffled_r2 = r2(randperm(length(r2)));
        shuffled_r1 = r1(randperm(length(r1)));

        % Control for lag-domain cross-correlation.
        [shuffle_corr, ~] = norm_xcorr(r1, shuffled_r2, corr_range);
        shuffle_correlograms(j, :) = same_conv(shuffle_corr, smooth_kernel);

        % Control for frequency-domain coupling.
        [~, ~, ~, shuffle_coh, ~] = compute_pair_spectra( ...
            r1, shuffled_r2, sample_rate, freqs, psd_window_sec, psd_overlap_frac);
        shuffle_coherences(j, :) = smooth_spectrum_for_plot(shuffle_coh, spec_smooth_window);

        % Controls for individual signal PSDs.
        [shuffle_pxx1, ~] = compute_signal_psd( ...
            shuffled_r1, sample_rate, freqs, psd_window_sec, psd_overlap_frac);
        [shuffle_pxx2, ~] = compute_signal_psd( ...
            shuffled_r2, sample_rate, freqs, psd_window_sec, psd_overlap_frac);
        shuffle_psd_r1(j, :) = smooth_spectrum_for_plot(shuffle_pxx1, spec_smooth_window);
        shuffle_psd_r2(j, :) = smooth_spectrum_for_plot(shuffle_pxx2, spec_smooth_window);

        fprintf('%d/%d shuffles finished\n', j, shuffle_N);
    end

    % Shuffled controls.
    shuffle_mean = mean_omitnan(shuffle_correlograms);
    shuffle_std = std_omitnan(shuffle_correlograms);

    coh_shuffle_mean = mean_omitnan(shuffle_coherences);
    coh_shuffle_std = std_omitnan(shuffle_coherences);

    psd_r1_shuffle_mean = mean_omitnan(shuffle_psd_r1);
    psd_r1_shuffle_std = std_omitnan(shuffle_psd_r1);

    psd_r2_shuffle_mean = mean_omitnan(shuffle_psd_r2);
    psd_r2_shuffle_std = std_omitnan(shuffle_psd_r2);

    %% Row 2: Cross-correlogram with shuffled control and significant segments.
    row = 2;
    tile = nexttile(i+n_state*(row-1));
    shuffle_upper = shuffle_mean + std_multiplier * shuffle_std;
    shuffle_lower = shuffle_mean - std_multiplier * shuffle_std;

    fill_control_band(tile, lags, shuffle_mean, shuffle_std, std_multiplier, ...
        [0.8, 0.8, 0.8], 'Shuffled mean ± 2SD');
    hold(tile, 'on');
    xline(tile, 0, 'k--', 'HandleVisibility', 'off');
    yline(tile, 0, 'k--', 'HandleVisibility', 'off');
    plot(tile, lags, smooth_correlogram, 'm-', 'LineWidth', 1, ...
        'DisplayName', 'Cross-corr');

    ccg_sig_mask = (smooth_correlogram > shuffle_upper) | ...
                   (smooth_correlogram < shuffle_lower);
    plot_significant_segments(tile, lags, smooth_correlogram, ccg_sig_mask, ...
        sig_min_run_bins, 'k', 'Significant');
    hold(tile, 'off');

    title(tile, sprintf('Correlogram: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Lag (ms)');
    ylabel(tile, 'Normalized correlation');
    legend(tile, 'show', 'Location', 'southeast');
    ylim(tile, [-0.01, 0.01]);

    %% Row 3: Auto-correlograms.
    row = 3;
    tile = nexttile(i+n_state*(row-1));
    plot(tile, lags, smooth_correlogram, 'm-', 'LineWidth', 2, ...
        'DisplayName', 'Cross-corr', 'Color', [1, 0, 1, 0.2]);
    hold(tile, 'on');
    xline(tile, 0, 'k--', 'HandleVisibility', 'off');
    yline(tile, 0, 'k--', 'HandleVisibility', 'off');
    plot(tile, lags, smooth_auto1, 'r-', 'LineWidth', 1, ...
        'DisplayName', 'Auto-corr 1');
    plot(tile, lags, smooth_auto2, 'b-', 'LineWidth', 1, ...
        'DisplayName', 'Auto-corr 2');
    hold(tile, 'off');
    title(tile, sprintf('Auto-correlograms: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Lag (ms)');
    ylabel(tile, 'Normalized correlation');
    legend(tile, 'show', 'Location', 'southeast');
    ylim(tile, [-0.03, 0.03]);

    %% Row 4: Spectral coherence, estimated from the original signals.
    % This replaces FFT(normalized cross-correlogram). Coherence is the
    % normalized cross-spectrum: |Sxy|^2 / (Sxx*Syy), so it is comparable
    % across frequencies and conditions.
    row = 4;
    tile = nexttile(i+n_state*(row-1));
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

    title(tile, sprintf('Spectral coherence: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Frequency (Hz)');
    ylabel(tile, 'Coherence');
    ylim(tile, [0, 1]);
    legend(tile, 'show', 'Location', 'northeast');

    %% Row 5: PSD of the original pair signals, with shuffled controls.
    row = 5;
    tile = nexttile(i+n_state*(row-1));

    fill_control_band(tile, spec_freqs, psd_r1_shuffle_mean, psd_r1_shuffle_std, ...
        std_multiplier, [1.0, 0.85, 0.85], 'Signal 1 shuffled ± 2SD');
    hold(tile, 'on');
    fill_control_band(tile, spec_freqs, psd_r2_shuffle_mean, psd_r2_shuffle_std, ...
        std_multiplier, [0.85, 0.85, 1.0], 'Signal 2 shuffled ± 2SD');

    plot(tile, spec_freqs, psd_r1_plot, 'r-', 'LineWidth', 1, ...
        'DisplayName', 'Signal 1 PSD');
    plot(tile, spec_freqs, psd_r2_plot, 'b-', 'LineWidth', 1, ...
        'DisplayName', 'Signal 2 PSD');

    psd_r1_upper = psd_r1_shuffle_mean + std_multiplier * psd_r1_shuffle_std;
    psd_r2_upper = psd_r2_shuffle_mean + std_multiplier * psd_r2_shuffle_std;
    psd_r1_sig_mask = psd_r1_plot > psd_r1_upper;
    psd_r2_sig_mask = psd_r2_plot > psd_r2_upper;

    plot_significant_segments(tile, spec_freqs, psd_r1_plot, psd_r1_sig_mask, ...
        sig_min_run_bins, [0.5, 0, 0], 'Signal 1 significant');
    plot_significant_segments(tile, spec_freqs, psd_r2_plot, psd_r2_sig_mask, ...
        sig_min_run_bins, [0, 0, 0.5], 'Signal 2 significant');
    hold(tile, 'off');

    title(tile, sprintf('Signal PSD: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Frequency (Hz)');
    ylabel(tile, 'PSD');

    y_candidates = [psd_r1_plot(:); psd_r2_plot(:); ...
        psd_r1_upper(:); psd_r2_upper(:)];
    y_candidates = y_candidates(isfinite(y_candidates));
    if ~isempty(y_candidates) && max(y_candidates) > 0
        ylim(tile, [0, 1.1 * max(y_candidates)]);
    end
    legend(tile, 'show', 'Location', 'northeast');

end

%% Export to pdf
% Save current figure as vector PDF
fig = gcf;

% ----- Figure settings -----
save_folder = fullfile(root, 'Figures', 'Paper');
check_path(save_folder);
filename = fullfile(save_folder, 'Figure1.pdf');

figWidth  = 16.0;   % inches
figHeight = 28.0;   % inches

resolution = 300;  % dpi; mainly affects rasterized components
% -------------------------

% Set figure size on screen
set(fig, 'Units', 'inches');
fig.Position(3:4) = [figWidth, figHeight];

% Set paper size for PDF export
set(fig, 'PaperUnits', 'inches');
set(fig, 'PaperSize', [figWidth, figHeight]);
set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);

% Make background white
set(fig, 'Color', 'w');

% Export as vector PDF
exportgraphics(fig, filename, ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);

% Export as jpg for quick preview
filename = fullfile(save_folder, 'Figure1_preview.jpg');
exportgraphics(fig, filename, ...
    'ContentType', 'image', ...
    'BackgroundColor', 'white', ...
    'Resolution', resolution);

%%
close(fig);


% %% Fig 2: J counting bar plot
% % counting
% pos_counts = zeros(1, length(areas));
% neg_counts = zeros(1, length(areas));
% max_counts = zeros(1, length(areas));
% for i = 1:length(areas)
%     % Load data
%     meta.area = areas{i};
%     meta.prepost = preposts{i};
%     meta.state = states{i};
%     meta.shuffle_idx = 0;
%     meta.kernel_name = 'DeltaPure';
%     meta.reg_name = 'L2=0_2';
%     meta.epoch = 3000;
%     meta.fold_idx = 0;

%     % Load raster data for neuron area info
%     meta.filename = generate_filename('raster', meta);
%     raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
%     fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
%     fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
%     fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

%     % split neurons by area
%     cell_area = raster_data.data.cell_area;
%     filter1 = ismember(cell_area, {'ACC'}); % filter for ACC neurons
%     filter2 = ismember(cell_area, {'VLPFC'}); % filter for VLPFC neurons

%     % Load GLM data for connectivity info
%     meta.filename = generate_filename('GLM', meta);
%     GLM_data = load(fullfile(root, 'Data', 'Working', 'GLM', meta.filename));
%     fprintf('Loaded GLM data for %s %s %s\n', meta.prepost, meta.state, meta.area);
%     N = GLM_data.meta.N;
%     J = GLM_data.data.model_par(:, (2:N+1)); % kernel 1 weights
%     err = GLM_data.data.model_err.total(:, (2:N+1));

%     % Count significant J
%     pos_count = sum(J(filter1, filter2) > err_multi*err(filter1, filter2), 'all') + ...
%     sum(J(filter2, filter1) > err_multi*err(filter2, filter1), 'all');
%     neg_count = sum(J(filter1, filter2) < -err_multi*err(filter1, filter2), 'all') + ...
%     sum(J(filter2, filter1) < -err_multi*err(filter2, filter1), 'all');
%     max_count = numel(J(filter1, filter2)) + numel(J(filter2, filter1));
%     pos_counts(i) = pos_count;
%     neg_counts(i) = neg_count;
%     max_counts(i) = max_count;
% end

% % Plot
% f = figure();
% % x axis: states. Groups: positive vs negative connections.
% pos_ratios = pos_counts ./ max_counts;
% neg_ratios = neg_counts ./ max_counts;




%% functions
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