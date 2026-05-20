%% Figure 2: Example session population-level ACC-VLPFC analyses.
% Row 1: network plot.
% Rows 2-5 correspond to Figure 1 rows 2-5, but use population-average
% firing rates from all ACC and VLPFC neurons in this session.

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
% Figure 1 row 1 = raster; Figure 2 row 1 = network plot.
f = figure('Color', 'w');
tiles = tiledlayout(5, n_state, "TileSpacing", "Compact", "Padding", "Compact");

%% parameters
shuffle_N = 10;       % Increase for final statistical bands, e.g. 100.
err_multi = 2;       % threshold for significant J, in multiples of GLM error estimate.
t_range = 1:60000;
corr_range = 150;   % ms, assuming 1 kHz sampling.
smooth_window = 3;  % ms.
spec_smooth_window = 5; % frequency-bin smoothing for FFT/PSD visualization only.
psd_window_sec = 10; % Welch PSD window length for row 5 signal spectra.
psd_overlap_frac = 0.5; % Welch PSD fractional overlap.
sample_rate = 1000; % Hz.
freqs = linspace(0, 100, 201); % Hz.

% Same triangular smoothing kernel as Figure 1.
smooth_kernel = 1 - abs(-smooth_window:smooth_window) / smooth_window;
smooth_kernel = smooth_kernel / sum(smooth_kernel);

% Optional: highlight the same example pair used in Figure 1 if plot_network
% supports highlighted node indices. These are global neuron indices.
highlight_global_pair = [2, 42];

%% Row 1: network plot
% for i = 1:n_state
%     % Load metadata and raster data for neuron area information.
%     meta.area = areas{i};
%     meta.prepost = preposts{i};
%     meta.state = states{i};
%     meta.shuffle_idx = 0;
%     meta.kernel_name = 'DeltaPure';
%     meta.reg_name = 'L2=0_2';
%     meta.epoch = 3000;
%     meta.fold_idx = 0;
% 
%     meta.filename = generate_filename('raster', meta);
%     raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
%     fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
%     fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
%     fprintf('Trial_num: %d\n', raster_data.meta.trial_num);
% 
%     cell_area = raster_data.data.cell_area;
%     [filter_acc, filter_vlpfc] = get_area_filters(cell_area);
%     check_area_filters(filter_acc, filter_vlpfc, meta);
% 
%     % Convert global neuron indices to local ACC/VLPFC indices for plotting.
%     acc_global_idx = find(filter_acc);
%     vlpfc_global_idx = find(filter_vlpfc);
%     highlight_acc_idx = find(acc_global_idx == highlight_global_pair(1), 1);
%     highlight_vlpfc_idx = find(vlpfc_global_idx == highlight_global_pair(2), 1);
% 
%     % Load GLM data for connectivity information.
%     meta.filename = generate_filename('GLM', meta);
%     GLM_data = load(fullfile(root, 'Data', 'Working', 'GLM', meta.filename));
%     fprintf('Loaded GLM data for %s %s %s\n', meta.prepost, meta.state, meta.area);
% 
%     N = GLM_data.meta.N;
%     J = GLM_data.data.model_par(:, 2:N+1);          % kernel 1 weights.
%     err = GLM_data.data.model_err.total(:, 2:N+1);  % error estimates.
% 
%     tile = nexttile(i);  % row 1.
%     call_plot_network(tile, ...
%         J(filter_acc, filter_vlpfc), J(filter_vlpfc, filter_acc), ...
%         err(filter_acc, filter_vlpfc), err(filter_vlpfc, filter_acc), ...
%         err_multi, highlight_acc_idx, highlight_vlpfc_idx);
% 
%     title(tile, sprintf('Network: %s, %s', meta.prepost, meta.state));
% end

%% Rows 2-5: population firing correlation and spectra
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
    t_idx = safe_time_range(t_range, size(raster, 2), meta);

    cell_area = raster_data.data.cell_area;
    [filter_acc, filter_vlpfc] = get_area_filters(cell_area);
    check_area_filters(filter_acc, filter_vlpfc, meta);

    r_acc = raster(filter_acc, t_idx);
    r_vlpfc = raster(filter_vlpfc, t_idx);

    N_acc = sum(filter_acc);
    N_vlpfc = sum(filter_vlpfc);

    % Population-average firing activity. Because the raster is sampled at
    % 1 kHz, mean(spikes per ms) * 1000 gives Hz for optional time-domain plots.
    pop_acc = mean(r_acc, 1);
    pop_vlpfc = mean(r_vlpfc, 1);

    % Correlation and auto-correlation of population-average activity.
    [correlogram, lags] = norm_xcorr(pop_acc, pop_vlpfc, corr_range);
    [auto_acc, ~] = norm_xcorr(pop_acc, pop_acc, corr_range);
    [auto_vlpfc, ~] = norm_xcorr(pop_vlpfc, pop_vlpfc, corr_range);

    smooth_correlogram = same_conv(correlogram, smooth_kernel);
    smooth_auto_acc = same_conv(auto_acc, smooth_kernel);
    smooth_auto_vlpfc = same_conv(auto_vlpfc, smooth_kernel);

    % Shuffled control. Circular shift preserves the temporal structure of
    % the VLPFC population signal while disrupting alignment with ACC.
    shuffle_correlograms = zeros(shuffle_N, length(correlogram));
    for j = 1:shuffle_N
        fprintf('%d/%d population shuffles started\n', j, shuffle_N);
        shift_val = random_circular_shift(length(pop_vlpfc), corr_range);
        shuffled_vlpfc = circshift(pop_vlpfc, shift_val);
        [shuffle_corr, ~] = norm_xcorr(pop_acc, shuffled_vlpfc, corr_range);
        shuffle_correlograms(j, :) = same_conv(shuffle_corr, smooth_kernel);
        fprintf('%d/%d population shuffles finished\n', j, shuffle_N);
    end
    shuffle_mean = mean(shuffle_correlograms, 1);
    shuffle_std = std(shuffle_correlograms, [], 1);

    % Row 2: population cross-correlogram.
    row = 2;
    tile = nexttile(i + n_state*(row-1));
    std_multiplier = 2;
    shuffle_upper = shuffle_mean + std_multiplier * shuffle_std;
    shuffle_lower = shuffle_mean - std_multiplier * shuffle_std;

    fill(tile, [lags, fliplr(lags)], [shuffle_upper, fliplr(shuffle_lower)], ...
        [0.8, 0.8, 0.8], 'EdgeColor', 'none', ...
        'DisplayName', 'Shifted mean ± 2SD');
    hold(tile, 'on');
    xline(tile, 0, 'k--', 'HandleVisibility', 'off');
    yline(tile, 0, 'k--', 'HandleVisibility', 'off');
    plot(tile, lags, smooth_correlogram, 'm-', 'LineWidth', 1, ...
        'DisplayName', 'ACC-VLPFC cross-corr');
    hold(tile, 'off');

    title(tile, sprintf('Population correlogram: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Lag (ms)');
    ylabel(tile, 'Normalized correlation');
    legend(tile, 'show', 'Location', 'southeast');
    ylim(tile, [-0.05, 0.05]);

    % Row 3: population auto-correlograms, with cross-correlation in gray/magenta.
    row = 3;
    tile = nexttile(i + n_state*(row-1));
    plot(tile, lags, smooth_correlogram, 'm-', 'LineWidth', 1, ...
        'DisplayName', 'Cross-corr');
    hold(tile, 'on');
    xline(tile, 0, 'k--', 'HandleVisibility', 'off');
    yline(tile, 0, 'k--', 'HandleVisibility', 'off');
    plot(tile, lags, smooth_auto_acc, 'r-', 'LineWidth', 1, ...
        'DisplayName', sprintf('ACC auto-corr, N=%d', N_acc));
    plot(tile, lags, smooth_auto_vlpfc, 'b-', 'LineWidth', 1, ...
        'DisplayName', sprintf('VLPFC auto-corr, N=%d', N_vlpfc));
    hold(tile, 'off');

    title(tile, sprintf('Population auto-correlograms: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Lag (ms)');
    ylabel(tile, 'Normalized correlation');
    legend(tile, 'show', 'Location', 'southeast');
    ylim(tile, [-0.05, 0.05]);

    % Row 4: spectrum of the population cross-correlogram.
    % Remove DC before FFT. Do not smooth before FFT; smooth only the
    % plotted spectrum for visualization.
    spec_Xcorr = myFT(remove_dc(correlogram), sample_rate, freqs);
    spec_Xcorr_plot = smooth_spectrum_for_plot(spec_Xcorr, spec_smooth_window);

    row = 4;
    tile = nexttile(i + n_state*(row-1));
    plot(tile, freqs, spec_Xcorr_plot, 'm-', 'LineWidth', 1, ...
        'DisplayName', 'Cross-corr spectrum');
    title(tile, sprintf('FFT of population cross-corr: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Frequency (Hz)');
    ylabel(tile, 'Amplitude');
    legend(tile, 'show', 'Location', 'northeast');

    % Row 5: PSD of the original population signals only.
    % Do not mix population-signal PSD with FFT(auto-corr), because those have
    % different normalization and scale in this script.
    [psd_acc, psd_freqs] = compute_signal_psd(remove_dc(pop_acc), sample_rate, freqs, ...
        psd_window_sec, psd_overlap_frac);
    [psd_vlpfc, ~] = compute_signal_psd(remove_dc(pop_vlpfc), sample_rate, freqs, ...
        psd_window_sec, psd_overlap_frac);

    psd_acc_plot = smooth_spectrum_for_plot(psd_acc, spec_smooth_window);
    psd_vlpfc_plot = smooth_spectrum_for_plot(psd_vlpfc, spec_smooth_window);

    row = 5;
    tile = nexttile(i + n_state*(row-1));
    plot(tile, psd_freqs, psd_acc_plot, 'r-', 'LineWidth', 1, ...
        'DisplayName', 'ACC population signal PSD');
    hold(tile, 'on');
    plot(tile, psd_freqs, psd_vlpfc_plot, 'b-', 'LineWidth', 1, ...
        'DisplayName', 'VLPFC population signal PSD');
    hold(tile, 'off');

    title(tile, sprintf('PSD of population signals: %s, %s', meta.prepost, meta.state));
    xlabel(tile, 'Frequency (Hz)');
    ylabel(tile, 'PSD');
    legend(tile, 'show', 'Location', 'northeast');
end

%% Export to pdf and preview image
fig = gcf;

save_folder = fullfile(root, 'Figures', 'Paper');
check_path(save_folder);

figWidth  = 16.0;   % inches.
figHeight = 28.0;   % inches.
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

function t_idx = safe_time_range(t_range, raster_len, meta)
    t_idx = t_range(t_range >= 1 & t_range <= raster_len);
    if isempty(t_idx)
        error('Requested t_range is outside raster length for %s %s %s.', ...
            meta.prepost, meta.state, meta.area);
    end
    if numel(t_idx) < numel(t_range)
        warning('Truncated t_range for %s %s %s to fit raster length %d.', ...
            meta.prepost, meta.state, meta.area, raster_len);
    end
end

function shift_val = random_circular_shift(signal_len, corr_range)
    % Avoid tiny shifts that leave near-zero-lag structure mostly unchanged.
    min_shift = min(corr_range + 1, signal_len);
    max_shift = max(signal_len - corr_range - 1, 1);

    if max_shift > min_shift
        shift_val = randi([min_shift, max_shift]);
    else
        shift_val = randi(signal_len);
    end
end

function call_plot_network(ax, J12, J21, err12, err21, err_multi, highlight_i, highlight_j)
    % plot_network appears to support highlighted ACC/VLPFC local indices in
    % the current project. This wrapper keeps the script runnable if that
    % signature differs across versions.
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

function [correlogram, lags] = norm_xcorr(r1, r2, max_lag)
%NORM_XCORR Pearson-normalized cross-correlation between two signals.
%
%   Positive lag means r1 is delayed relative to r2.

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

        x = x - mean(x);
        y = y - mean(y);

        denom = sqrt(sum(x.^2) * sum(y.^2));

        if denom > 0
            correlogram(i) = sum(x .* y) / denom;
        else
            correlogram(i) = NaN;
        end
    end
end

function convolved = same_conv(signal, kernel)
    % Convolve signal and kernel, return the central part. Normalize borders
    % by the overlapping kernel weight.
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

function [pxx, f] = compute_signal_psd(signal, sample_rate, freqs, window_sec, overlap_frac)
    % Estimate signal PSD with Welch averaging. The input should already be
    % DC-removed. If pwelch is unavailable, fall back to a power-like direct
    % spectrum so the plotting script still runs.
    signal = signal(:);
    signal = remove_dc(signal);
    freqs = freqs(:);
    N = length(signal);

    if N < 4
        f = freqs.';
        pxx = nan(size(f));
        return;
    end

    win_len = min(max(round(window_sec * sample_rate), 4), N);
    if win_len == 1
        win = 1;
    else
        win = 0.5 - 0.5*cos(2*pi*(0:win_len-1).'/(win_len-1)); % Hann window without toolbox dependency.
    end
    noverlap = min(floor(overlap_frac * win_len), win_len - 1);

    try
        [pxx_col, f_col] = pwelch(signal, win, noverlap, freqs, sample_rate);
        pxx = pxx_col(:).';
        f = f_col(:).';
    catch ME
        warning('pwelch failed (%s). Falling back to direct power-like spectrum.', ME.message);
        amp = myFT(signal, sample_rate, freqs);
        pxx = amp(:).'.^2;
        f = freqs(:).';
    end
end

function amp = myFT(signal, sample_rate, freq)
%MYFT Real amplitude spectrum at specified frequencies.

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

    dc_idx = (freq == 0);
    amp(dc_idx) = abs(spec(dc_idx)) / N;

    amp = reshape(amp, freq_shape);
end
