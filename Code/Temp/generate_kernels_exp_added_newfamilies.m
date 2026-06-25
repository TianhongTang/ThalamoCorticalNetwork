%% generate_kernels_exp_added_newfamilies.m - Generate long-timescale GLM kernels
% New kernel families only. All kernels in this script use the same length.
%
% Output format follows existing kernel files:
%   data.conn_kernels, data.PS_kernels
%   meta.kernel_name, meta.n_conn_kernel, meta.n_PS_kernel, meta.kernel_len
%
% Families generated:
%   1. LongExp*       : exponential decay kernels, tau = 5 ... 1000 ms
%   2. LongGaussC*   : Gaussian kernels, center = 5 ... 1000 ms
%   3. LongGDeriv*   : positive Gaussian-derivative-like kernels
%   4. LongStepB*    : non-overlapping adaptive box kernels
%
% For every family:
%   - single-kernel files are generated with one connection kernel per file
%   - one composite file is generated with all kernels in that family
%   - composite kernel shapes are plotted
%
% All newly generated kernels are connection-only: PS_kernels = {}.

clear;
set(0, 'DefaultFigureVisible', 'off');

%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
if ~isfolder(fullfile(root, 'Data')) && isfolder(fullfile(pwd, 'Data'))
    root = pwd;
end
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Shared settings
save_folder = fullfile(root, 'Data', 'Working', 'kernel');
check_path(save_folder);

T = 3000;              % ms
sample_dt = 1;         % ms, current code assumes 1 ms bins
t = 0:sample_dt:T;
kernel_len = numel(t);

% Log-like spacing from fast to slow timescales.
time_scales_ms = [5, 10, 20, 40, 80, 160, 320, 640, 1000];

%% 1. Long exponential decay kernels
exp_kernels = cell(1, numel(time_scales_ms));
exp_labels = cell(1, numel(time_scales_ms));

for idx = 1:numel(time_scales_ms)
    tau = time_scales_ms(idx);
    k = exp(-t ./ tau);
    k = normalize_kernel(k);

    kernel_name = sprintf('LongExp%d', tau);
    save_kernel_set(save_folder, kernel_name, {k}, {}, kernel_len);

    exp_kernels{idx} = k;
    exp_labels{idx} = sprintf('Exp tau=%d ms', tau);
end
save_kernel_set(save_folder, 'LongExpMix', exp_kernels, {}, kernel_len);
plot_kernel_family(save_folder, 'LongExpMix', t, exp_kernels, exp_labels);

%% 2. Long Gaussian kernels
% Centers follow the same timescale grid. Width is adaptive: narrow early,
% broader late, with a lower bound to avoid near-delta kernels.
gauss_kernels = cell(1, numel(time_scales_ms));
gauss_labels = cell(1, numel(time_scales_ms));

for idx = 1:numel(time_scales_ms)
    center = time_scales_ms(idx);
    sigma = adaptive_gaussian_sigma(center);
    k = exp(-((t - center).^2) ./ (2 * sigma^2));
    k = normalize_kernel(k);

    kernel_name = sprintf('LongGaussC%d', center);
    save_kernel_set(save_folder, kernel_name, {k}, {}, kernel_len);

    gauss_kernels{idx} = k;
    gauss_labels{idx} = sprintf('Gauss c=%d, sigma=%.1f ms', center, sigma);
end
save_kernel_set(save_folder, 'LongGaussMix', gauss_kernels, {}, kernel_len);
plot_kernel_family(save_folder, 'LongGaussMix', t, gauss_kernels, gauss_labels);

%% 3. Positive Gaussian-derivative-like kernels
% Shape: t * exp(-t^2 / (2 tau^2)). This is the positive-time lobe of a
% Gaussian derivative up to scaling/sign, with one timescale parameter tau.
gderiv_kernels = cell(1, numel(time_scales_ms));
gderiv_labels = cell(1, numel(time_scales_ms));

for idx = 1:numel(time_scales_ms)
    tau = time_scales_ms(idx);
    k = t .* exp(-(t.^2) ./ (2 * tau^2));
    k = normalize_kernel(k);

    kernel_name = sprintf('LongGDeriv%d', tau);
    save_kernel_set(save_folder, kernel_name, {k}, {}, kernel_len);

    gderiv_kernels{idx} = k;
    gderiv_labels{idx} = sprintf('t exp(-t^2/2tau^2), tau=%d ms', tau);
end
save_kernel_set(save_folder, 'LongGDerivMix', gderiv_kernels, {}, kernel_len);
plot_kernel_family(save_folder, 'LongGDerivMix', t, gderiv_kernels, gderiv_labels);

%% 4. Adaptive non-overlapping box kernels
% Edges are anchored to the same timescale grid and then extended to T.
% These boxes are non-overlapping and cover the full [0, T] time axis.
step_edges = unique([0, time_scales_ms, T]);
if step_edges(end) ~= T
    step_edges = [step_edges, T];
end
step_edges = sort(step_edges);
step_kernels = cell(1, numel(step_edges) - 1);
step_labels = cell(1, numel(step_edges) - 1);

for idx = 1:(numel(step_edges) - 1)
    start_ms = step_edges(idx);
    end_ms = step_edges(idx + 1);

    k = zeros(size(t));
    if idx < numel(step_edges) - 1
        mask = (t >= start_ms) & (t < end_ms);
    else
        mask = (t >= start_ms) & (t <= end_ms);
    end
    k(mask) = 1;
    k = normalize_kernel(k);

    kernel_name = sprintf('LongStepB%d_%d', start_ms, end_ms);
    save_kernel_set(save_folder, kernel_name, {k}, {}, kernel_len);

    step_kernels{idx} = k;
    step_labels{idx} = sprintf('Box %d-%d ms', start_ms, end_ms);
end
save_kernel_set(save_folder, 'LongStepAdaptive', step_kernels, {}, kernel_len);
plot_kernel_family(save_folder, 'LongStepAdaptive', t, step_kernels, step_labels);

fprintf('Generated new long kernel families in:\n  %s\n', save_folder);
fprintf('Shared kernel length: %d samples (%d ms at %d ms/bin).\n', kernel_len, T, sample_dt);

%% Local helper functions
function k = normalize_kernel(k)
    k = double(k(:))';
    k(~isfinite(k)) = 0;
    s = sum(abs(k));
    if s == 0
        error('Kernel has zero total mass and cannot be normalized.');
    end
    k = k ./ s;
end

function sigma = adaptive_gaussian_sigma(center)
    % Keep very early kernels reasonably narrow, and broaden later kernels.
    sigma = max(3, 0.25 * center);
end

function save_kernel_set(save_folder, kernel_name, conn_kernels, PS_kernels, kernel_len)
    n_conn_kernel = numel(conn_kernels);
    n_PS_kernel = numel(PS_kernels);

    for i = 1:n_conn_kernel
        assert(isrow(conn_kernels{i}), 'conn_kernels{%d} must be a row vector.', i);
        assert(numel(conn_kernels{i}) == kernel_len, ...
            'conn_kernels{%d} has length %d, expected %d.', i, numel(conn_kernels{i}), kernel_len);
    end
    for i = 1:n_PS_kernel
        assert(isrow(PS_kernels{i}), 'PS_kernels{%d} must be a row vector.', i);
        assert(numel(PS_kernels{i}) == kernel_len, ...
            'PS_kernels{%d} has length %d, expected %d.', i, numel(PS_kernels{i}), kernel_len);
    end

    meta = struct('kernel_name', kernel_name, ...
                  'n_conn_kernel', n_conn_kernel, ...
                  'n_PS_kernel', n_PS_kernel, ...
                  'kernel_len', kernel_len);
    data = struct('conn_kernels', {conn_kernels}, 'PS_kernels', {PS_kernels});
    save_name = generate_filename('kernel', meta);
    meta.file_name = save_name;
    save_path = fullfile(save_folder, save_name);
    save(save_path, 'meta', 'data');
end

function plot_kernel_family(save_folder, kernel_name, t, conn_kernels, legend_names)
    f = figure('Visible', 'off', 'Color', 'w');
    hold on;
    for i = 1:numel(conn_kernels)
        plot(t, conn_kernels{i}, 'LineWidth', 1);
    end
    hold off;
    xlabel('Time (ms)');
    ylabel('Normalized amplitude');
    title(['Kernel: ', kernel_name], 'Interpreter', 'none');
    legend(legend_names, 'Location', 'northeastoutside', 'Interpreter', 'none');
    xlim([0, max(t)]);
    saveas(f, fullfile(save_folder, ['kernel_', kernel_name, '.png']));
    close(f);
end
