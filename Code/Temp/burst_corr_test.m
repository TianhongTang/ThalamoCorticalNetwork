%% Simulate two binary signals driven by noisy phase oscillators

clear; clc; close all;

%% Parameters
dt = 0.001;              % 1 ms
T_total = 20;            % seconds
t = 0:dt:T_total-dt;
N = length(t);

period = 0.200;          % 200 ms
freq = 1 / period;       % 5 Hz

mean_dphi = 2*pi*dt/period;   % expected phase increase per time step
std_dphi  = 5 * mean_dphi; % phase increment noise level

phase_diff = 0;            % initial / fixed phase offset

base_prob = 0.05;
mod_amp = 0.04;

max_lag_ms = 1000;
max_lag = round(max_lag_ms / 1000 / dt);

%% Phase simulation
% Option 1: same phase noise, fixed phase difference
% Option 2: independent phase noise, phase difference slowly drifts

use_same_phase_noise = true;

dphi1 = mean_dphi + std_dphi * randn(1, N);

if use_same_phase_noise
    dphi2 = dphi1;   % phase difference remains fixed
else
    dphi2 = mean_dphi + std_dphi * randn(1, N); % phase difference drifts
end

phi1 = cumsum(dphi1);
phi2 = cumsum(dphi2) + phase_diff;

%% Hidden variables
h1 = sin(phi1);
h2 = sin(phi2);

%% Convert hidden variables to probabilities, with noise
sigmoid = @(x) 1 ./ (1 + exp(-x));
logit = @(p) log(p ./ (1 - p));

base_prob = 0.05;        % desired baseline probability
mod_gain = 1.2;          % modulation strength in logit space
prob_noise_std = 0.4;    % noise level in logit space

eta1_clean = logit(base_prob) + mod_gain * h1;
eta2_clean = logit(base_prob) + mod_gain * h2;

eta1 = eta1_clean + prob_noise_std * randn(size(t));
eta2 = eta2_clean + prob_noise_std * randn(size(t));

p1_clean = sigmoid(eta1_clean);
p2_clean = sigmoid(eta2_clean);

p1 = sigmoid(eta1);
p2 = sigmoid(eta2);

%% Generate binary signals
sig1 = rand(size(t)) < p1;
sig2 = rand(size(t)) < p2;

%% Correlograms
x1 = double(sig1) - mean(sig1);
x2 = double(sig2) - mean(sig2);

[acg1, lags] = xcorr(x1, x1, max_lag, 'coeff');
[acg2, ~]    = xcorr(x2, x2, max_lag, 'coeff');
[ccg, ~]     = xcorr(x1, x2, max_lag, 'coeff');

lags_ms = lags * dt * 1000;

%% Plot hidden variables and phase increments
figure;

subplot(5,1,1);
plot(t, dphi1, 'LineWidth', 1); hold on;
plot(t, dphi2, 'LineWidth', 1);
xlim([0 2]);
ylabel('\Delta phase');
legend('dphi1', 'dphi2');
title('Noisy phase increments');

subplot(5,1,2);
plot(t, phi1, 'LineWidth', 1); hold on;
plot(t, phi2, 'LineWidth', 1);
xlim([0 2]);
ylabel('Phase');
legend('phi1', 'phi2');
title('Accumulated phase');

subplot(5,1,3);
plot(t, h1, 'LineWidth', 1.2); hold on;
plot(t, h2, 'LineWidth', 1.2);
xlim([0 2]);
ylabel('Hidden variable');
legend('h1', 'h2');
title('Noisy oscillatory hidden variables');

subplot(5,1,4);
plot(t, p1, 'LineWidth', 1.2); hold on;
plot(t, p2, 'LineWidth', 1.2);
xlim([0 2]);
ylabel('P(signal = 1)');
legend('p1', 'p2');
title('Time-varying probability');

subplot(5,1,5);
plot(t, sig1, '.'); hold on;
plot(t, sig2 + 1.2, '.');
xlim([0 2]);
ylim([-0.2 2.4]);
xlabel('Time (s)');
ylabel('Binary signals');
yticks([0 1.2]);
yticklabels({'sig1', 'sig2'});
title('Generated binary signals');

%% Plot correlograms
figure;

subplot(3,1,1);
plot(lags_ms, acg1, 'LineWidth', 1.5);
xlabel('Lag (ms)');
ylabel('Correlation');
title('Auto correlogram: signal 1');
xline(0, '--');

subplot(3,1,2);
plot(lags_ms, acg2, 'LineWidth', 1.5);
xlabel('Lag (ms)');
ylabel('Correlation');
title('Auto correlogram: signal 2');
xline(0, '--');

subplot(3,1,3);
plot(lags_ms, ccg, 'LineWidth', 1.5);
xlabel('Lag (ms)');
ylabel('Correlation');
title('Cross correlogram: signal 1 vs signal 2');
xline(0, '--');