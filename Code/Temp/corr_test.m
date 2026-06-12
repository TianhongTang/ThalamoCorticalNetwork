%% Generate two independent normal-distributed variables
clear; clc; close all;

% rng(24);          % for reproducibility
N = 1000;        % number of points

x = randn(N, 1); % x coordinates ~ N(0,1)
y = randn(N, 1); % y coordinates ~ N(0,1), independent from x

mix_ratio = 0.03;
xx = (1-mix_ratio)*x+mix_ratio*y;
yy = (1-mix_ratio)*y+mix_ratio*x;
x=xx; y=yy;

%% Scatter plot
figure;
scatter(x, y, 20, 'filled');
xlabel('x');
ylabel('y');
title('Scatter plot of two independent normal variables');
axis equal;
grid on;

%% Pearson correlation
[r_pearson, p_pearson] = corr(x, y, 'Type', 'Pearson');

%% Spearman correlation
[r_spearman, p_spearman] = corr(x, y, 'Type', 'Spearman');

%% Print results
fprintf('Pearson correlation:\n');
fprintf('r = %.4f, p = %.4g\n\n', r_pearson, p_pearson);

fprintf('Spearman correlation:\n');
fprintf('rho = %.4f, p = %.4g\n', r_spearman, p_spearman);