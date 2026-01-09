%% basic_statistics_main.m -  Basic statistics of each session/state: Synchrony index, firing rate, fano factor, etc.


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

%% Main
states = {'RestOpen', 'RestClose'};
state_names = {'Eye Open', 'Eye Closed'};
n_states = length(states);
% histogram of fano factor and synchrony chi
data_folder = fullfile(root, 'Data', 'Working', 'basic');
check_path(data_folder);
load(fullfile(data_folder, 'basic_statistics.mat'), 'fano_factor_all', 'sync_chi_all');

fano_diff = fano_factor_all(:, 2)-fano_factor_all(:, 1);
sync_diff = sync_chi_all(:, 2)-sync_chi_all(:, 1);

gap_idx = [241, 608, 630, 644];
gap_filter = ismember(1:length(fano_diff), gap_idx);
fano_diff_gap = fano_diff(gap_filter);
sync_diff_gap = sync_diff(gap_filter);
fano_diff = fano_diff(~gap_filter);
sync_diff = sync_diff(~gap_filter);

nan_filter = isnan(fano_diff) | isnan(sync_diff);
fano_diff = fano_diff(~nan_filter);
sync_diff = sync_diff(~nan_filter);
nan_filter_gap = isnan(fano_diff_gap) | isnan(sync_diff_gap);
fano_diff_gap = fano_diff_gap(~nan_filter_gap);
sync_diff_gap = sync_diff_gap(~nan_filter_gap);

figure();
t = tiledlayout(3, 1);

nexttile;
% hold on;
% for state_idx = 1:2
%     histogram(fano_factor_all(:, state_idx), 'Normalization', 'pdf', 'FaceAlpha', 0.5);
%     title(sprintf('Fano Factor'));
%     xlabel('Fano Factor');
%     ylabel('Probability Density');
% end
% hold off;
histogram(fano_diff, 'Normalization', 'pdf', 'FaceAlpha', 0.5,'BinWidth',0.2);
title(sprintf('Fano Factor Close - open'));
xlabel('Fano Factor');
ylabel('Probability Density');
% legend(state_names);

nexttile;
% hold on;
% for state_idx = 1:2
%     histogram(sync_chi_all(:, state_idx), 'Normalization', 'pdf', 'FaceAlpha', 0.5);
%     title(sprintf('Synchrony Chi - %s', state_names{state_idx}));
%     xlabel('Synchrony Chi');
%     ylabel('Probability Density');
% end
% hold off;
histogram(sync_diff, 'Normalization', 'pdf', 'FaceAlpha', 0.5,'BinWidth',0.005);
title(sprintf('Synchrony Chi Close - open'));
xlabel('Synchrony Chi');
ylabel('Probability Density');
% legend(state_names);

%% pearson correlation and p value
R = corrcoef(fano_diff, sync_diff);
[~, p] = corr(fano_diff, sync_diff);

nexttile;
hold on;
scatter(fano_diff, sync_diff, 'b', 'filled');
scatter(fano_diff_gap, sync_diff_gap, 'r', 'filled');
plot([0, 0], [-0.1, 0.06], 'k-');
plot([-10, 15], [0, 0], 'k-');

% linear fit
p_fit = polyfit(fano_diff, sync_diff, 1);
x_fit = linspace(min(fano_diff), max(fano_diff), 100);
y_fit = polyval(p_fit, x_fit);
plot(x_fit, y_fit, 'r--', 'LineWidth', 1);

xlabel('fano difference (Close - Open)');
ylabel('sync difference (Close - Open)');
title(sprintf('Fano vs Sync diff, Rho=%.2f, p=%.3f', R(1,2), p));

% also put text of R and p
text_loc_x = min(fano_diff) + 0.1 * (max(fano_diff) - min(fano_diff));
text_loc_y = max(sync_diff) - 0.1 * (max(sync_diff) - min(sync_diff));
text(text_loc_x, text_loc_y, sprintf('Rho=%.2f\np=%.3f', R(1,2), p), 'FontSize', 10, 'BackgroundColor', 'w');
    
hold off;

figure_folder = fullfile(root, 'Figures', 'basics');
check_path(figure_folder);
figure_name = 'basic_statistics_histograms.png';
saveas(gcf, fullfile(figure_folder, figure_name));