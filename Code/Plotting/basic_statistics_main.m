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

figure();
% t = tiledlayout(3, 1);

% nexttile;
% % hold on;
% % for state_idx = 1:2
% %     histogram(fano_factor_all(:, state_idx), 'Normalization', 'pdf', 'FaceAlpha', 0.5);
% %     title(sprintf('Fano Factor'));
% %     xlabel('Fano Factor');
% %     ylabel('Probability Density');
% % end
% % hold off;
% histogram(fano_factor_all(:, 2)-fano_factor_all(:, 1), 'Normalization', 'pdf', 'FaceAlpha', 0.5,'BinWidth',0.1);
% title(sprintf('Fano Factor Close - open'));
% xlabel('Fano Factor');
% ylabel('Probability Density');
% % legend(state_names);
% 
% nexttile;
% % hold on;
% % for state_idx = 1:2
% %     histogram(sync_chi_all(:, state_idx), 'Normalization', 'pdf', 'FaceAlpha', 0.5);
% %     title(sprintf('Synchrony Chi - %s', state_names{state_idx}));
% %     xlabel('Synchrony Chi');
% %     ylabel('Probability Density');
% % end
% % hold off;
% histogram(sync_chi_all(:, 2)-sync_chi_all(:, 1), 'Normalization', 'pdf', 'FaceAlpha', 0.5,'BinWidth',0.001);
% title(sprintf('Synchrony Chi Close - open'));
% xlabel('Synchrony Chi');
% ylabel('Probability Density');
% % legend(state_names);

nexttile;
scatter(fano_factor_all(:, 2)-fano_factor_all(:, 1), sync_chi_all(:, 2)-sync_chi_all(:, 1));
xlabel('fano difference');
ylabel('sync difference');

figure_folder = fullfile(root, 'Results', 'Figures', 'basics');
check_path(figure_folder);
figure_name = 'basic_statistics_histograms.png';
saveas(gcf, fullfile(figure_folder, figure_name));