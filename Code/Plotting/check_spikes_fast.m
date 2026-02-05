%% check_spikes_fast.m - Quick visualization of rasters on screen

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
raster_file_folder = fullfile(root, 'Data', 'Working', 'raster');
sort_idx_folder = fullfile(root, 'Data', 'Working', 'sort_idx');

files = {'raster_EmperorMusPreRestOpenCortexAlignLast_2.mat';
         'raster_EmperorMusPreRestOpenCortexAlignLast_2.mat';
         'raster_EmperorMusPreRestOpenCortexAlignLast_2.mat';
         'raster_EmperorMusPostRestOpenCortexAlignLast_2.mat';
         'raster_EmperorMusPostRestOpenCortexAlignLast_2.mat';
         'raster_EmperorMusPostRestOpenCortexAlignLast_2.mat';
};
idx_files = {'sortidx_EmperorMusPreCortex_2.mat';
              'sortidx_EmperorMusPreCortex_2.mat';
              'sortidx_EmperorMusPreCortex_2.mat';
              'sortidx_EmperorMusPostCortex_2.mat';
              'sortidx_EmperorMusPostCortex_2.mat';
              'sortidx_EmperorMusPostCortex_2.mat';
};

neuron_idxs = {8, 34, 36, 8, 34, 36};
neuron_num = numel(neuron_idxs);

plotted_rasters = cell(neuron_num, 1);

figure;
hold on;

for i = 1:neuron_num
    file_path = fullfile(raster_file_folder, files{i});
    load(file_path, 'rasters');
    sort_idx_path = fullfile(sort_idx_folder, idx_files{i});
    load(sort_idx_path, 'sort_idx');
    neuron_idx = neuron_idxs{i};

    current_color = rand(1, 3);
    trial_num = size(rasters, 2);
    pointer = 1;
    for trial_idx = 1:trial_num
        raster = rasters{trial_idx};
        raster = raster(sort_idx, :);
        raster = raster(neuron_idx, :);
        
        trial_len = size(raster, 2);
        plot(pointer:(pointer + trial_len - 1), raster + i - 0.5, 'Color', current_color);
        plotted_rasters{i} = [plotted_rasters{i}, raster];
        % vertical black line between trials
        plot([pointer + trial_len - 0.5, pointer + trial_len - 0.5], [i-0.5, i+0.5], 'k--', 'LineWidth', 1);
        pointer = pointer + trial_len;
    end
end
% reverse y-axis
set(gca, 'YDir','reverse');
xlabel('Time (ms)');
ylabel('Neuron No.');
hold off;

smooth_kernel = ones(1, 50) / 50;
% smooth, then plot pearson correlation matrix
corr_mat = zeros(neuron_num, neuron_num);
for i = 1:neuron_num
    raster_i = conv(plotted_rasters{i}, smooth_kernel, 'same');
    for j = 1:neuron_num
        raster_j = conv(plotted_rasters{j}, smooth_kernel, 'same');
        corr_mat(i, j) = corr(raster_i', raster_j');
    end
end


figure;
hold on;
imagesc(corr_mat);
colorbar;
title('Pearson Correlation Matrix of Smoothed Rasters');
xlabel('Neuron Index');
ylabel('Neuron Index');
hold off;