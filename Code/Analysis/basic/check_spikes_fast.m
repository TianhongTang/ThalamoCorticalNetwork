%% check_spikes_fast.m - Quick visualization of spikes/rasters

clear;
%% Get root folder
code_depth = 4;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
% include code folder and utils
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Main
% mode = 'spikes';
% mode = 'raster';
state = 'SlayerMusPreRestOpenFull';
session_idx = 6;

% load data
data_folder = fullfile(root, 'Data', 'Working', 'raster');
data_name = sprintf('raster_%s_%d.mat', state, session_idx);
data_path = fullfile(data_folder, data_name);
load(data_path, "rasters");

% Smooth kernel
myGaussian = @(x, mu, sigma) exp(-((x - mu).^2) / (2*sigma^2)) / (sigma*sqrt(2*pi));
smooth_kernel = myGaussian(-200:200, 0, 50);
smoothed_rasters = cellfun(@(x) conv2(x, smooth_kernel, 'valid'), rasters, 'UniformOutput', false);

% load session_info
session_info_folder = fullfile(root, 'Data', 'Working', 'Meta');
data_name = sprintf('all_session_info_KZ.mat');
data_path = fullfile(session_info_folder, data_name);
load(data_path, 'all_session_info');
session_info = all_session_info(session_idx);
neuron_info = session_info.neuronList;
thal_filter = cellfun(@(x) strcmp(x, 'Thalamus'), {neuron_info.NeuralTargetsAnatomy});

%% concatenate rasters
concatenated_rasters = cell2mat(rasters);
concatenated_smoothed = cell2mat(smoothed_rasters);
trial_borders = cumsum(cellfun(@(x) size(x, 2), rasters));

%% plot raster
% raster_visualization(rasters{1}(thal_filter, :));
% raster_visualization(rasters{1}(:, :), trial_borders);
raster_visualization(concatenated_rasters, trial_borders);