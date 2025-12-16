%% basic_statistics_main.m -  Basic statistics of each session/state: Synchrony index, firing rate, fano factor, etc.


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

% Load session metadata
meta_folder = fullfile(root, 'Data', 'Working', 'Meta');
check_path(meta_folder);
meta_file_name = 'all_session_info_KZ.mat';
meta_file_path = fullfile(meta_folder, meta_file_name);
load(meta_file_path, 'all_session_info', 'segmentNames');
session_num = size(all_session_info, 2);

myGaussian = @(x, mu, sigma) exp(-((x - mu).^2) / (2*sigma^2)) / (sigma*sqrt(2*pi));
smooth_kernel = myGaussian(-200:200, 0, 50); 
smooth_kernel = smooth_kernel / sum(smooth_kernel);

states = {'RestOpen', 'RestClose'};
state_names = {'Eye Open', 'Eye Closed'};
n_states = length(states);
fano_factor_all = NaN(session_num, n_states); % (session, state)
sync_chi_all = NaN(session_num, n_states);

for session_idx = 1:session_num
    session_info = all_session_info(session_idx);
    fprintf('Processing session %d: %s\n', session_idx, session_info.sessionname);
    neuron_info = session_info.neuronList;
    thalamus_filter = cellfun(@(x) strcmp(x, 'Thalamus'), {neuron_info.NeuralTargetsAnatomy});
    thalamus_count = sum(thalamus_filter);
    if thalamus_count < 5
        fprintf('Session %d: Not enough thalamus neurons(%d). Skipping...\n', session_idx, thalamus_count);
        continue;
    else
        fprintf('Session %d: %d thalamus neurons found.\n', session_idx, thalamus_count);
    end
    
    for state_idx = 1:n_states
        state = states{state_idx};
        session_type = ['KZ', state];

        % load raster data
        raster_folder = fullfile(root, 'Data', 'Working', 'raster');
        raster_name = sprintf('raster_%s_%d.mat', session_type, session_idx);
        raster_path = fullfile(raster_folder, raster_name);
        load(raster_path, 'rasters', 'N', 'trial_num');
        if trial_num == 0
            fprintf('Session %d, State %s: No trials found. Skipping...\n', session_idx, state);
            continue;
        else
            fprintf('Session %d, State %s: %d trials found.\n', session_idx, state, trial_num);
        end

        % concatenate trials
        smoothed_rasters = cellfun(@(x) conv2(x, smooth_kernel, 'valid'), rasters, 'UniformOutput', false);
        total_length = sum(cellfun(@(x) size(x, 2), rasters));
        total_length_smooth = sum(cellfun(@(x) size(x, 2), smoothed_rasters));
        raster_concat = zeros(N, total_length);
        raster_smooth_concat = zeros(N, total_length_smooth);
        pointer = 1;
        pointer_smooth = 1;
        for t = 1:trial_num
            t_len = size(rasters{t}, 2);
            raster_concat(:, pointer:(pointer + t_len - 1)) = rasters{t};
            pointer = pointer + t_len;

            t_len_smooth = size(smoothed_rasters{t}, 2);
            raster_smooth_concat(:, pointer_smooth:(pointer_smooth + t_len_smooth - 1)) = smoothed_rasters{t};
            pointer_smooth = pointer_smooth + t_len_smooth;
        end

        % use thalamus neurons only
        raster_concat = raster_concat(thalamus_filter, :);

        % calculate fano factor and synchrony chi
        fano = mean(fano_factor(raster_concat), 'omitnan');
        chi = calc_synchrony(raster_smooth_concat);
        fano_factor_all(session_idx, state_idx) = fano;
        sync_chi_all(session_idx, state_idx) = chi;
        fprintf('Fano Factor = %.4f, Synchrony Chi = %.4f\n', fano, chi);
    end
end

fprintf('Summary of statistics:\n');
for state_idx = 1:n_states
    valid_session_num = sum(~isnan(fano_factor_all(:, state_idx)));
    ave_fano = mean(fano_factor_all(:, state_idx), 'omitnan');
    se_fano = std(fano_factor_all(:, state_idx), 0, 'omitnan') / sqrt(valid_session_num);
    valid_session_num = sum(~isnan(sync_chi_all(:, state_idx)));
    ave_chi = mean(sync_chi_all(:, state_idx), 'omitnan');
    se_chi = std(sync_chi_all(:, state_idx), 0, 'omitnan') / sqrt(valid_session_num);
    fprintf('State: %s\n', state_names{state_idx});
    fprintf('Average Fano Factor: %.4f ± %.4f\n', ave_fano, se_fano);
    fprintf('Average Synchrony Chi: %.4f ± %.4f\n', ave_chi, se_chi);
end

save_folder = fullfile(root, 'Data', 'Working', 'basic');
check_path(save_folder);
save(fullfile(save_folder, 'basic_statistics.mat'), 'fano_factor_all', 'sync_chi_all');