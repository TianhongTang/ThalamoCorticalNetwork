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
REPLOT = false;

metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums');
states = {'RestOpen', 'RestClose'};
prepost = {'Pre', 'Post'};
area_types = {'Full', 'Cortex'};
aligns = { '', 'AlignFirst', 'AlignLast', 'AlignLongest'};
for dataset_idx = 1:dataset_num
    dataset_name = dataset_names{dataset_idx};
    session_num = session_nums(dataset_idx);
    for session_idx = 1:session_num
        fprintf('=========================\n');
        fprintf('Dataset: %s, Session: %d/%d\n', dataset_name, session_idx, session_num);
        for state_idx = 1:length(states)
            state = states{state_idx};
            for prepost_idx = 1:length(prepost)
                prepost_str = prepost{prepost_idx};
                for area_type_idx = 1:length(area_types)
                    area_type = area_types{area_type_idx};
                    for align_idx = 1:length(aligns)
                        align = aligns{align_idx};
                        fprintf('State: %s, %s, %s, Area: %s\n', state, prepost_str, align, area_type);

                        % check if figure exists
                        figure_folder = fullfile(root, 'Figures', 'Rasters_PDS');
                        check_path(figure_folder);
                        figure_name = sprintf('raster_%s%s%s%s%s_%d.png', dataset_name, prepost_str, state, area_type, align, session_idx);
                        figure_path = fullfile(figure_folder, figure_name);
                        
                        if isfile(figure_path) && ~REPLOT
                            fprintf('Figure already exists: %s\n', figure_path);
                            continue;
                        end
                        
                        % load data
                        data_folder = fullfile(root, 'Data', 'Working', 'raster');
                        data_name = sprintf('raster_%s%s%s%s%s_%d.mat', dataset_name, prepost_str, state, area_type, align, session_idx);
                        data_path = fullfile(data_folder, data_name);
                        if ~isfile(data_path)
                            fprintf('Data file not found: %s\n', data_path);
                            continue;
                        end
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
                        f = raster_visualization(concatenated_rasters, trial_borders, 'off');

                        % save figure
                        figure_folder = fullfile(root, 'Figures', 'Rasters_PDS');
                        check_path(figure_folder);
                        figure_name = sprintf('raster_%s%s%s%s%s_%d.png', dataset_name, prepost_str, state, area_type, align, session_idx);
                        figure_path = fullfile(figure_folder, figure_name);
                        saveas(f, figure_path);
                        close(f);
                    end
                end
            end
        end
    end
end


