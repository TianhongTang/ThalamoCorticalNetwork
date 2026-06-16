%% check_spikes_complete.m - Visualization of rasters

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
% mode = 'spikes';
% mode = 'raster';
REPLOT = false;

segment_len_ms = 10000; % Each output raster figure shows this many ms.
% load metadata
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums');

% select states to plot
states = {'RestOpen', 'RestClose'};
prepost = {'Pre', 'Post'};
area_types = {'Full', 'Cortex'};
aligns = {'None', 'Last', 'Longest'};

% main loop
for dataset_idx = 1:3
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
                        % Output folder. Segment-level existence is checked after loading data,
                        % because the number of segments depends on the concatenated raster length.
                        figure_folder = fullfile(root, 'Figures', 'Rasters_PDS');
                        check_path(figure_folder);

                        % construct meta for loading
                        meta = struct();

                        if contains(dataset_name, 'Slayer')
                            meta.animal_name = 'Slayer';
                        elseif contains(dataset_name, 'Zeppelin')
                            meta.animal_name = 'Zeppelin';
                        elseif contains(dataset_name, 'Emperor')
                            meta.animal_name = 'Emperor';
                        end
                        if contains(dataset_name, 'Mus')
                            meta.injection = 'Muscimol';
                        elseif contains(dataset_name, 'Sal')
                            meta.injection = 'Saline';
                        elseif contains(dataset_name, 'Noinj')
                            meta.injection = 'No injection';
                        end
                        meta.prepost     = prepost_str;
                        meta.state       = state;
                        meta.area        = area_type;
                        meta.align       = align;
                        meta.session_idx = session_idx;
                        meta.file_name   = generate_filename('raster', meta);

                        % load data
                        data_folder = fullfile(root, 'Data', 'Working', 'raster');
                        data_path = fullfile(data_folder, meta.file_name);
                        if ~isfile(data_path)
                            fprintf('Data file not found: %s\n', data_path);
                            continue;
                        end
                        d = load(data_path);
                        rasters = d.data.rasters;
                        cell_area = d.data.cell_area;
                        % channel = d.data.channel;
                        N = d.meta.N;

                        sortidx_folder = fullfile(root, 'Data', 'Working', 'sortidx');
                        meta.criterion = 'channel';
                        meta.file_name   = generate_filename('sortidx', meta);
                        sortidx_path = fullfile(sortidx_folder, meta.file_name);
                        if ~isfile(sortidx_path)
                            fprintf('Sort index file not found: %s\n', sortidx_path);
                            continue;
                        end
                        sort_idx = load(sortidx_path).data.sort_idx;

                        % apply sorting
                        cell_area = cell_area(sort_idx);
                        for r_idx = 1:length(rasters)
                            rasters{r_idx} = rasters{r_idx}(sort_idx, :);
                        end

                        % assign colors based on cell area
                        plot_colors = zeros(numel(cell_area), 3);
                        for i = 1:numel(cell_area)
                            switch cell_area{i}
                                case 'Thalamus'
                                    plot_colors(i, :) = [1, 0, 1];
                                case 'ACC'
                                    plot_colors(i, :) = [0, 0, 1]; % blue
                                case 'VLPFC'
                                    plot_colors(i, :) = [1, 0, 0]; % red
                                otherwise
                                    plot_colors(i, :) = [0, 0, 0]; % black
                            end
                        end

                        % Smooth kernel
                        myGaussian = @(x, mu, sigma) exp(-((x - mu).^2) / (2*sigma^2)) / (sigma*sqrt(2*pi));
                        smooth_kernel = myGaussian(-200:200, 0, 50);
                        smoothed_rasters = cellfun(@(x) conv2(x, smooth_kernel, 'valid'), rasters, 'UniformOutput', false);

                        % load session_info
                        session_info_folder = fullfile(root, 'Data', 'Working', 'Meta');
                        data_name = sprintf('all_session_info_KZ.mat');
                        data_path = fullfile(session_info_folder, data_name);
                        d = load(data_path, 'all_session_info');
                        all_session_info = d.all_session_info;
                        session_info = all_session_info(session_idx);
                        neuron_info = session_info.neuronList;
                        thal_filter = cellfun(@(x) strcmp(x, 'Thalamus'), {neuron_info.NeuralTargetsAnatomy});
                        %% concatenate rasters
                        concatenated_rasters = cell2mat(rasters);
                        concatenated_smoothed = cell2mat(smoothed_rasters); %#ok<NASGU>
                        trial_borders = cumsum(cellfun(@(x) size(x, 2), rasters));

                        %% plot raster in fixed-length segments
                        % Instead of compressing the full concatenated raster into one figure,
                        % split it into 10000 ms windows. The last segment is allowed to be
                        % shorter, but xlim is kept at [1, segment_len_ms] to preserve scale.
                        total_len = size(concatenated_rasters, 2);
                        segment_num = ceil(total_len / segment_len_ms);

                        for segment_idx = 1:segment_num
                            segment_start = (segment_idx - 1) * segment_len_ms + 1;
                            segment_end = min(segment_idx * segment_len_ms, total_len);
                            segment_range = segment_start:segment_end;
                            segment_raster = concatenated_rasters(:, segment_range);

                            segment_trial_borders = trial_borders( ...
                                trial_borders >= segment_start & trial_borders <= segment_end) ...
                                - segment_start + 1;

                            figure_name = sprintf('raster_%s%s%s%s%s_%d_part%03d_%06d-%06dms.png', ...
                                dataset_name, prepost_str, state, area_type, align, session_idx, ...
                                segment_idx, segment_start, segment_end);
                            figure_path = fullfile(figure_folder, figure_name);

                            if isfile(figure_path) && ~REPLOT
                                fprintf('Figure already exists: %s', figure_path);
                                continue;
                            end

                            f = raster_visualization(segment_raster, plot_colors, segment_trial_borders, 'off');
                            set_raster_xaxis_to_fixed_window(f, segment_len_ms);

                            sgtitle(sprintf('%s %s %s %s %s session %d | part %d/%d | %d-%d ms', ...
                                dataset_name, prepost_str, state, area_type, align, session_idx, ...
                                segment_idx, segment_num, segment_start, segment_end), ...
                                'Interpreter', 'none');

                            saveas(f, figure_path);
                            close(f);
                        end

                    end
                end
            end
        end
    end
end


function set_raster_xaxis_to_fixed_window(fig_handle, segment_len_ms)
    axes_handles = findall(fig_handle, 'Type', 'axes');
    for ax_idx = 1:numel(axes_handles)
        ax = axes_handles(ax_idx);
        try
            xlim(ax, [1, segment_len_ms]);
        catch
            % Some non-data axes may not support xlim. Ignore them.
        end
    end
end
