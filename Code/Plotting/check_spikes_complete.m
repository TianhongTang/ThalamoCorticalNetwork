%% check_spikes_complete.m - Complete visualization of spikes/rasters, plot and save as files.

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
state_colors.Resting = [0.25, 0.25, 1]; % blue for eye closed
state_colors.Arousal = [1, 0, 0]; % red for arousal
state_colors.Lowsaccade = [0, 1, 0]; % green for low saccade
state_colors.Default = [0, 0, 0]; % black for default
state_colors.EyeOpen = [1, 0.25, 1]; % magenta for eye open
window_size = 20000; % 20s
dt = 0.001; % 1 ms
mode = 'session';
% mode = 'state';

% load session metadata
meta_folder = fullfile(root, 'Data', 'Working', 'Meta');
meta_file_name = 'all_session_info_KZ.mat';
meta_file_path = fullfile(meta_folder, meta_file_name);
load(meta_file_path, 'all_session_info', 'segmentNames');
session_num = size(all_session_info, 2);

switch mode
case 'session'
    parfor session_idx = 1:session_num
        % Load spikes data
        session_name = all_session_info(session_idx).sessionname;
        fprintf('Processing session %d: %s\n', session_idx, session_name);

        spikes_folder = fullfile(root, 'Data', 'Working', 'Spikes');
        spikes_file_name = sprintf('spikes_%s.mat', session_name);
        spikes_file_path = fullfile(spikes_folder, spikes_file_name);
        loaded = load(spikes_file_path, 'all_spikes', 'neuron_info', 'N', 'state_windows', 'session_length', 'session_name');
        all_spikes = loaded.all_spikes;
        neuron_info = loaded.neuron_info;
        N = loaded.N;
        state_windows = loaded.state_windows;
        session_length = loaded.session_length;
        session_name = loaded.session_name;

        %% sort state windows: merge all state windows into a single timeline
        states = {'Resting', 'Arousal', 'Lowsaccade'};

        % handle window overlap
        for state_i = 1:length(states)
            state = states{state_i};
            for state_j = (state_i+1):length(states)
                other_state = states{state_j};
                state_ranges = state_windows.(['range', state]);
                other_state_ranges = state_windows.(['range', other_state]);
                for i = 1:size(state_ranges, 1)
                    state_range = state_ranges(i, :);
                    for j = 1:size(other_state_ranges, 1)
                        other_state_range = other_state_ranges(j, :);
                        % check overlap. If overlap, remove the overlapping part from the other state range
                        include_head = state_range(1) <= other_state_range(1) && other_state_range(1) <= state_range(2);
                        include_tail = state_range(1) <= other_state_range(2) && other_state_range(2) <= state_range(2);
                        if include_head && include_tail
                            % other state range fully included in state range: remove other state range
                            other_state_ranges(j, :) = [-1, -1]; % mark for removal
                        elseif include_head
                            % overlap at head: adjust other state range start
                            other_state_ranges(j, 1) = state_range(2) + 1;  
                        elseif include_tail
                            % overlap at tail: adjust other state range end
                            other_state_ranges(j, 2) = state_range(1) - 1;  
                        end
                    end
                end
                filter = other_state_ranges(:,1) ~= -1;
                other_state_ranges = other_state_ranges(filter, :);
                state_windows.(['range', other_state]) = other_state_ranges;
            end
        end

        % print windows after overlap handling
        fprintf('State windows after overlap handling for session %d:\n', session_idx);
        for state_i = 1:length(states)
            state = states{state_i};
            state_ranges = state_windows.(['range', state]);
            state_ranges = state_ranges(state_ranges(:,1) ~= -1, :);
            state_windows.(['range', state]) = state_ranges;
            fprintf('  State: %s\n', state);
            for i = 1:size(state_ranges, 1)
                fprintf('    Window: [%d, %d]\n', state_ranges(i,1), state_ranges(i,2));
            end
        end

        % states = {'Resting', 'Arousal'};
        state_num = length(states);
        current_pointer = 0; % time pointer in ms
        state_pointers = ones(state_num, 1); % state window idx for each state
        sorted_states = {};
        sorted_state_windows = [];
        while current_pointer < session_length
            % find next state window
            next_state_idx = -1;
            next_start = session_length + 1;
            for state_idx = 1:state_num
                if state_pointers(state_idx) > size(state_windows.(['range', states{state_idx}]), 1)
                    % no more windows for this state
                    continue;
                end
                state = states{state_idx};
                state_ranges = state_windows.(['range', state]);
                state_range = state_ranges(state_pointers(state_idx), :);
                if state_range(1) < next_start
                    next_start = state_range(1);
                    next_state_idx = state_idx;
                end
            end

            % no more state windows: add the final default window and break
            if next_state_idx == -1
                if current_pointer < session_length
                    sorted_states{end+1} = 'Default'; %#ok<SAGROW>
                    sorted_state_windows = [sorted_state_windows; current_pointer, session_length - 1]; 
                end
                break;
            end

            next_state = states{next_state_idx};

            % check for overlap
            if next_start < current_pointer
                error('Window overlap detected: session %d, current pointer %d, state %s start %d', session_idx, current_pointer, next_state, next_start);
            end

            % add default state window if any gap
            if next_start > current_pointer
                sorted_states{end+1} = 'Default'; %#ok<SAGROW>
                sorted_state_windows = [sorted_state_windows; current_pointer, next_start - 1]; 
                current_pointer = next_start;
            end

            % add next state window
            if next_state_idx ~= -1
                state_ranges = state_windows.(['range', next_state]);
                state_range = state_ranges(state_pointers(strcmp(states, next_state)), :); 
                sorted_states{end+1} = next_state; %#ok<SAGROW>
                sorted_state_windows = [sorted_state_windows; state_range]; 
                current_pointer = state_range(2) + 1;
                state_pointers(strcmp(states, next_state)) = state_pointers(strcmp(states, next_state)) + 1;
            end
        end

        % print sorted states for verification
        fprintf('Sorted state windows for session %d:\n', session_idx);
        for i = 1:length(sorted_states)
            fprintf('  State: %s, Window: [%d, %d]\n', sorted_states{i}, sorted_state_windows(i,1), sorted_state_windows(i,2));
        end
        drawnow('update');

        %% rasterize spikes
        raster_edges = 0:dt:(session_length*dt);
        raster = zeros(N, session_length);

        for neuron_idx = 1:N
            spikes = all_spikes{neuron_idx}; % in ms
            raster(neuron_idx, :) = histcounts(spikes, raster_edges);
        end
        raster(raster > 1) = 1; % binarize

        % plot raster
        window_num = ceil(session_length / window_size);
        state_pointer = 1; % pointer for sorted states

        for window_idx = 1:window_num
            f = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
            window_start = (window_idx - 1) * window_size + 1;
            window_end = min(window_idx * window_size, session_length);
            hold on;

            % plot states in the window
            while sorted_state_windows(state_pointer, 1) <= window_end
                state_start = max(sorted_state_windows(state_pointer, 1), window_start);
                state_end = min(sorted_state_windows(state_pointer, 2), window_end);
                state = sorted_states{state_pointer};
                state_color = state_colors.(state);
                for i = 1:N
                    plot(raster(i, state_start:state_end) + i - 0.5, 'Color', state_color, 'LineWidth', 1);
                end
                % if fully covered the state window, continue to next state
                if sorted_state_windows(state_pointer, 2) <= window_end
                    state_pointer = state_pointer + 1;
                end
            end
            hold off;
            xlabel("Time (ms)");
            ylabel("Neuron No.")
            xlim([1, window_size]);
            title(sprintf('Session %s, Window %d (%d-%d ms)', ...
                session_name, window_idx, window_start, window_end));
            % save figure
            save_folder = fullfile(root, 'Figures', 'rasters_overview_session');
            check_path(save_folder);
            save_name = sprintf('raster_session%d_window%d.png', session_idx, window_idx);
            saveas(f, fullfile(save_folder, save_name));
            close(f);
        end
    end

case 'state'
    parfor session_idx = 1:session_num
        fprintf('Processing session %d: %s\n', session_idx, all_session_info(session_idx).sessionname);
        neuron_info = all_session_info(session_idx).neuronList;
        areas = {neuron_info.NeuralTargetsAnatomy};
        % thalamus_filter = cellfun(@(x) strcmp(x, 'Thalamus'), {neuron_info.NeuralTargetsAnatomy});

        for state_idx = 1:length(states)
            state = states{state_idx};
            % load data
            data_folder = fullfile(root, 'Data', 'Working', 'raster');
            data_name = sprintf('raster_KZ%s_%d.mat', state, session_idx);
            data_path = fullfile(data_folder, data_name);
            rasters = load(data_path).rasters;
            trial_num = length(rasters);
            if trial_num == 0
                fprintf('Session %d, State %s: No trials found. Skipping...\n', session_idx, state);
                continue;
            end

            for trial_idx = 1:trial_num
                raster = rasters{trial_idx};
                [N, T] = size(raster);
                window_num = ceil(T / window_size);
                for window_idx = 1:window_num
                    start_t = (window_idx - 1) * window_size + 1;
                    end_t = min(window_idx * window_size, T);
                    raster_window = raster(:, start_t:end_t);
                    % plot raster
                    f = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
                    hold on;
                    for i = 1:N
                        if strcmp(state, 'RestOpen')
                            line_color = [1, 0.25, 1]; % magenta for eye open
                        else
                            line_color = [0.25, 0.25, 1]; % blue for eye closed
                        end

                        if ~strcmp(areas{i}, 'Thalamus')
                            line_color = 1 - ((1-line_color) * 0.3); % light color for non-thalamus neurons
                        end
                        plot(raster_window(i,:) + i - 0.5, 'Color', line_color, 'LineWidth', 1);
                    end
                    hold off;
                    xlabel("Time (ms)");
                    ylabel("Neuron No.")
                    xlim([1, window_size]);
                    title(sprintf('Session %d, State %s, Trial %d. Window %d (%d-%d ms)', ...
                        session_idx, state, trial_idx, window_idx, start_t, end_t));
                    % save figure
                    save_folder = fullfile(root, 'Figures', 'rasters_overview');
                    check_path(save_folder);
                    save_name = sprintf('raster_session%d_%s_trial%d_window%d.png', session_idx, state, trial_idx, window_idx);
                    saveas(f, fullfile(save_folder, save_name));
                    close(f);
                end
            end
        end
        fprintf('Completed session %d.\n', session_idx);
    end
end