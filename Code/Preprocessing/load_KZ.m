% load preprocessed KZ data (Resting state only)
% extract spike timings and rasters.

%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end

%% Main
dt=0.001;
cutoff = 250; % ms, removing beginning and end of states

% Load session metadata
meta_folder = fullfile(root, 'Data', 'Working', 'Meta');
check_path(meta_folder);
meta_file_name = 'all_session_info_KZ.mat';
meta_file_path = fullfile(meta_folder, meta_file_name);
load(meta_file_path, 'all_session_info', 'segmentNames');

session_num = length(all_session_info);
area_num = length(segmentNames);
fprintf('Total number of sessions: %d\n', session_num);

valid_session_count = 0;
for session_idx = 1:session_num
    fprintf('-------------------------\n');
    session_info = all_session_info(session_idx);
    monkey_name = session_info.monkeyName;
    session_name = session_info.sessionname;
    session_length = session_info.sessionlength;
    fprintf('Session(%d/%d): %s, Monkey: %s, Length: %d ms\n', session_idx, session_num, session_name, monkey_name, session_length);

    % Load spikes data
    spikes_folder = fullfile(root, 'Data', 'Working', 'Spikes');
    spikes_file_name = sprintf('spikes_%s.mat', session_name);
    spikes_file_path = fullfile(spikes_folder, spikes_file_name);
    load(spikes_file_path, 'all_spikes', 'neuron_info', 'N', 'state_windows', 'session_length', 'session_name');
    assert(strcmp(session_name, session_info.sessionname), 'Session name mismatch!');
    assert(session_length == session_info.sessionlength, 'Session length mismatch!');
    fprintf('Total Neurons: %d. \n', N);

    % State windows
    resting_ranges = state_windows.rangeResting;
    resting_durations = resting_ranges(:, 2) - resting_ranges(:, 1) + 1;
    resting_num = size(resting_ranges, 1);
    arousal_ranges = state_windows.rangeArousal;
    arousal_durations = arousal_ranges(:, 2) - arousal_ranges(:, 1) + 1;
    arousal_num = size(arousal_ranges, 1);
    lowsaccade_ranges = state_windows.rangeLowsaccade;
    lowsaccade_durations = lowsaccade_ranges(:, 2) - lowsaccade_ranges(:, 1) + 1;
    lowsaccade_num = size(lowsaccade_ranges, 1);
    non_resting_num = resting_num + 1;
    non_resting_ranges = zeros(non_resting_num, 2);
    if resting_num == 0
        non_resting_ranges(1, :) = [0, session_length - 1];
    else
        % before first resting
        non_resting_ranges(1, :) = [0, resting_ranges(1, 1) - 1];
        for i = 1:resting_num - 1
            non_resting_ranges(i + 1, :) = [resting_ranges(i, 2) + 1, resting_ranges(i + 1, 1) - 1];
        end
        % after last resting
        non_resting_ranges(non_resting_num, :) = [resting_ranges(resting_num, 2) + 1, session_length - 1];
    end
    non_resting_durations = non_resting_ranges(:, 2) - non_resting_ranges(:, 1) + 1;

    % Process spikes
    states = {'RestOpen', 'RestClose'};
    state_num = length(states);
    for state_idx = 1:state_num
        state = states{state_idx};
        switch state
        case 'RestOpen'
            ranges = non_resting_ranges;
            durations = non_resting_durations;
        case 'RestClose'
            ranges = resting_ranges;
            durations = resting_durations;
        end

        % Remove short trials
        original_trial_num = size(ranges, 1);
        trial_filter = durations > 2 * cutoff;
        selected_ranges = ranges(trial_filter, :);
        trial_num = size(selected_ranges, 1);

        fprintf('  State: %s, Valid Trials: %d / %d\n', state, trial_num, original_trial_num);
        for t = 1:original_trial_num
            if trial_filter(t)
                fprintf('    Include trial %d: [%d, %d], duration: %d ms\n', t, ranges(t, 1), ranges(t, 2), durations(t));
            else
                fprintf('    Exclude trial %d: [%d, %d], duration: %d ms\n', t, ranges(t, 1), ranges(t, 2), durations(t));
            end
        end
        spikes = cell(N, trial_num);
        rasters = cell(1, trial_num);
        firing_rates = cell(1, trial_num);
        channel = [neuron_info.electrodeID];
        cell_id = {neuron_info.name};
        trial_len = zeros(1, trial_num);
        for trial_idx=1:trial_num
            trial_start = selected_ranges(trial_idx, 1) + cutoff;
            trial_end = selected_ranges(trial_idx, 2) - cutoff;
            trial_duration = trial_end - trial_start + 1;
            trial_edges = 0:dt:(trial_duration/1000);
            trial_B = length(trial_edges) - 1;
            trial_len(trial_idx) = trial_B;
            rasters{trial_idx} = zeros(N, trial_B);

            for neuron_idx = 1:N
                neuron_spikes = all_spikes{neuron_idx};
                % only keep spikes from trial_start to trial_end,
                % align to trial_start.
                spike_trial = neuron_spikes(neuron_spikes>trial_start/1000 & neuron_spikes<trial_end/1000);
                spike_trial = spike_trial - trial_start/1000;

                % Convert spike times to raster
                spikes{neuron_idx, trial_idx} = spike_trial;
                raster = histcounts(spike_trial, trial_edges);
                raster(raster>1) = 1;
                rasters{trial_idx}(neuron_idx, :) = raster;
                firing_rates{trial_idx}(neuron_idx) = mean(raster);
            end
        end

        session_name_save = ['KZ', state];
        save_folder = fullfile(root, 'Data', 'Working', 'raster');
        save_name = sprintf('raster_%s_%d.mat', session_name_save, session_idx);
        check_path(save_folder);
        save_path = fullfile(save_folder, save_name);
        save(save_path,...
            "rasters", "spikes", "firing_rates", "trial_num", "trial_len", ...
            "session_name", "N", "cell_id", "channel", 'dt');
    end
end