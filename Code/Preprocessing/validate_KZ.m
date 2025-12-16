%% View and validate dataset from Kaining.
% Uses preprocessed data from load_raw_KZ.m

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
tic;
% load session metadata
meta_folder = fullfile(root, 'Data', 'Working', 'Meta');
check_path(meta_folder);
meta_file_name = 'all_session_info_KZ.mat';
meta_file_path = fullfile(meta_folder, meta_file_name);
load(meta_file_path, 'all_session_info', 'segmentNames');

session_num = length(all_session_info);
area_num = length(segmentNames);
fprintf('Total number of sessions: %d\n', session_num);

% initialize stats
arousal_total_count = 0;
resting_total_count = 0;
lowsaccade_total_count = 0;
arousal_total_duration = 0;
resting_total_duration = 0;
lowsaccade_total_duration = 0;
total_aligned_duration = 0;

total_neuron_count = 0;
max_neurons_per_session = 0;
min_neurons_per_session = Inf;
total_area_counts = zeros(area_num, 1);
total_area_session_counts = zeros(area_num, 1);

% special stats
special.thalamus_sessions = {};
thalamus_idx = find(strcmp(segmentNames, 'Thalamus'), 1);
fprintf('Thalamus area index: %d\n', thalamus_idx);

for session_idx = 1:session_num
    fprintf('-------------------------\n');
    session_info = all_session_info(session_idx);
    monkey_name = session_info.monkeyName;
    session_name = session_info.sessionname;
    session_length = session_info.sessionlength;
    fprintf('Session(%d/%d): %s, Monkey: %s, Length: %d ms\n', session_idx, session_num, session_name, monkey_name, session_length);

    % Load state data
    stateWindows = all_session_info(session_idx).stateWindows;
    resting_ranges = stateWindows.rangeResting;
    resting_durations = resting_ranges(:, 2) - resting_ranges(:, 1) + 1;
    arousal_ranges = stateWindows.rangeArousal;
    arousal_durations = arousal_ranges(:, 2) - arousal_ranges(:, 1) + 1;
    lowsaccade_ranges = stateWindows.rangeLowsaccade;
    lowsaccade_durations = lowsaccade_ranges(:, 2) - lowsaccade_ranges(:, 1) + 1;
    fprintf('  Arousal: %d, total dur: %d, mean: %.2f, min: %d, max: %d\n', ...
        size(arousal_ranges, 1), sum(arousal_durations), mean(arousal_durations), ...
        min([arousal_durations; Inf]), max(arousal_durations));
    fprintf('  Resting: %d, total dur: %d, mean: %.2f, min: %d, max: %d\n', ...
        size(resting_ranges, 1), sum(resting_durations), mean(resting_durations), ...
        min([resting_durations; Inf]), max(resting_durations));
    fprintf('  Lowsaccade: %d, total dur: %d, mean: %.2f, min: %d, max: %d\n', ...
        size(lowsaccade_ranges, 1), sum(lowsaccade_durations), mean(lowsaccade_durations), ...
        min([lowsaccade_durations; Inf]), max(lowsaccade_durations));
    arousal_total_count = arousal_total_count + size(arousal_ranges, 1);
    resting_total_count = resting_total_count + size(resting_ranges, 1);
    lowsaccade_total_count = lowsaccade_total_count + size(lowsaccade_ranges, 1);
    session_arousal_duration = sum(arousal_durations(~isnan(arousal_durations)));
    session_resting_duration = sum(resting_durations(~isnan(resting_durations)));
    session_lowsaccade_duration = sum(lowsaccade_durations(~isnan(lowsaccade_durations)));
    arousal_total_duration = arousal_total_duration + session_arousal_duration;
    resting_total_duration = resting_total_duration + session_resting_duration;
    lowsaccade_total_duration = lowsaccade_total_duration + session_lowsaccade_duration;
    aligned_duration = min(session_arousal_duration, session_resting_duration);
    total_aligned_duration = total_aligned_duration + aligned_duration;
    fprintf('  Aligned Duration (min of Arousal and Resting): %d ms\n', aligned_duration);

    % Load anatomy data
    neuron_list = all_session_info(session_idx).neuronList;
    neuron_num = all_session_info(session_idx).neuronNum;

    % Load sorted cell data
    spike_folder = fullfile(root, 'Data', 'Working', 'Spikes');
    check_path(spike_folder);
    spike_file_name = sprintf('spikes_%s.mat', session_name);
    spike_file_path = fullfile(spike_folder, spike_file_name);
    load(spike_file_path, 'spikes', 'neuron_info', 'N', 'state_windows', 'session_length', 'session_name');

    % Process each unit
    SPK_count = 0;
    MUA_count = 0;
    unknown_count = 0;
    session_max_spike_time = 0;
    session_min_spike_time = Inf;
    anatomy_count = 0;
    area_counts = zeros(area_num, 1);

    for neuron_idx = 1:neuron_num
        unit_name = neuron_info(neuron_idx).name;
        if strcmp(unit_name(1:3), 'SPK')
            SPK_count = SPK_count + 1;
        elseif strcmp(unit_name(1:3), 'MUA')
            MUA_count = MUA_count + 1;
            continue; % skip MUA for anatomy stats
        else
            unknown_count = unknown_count + 1;
        end

        % area count
        neuron_area = neuron_list(neuron_idx).NeuralTargetsAnatomy;
        area_idx = find(strcmp(segmentNames, neuron_area), 1);
        if ~isempty(area_idx)
            area_counts(area_idx) = area_counts(area_idx) + 1;
            anatomy_count = anatomy_count + 1;
        else
            fprintf('  Warning: Neuron %s has unknown area: %s\n', unit_name, neuron_area);
        end

        % spikes
        spike_times = spikes{neuron_idx};
        neuron_max_spike_time = max(spike_times);
        neuron_min_spike_time = min(spike_times);
        session_max_spike_time = max(session_max_spike_time, neuron_max_spike_time);
        session_min_spike_time = min(session_min_spike_time, neuron_min_spike_time);

    end
    non_empty_areas = segmentNames(area_counts > 0);
    non_empty_counts = area_counts(area_counts > 0);
    total_neuron_count = total_neuron_count + SPK_count; % only count SPK
    max_neurons_per_session = max(max_neurons_per_session, SPK_count);
    min_neurons_per_session = min(min_neurons_per_session, SPK_count);
    total_area_counts = total_area_counts + area_counts;
    total_area_session_counts(area_counts > 0) = total_area_session_counts(area_counts > 0) + 1; 

    % Special stats
    if area_counts(thalamus_idx) >= 5
        special.thalamus_sessions{end+1} = session_name; 
    end

    fprintf('Session(%d/%d): %s, Monkey: %s. \n', session_idx, session_num, session_name, monkey_name);
    fprintf('Total Neurons: %d, SPK: %d, MUA: %d, Unknown: %d. \n', ...
        neuron_num, SPK_count, MUA_count, unknown_count);
    fprintf('Min Spike Time: %.4f s, Max Spike Time: %.4f s\n', ...
        session_min_spike_time, session_max_spike_time);
    fprintf('Non-empty Areas:%d, mean units per area: %.2f, max units per area: %.2f\n', ...
        length(non_empty_areas), mean(non_empty_counts), max(non_empty_counts));
    fprintf('Areas: ');
    for area_idx = 1:length(non_empty_areas)
        fprintf('%s(%d), ', non_empty_areas{area_idx}, non_empty_counts(area_idx));
    end
    fprintf('\n');
end
toc;

fprintf('-------------------------\n');
fprintf('Summary across all sessions:\n');
fprintf('Session num: %d\n', session_num);

fprintf('Arousal total count: %d\n', arousal_total_count);
fprintf('Arousal total duration: %d ms\n', arousal_total_duration);
fprintf('Resting total count: %d\n', resting_total_count);
fprintf('Resting total duration: %d ms\n', resting_total_duration);
fprintf('Lowsaccade total count: %d\n', lowsaccade_total_count);
fprintf('Lowsaccade total duration: %d ms\n', lowsaccade_total_duration);
fprintf('Total Aligned Duration across all sessions: %d ms\n', total_aligned_duration);

fprintf('Total SPK Neuron Count across all sessions: %d\n', total_neuron_count);
fprintf('Max Neurons in a session: %d\n', max_neurons_per_session);
fprintf('Min Neurons in a session: %d\n', min_neurons_per_session);
fprintf('Total Neuron Count by Area across all sessions:\n');
for area_idx = 1:area_num
    fprintf('  %s: %d sessions, %d neurons\n', segmentNames{area_idx}, total_area_session_counts(area_idx), total_area_counts(area_idx));
end

fprintf('Sessions with at least 5 neurons in Thalamus area:\n');
for i = 1:length(special.thalamus_sessions)
    fprintf('  %d:%s\n', i, special.thalamus_sessions{i});
end