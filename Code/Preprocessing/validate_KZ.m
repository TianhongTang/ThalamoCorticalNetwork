%% View and validate dataset from Kaining.

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
% load session metadata
session_data_folder = fullfile(root, 'Data', 'Experimental', '15SecondsThreshold');
session_data_path = fullfile(session_data_folder, 'StateInfoSummary.mat');
load(session_data_path, 'savestruct');
session_num = length(savestruct);
fprintf('Total number of sessions: %d\n', session_num);

% load brain area names
segment_data_path = fullfile(root, 'Data', 'Experimental', 'AnatomyV3', 'segmentNames.mat');
load(segment_data_path, 'segmentNames'); % segmentNames
area_num = length(segmentNames);
total_area_counts = zeros(area_num, 1);
total_area_session_counts = zeros(area_num, 1);

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

for session_idx = 1:session_num
    fprintf('-------------------------\n');
    session_info = savestruct(session_idx);
    monkey_name = session_info.monkeyName;
    session_name = session_info.sessionname;
    session_length = session_info.sessionlength;
    fprintf('Session(%d/%d): %s, Monkey: %s, Length: %d ms\n', session_idx, session_num, session_name, monkey_name, session_length);

    % Load state data
    state_data_folder = fullfile(session_data_folder, 'StateInfo', monkey_name);
    state_data_name = sprintf('%s.mat', session_name);
    state_data_path = fullfile(state_data_folder, state_data_name);
    load(state_data_path, 'stateInfo'); % stateInfo.stateWindows.range[Resting/Arousal/Lowsaccade]
    resting_ranges = stateInfo.stateWindows.rangeResting;
    resting_durations = resting_ranges(:, 2) - resting_ranges(:, 1) + 1;
    arousal_ranges = stateInfo.stateWindows.rangeArousal;
    arousal_durations = arousal_ranges(:, 2) - arousal_ranges(:, 1) + 1;
    lowsaccade_ranges = stateInfo.stateWindows.rangeLowsaccade;
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
    anatomy_data_folder = fullfile(root, 'Data', 'Experimental', 'AnatomyV3', monkey_name);
    anatomy_data_name = sprintf('%s_tasksession_mega_file.mat', session_name);
    anatomy_data_path = fullfile(anatomy_data_folder, anatomy_data_name);
    load(anatomy_data_path, 'Neuronlist');

    % Load sorted cell data
    spike_data_folder = fullfile(root, 'Data', 'Experimental', 'SortedData', monkey_name);
    spike_data_name = sprintf('%s_tasksession_mega_file.mat', session_name);
    spike_data_path = fullfile(spike_data_folder, spike_data_name);
    spike_data = load(spike_data_path); 
    unit_num = length(spike_data.Neuronlist); % SPK and MUA

    % Process each unit
    SPK_count = 0;
    MUA_count = 0;
    unknown_count = 0;
    session_max_spike_time = 0;
    anatomy_count = 0;
    area_counts = zeros(area_num, 1);

    for unit_idx = 1:unit_num
        unit_name = spike_data.Neuronlist(unit_idx).name;
        if strcmp(unit_name(1:3), 'SPK')
            SPK_count = SPK_count + 1;
        elseif strcmp(unit_name(1:3), 'MUA')
            MUA_count = MUA_count + 1;
            continue; % skip MUA for anatomy stats
        else
            unknown_count = unknown_count + 1;
        end

        % spikes
        unit_data = spike_data.(unit_name);
        spike_times = unit_data.sptimes{1}; % Spike times in seconds
        unit_max_spike_time = max(spike_times);
        session_max_spike_time = max(session_max_spike_time, unit_max_spike_time);

        % anatomy info
        anatomy_unit_idx = find(strcmp({Neuronlist.name}, unit_name), 1);
        if ~isempty(anatomy_unit_idx)
            anatomy_count = anatomy_count + 1;
            anatomy_info = Neuronlist(anatomy_unit_idx);
            electrode = anatomy_info.electrodeID;
            depth = anatomy_info.electrodeDepth;
            area = anatomy_info.NeuralTargetsAnatomy;
            area_idx = find(strcmp(segmentNames, area), 1);
            if ~isempty(area_idx)
                area_counts(area_idx) = area_counts(area_idx) + 1;
            else
                fprintf('  Warning: Unit %s has unknown area: %s\n', unit_name, area);
            end

        else
            electrode = -1;
            depth = -1;
            area = 'Unknown';
        end
    end
    non_empty_areas = segmentNames(area_counts > 0);
    non_empty_counts = area_counts(area_counts > 0);
    total_neuron_count = total_neuron_count + SPK_count; % only count SPK
    max_neurons_per_session = max(max_neurons_per_session, SPK_count);
    min_neurons_per_session = min(min_neurons_per_session, SPK_count);
    total_area_counts = total_area_counts + area_counts;
    total_area_session_counts(area_counts > 0) = total_area_session_counts(area_counts > 0) + 1;

    fprintf('Session(%d/%d): %s, Monkey: %s. \n', session_idx, session_num, session_name, monkey_name);
    fprintf('Total Units: %d, SPK: %d, MUA: %d, Unknown: %d. \n', ...
        unit_num, SPK_count, MUA_count, unknown_count);
    fprintf('Max Spike Time: %.4f s. \n', session_max_spike_time);
    fprintf('Non-empty Areas:%d, mean units per area: %.2f, max units per area: %.2f\n', ...
        length(non_empty_areas), mean(non_empty_counts), max(non_empty_counts));
    fprintf('Areas: ');
    for area_idx = 1:length(non_empty_areas)
        fprintf('%s(%d), ', non_empty_areas{area_idx}, non_empty_counts(area_idx));
    end
    fprintf('\n');
end

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
    fprintf('  %s: %d\n', segmentNames{area_idx}, total_area_counts(area_idx));
end