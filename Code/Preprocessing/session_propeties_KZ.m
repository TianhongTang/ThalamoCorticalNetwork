% load preprocessed KZ data, extract session properties and save.
% general analysis of session features. Save filters to Meta folder.

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

% filters: more than 5 thalamus neurons, have resting state data, (structure)
session_properties = struct();
property_filters = struct();

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
    load(spikes_file_path, 'neuron_info', 'N', 'state_windows', 'session_length', 'session_name');
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

    total_resting_dur = sum(resting_durations);
    total_arousal_dur = sum(arousal_durations);
    total_lowsaccade_dur = sum(lowsaccade_durations);

    if isnan(resting_ranges(1, 1))
        resting_num = 0;
    end

    % save session properties
    thalamus_filter = cellfun(@(x) strcmp(x, 'Thalamus'), {neuron_info.NeuralTargetsAnatomy});
    session_properties(session_idx).thalamus_count = sum(thalamus_filter);
    session_properties(session_idx).resting_count = resting_num;
    session_properties(session_idx).session_length = session_length;
    session_properties(session_idx).resting_dur = total_resting_dur;
    session_properties(session_idx).arousal_dur = total_arousal_dur;
    session_properties(session_idx).lowsaccade_dur = total_lowsaccade_dur;
    session_properties(session_idx).default_dur = session_length - total_resting_dur - total_arousal_dur - total_lowsaccade_dur;
    session_properties(session_idx).eyeClosed_dur = total_resting_dur;
    session_properties(session_idx).eyeOpen_dur = session_length - total_resting_dur;

    % save filters
    property_filters(session_idx).thalamus = session_properties(session_idx).thalamus_count >= 5;
    property_filters(session_idx).resting = session_properties(session_idx).resting_count > 0;
    property_filters(session_idx).combined = property_filters(session_idx).thalamus & property_filters(session_idx).resting;
end

session_filters = repmat(struct('thalamus', false, 'resting', false, 'combined', false), 1, session_num);

% save filters
filter_save_path = fullfile(meta_folder, 'session_filters_KZ.mat');
save(filter_save_path, 'property_filters', 'session_properties');