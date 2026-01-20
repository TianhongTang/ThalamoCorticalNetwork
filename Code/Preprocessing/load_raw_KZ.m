%% Dataset from Kaining. Load raw data and extract useful information.
% raw data is too large.
% for each session:

% After this, processed data:
% Data/Working/Meta/all_session_info_KZ.mat (all_session_info, segmentNames)
% Data/Working/Spikes/spikes_<session_name>.mat (spikes, neuron_info, N, state_windows, session_length, session_name)


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

% save meta session info:
all_session_info = savestruct;

for session_idx = 1:session_num
    % fprintf('-------------------------\n');
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
    all_session_info(session_idx).stateWindows = stateInfo.stateWindows;

    % Load anatomy data
    anatomy_data_folder = fullfile(root, 'Data', 'Experimental', 'AnatomyV3', monkey_name);
    anatomy_data_name = sprintf('%s_tasksession_mega_file.mat', session_name);
    anatomy_data_path = fullfile(anatomy_data_folder, anatomy_data_name);
    load(anatomy_data_path, 'Neuronlist');
    neuron_filter = arrayfun(@(x) strcmp(x.name(1:3), 'SPK'), Neuronlist);
    Neuronlist = Neuronlist(neuron_filter); % only keep SPK
    neuron_num = length(Neuronlist);
    all_session_info(session_idx).neuronList = Neuronlist;
    all_session_info(session_idx).neuronNum = neuron_num;

    % Load sorted cell data
    spike_data_folder = fullfile(root, 'Data', 'Experimental', 'SortedData', monkey_name);
    spike_data_name = sprintf('%s_tasksession_mega_file.mat', session_name);
    spike_data_path = fullfile(spike_data_folder, spike_data_name);
    spike_data = load(spike_data_path); 
    
    all_spikes = cell(neuron_num, 1);
    % Process each unit
    for neuron_idx = 1:neuron_num   
        neuron_name = Neuronlist(neuron_idx).name;

        % spikes
        neuron_spikes = spike_data.(neuron_name);
        spike_times = neuron_spikes.sptimes{1}; % Spike times in seconds
        all_spikes{neuron_idx} = spike_times;
        max_spike_time = max(spike_times);
        % session_max_spike_time = max(session_max_spike_time, max_spike_time);
    end

    % save spikes
    neuron_info = Neuronlist;
    N = neuron_num;
    state_windows = stateInfo.stateWindows;

    spike_folder = fullfile(root, 'Data', 'Working', 'Spikes');
    check_path(spike_folder);
    spike_file_name = sprintf('spikes_%s.mat', session_name);
    spike_file_path = fullfile(spike_folder, spike_file_name);
    save(spike_file_path, 'all_spikes', 'neuron_info', 'N', 'state_windows', 'session_length', 'session_name', '-v7.3');
    fprintf('Saved spikes to %s\n', spike_file_path);
end
fprintf('--------------------------\n');
fprintf('Processed %d sessions.\n', session_num);
% save session meta info
meta_folder = fullfile(root, 'Data', 'Working', 'Meta');
check_path(meta_folder);
meta_file_name = 'all_session_info_KZ.mat';
meta_file_path = fullfile(meta_folder, meta_file_name);
save(meta_file_path, 'all_session_info', 'segmentNames', '-v7.3');
fprintf('Saved all session info to %s\n', meta_file_path);