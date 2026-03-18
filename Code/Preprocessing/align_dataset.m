%% Align raster files of all states to the same duration
% Also discard short resting epochs.

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
STATUS_LOG = false;
SKIP_EXISTING = true;
resting_dur_threshold = 10; % minimum duration in seconds for resting epochs to be included in the analysis

% load data
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums', 'cortex_files', 'thalamus_files', 'eyeID_files');

% Load data
area_types = {'Full', 'Cortex'};
prepost_types = {'Pre', 'Post'};
states = {'RestOpen', 'RestClose'};
kernel = 'DeltaPure';

% load kernel length
folder_name = fullfile(root, 'Data', 'Working', 'kernel');
file_name = sprintf('kernel_%s.mat', kernel);
file_path = fullfile(folder_name, file_name);
load(file_path, 'meta');
kernel_len = meta.kernel_len;
align_kernel_name = kernel;
align_kernel_len = kernel_len;

tasks = {};
% register tasks
for dataset_idx = 1:dataset_num
    dataset_name = dataset_names{dataset_idx};
    session_num = session_nums(dataset_idx);
    for session_idx = 1:session_num
        for area_idx = 1:length(area_types)
            area_type = area_types{area_idx};
            task = struct();
            task.dataset_name        = dataset_name;
            task.area_type           = area_type;
            task.session_idx         = session_idx;
            [animal_name, injection] = split_dataset_name(dataset_name);
            task.animal_name         = animal_name;
            task.injection           = injection;

            prepost_states = {};
            for prepost_idx = 1:length(prepost_types)
                prepost = prepost_types{prepost_idx};
                if strcmp(area_type, 'Full') && strcmp(prepost, 'Post')
                    % All thalamus data in post sessions are not available
                    continue;
                end
                if contains(dataset_name, 'Noinj') && strcmp(prepost, 'Post')
                    % No injection session does not have post sessions 
                    continue;
                end

                for state_idx = 1:length(states)
                    state = states{state_idx};
                    if strcmp(dataset_name, 'SlayerNoinj') && contains(state, 'Rest')
                        % SlayerNoinj does not have RestOpen and RestClose sessions
                        continue;
                    end
                    prepost_states{end+1} = struct('prepost', prepost, 'state', state); %#ok<SAGROW>
                end
            end
            task.prepost_states = prepost_states;
            tasks{end+1} = task; %#ok<SAGROW>
        end
    end
end

% run tasks
task_num = length(tasks);
for task_idx = 1:task_num
    task = tasks{task_idx};
    dataset_name   = task.dataset_name;
    animal_name    = task.animal_name;
    injection      = task.injection;
    session_idx    = task.session_idx;
    area_type      = task.area_type;
    prepost_states = task.prepost_states;
    fprintf('----------------------------------------\n');
    fprintf('Task %d/%d: %s Session %d Area %s with %d prepost_states\n', ...
        task_idx, task_num, dataset_name, session_idx, area_type, numel(prepost_states));
    
    % find min duration among prepost_states
    min_duration = 10000000;
    min_longest = 10000000;
    for ps_idx = 1:numel(prepost_states)
        prepost = prepost_states{ps_idx}.prepost;
        state = prepost_states{ps_idx}.state;

        % load length info
        raster_info             = struct();
        raster_info.animal_name = animal_name;
        raster_info.injection   = injection;
        raster_info.prepost     = prepost;
        raster_info.state       = state;
        raster_info.area        = area_type;
        raster_info.align       = 'None';
        raster_info.session_idx = session_idx;

        folder_name = fullfile(root, 'Data', 'Working', 'raster');
        file_name = generate_filename('raster', raster_info);
        raster_path = fullfile(folder_name, file_name);
        fprintf('Loading %s ... \n', raster_path);
        load(raster_path, 'meta', 'data');
        trial_len = data.trial_len;
        N = meta.N;

        trial_filter = (trial_len > align_kernel_len) & (trial_len > resting_dur_threshold*1000); % only include trials longer than kernel length and threshold
        trial_len = trial_len(trial_filter);
        trial_num = sum(trial_filter);

        B = sum(trial_len - align_kernel_len + 1); % effective length
        longest = max(trial_len); % longest trial length

        min_duration = min(min_duration, B);
        if ~isempty(longest)
            min_longest = min(min_longest, longest);
        end
        fprintf('%s%s, N = %d, B = %d, trial_num = %d, ave = %d, max = %d.\n', prepost, state, N, B, trial_num, floor(B/trial_num), max(trial_len));
    end
    seconds = floor(min_duration/1000);
    seconds_longest = floor(min_longest/1000);
    fprintf('%s%s Session %d, min_duration = %d (%d:%d), min_longest = %d (%d:%d)\n',...
         dataset_name, area_type, session_idx, min_duration, floor(seconds/60),...
         mod(seconds, 60), min_longest, floor(seconds_longest/60), mod(seconds_longest, 60));

    % align data to min_duration
    modes = {'First', 'Last', 'Longest'};
    for ps_idx = 1:numel(prepost_states)
        prepost = prepost_states{ps_idx}.prepost;
        state = prepost_states{ps_idx}.state;
        if STATUS_LOG, fprintf('Aligning %s ...\n', [prepost, state]); end %#ok<UNRCH>

        for mode_idx = 1:numel(modes)
            % skip task state for Longest mode
            if strcmp(modes{mode_idx}, 'Longest') && strcmp(state, 'Task')
                continue;
            end
            % raster file and border file
            if STATUS_LOG, fprintf('Loading...'); end %#ok<UNRCH>
            
            % if already exists, skip
            aligned_meta = struct();
            aligned_meta.animal_name = animal_name;
            aligned_meta.injection = injection;
            aligned_meta.prepost = prepost;
            aligned_meta.state = state;
            aligned_meta.area = area_type;
            aligned_meta.align = modes{mode_idx};
            aligned_meta.session_idx = session_idx;

            save_folder = fullfile(root, 'Data', 'Working', 'raster');
            check_path(save_folder);
            file_name = generate_filename('raster', aligned_meta);
            file_path = fullfile(save_folder, file_name);
            if isfile(file_path) && SKIP_EXISTING
                fprintf('   - Already exists for %s %s %s session%d. Skipping...\n\n', dataset_name, area_type, modes{mode_idx}, session_idx);
                continue;
            end
            
            % load unaligned data
            unaligned_meta = aligned_meta;
            unaligned_meta.align = 'None';
            
            folder_name = fullfile(root, 'Data', 'Working', 'raster');
            file_name = generate_filename('raster', unaligned_meta);
            raster_path = fullfile(folder_name, file_name);
            load(raster_path, 'meta', 'data');
            unaligned_meta = meta;
            unaligned_data = data;
            trial_len      = data.trial_len;
            rasters        = data.rasters;

            % Align rasters
            if STATUS_LOG, fprintf('Done\n'); end %#ok<UNRCH>

            mode = modes{mode_idx};
            if STATUS_LOG, fprintf('Mode: %s\n', mode); end %#ok<UNRCH>

            raster_aligned = cell(1, 0);
            trial_len_aligned = zeros(1, 0);
            current_eff_dur = 0;
            
            % Exclude short trials
            trial_filter = (trial_len > align_kernel_len) & (trial_len > resting_dur_threshold*1000); % only include trials longer than kernel length and threshold
            trial_len = trial_len(trial_filter);
            rasters = rasters(trial_filter);
            trial_num = sum(trial_filter);

            switch mode
                case 'First'
                    for trial = 1:trial_num
                        trial_eff_dur = trial_len(trial) - kernel_len + 1;
                        if current_eff_dur + trial_eff_dur >= min_duration
                            cut_len = min_duration - current_eff_dur + kernel_len - 1;
                            raster_aligned{end+1} = rasters{trial}(:, 1:cut_len); %#ok<SAGROW>
                            trial_len_aligned(end+1) = cut_len; %#ok<SAGROW>
                            break;
                        else
                            raster_aligned{end+1} = rasters{trial}; %#ok<SAGROW>
                            current_eff_dur = current_eff_dur + trial_eff_dur;
                            trial_len_aligned(end+1) = trial_eff_dur + kernel_len - 1; %#ok<SAGROW>
                        end
                    end
                case 'Last'
                    for trial = trial_num:-1:1
                        trial_eff_dur = trial_len(trial) - kernel_len + 1;
                        if current_eff_dur + trial_eff_dur >= min_duration
                            cut_len = min_duration - current_eff_dur + kernel_len - 1;
                            raster_aligned{end+1} = rasters{trial}(:, end-cut_len+1:end); %#ok<SAGROW>
                            trial_len_aligned(end+1) = cut_len; %#ok<SAGROW>
                            break;
                        else
                            raster_aligned{end+1} = rasters{trial}; %#ok<SAGROW>
                            current_eff_dur = current_eff_dur + trial_eff_dur;
                            trial_len_aligned(end+1) = trial_eff_dur + kernel_len - 1; %#ok<SAGROW>
                        end
                    end
                case 'Longest'
                    % Center part of the longest trial
                    longest_idx = find(trial_len == max(trial_len), 1);
                    if ~isempty(longest_idx)
                        longest_len = trial_len(longest_idx);
                        t_start = floor((longest_len - min_longest)/2) + 1;
                        t_end = t_start + min_longest - 1;
                        raster_aligned{1} = rasters{longest_idx}(:, t_start:t_end);
                        trial_len_aligned(1) = min_longest; 
                    end
            
                otherwise
                    error('Unknown mode: %s', mode);

            end

            if strcmp(mode, 'Last')
                raster_aligned = fliplr(raster_aligned);
                trial_len_aligned = fliplr(trial_len_aligned);
            end

            firing_rates = cell(1, trial_num);
            for i = 1:trial_num
                firing_rates{i} = mean(rasters{i}, 2);
            end

            % Save aligned data
            % if STATUS_LOG, fprintf('Saving %s...\n', full_name); end %#ok<UNRCH>
            
            meta = unaligned_meta;
            data = unaligned_data;

            meta.align            = modes{mode_idx};
            meta.align_kernel     = align_kernel_name;
            meta.align_kernel_len = align_kernel_len;
            meta.trial_num        = numel(raster_aligned);

            data.rasters      = raster_aligned;
            data.trial_len    = trial_len_aligned;
            data.firing_rates = firing_rates;

            save_folder = fullfile(root, 'Data', 'Working', 'raster');
            check_path(save_folder);
            file_name = generate_filename('raster', meta);
            file_path = fullfile(save_folder, file_name);
            save(file_path, 'meta', 'data', '-v7.3');

            if STATUS_LOG, fprintf('Done.\n'); end %#ok<UNRCH>
        end
    end
end

