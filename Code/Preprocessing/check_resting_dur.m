% load cell-based PDS file, print detailed session info.
% extract trial-based spike timings and rasters.

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
dt=0.001;

% load dataset metadata
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums', 'cortex_files', 'thalamus_files', 'eyeID_files');

areas = {'ACC', 'VLPFC', 'Thalamus'};
area_num = length(areas);
cell_session_names = cell(1, area_num);
cell_ids = cell(1, area_num);
cell_filenames = cell(1, area_num);

for area_idx = 1:area_num
    % Load all cell IDs and session names for this area
    area = areas{area_idx};
    folder_name = fullfile(root, 'Data', 'Experimental', 'PDS', area);
    file_names = {dir(folder_name).name}.'; % format: LemmyKim-#date-00#-MYInfoPavChoice_cl##_PDS.mat
    splited = cellfun(@(name) split(name, ["-", "_"]), file_names, 'UniformOutput', false);
    splited = splited(3:end); % remove . and ..
    cell_filter = cellfun(@(x) length(x)>=5, splited, 'UniformOutput', true);
    splited = splited(cell_filter);
    cell_file_num = length(splited);
    fprintf('Area: %s, Total files: %d, Valid cells: %d\n', ...
        area, length(file_names)-2, cell_file_num);
    cell_session_names{area_idx} = cellfun(@(x) [x{2}, '-', x{3}], splited, 'UniformOutput', false);
    file_names = file_names(3:end); % remove . and ..
    cell_filenames{area_idx} = file_names(cell_filter);
    % simrec naming: LemmyKim-#date-00#-MYInfoPavChoice_NovelReward_cl##_PDS.mat
    simrec_filter = cellfun(@(x) strcmp(x{5}, 'NovelReward'), splited, 'UniformOutput', true);
    cell_ids{area_idx} = cell(1, cell_file_num);
    for i = 1:cell_file_num
        if simrec_filter(i)
            cell_ids{area_idx}{i} = splited{i}{6};
        else
            cell_ids{area_idx}{i} = splited{i}{5};
        end
    end
end

% disp(cell_filenames);
% disp(cell_session_names);
all_valid_times = cell(dataset_num, 1);

for dataset_idx = 1:dataset_num
    dataset_name = dataset_names{dataset_idx};
    session_num = session_nums(dataset_idx);
    fprintf('=========================\n');
    valid_times = zeros(session_num, 2); % session x (pre/post)

    for session_idx = 1:session_num
        cortex_file = cortex_files{dataset_idx}{session_idx};
        thalamus_file = thalamus_files{dataset_idx}{session_idx};
        eyeID_file = eyeID_files{dataset_idx}{session_idx};
        has_cortex = ~isempty(cortex_file);
        has_thalamus = ~isempty(thalamus_file);
        has_eyeID = ~isempty(eyeID_file);

        fprintf('-------------------------\n');
        fprintf('Dataset(%d/%d): %s, Session(%d/%d):\n - Cortex file: %s\n - Thalamus file: %s\n - EyeID file: %s\n\n', ...
            dataset_idx, dataset_num, dataset_name, session_idx, session_num, cortex_file, thalamus_file, eyeID_file);

        if ~has_eyeID
            fprintf(' - No EyeID file\n\n');
        else
            % load eyeID data
            folder_name_eye = fullfile(root, 'Data', 'Experimental', 'eyeID');
            file_path_eye = fullfile(folder_name_eye, sprintf('eyeID-%s.mat', eyeID_file));
            load(file_path_eye, "r");
            % if r has pre and post fields, use them
            if isfield(r, 'pre') && isfield(r, 'post')
                eyeID_pre = r.pre;
                eyeID_post = r.post;
                has_post = true;
            else
                eyeID_pre = r;
                eyeID_post = [];
                has_post = false;
            end
            fprintf(' - EyeID loaded.');
            if has_post
                fprintf(' (Pre and Post)\n\n');
            else
                fprintf(' (Pre only)\n\n');
            end
        end

        states = {'Task', 'RestOpen', 'RestClose'};
        for state_idx = 2:3
            state = states{state_idx};

            if strcmp(state, 'RestOpen') || strcmp(state, 'RestClose')
                if ~has_eyeID
                    fprintf('   - No EyeID file, skip %s %s %s session%d\n\n', dataset_name, area, state, session_idx);
                    continue;
                end
            end

            % task trials
            subsession_types = {'Pre', 'Post'};
            
            % simrec have different taskIDs than other two
            if contains(dataset_name, 'Noinj')
                taskIDs = [1, 2];
            else
                taskIDs = [1.1, 1.2];
            end
            
            for subsession_idx = 1:2
                subsession = subsession_types{subsession_idx};
                taskID = taskIDs(subsession_idx);

                if ~has_post && strcmp(subsession, 'Post')
                    fprintf('   - No Post session for EyeID, skip %s %s %s session%d\n\n', dataset_name, area, subsession, session_idx);
                    continue;
                end
                
                if strcmp(subsession, 'Pre')
                    eyeID = eyeID_pre;
                else
                    eyeID = eyeID_post;
                end
                    
                fprintf(" - Loading: %s, %s...\n", subsession, state);

                state_last_spike = Inf;

                for area_idx = 1:area_num
                    area = areas{area_idx};
                    folder_name = fullfile(root, 'Data', 'Experimental', 'PDS', area);

                    if strcmp(area, 'Thalamus')
                        session_filename = thalamus_file;
                        if ~has_thalamus
                            fprintf(' - No Thalamus file, skip %s %s session%d\n\n', dataset_name, area, session_idx);
                            continue;
                        end
                        if strcmp(subsession, 'Post')
                            continue;
                        end
                    else
                        session_filename = cortex_file;
                        if ~has_cortex
                            fprintf(' - No Cortex file, skip %s %s session%d\n\n', dataset_name, area, session_idx);
                            continue;
                        end
                    end
                    session_file_idx = strcmp(cell_session_names{area_idx}, session_filename);
                    N = sum(session_file_idx);
                    cell_id = cell_ids{area_idx}(session_file_idx);
                    file_names = cell_filenames{area_idx}(session_file_idx);

                    % area_cellnames = cell_filenames{area_idx};
                    % disp(area_cellnames);
                    % disp(file_names);
                    
                    fprintf('   - Area: %s, Cell num: %d\n', area, N);
                    % Loading a single session&state starts here

                    trial_num    = NaN;
                    spikes       = NaN;
                    rasters      = NaN;
                    firing_rates = NaN;
                    cuetype      = NaN;
                    channel      = NaN;
                    trial_len    = NaN;
                    max_duration  = 0;
                    min_duration  = 9999;
                    max_trial_len = 0;
                    min_trial_len = 999999;
                    
                    area_last_spike = 0;

                    % loop over neurons to check consistency
                    first_neuron = true;
                    for i = 1:N
                        % fprintf("Loading: %s, %s, %s, session%d, #%d...\n", control, area, subsession, session_idx, i);
                        file_name = file_names{i};
                        filepath = fullfile(folder_name, file_name);
                        load(filepath, "PDS");
                        
                        % ---choose trials
                        if strcmp(state, 'Task')
                            if first_neuron
                                % % taskID counting
                                % fprintf('   - TaskID counts:\n');
                                % all_taskIDs = unique(PDS.taskID);
                                % taskID_counts = zeros(length(all_taskIDs), 1);
                                % for t = 1:length(all_taskIDs)
                                %     taskID_counts(t) = sum(PDS.taskID==all_taskIDs(t));
                                %     fprintf('     - TaskID %.1f: %d trials\n', all_taskIDs(t), taskID_counts(t));
                                % end
                                % fprintf('\n');
                            end
                            
                            % taskID fix: NaN, 1000, 1001 to 1.1, 1.1, 1.2, due to old data saving issue.
                            if any(isnan(PDS.taskID))
                                % if first_neuron
                                %     fprintf("     - %d NaN found in task ID, file: %s \n", sum(isnan(PDS.taskID)), file_name);
                                % end
                                PDS.taskID(isnan(PDS.taskID)) = 1.1;
                            end
                            if any(PDS.taskID==1000)
                                % if first_neuron
                                %     fprintf("     - %d 1000 found in task ID, file: %s \n", sum(PDS.taskID==1000), file_name);
                                % end
                                PDS.taskID(PDS.taskID==1000) = 1.1;
                            end
                            if any(PDS.taskID==1001)
                                % if first_neuron
                                %     fprintf("     - %d 1001 found in task ID, file: %s \n", sum(PDS.taskID==1001), file_name);
                                % end
                                PDS.taskID(PDS.taskID==1001) = 1.2;
                            end

                            % -choice and pav trials
                            valid_trials = (PDS.taskID==taskID)&(PDS.goodtrial);
                            % -pav trials only
                            pav_trials = (PDS.taskID==taskID)&(PDS.goodtrial)&isnan(PDS.timeoffer1cho)&isnan(PDS.timeoffer2cho);
                            % -choice trials only
                            choice_trials = (PDS.taskID==taskID)&(PDS.goodtrial)&(~isnan(PDS.timeoffer1cho)|~isnan(PDS.timeoffer2cho));
                            
                            % Use choice trials
                            selected_trial = choice_trials;
                            selected_spikes = PDS.sptimes(selected_trial);

                            % if first_neuron
                            % fprintf('   - Total trials: %d.\n', length(PDS.taskID));
                            % fprintf('     - Valid trials (good, taskID %.1f): %d.\n', taskID, sum(valid_trials));
                            % fprintf('     - Pav trials: %d.\n', sum(pav_trials));
                            % fprintf('     - Choice trials: %d.\n\n', sum(choice_trials));
                            % end

                        end

                        % ---choose phases in trial
                        switch state
                        case 'All'
                            % -all
                            selected_start = PDS.timetargeton(selected_trial);
                            selected_end   = PDS.timeInfotargetoff(selected_trial);

                        case 'Offer1'
                            % -after offer1
                            selected_start = PDS.timetargeton(selected_trial);
                            selected_end   = selected_start + 1;

                        case 'Offer2'
                            % -after offer2
                            selected_start = PDS.timetargeton1(selected_trial);
                            selected_end   = selected_start + 1;

                        case 'Decision'
                            % -before decision
                            decision_time = [PDS.timeoffer1cho;PDS.timeoffer2cho];
                            decision_time(isnan(decision_time)) = 0;
                            decision_time = sum(decision_time, 1);
                            selected_end = decision_time(selected_trial);
                            selected_start = selected_end - 1;

                        case 'InfoAnti'
                            % -before info cue
                            selected_end = PDS.timeInfotargeton(selected_trial);
                            selected_start = selected_end - 1;

                        case 'InfoResp'
                            % -after info cue
                            selected_start = PDS.timeInfotargeton(selected_trial);
                            selected_end = selected_start + 1;

                        case 'Reward'
                            % -after reward
                            selected_start = PDS.timereward(selected_trial);
                            selected_end = selected_start + 1;

                        case 'RandomA'
                            % -random 1s in the trial
                            if ~isa(randomA_start, "cell")
                                randomA_start = cell(session_num, 2);
                            end
                            if isempty(randomA_start{session_idx, subsession_idx})
                                randomA_start{session_idx, subsession_idx} = rand(1, sum(selected_trial));
                            end
                            trial_range = PDS.trialendtime(selected_trial) - PDS.trialstarttime(selected_trial) - 1;
                            selected_start = randomA_start{session_idx, subsession_idx} .* trial_range;
                            selected_end = selected_start + 1;

                        case 'RandomB'
                            % -random 1s in the trial
                            if ~isa(randomB_start, "cell")
                                randomB_start = cell(session_num, 2);
                            end
                            if isempty(randomB_start{session_idx, subsession_idx})
                                randomB_start{session_idx, subsession_idx} = rand(1, sum(selected_trial));
                            end
                            trial_range = PDS.trialendtime(selected_trial) - PDS.trialstarttime(selected_trial) - 1;
                            selected_start = randomB_start{session_idx, subsession_idx} .* trial_range;
                            selected_end = selected_start + 1;

                        case 'RandomShort'
                            % -random 2s in the trial
                            if ~isa(randomShort_start, "cell")
                                randomShort_start = cell(session_num, 2);
                            end
                            if isempty(randomShort_start{session_idx, subsession_idx})
                                randomShort_start{session_idx, subsession_idx} = rand(1, sum(selected_trial));
                            end
                            trial_range = PDS.trialendtime(selected_trial) - PDS.trialstarttime(selected_trial) - 2.00005;
                            selected_start = randomShort_start{session_idx, subsession_idx} .* trial_range;
                            selected_end = selected_start + 2.00005;

                        case 'RandomLong'
                            % -random 5s in the trial
                            if ~isa(randomLong_start, "cell")
                                randomLong_start = cell(session_num, 2);
                            end
                            if isempty(randomLong_start{session_idx, subsession_idx})
                                randomLong_start{session_idx, subsession_idx} = rand(1, sum(selected_trial));
                            end
                            trial_range = PDS.trialendtime(selected_trial) - PDS.trialstarttime(selected_trial) - 5.00005;
                            selected_start = randomLong_start{session_idx, subsession_idx} .* trial_range;
                            selected_end = selected_start + 5.00005;

                        case 'Task'
                            % -full task
                            selected_start = zeros(1, sum(selected_trial));
                            selected_end = PDS.trialendtime(selected_trial) - PDS.trialstarttime(selected_trial);
                                                    
                        case 'RestOpen'
                            % -resting state eye open
                            selected_start = eyeID.eCloseOnset_t - eyeID.eOpenDur;
                            selected_end = eyeID.eCloseOnset_t;

                        case 'RestClose'
                            % -resting state eye close
                            selected_start = eyeID.eCloseOnset_t;
                            selected_end = eyeID.eCloseOnset_t + eyeID.eCloseDur;
                            
                        end

                        max_dur_cell = max(selected_end - selected_start);
                        min_dur_cell = min(selected_end - selected_start);
                        if max_dur_cell > max_duration
                            max_duration = max_dur_cell;
                        end
                        if min_dur_cell < min_duration
                            min_duration = min_dur_cell;
                        end

                        %% Process spikes for this neuron
                        if first_neuron
                            % if this is the first neuron, initialize output variables
                            if strcmp(state, 'Task')
                                cuetype = PDS.Cuetype(selected_trial);
                                trial_num = sum(selected_trial);
                            else
                                trial_num = length(selected_start);
                            end
                            spikes = cell(N, trial_num);
                            rasters = cell(1, trial_num);
                            firing_rates = cell(1, trial_num);
                            channel = zeros(1, N);
                            trial_len = ones(1, trial_num) * NaN;
                            for j=1:trial_num
                                rasters{j} = NaN;
                                firing_rates{j} = ones(N, 1) * NaN;
                            end
                            first_neuron = false;
                        else
                            % check consistency
                            if strcmp(state, 'Task')
                                assert(trial_num == sum(selected_trial), 'Inconsistent trial num in %s, neuron %d', session_filename, i);
                                assert(all(cuetype == PDS.Cuetype(selected_trial)), 'Inconsistent cuetype in %s, neuron %d', session_filename, i);
                            else
                                assert(trial_num == length(selected_start), 'Inconsistent trial num in %s, neuron %d', session_filename, i);
                            end
                        end

                        channel(i) = PDS.channel;
                        % spikes & rasters for each trial
                        last_spike_time = 0;
                        for j=1:trial_num
                            if strcmp(state, 'Task')
                                spike_trial = selected_spikes{j};
                            else
                                switch subsession
                                case 'Pre'
                                    % -resting state eye open
                                    % % skip Muscimol session 2
                                    % if strcmp(session_name, '11012023') 
                                    %     spike_trial = [];
                                    % else
                                    %     if strcmp(area, 'Thalamus')
                                    %         spike_trial = PDS.sptimes_aftertask;
                                    %     else
                                    %         spike_trial = PDS.sptimes_betweentask{1};
                                    %     end
                                    % end   
                                    
                                    if strcmp(area, 'Thalamus')
                                        spike_trial = PDS.sptimes_aftertask;
                                    else
                                        spike_trial = PDS.sptimes_betweentask{1};
                                    end

                                case 'Post'
                                    % -resting state eye close
                                    spike_trial = PDS.sptimes_aftertask;
                                end
                            end

                            last_spike_time = max([last_spike_time; spike_trial(:)]);

                            % % only keep spikes from cue to reward,
                            % % align to cue.
                            % spike_trial = spike_trial(spike_trial>selected_start(j) & spike_trial<selected_end(j));
                            % spike_trial = spike_trial - selected_start(j);
                            
                            % % auto fit start and end
                            % trial_duration = selected_end(j) - selected_start(j);

                            % % fixed trial duration
                            % % trial_duration = 1;

                            % trial_edges = 0:dt:trial_duration;
                            % trial_B = length(trial_edges) - 1;
                            % trial_len(j) = trial_B;

                            % if trial_B > max_trial_len
                            %     max_trial_len = trial_B;
                            % end
                            % if trial_B < min_trial_len
                            %     min_trial_len = trial_B;
                            % end

                            % % Convert spike times to raster
                            % spikes{i, j} = spike_trial;
                            % raster = histcounts(spike_trial, trial_edges);
                            % raster(raster>1) = 1;
                            % if isnan(rasters{j})
                            %     rasters{j} = zeros(N, trial_B);
                            % end
                            % rasters{j}(i, :) = raster;
                            % firing_rates{j}(i) = mean(raster);
                        end % trial loop
                        fprintf('     - Neuron %d/%d processed. Area: %s, last spike time: %.3f s\n', i, N, area, last_spike_time);
                        
                        area_last_spike = max(last_spike_time, area_last_spike);
                    end % neuron loop
                    fprintf('   - Area %s last spike time: %.3f s\n', area, area_last_spike);
                    if N>0
                        state_last_spike = min(area_last_spike, state_last_spike); % all three areas have spikes up to this time
                    end
                end % area loop
                fprintf(' - State %s last spike time: %.3f s\n\n', state, state_last_spike);
                max_eyeID = max(eyeID.eCloseOffset_t);
                if isempty(max_eyeID)
                    max_eyeID = 0;
                end
                fprintf(' - EyeID max time: %.3f s\n\n', max_eyeID);
                valid_time = min(state_last_spike, max_eyeID);
                fprintf(' - Valid resting time duration: %.3f s\n\n', valid_time);
                valid_times(session_idx, subsession_idx) = valid_time;
            end % subsession loop
        end % state loop
    end % session loop
    all_valid_times{dataset_idx} = valid_times;
end % dataset loop

%% Save results
save_folder = fullfile(root, 'Data', 'Working', 'Meta');
save_path = fullfile(save_folder, 'resting_valid_durations.mat');
save(save_path, 'all_valid_times', '-v7.3');