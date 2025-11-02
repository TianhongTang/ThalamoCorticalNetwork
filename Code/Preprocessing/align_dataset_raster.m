% Align task, eyeclose and eyeopen datasets to have the same duration

%% Load data
control = 'Saline';
% control = 'Muscimol';
kernel = 'DeltaPure';

tasks = {...
    'PreTask_full', 'PreTask_cortex', 'PostTask_cortex',...
    'PreRestClose_full', 'PreRestClose_cortex', 'PostRestClose_cortex',...
    'PreRestOpen_full', 'PreRestOpen_cortex', 'PostRestOpen_cortex'};

% check data size
for session = 1:10
    min_duration = 10000000;
    for task_idx = 1:length(tasks)
        task = tasks{task_idx};

        % file_path = ['../GLM_data/', control, task, '/GLMdata_', control, task, '_', int2str(session), '_', kernel, '_0.mat'];
        % load(file_path, 'N', 'B');
        file_path = ['../GLM_data/', control, task, '/raster_', control, task, '_', int2str(session), '_0.mat'];
        load(file_path, 'n_trial', 'trial_len', 'N');
        file_path = ['../GLM_data/kernel_', kernel, '.mat'];
        load(file_path, 'kernel_len');

        effective_len = sum(trial_len - kernel_len + 1);
        B = effective_len;

        min_duration = min(min_duration, B);
        fprintf('Session %d, task %s, N = %d, B = %d, n_trial = %d\n', session, task, N, B, n_trial);
        fprintf('Session %d, task %s, effective_len = %d\n', session, task, effective_len);
    end
    seconds = floor(min_duration/1000);
    fprintf('Session %d, min_duration = %d, %d:%d\n', session, min_duration, floor(seconds/60), mod(seconds, 60));

    % align data to min_duration
    modes = {'First', 'Last'};
    for task_idx = 1:length(tasks)
        task = tasks{task_idx};
        fprintf('Aligning %s %s\n', control, task);
        fprintf('Loading...');



        for mode_idx = 1:length(modes)
            % raster file and border file
            file_path = ['../GLM_data/', control, task, '/raster_', control, task, '_', int2str(session), '_0.mat'];
            load(file_path, 'n_trial', 'cell_id', 'cell_area', 'channel', 'session_name_full', 'trial_len', 'rasters');
            file_path = ['../GLM_data/', control, task, '/borders_', control, task, '_', int2str(session), '.mat'];
            load(file_path, 'borders');
            fprintf('Done\n');

            mode = modes{mode_idx};
            fprintf('Mode: %s\n', mode);

            raster_aligned = cell(1, 0);
            trial_len_aligned = zeros(1, 0);
            current_eff_dur = 0;

            switch mode
                case 'First'
                    for trial = 1:n_trial
                        trial_eff_dur = trial_len(trial) - kernel_len + 1;
                        if current_eff_dur + trial_eff_dur >= min_duration
                            cut_len = min_duration - current_eff_dur + kernel_len - 1;
                            raster_aligned{end+1} = rasters{trial}(:, 1:cut_len);
                            trial_len_aligned(end+1) = cut_len;
                            break;
                        else
                            raster_aligned{end+1} = rasters{trial};
                            current_eff_dur = current_eff_dur + trial_eff_dur;
                            trial_len_aligned(end+1) = trial_eff_dur + kernel_len - 1;
                        end
                    end
                case 'Last'
                    for trial = n_trial:-1:1
                        trial_eff_dur = trial_len(trial) - kernel_len + 1;
                        if current_eff_dur + trial_eff_dur >= min_duration
                            cut_len = min_duration - current_eff_dur + kernel_len - 1;
                            raster_aligned{end+1} = rasters{trial}(:, end-cut_len+1:end);
                            trial_len_aligned(end+1) = cut_len;
                            break;
                        else
                            raster_aligned{end+1} = rasters{trial};
                            current_eff_dur = current_eff_dur + trial_eff_dur;
                            trial_len_aligned(end+1) = trial_eff_dur + kernel_len - 1;
                        end
                    end
            end

            rasters = raster_aligned;
            trial_len = trial_len_aligned;
            n_trial = length(rasters);

            firing_rates = cell(1, n_trial);
            for i = 1:n_trial
                firing_rates{i} = mean(rasters{i}, 2);
            end

            full_name = [control, task, '_Align', mode];
            fprintf('Saving...');
            check_path(['../GLM_data/', full_name]);

            % raster file and border file
            file_path = ['../GLM_data/', full_name, '/raster_', full_name, '_', int2str(session), '_0.mat'];
            save(file_path, 'n_trial', 'cell_id', 'cell_area', 'channel', 'session_name_full', 'trial_len', 'rasters', 'firing_rates');
            file_path = ['../GLM_data/', full_name, '/borders_', full_name, '_', int2str(session), '.mat'];
            save(file_path, 'borders');

            fprintf('Done.\n');
        end
    end
end
