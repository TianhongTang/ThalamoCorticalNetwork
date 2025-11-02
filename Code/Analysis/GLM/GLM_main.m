%% GLM raster fitting
% Input: raw raster file "../GLM_data/[dataset_name]/raster_[dataset_name]_[session]_0.mat"

% Required variables in the input file:
% "n_trial": integer, trial number.
% "trial_len": int(1, n_trial), time bin number of each trial.
% "rasters": cell(1, n_trial), each element is a trial.
% Each raster is a (N, trial_len(i)) binary matrix, N is number of neurons.

%% Get root folder
code_depth = 4;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
% include code folder and utils
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Main
clear
% task_names = {'MuscimolPre_cortex', 'MuscimolPost_cortex', 'MuscimolPre_full'};
task_names = {...

    };
force_rebuild = false;
force_retrain = false;
debug = true;
total_training = 0;
skipped = 0;
failed= 0;
success = 0;

failed_list = {};

% task_names = {'MuscimolPreDecision_full', 'SalinePreDecision_full', 'SimRecPreDecision_full'};
% task_names = {'MuscimolPostDecision_full', 'SalinePostDecision_full'};
% trial_names = {'100B', '50BI', '50BN', '100S', '0'};
total_time = 0;
for task_idx=1:length(task_names)
    % for cuetype=1:5
    % % compare if is Muscimol sessions
    % if strcmp(task_names{task_idx}(1:8), 'Muscimol')
    %     session_idxs = [6,7,8,9,10,1,4,5,2,3];
    % else
    %     session_idxs = [1,4,5,2,3];
    % end
    session_idxs = 1:10;

    for session_idx = session_idxs
        for trial_idx = 1:1
            try
                % if task_idx==1 && session_idx<10
                %     continue;
                % end
                % session_idx = session*100+cuetype;
                % fprintf("Main: session%d, cue%d\n", session, cuetype);

                % dataset_name = [task_names{task_idx}, '_', trial_names{trial_idx}];
                tick_session = tic;
                dataset_name = task_names{task_idx};
                fprintf("Main: %s, session%d\n", dataset_name, session_idx);
                total_training = total_training + 1;
                skip_flag = true;

                %% parameters
                % dataset_name = 'generated';
                % % session = 0;
                % kernel_name = 'expDecay10';
                
                % dataset_name = 'taskCue';
                % % session = 1;
                % kernel_name = 'expDecay10';
                
                % task
                % session = 1;
                % kernel_name = 'expDecay10';
                % kernel_name = 'expMulti200';
                kernel_name = 'DeltaPure';
                % kernel_name = 'DoublePure';
                
                % reg.l1=5;
                % reg.l2=0;
                % reg.name='L1=5';    
                
                reg.l1=0;
                reg.l2=2;
                reg.name='L2=2';
                
                % reg.l1=2;
                % reg.l2=0;
                % reg.name='L1=2';
                
                % reg.l1=0;
                % reg.l2=0;
                % reg.name='NoReg';
                
                shuffle_size=0;
                max_epoch=2500;
                
                
                %% generate shuffled raster
                fprintf("Shuffle rasters\n");
                tic;
                for shuffle_seed=1:shuffle_size
                    % skip if already exists
                    target_path = ['../GLM_data/', dataset_name, '/raster_', ...
                        dataset_name, '_', int2str(session_idx),  '_', int2str(shuffle_seed), '.mat'];
                    if isfile(target_path) && ~force_rebuild
                        fprintf("Skip %d. \n", shuffle_seed);
                        continue;
                    end

                    skip_flag = false;
                    shuffle_across_trial=(shuffle_seed<2);
                    shuffle(dataset_name, session_idx, shuffle_seed, shuffle_across_trial);
                end
                toc;

                %% convolve predj and combine trials
                fprintf("Convolution\n");
                tic;
                for shuffle_seed=0:shuffle_size % seed=0: original data (no shuffle)
                    % skip if already exists
                    target_path = ['../GLM_data/', dataset_name, '/GLMdata_', dataset_name,...
                        '_', int2str(session_idx), '_', kernel_name,  '_', int2str(shuffle_seed), '.mat'];
                    if isfile(target_path) && ~force_rebuild
                        fprintf("Skip %d. \n", shuffle_seed);
                        continue;
                    end

                    skip_flag = false;
                    convolution(dataset_name, session_idx, kernel_name, shuffle_seed);
                end
                toc;
                %% GLM inference
                for shuffle_seed=0:shuffle_size
                    fprintf("Training %d\n", shuffle_seed);

                    % skip if already exists
                    foldername = ['../GLM_model/', dataset_name];
                    target_path = [foldername, '/GLM_', dataset_name, '_', ...
                    int2str(session_idx), '_', kernel_name, '_', int2str(shuffle_seed), '_', ...
                    reg.name, '_', int2str(max_epoch)];
                    if isfile([target_path, '.mat']) && ~force_retrain
                        fprintf("Skip. \n");
                        continue;
                    end

                    skip_flag = false;
                    tic;
                    GLM_multi_kernel_err(dataset_name, session_idx, kernel_name, shuffle_seed, max_epoch, reg, 1, 5e-3);
                    toc;
                end

                %% plot
                % plot_GLM(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size);
                % type_file = ['../GLM_data/', dataset_name, '/celltype_', dataset_name, '_', int2str(session_idx), ...
                % '.mat'];
                % load(type_file, "cell_type");

                fprintf("Plotting\n");
                tic;
                channel_file = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', int2str(session_idx), ...
                '_0.mat'];
                load(channel_file, "channel");
                plot_GLM_sorted(dataset_name, session_idx, kernel_name, max_epoch, reg, shuffle_size, "idx", channel);
                toc;
                % 
                % %% plot gen
                % % plot_generated(dataset_name)
                % session_time = toc(tick_session);
                % fprintf("Session time: %f\n", session_time);
                % total_time = total_time + session_time;
                % fprintf("Total time: %f\n", total_time);

                if skip_flag
                    skipped = skipped + 1;
                else
                    success = success + 1;
                end

            catch ME
                fprintf("Failed: %s\n", ME.message);
                failed = failed + 1;
                failed_list{end+1} = {dataset_name, int2str(session_idx), ME.message};
                if debug
                    throw(ME);
                end
            end
        end
    end
    % plot_rest(session_idx, kernel_name, max_epoch, reg, shuffle_size);
end

fprintf("Total: %d, Success: %d, Skipped: %d, Failed: %d\n", total_training, success, skipped, failed);
% save failed_list
save('../GLM_data/failed_list.mat', 'failed_list');
if failed>0
    for i=1:length(failed_list)
        fprintf("Failed: %s, %s, %s\n", failed_list{i}{1}, failed_list{i}{2}, failed_list{i}{3});
    end
end