%% J_distance.m - Scatter plot of J vs channel distance
% 

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

training_tasks = {'Slayer', 'Zeppelin', 'Emperor'};
debug = false;

total_time = 0;
total_task_num = 0;

failed_list = {};
skipped = 0;
failed = 0;
success = 0;

for training_idx = 1:length(training_tasks)
    training_task = training_tasks{training_idx};

    % load training task
    training_task_folder = fullfile(root, 'Data', 'Working', 'Meta', 'training_tasks');
    training_task_name = sprintf('training_task_%s.mat', training_task);
    training_task_path = fullfile(training_task_folder, training_task_name);
    if ~isfile(training_task_path)
        fprintf('Training task file not found: %s\n', training_task_path);
        continue;
    end
    load(training_task_path, 'tasks');
    task_num = length(tasks);
    fprintf('Loaded training task: %s, with %d sessions.\n', training_task, task_num);
    total_task_num = total_task_num + task_num;

    for kernel_idx = 1:3
        all_J = [];
        all_dist = [];

        % run each session
        for task_idx = 1:task_num
            task = tasks{task_idx};
            try

                tick_session = tic;
                fprintf("Task %d/%d: %s, session%d\n", task_idx, task_num, task.dataset_name, task.session_idx);
                skip_flag = true;
                dataset_name = task.dataset_name;
                border_name = task.border_name;
                session_idx = task.session_idx;

                config = task.config;
                kernel_name = config.kernel;
                reg = config.reg;
                shuffle_size=config.shuffle_size;
                max_epoch=config.max_epochs;
                fold_num = config.crossval_fold_num;
                n_PS_kernel = 0;
                
                % load model parameters
                model_folder = fullfile(root, 'Data', 'Working', 'GLM_models');
                model_file = sprintf('GLM_%s_s%d_shuffle%d_%s_%s_epoch%d_fold%d.mat', dataset_name, session_idx, 0, kernel_name, ...
                    reg.name, max_epoch, 0);
                target_path = fullfile(model_folder, model_file);

                load(target_path, 'model_par', 'model_err', 'N');
                J_mat = model_par(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));
                J_err = model_err.total(:, (N*(kernel_idx-1) + n_PS_kernel + 2):(N*kernel_idx + n_PS_kernel + 1));

                % load channel
                channel_folder = fullfile(root, 'Data', 'Working', 'raster');
                channel_file = sprintf('raster_%s_%d.mat', dataset_name, session_idx);
                channel_path = fullfile(channel_folder, channel_file);
                load(channel_path, "channel", 'cell_area');

                assert(numel(channel) == N, "Channel number does not match model parameter dimension");
                assert(numel(cell_area) == N, "Cell area number does not match model parameter dimension");

                % filter within area connections
                areas = {'ACC', 'VLPFC', 'Thalamus'};
                area_num = numel(areas);

                J_areas = cell(area_num, 1);
                dist_areas = cell(area_num, 1);

                for area_idx = 1:area_num
                    area = areas{area_idx};
                    area_channel_idx = ismember(cell_area, area);
                    area_channel = channel(area_channel_idx);
                    area_mat = J_mat(area_channel_idx, area_channel_idx);
                    area_err = J_err(area_channel_idx, area_channel_idx);
                    area_N = sum(area_channel_idx);

                    for i = 1:area_N
                        for j = 1:area_N
                            ch_i = area_channel(i);
                            ch_j = area_channel(j);
                            dist_ij = abs(ch_i - ch_j);
                            J_ij = area_mat(i, j);
                            if isnan(J_ij) || i==j
                                continue;
                            end

                            all_J(end+1) = J_ij; %#ok<SAGROW>
                            all_dist(end+1) = dist_ij; %#ok<SAGROW>
                            J_areas{area_idx}(end+1) = J_ij;
                            dist_areas{area_idx}(end+1) = dist_ij;
                        end
                    end
                end

                % scatter plot
                f = figure('Position', [100, 100, 1800, 1200], 'Visible', 'off');
                tiledlayout(2, 2);
                for area_idx = 1:area_num
                    nexttile;
                    scatter(dist_areas{area_idx}, J_areas{area_idx}, 3, 'filled', 'MarkerFaceAlpha', 0.5);
                    xlabel('Channel Distance');
                    ylabel('J Value');
                    title(sprintf('%s - %s', dataset_name, areas{area_idx}));
                end

                % last tile: all areas
                nexttile;
                all_J_task = [J_areas{1}, J_areas{2}, J_areas{3}];
                all_dist_task = [dist_areas{1}, dist_areas{2}, dist_areas{3}];
                scatter(all_dist_task, all_J_task, 3, 'filled', 'MarkerFaceAlpha', 0.5);
                xlabel('Channel Distance');
                ylabel('J Value');
                title(sprintf('%s - All Areas', dataset_name));

                % save figure
                figure_folder = fullfile(root, 'Figures', 'J_distance');
                check_path(figure_folder);
                figure_file = sprintf('J_distance_%s_kernel%d.png', dataset_name, kernel_idx);
                save_path = fullfile(figure_folder, figure_file);
                saveas(f, save_path);

                % time
                elapsed_time = toc(tick_session);
                total_time = total_time + elapsed_time;
            
            catch ME
                fprintf("Failed: %s\n", ME.message);
                failed = failed + 1;
                failed_list{end+1} = {dataset_name, int2str(session_idx), ME.message}; %#ok<SAGROW> 
                if debug
                    rethrow(ME); %#ok<*UNRCH>
                end
            end
        end
        % all sessions scatter
        f = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
        scatter(dist_areas{area_idx}, J_areas{area_idx}, 3, 'filled', 'MarkerFaceAlpha', 0.5);
        xlabel('Channel Distance');
        ylabel('J Value');
        title(sprintf('%s - All Sessions', training_task));
        % save figure
        figure_folder = fullfile(root, 'Figures', 'J_distance');
        check_path(figure_folder);
        figure_file = sprintf('J_distance_%s_kernel%d_all_sessions.png', training_task, kernel_idx);
        save_path = fullfile(figure_folder, figure_file);
        saveas(f, save_path);    
    end
end

% fprintf("Total: %d, Success: %d, Skipped: %d, Failed: %d\n", total_task_num, success, skipped, failed);
% save failed_list
% folder = fullfile(root, 'Data', 'Working', 'log');
% check_path(folder);
% save(fullfile(folder, 'failed_list.mat'), 'failed_list');
for i=1:numel(failed_list)
    fprintf("Failed: %s, %s, %s\n", failed_list{i}{1}, failed_list{i}{2}, failed_list{i}{3});
end