% Align task, eyeclose and eyeopen datasets, keep data length (B) the same.

%% Load data
% control = 'Muscimol';
control = 'Saline';
kernel = 'DeltaPure';

tasks = {'PreTask_full', 'PreTask_cortex', 'PostTask_cortex',...
    'PreRestClose_full', 'PreRestClose_cortex', 'PostRestClose_cortex',...
    'PreRestOpen_full', 'PreRestOpen_cortex', 'PostRestOpen_cortex'};

if strcmp(control, 'Muscimol')
    sessions = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
elseif strcmp(control, 'Saline')
    sessions = [1, 2, 3, 4, 5];
end

% check data size
for session = sessions
    min_duration = 10000000;
    for task_idx = 1:length(tasks)
        task = tasks{task_idx};

        file_path = ['../GLM_data/', control, task, '/GLMdata_', control, task, '_', int2str(session), '_', kernel, '_0.mat'];
        load(file_path, 'N', 'B');
        file_path = ['../GLM_data/', control, task, '/raster_', control, task, '_', int2str(session), '_0.mat'];
        load(file_path, 'n_trial');

        min_duration = min(min_duration, B);
        
        fprintf('Session %d, task %s, N = %d, B = %d, n_trial = %d\n', session, task, N, B, n_trial);
    end
    seconds = floor(min_duration/1000);
    fprintf('Session %d, min_duration = %d, %d:%d\n', session, min_duration, floor(seconds/60), mod(seconds, 60));

    % align data to min_duration
    % modes = {'First', 'Last', 'Random'};
    modes = {'Last'};
    for task_idx = 1:length(tasks)
        task = tasks{task_idx};
        fprintf('Aligning %s %s\n', control, task);
        fprintf('Loading...');
        file_path = ['../GLM_data/', control, task, '/GLMdata_', control, task, '_', int2str(session), '_', kernel, '_0.mat'];
        load(file_path, 'N', 'B', 'PS_kernels', 'conn_kernels', 'kernel_len', 'n_PS_kernel',...
            'n_conn_kernel', 'predjs_PS', 'predjs_conn', 'raster');
        raster_full = raster;
        predjs_PS_full = predjs_PS;
        predjs_conn_full = predjs_conn;
        B_full = B;

        % raster and border file 
        file_path = ['../GLM_data/', control, task, '/raster_', control, task, '_', int2str(session), '_0.mat'];
        load(file_path, 'n_trial', 'cell_id', 'cell_area', 'channel', 'session_name_full');
        file_path = ['../GLM_data/', control, task, '/borders_', control, task, '_', int2str(session), '.mat'];
        load(file_path, 'borders');
        fprintf('Done\n');

        for mode_idx = 1:length(modes)
            mode = modes{mode_idx};
            fprintf('Mode: %s\n', mode);
            switch mode
                case 'First'
                    selected = 1:min_duration;
                    % fprintf('First\n');
                case 'Last'
                    selected = (B_full-min_duration+1):B_full;
                    % fprintf('Last\n');
                case 'Random'
                    selected = randperm(B_full, min_duration);
                    % fprintf('Random\n');
            end
            % if strcmp(mode, 'First')
            %     selected = 1:min_duration;
            %     fprintf('First\n');
            % elseif strcmp(mode, 'Last')
            %     selected = (B-min_duration+1):B;
            %     fprintf('Last\n');
            % elseif strcmp(mode, 'Random')
            %     selected = randperm(B, min_duration);
            %     fprintf('Random\n');
            % end

            raster = raster_full(:, selected);
            predjs_PS = predjs_PS_full(:, selected, :);
            predjs_conn = predjs_conn_full(:, selected, :);
            % if n_PS_kernel > 0
            %     predjs_PS = predjs_PS_full(:, selected, :);
            % else
            %     predjs_PS = [];
            % end
            % if n_conn_kernel > 0
            %     predjs_conn = predjs_conn_full(:, selected, :);
            % else
            %     predjs_conn = [];
            % end
            B = min_duration;

            full_name = [control, task, '_Align', mode];
            fprintf('Saving...');
            file_path = ['../GLM_data/', full_name, '/GLMdata_', full_name, '_', int2str(session), '_', kernel, '_0.mat'];
            check_path(['../GLM_data/', full_name]);
            save(file_path, 'N', 'B', 'PS_kernels', 'conn_kernels', 'kernel_len', 'n_PS_kernel',...
                'n_conn_kernel', 'predjs_PS', 'predjs_conn', 'raster');

            % raster file and border file
            file_path = ['../GLM_data/', full_name, '/raster_', full_name, '_', int2str(session), '_0.mat'];
            save(file_path, 'n_trial', 'cell_id', 'cell_area', 'channel', 'session_name_full');
            file_path = ['../GLM_data/', full_name, '/borders_', full_name, '_', int2str(session), '.mat'];
            save(file_path, 'borders');

            fprintf('Done.\n');
            fprintf('CheckSum: %f\n', sum(raster, 'all'));
        end
    end
end
