%% generate_test.m - Generate test data for J counting analysis
% 

clear;
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

prepost_strs = {'Pre', 'Post'};
state_strs = {'RestOpen', 'RestClose'};

%% Main
N = 30;
borders = [1, 11, 21];
session_type = 'Test';
area_type = 'Test';
align = 'AlignLast15';
kernel = 'DeltaPure';
reg = 'L2=0_2';
epoch = '3000';


for prepost_idx = 1:length(prepost_strs)
    prepost_str = prepost_strs{prepost_idx};
    for state_idx = 1:length(state_strs)
        state_str = state_strs{state_idx};
        J = zeros(N, N);

        J(1:2, 6:10) = 1; % 10 pos in area 1
        J(11:12, 16:20) = 1; % 10 pos in area 2
        J(21:23, 26:30) = 1; % 15 pos in area 3
        J(26:30, 21:23) = 1; % 15 neg in area 3

        J(1:4, 11:15) = 1; % 20 pos in area 1-2
        J(5:6, 11:15) = -1; % 10 neg in area 1-2
        J(11:15, 1:4) = 1; % 20 pos in area 2-1
        J(11:15, 5:6) = -1; % 10 neg in area 2-1

        J(1:3, 21:25) = 1; % 15 pos in area 1-3
        J(4:6, 21:25) = -1; % 15 neg in area 1-3
        J(11:11, 21:25) = 1; % 5 pos in area 2-3
        J(12:12, 21:25) = -1; % 5 neg in area 2-3
        J(21:25, 11:11) = 1; % 5 pos in area 3-2
        J(21:25, 12:12) = -1; % 5 neg in area 3-2

        if strcmp(state_str, 'RestClose')
            J(7:8, 16:20) = 1; % 10 pos in area 1-2
        end

        if strcmp(prepost_str, 'Post')
            J(9:10, 11:15) = 1; % 10 pos in area 1-2
            J(9:10, 16:20) = -1; % 10 neg in area 1-2
        end

        model_par = zeros(N, 3*N+1);
        model_par(:, 2:N+1) = J;
        model_par(:, N+2:2*N+1) = J;
        model_par(:, 2*N+2:3*N+1) = J;
        model_err.total = ones(N, 3*N+1) * 0.5;
    
        for session_idx = 1:2
            border_folder = fullfile(root, 'Data', 'Working', 'border');
            check_path(border_folder);
            border_name = sprintf('borders_%s%s_%d.mat', session_type, area_type, session_idx);
            border_path = fullfile(border_folder, border_name);
            save(border_path, 'borders');
    
            folder_name = fullfile(root, 'Data', 'Working', 'GLM_models');
            check_path(folder_name);
            file_name = sprintf('GLM_%s%s%s%s%s_s%d_shuffle0_%s_%s_epoch%s_fold0.mat', ...
                session_type, prepost_str, state_str, area_type, align, session_idx, kernel, reg, epoch); 
            file_path = fullfile(folder_name, file_name);
            save(file_path, 'model_par', 'model_err', 'N');
        end

        % plot and save J matrix. clim = [-1, 1]
        figure;
        imagesc(J, [-1, 1]);
        title(sprintf('J matrix - %s %s', prepost_str, state_str));
        colorbar;
        plot_folder = fullfile(root, 'Data', 'Working', 'J_counting');
        check_path(plot_folder);
        plot_name = sprintf('J_test_%s%s_%s.png', prepost_str, state_str);
        plot_path = fullfile(plot_folder, plot_name);
        saveas(gcf, plot_path);
    end
end

%% generate_test.m - Generate test data for J counting analysis
% 

clear;
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

prepost_strs = {'Pre', 'Post'};
state_strs = {'RestOpen', 'RestClose'};

%% Main
N = 30;
borders = [1, 11, 21];
session_type = 'Test';
area_type = 'Test';
align = 'AlignLast15';
kernel = 'DeltaPure';
reg = 'L2=0_2';
epoch = '3000';


for prepost_idx = 1:length(prepost_strs)
    prepost_str = prepost_strs{prepost_idx};
    for state_idx = 1:length(state_strs)
        state_str = state_strs{state_idx};
        J = zeros(N, N);

        J(1:2, 6:10) = 1; % 10 pos in area 1
        J(11:12, 16:20) = 1; % 10 pos in area 2
        J(21:23, 26:30) = 1; % 15 pos in area 3
        J(26:30, 21:23) = 1; % 15 neg in area 3

        J(1:4, 11:15) = 1; % 20 pos in area 1-2
        J(5:6, 11:15) = -1; % 10 neg in area 1-2
        J(11:15, 1:4) = 1; % 20 pos in area 2-1
        J(11:15, 5:6) = -1; % 10 neg in area 2-1

        J(1:3, 21:25) = 1; % 15 pos in area 1-3
        J(4:6, 21:25) = -1; % 15 neg in area 1-3
        J(11:11, 21:25) = 1; % 5 pos in area 2-3
        J(12:12, 21:25) = -1; % 5 neg in area 2-3
        J(21:25, 11:11) = 1; % 5 pos in area 3-2
        J(21:25, 12:12) = -1; % 5 neg in area 3-2

        if strcmp(state_str, 'RestClose')
            J(7:8, 16:20) = 1; % 10 pos in area 1-2
        end

        if strcmp(prepost_str, 'Post')
            J(9:10, 11:15) = 1; % 10 pos in area 1-2
            J(9:10, 16:20) = -1; % 10 neg in area 1-2
        end

        model_par = zeros(N, 3*N+1);
        model_par(:, 2:N+1) = J;
        model_par(:, N+2:2*N+1) = J;
        model_par(:, 2*N+2:3*N+1) = J;
        model_err.total = ones(N, 3*N+1) * 0.5;
    
        for session_idx = 1:2
            border_folder = fullfile(root, 'Data', 'Working', 'border');
            check_path(border_folder);
            border_name = sprintf('borders_%s%s_%d.mat', session_type, area_type, session_idx);
            border_path = fullfile(border_folder, border_name);
            save(border_path, 'borders');
    
            folder_name = fullfile(root, 'Data', 'Working', 'GLM_models');
            check_path(folder_name);
            file_name = sprintf('GLM_%s%s%s%s%s_s%d_shuffle0_%s_%s_epoch%s_fold0.mat', ...
                session_type, prepost_str, state_str, area_type, align, session_idx, kernel, reg, epoch); 
            file_path = fullfile(folder_name, file_name);
            save(file_path, 'model_par', 'model_err', 'N');
        end

        % plot and save J matrix. clim = [-1, 1]
        figure;
        imagesc(J, [-1, 1]);
        title(sprintf('J matrix - %s %s', prepost_str, state_str));
        colorbar;
        plot_folder = fullfile(root, 'Data', 'Working', 'J_counting');
        check_path(plot_folder);
        plot_name = sprintf('J_test_%s%s_%s.png', prepost_str, state_str);
        plot_path = fullfile(plot_folder, plot_name);
        saveas(gcf, plot_path);
    end
end
