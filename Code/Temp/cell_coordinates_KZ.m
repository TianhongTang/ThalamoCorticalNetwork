%% cell_coordinates_KZ.m - Visualize KZ dataset cells in 3D space

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

target_root = fullfile(root, 'Copy');
check_path(target_root);

%% Main
% load cell coordinates
meta_folder = fullfile(root, 'Data', 'Working', 'Meta');
check_path(meta_folder);
meta_file_name = 'all_session_info_KZ.mat';
meta_file_path = fullfile(meta_folder, meta_file_name); 
load(meta_file_path, 'all_session_info', 'segmentNames');

thal_all = zeros(0, 3); % x, y, z
thal_ant = zeros(0, 3);
thal_post = zeros(0, 3);

%% 
session_num = length(all_session_info);
area_num = length(segmentNames);
fprintf('Total number of sessions: %d\n', session_num);
for session_idx = 1:session_num
    session_info = all_session_info(session_idx);
    session_name = session_info.sessionname;
    monkey_name = session_info.monkeyName;
    if ~strcmp(monkey_name, 'Lemmy')
        continue;
    end

    neuron_info = session_info.neuronList;
    
    cell_area = {neuron_info.NeuralTargetsAnatomy};
    cell_coordinates = cat(1, neuron_info.fcsvCoordinates); % N x 3

    % select thalamus neurons
    thal_filter_all = strcmp(cell_area, 'Thalamus');
    thal_ant_filter = thal_filter_all;
    thal_post_filter = thal_filter_all;
    for i = 1:length(cell_area)
        if thal_filter_all(i)
            coord = cell_coordinates(i, :);
            if coord(3) < -18 || coord(2) < 37
                thal_ant_filter(i) = false;
            else
                thal_post_filter(i) = false;
            end
        end
    end

    thal_all = [thal_all; cell_coordinates(thal_filter_all, :)];
    thal_ant = [thal_ant; cell_coordinates(thal_ant_filter, :)];
    thal_post = [thal_post; cell_coordinates(thal_post_filter, :)];
end

% plot
figure;
scatter3(thal_ant(:, 1), thal_ant(:, 2), thal_ant(:, 3), 10, 'r', 'filled');
hold on;    
scatter3(thal_post(:, 1), thal_post(:, 2), thal_post(:, 3), 10, 'b', 'filled');
xlabel('X');
ylabel('P-A');
zlabel('Z');
title('Thalamus Neuron Coordinates (KZ Dataset)');
legend('Thalamus Anterior', 'Thalamus Posterior');
grid on;