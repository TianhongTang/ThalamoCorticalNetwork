function shuffle(dataset_name, session, shuffle_id, shuffle_seed, shuffle_type)
%% shuffle spike trains and keep the firing rate.
% shuffle_type: "None", "Within trial", "Across trial"

%% Get root folder
code_depth = 4;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
% include code folder and utils
% addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Main
% load original raster
file_folder = fullfile(root, 'Data', 'Working', 'raster');
file_name = sprintf('raster_%s_%d.mat', dataset_name, session);
raster_file = fullfile(file_folder, file_name);
load(raster_file, "rasters", "n_trial", "trial_len", "firing_rates");
n_raster = length(rasters);

rasters_shuffle = cell(1, n_raster);
rng(shuffle_seed);

switch shuffle_type
    case "None"
        rasters_shuffle = rasters;
    case "Within trial"
        % shuffle within trial
        for i =1:n_raster
            raster = rasters{i};
            [N, B] = size(raster);
            raster_shuffle = zeros(N, B);
            for j = 1:N
                shuffle_idx = randperm(B);
                raster_shuffle(j, :) = raster(j, shuffle_idx);
            end
            rasters_shuffle{i} = raster_shuffle;
        end
    case "Across trial"
        % shuffle across trial
        [N, ~] = size(rasters{1});
        full_raster = zeros(N, 0);
        for i=1:n_raster
            full_raster = [full_raster, rasters{i}];
        end
        [~, B] = size(full_raster);
        raster_shuffle = zeros(N, B);
        for j = 1:N
            shuffle_idx = randperm(B);
            raster_shuffle(j, :) = raster_shuffle(j, shuffle_idx);
        end
        pointer=0;
        for i=1:n_raster
            Bi = size(rasters{i}, 2);
            rasters_shuffle{i} = raster_shuffle(:, pointer+1:pointer+Bi);
            pointer = pointer+Bi;
        end
end

save_folder = fullfile(root, 'Data', 'Working', 'raster');
save_name = sprintf('shuffled_%s_%d_%d.mat', dataset_name, session, shuffle_id);
raster_file_shuffle = fullfile(save_folder, save_name);
rasters = rasters_shuffle;
save(raster_file_shuffle, "N", "n_trial", "trial_len", "rasters", "firing_rates", "shuffle_type", "shuffle_seed", '-v7.3');
end