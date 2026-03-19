function shuffle(raster_meta, shuffle_id, shuffle_seed, shuffle_type)
%% shuffle spike trains and keep the firing rate.
% shuffle_type: "None", "Within trial", "Across trial"

% required input files:
% - raster file: 'raster_%s_%d.mat' (dataset_name, session)
% output files:
% - shuffled raster file: 'shuffled_%s_%d_%d.mat' (dataset_name, session, shuffle_id)

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
file_name = generate_filename('raster', raster_meta);
raster_file = fullfile(file_folder, file_name);
% load(raster_file, "N", "rasters", "trial_num", "trial_len", "firing_rates");
raster_file = load(raster_file, "meta", "data");

N            = raster_file.meta.N;
trial_num    = raster_file.meta.trial_num;
rasters      = raster_file.data.rasters;
trial_len    = raster_file.data.trial_len;
n_raster     = length(rasters);

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

% recalculate firing rates
firing_rates = cell(1, n_raster);
for i = 1:n_raster
    firing_rates{i} = mean(rasters_shuffle{i}, 2)/raster_file.meta.dt; % in Hz
end

% construct meta and data for saving
meta = struct();
meta.animal_name  = raster_file.meta.animal_name;
meta.injection    = raster_file.meta.injection;
meta.prepost      = raster_file.meta.prepost;
meta.state        = raster_file.meta.state;
meta.area         = raster_file.meta.area;
meta.align        = raster_file.meta.align;
meta.session_idx  = raster_file.meta.session_idx;
meta.shuffle_mode = shuffle_type;
meta.shuffle_idx   = shuffle_id;  % not a typo.
meta.shuffle_seed = shuffle_seed;
meta.file_name    = generate_filename('shuffled', meta);
meta.N            = N;
meta.dt           = raster_file.meta.dt;
meta.trial_num    = trial_num;

data = struct();
data.rasters = rasters_shuffle;
data.trial_len = trial_len;
data.firing_rates = firing_rates;

save_folder = fullfile(root, 'Data', 'Working', 'shuffled');
save_name = meta.file_name;
raster_file_shuffle = fullfile(save_folder, save_name);
save(raster_file_shuffle, 'meta', 'data');
fprintf('Shuffled raster saved to: %s\n', raster_file_shuffle);
end