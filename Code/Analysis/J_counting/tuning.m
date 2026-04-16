%% tuning.m - J difference based on tuning

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
% model config
config.align                 = 'Last';
config.shuffle_idx           = 0;
config.kernel_name           = 'DeltaPure';
config.reg_name              = 'L2=0_2';
config.epoch                 = 3000;
config.fold_idx              = 1;
config.resting_dur_threshold = 15;

kernel_num = 3;

% Load metadata
metadata = load_meta(root);
tuning_metas = metadata.tuning;
tuning_num = numel(tuning_metas);

for i = 1:tuning_num
    % load tuning
    meta = tuning_metas{i};
    file_folder = fullfile(root, 'Data', 'Working', 'tuning');
    file_name = meta.file_name;
    file_path = fullfile(file_folder, file_path);
    load(file_path, 'meta', 'data');

    N = meta.N;

    J_all = zeros(N, N, kernel_num, 2, 2); % (N, N, n_kernels, pre/post, state);
    prepost_str = {'pre', 'post'};
    state_str = {'RestOpen', 'RestClosed'};
    for prepost_idx = 1:numel(prepost_str)
        for state_idx = 1:numel(state_str)
            % construct model meta
            model_meta = meta;
            model_meta.prepost = prepost_str{prepost_idx};
            model_meta.state = state_str{state_idx};
            for field = fieldnames(config)'
                model_meta.(field{1}) = config.(field{1});
            end

            % get model filename
            model_file_name = generate_filename(model_meta, 'GLM');
            model_file_path = fullfile(root, 'Data', 'Working', 'GLM', model_file_name);

            % load model
            model_file = load(model_file_path);
            model_par = model_file.data.model_par;
            model_N = model_file.meta.N;
            n_PS_kernel = model_file.data.kernel.meta.n_PS_kernel;
            assert(model_N == N, 'Model N does not match meta N');

            % extract J matrix
            for kernel_idx = 1:kernel_num
                start_idx = n_PS_kernel + 2 + (kernel_idx-1)*N;
                end_idx = start_idx + N - 1;
                J_mat = model_par(:, start_idx:end_idx);
                J_all(:, :, kernel_idx, prepost_idx, state_idx) = J_mat;
            end
        end
    end


    % define comparason groups
    % tuning format: (N, pre/post, offer1/offer2/BeforeCho/Cho2Info/Inof2Rew/AfterRew, ER/IxU/Unc/Sub)
    % 1. choice to info +/-, IxU
    filter = data.tuning(:, :, 4, 2); 



