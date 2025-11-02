function convolution(dataset_name, session, shuffle_id, kernel_name)
%% convolve raster folds and kernel(s) to predj.

% Prerequisite files: 
% split raster folds: from crossval_split.m 
% kernel: from generate_kernels.m

%% default parameters

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
% load data
folder_name = fullfile(root, 'Data', 'Working', 'crossval_split');
file_name = sprintf('crossval_%s_%d_%d.mat', dataset_name, session, shuffle_id);
raster_path = fullfile(folder_name, file_name);
folder_name = fullfile(root, 'Data', 'Working', 'kernel');
file_name = sprintf('kernel_%s.mat', kernel_name);
kernel_path = fullfile(folder_name, file_name);
load(raster_path, "N", "fold_num", "fold_rasters", "fold_trial_lens");
load(kernel_path, "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");

folds = cell(1, fold_num);
for fold_id = 1:fold_num
    fold = struct();
    rasters = fold_rasters{fold_id};
    trial_num = length(rasters);
    trial_len = fold_trial_lens{fold_id};

    valid_trial = trial_len>kernel_len;
    B = sum(trial_len(valid_trial)-(kernel_len-1));

    % clip valid raster parts and concatenate
    raster = zeros(N, B);
    pointer = 1;
    for t = 1:trial_num
        if ~valid_trial(t)
            continue
        end
        t_start = pointer;
        t_end = pointer + trial_len(t)-kernel_len;
        pointer = t_end+1;
        for i = 1:N
            raster(i, t_start:t_end) = rasters{t}(i, kernel_len:end);
        end
    end

    % calculate predj
    predjs_conn = zeros(N, B, n_conn_kernel);
    % connection kernels 
    for k = 1:n_conn_kernel
        kernel = conn_kernels{k};
        predj = zeros(N, B);
        pointer = 1;
        for t = 1:trial_num
            if ~valid_trial(t)
                continue
            end
            t_start = pointer;
            t_end = pointer + trial_len(t)-kernel_len;
            pointer = t_end+1;
            for i = 1:N
                predj(i, t_start:t_end) = conv(rasters{t}(i, :), kernel, "valid");
            end
        end
        predjs_conn(:, :, k) = predj;
    end

    % post-spike kernels (same as conn kernels)
    predjs_PS = zeros(N, B, n_PS_kernel);
    for k = 1:n_PS_kernel
        kernel = PS_kernels{k};
        predj = zeros(N, B);
        pointer = 1;
        for t = 1:trial_num
            if ~valid_trial(t)
                continue
            end
            t_start = pointer;
            t_end = pointer + trial_len(t)-kernel_len;
            pointer = t_end+1;
            for i = 1:N
                predj(i, t_start:t_end) = conv(rasters{t}(i, :), kernel, "valid");
            end
        end
        predjs_PS(:, :, k) = predj;
    end

    fold.raster = raster;
    fold.B = B;
    fold.predjs_conn = predjs_conn;
    fold.predjs_PS = predjs_PS;
    folds{fold_id} = fold;  
end

% pack kernel info
kernel = struct();
kernel.name = kernel_name;
kernel.conn_kernels = conn_kernels;
kernel.PS_kernels = PS_kernels;
kernel.n_conn_kernel = n_conn_kernel;
kernel.n_PS_kernel = n_PS_kernel;
kernel.kernel_len = kernel_len;

% save data
save_folder = fullfile(root, 'Data', 'Working', 'GLM_data');
check_path(save_folder);
save_name = sprintf('GLMdata_%s_%d_%d_%s.mat', dataset_name, session, shuffle_id, kernel_name);
save_path = fullfile(save_folder, save_name);
save(save_path, "N", "fold_num", "folds", "kernel", '-v7.3');

end