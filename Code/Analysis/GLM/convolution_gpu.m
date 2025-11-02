function convolution_gpu(dataset_name, session, kernel_name, shuffle_seed)
%% convolve raster and kernel(s) to predj.
% load data
raster_file = ['../GLM_data/', dataset_name, '/raster_', dataset_name, '_', ...
    int2str(session), '_', int2str(shuffle_seed), '.mat'];
kernel_file = ['../GLM_data/kernel_', kernel_name, '.mat'];
load(raster_file, "rasters", "n_trial", "trial_len");
load(kernel_file, "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len");

N = size(rasters{1}, 1);

% convolution
valid_trial = trial_len>kernel_len;
B = sum(trial_len(valid_trial)-(kernel_len-1));

% clip valid raster parts and concatenate
raster = zeros(N, B);
pointer = 1;
for t = 1:n_trial
    if ~valid_trial(t)
        continue
    end
    t_start = pointer;
    t_end = pointer + trial_len(t)-kernel_len;
    pointer = t_end+1;
    for i = 1:N
        raster(i, t_start:t_end) = rasters{t}(i, (kernel_len):end);
    end
end

% connection kernels 
predjs_conn = zeros(N, B, n_conn_kernel);
for k = 1:n_conn_kernel
    kernel = conn_kernels{k};
    kernel_gpu = gpuArray(kernel);
    predj = zeros(N, B);
    pointer = 1;
    for t = 1:n_trial
        if ~valid_trial(t)
            continue
        end
        t_start = pointer;
        t_end = pointer + trial_len(t)-kernel_len;
        pointer = t_end+1;
        for i = 1:N
            raster_gpu = gpuArray(rasters{t}(i, :));
            predj(i, t_start:t_end) = gather(conv(raster_gpu, kernel_gpu, "valid"));
        end
    end
    predjs_conn(:, :, k) = predj;
end

% post-spike kernels (same as conn kernels)
predjs_PS = zeros(N, B, n_PS_kernel);
for k = 1:n_PS_kernel
    kernel = PS_kernels{k};
    kernel_gpu = gpuArray(kernel);
    predj = zeros(N, B);
    pointer = 1;
    for t = 1:n_trial
        if ~valid_trial(t)
            continue
        end
        t_start = pointer;
        t_end = pointer + trial_len(t)-kernel_len;
        pointer = t_end+1;
        for i = 1:N
            raster_gpu = gpuArray(rasters{t}(i, :));
            predj(i, t_start:t_end) = gather(conv(raster_gpu, kernel_gpu, "valid"));
        end
    end
    predjs_PS(:, :, k) = predj;
end

% save data
save_file = ['../GLM_data/', dataset_name, '/GLMdata_', dataset_name,...
        '_', int2str(session), '_', kernel_name,  '_', int2str(shuffle_seed), '.mat'];
save(save_file, "raster", "predjs_conn", "predjs_PS", "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len", "N", "B");


end