function GLM_multi_kernel_crossval(dataset_name, session, kernel_name, shuffle_id, max_epoch, reg, log_level, lr, test_fold)
%% GLM inference with cross-validation
% test_fold: which fold to use as test set. If 0, use all data for training and testing.
%%%% required input : (from convolution)
%%%% dataset: "../GLM_data/[dataset_name]/
%%%%   GLMdata_[dataset_name]_[session]_[shuffle_id]_[kernel_name].mat"

%% default parameters
if nargin < 7
    log_level=2;
end
if nargin < 6
    reg.l1=0;
    reg.l2=0;
    reg.name='None';
end
if nargin < 5
    max_epoch=1000;
end
if nargin < 8
    lr=1e-3;
end
if nargin < 9
    test_fold=1;
end

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
folder_name= fullfile(root, 'Data', 'Working', 'GLM_data');
file_name = sprintf('GLMdata_%s_%d_%d_%s.mat', dataset_name, session, shuffle_id, kernel_name);
file_path = fullfile(folder_name, file_name);
load(file_path, "N", "fold_num", "folds", "kernel");
conn_kernels = kernel.conn_kernels;
PS_kernels = kernel.PS_kernels;
n_conn_kernel = kernel.n_conn_kernel;
n_PS_kernel = kernel.n_PS_kernel;
kernel_len = kernel.kernel_len;

% concatenate training folds
raster = [];
predjs_conn = [];
predjs_PS = [];
B_train = 0;
for fold_id = 1:fold_num
    if fold_id == test_fold
        continue
    end
    raster = cat(2, raster, folds{fold_id}.raster);
    predjs_conn = cat(2, predjs_conn, folds{fold_id}.predjs_conn);
    predjs_PS = cat(2, predjs_PS, folds{fold_id}.predjs_PS);
    B_train = B_train + folds{fold_id}.B;
end

% test set
if test_fold == 0
    raster_test = raster;
    predjs_conn_test = predjs_conn;
    predjs_PS_test = predjs_PS;
    B_test = B_train;
else
    raster_test = folds{test_fold}.raster;
    predjs_conn_test = folds{test_fold}.predjs_conn;
    predjs_PS_test = folds{test_fold}.predjs_PS;
    B_test = folds{test_fold}.B;
end 

% filter: ignore 0 firing rate neurons in inference
raster_filter = sum(raster, 2)>0;
% N_original = N;
% raster_original = raster;
N_filtered = sum(raster_filter);
raster = raster(raster_filter, :);
predjs_PS = predjs_PS(raster_filter, :, :);
predjs_conn = predjs_conn(raster_filter, :, :);
raster_test = raster_test(raster_filter, :);
predjs_PS_test = predjs_PS_test(raster_filter, :, :);
predjs_conn_test = predjs_conn_test(raster_filter, :, :);

% Initial condition (can be better):
par0=zeros(N_filtered, 1 + n_PS_kernel + N_filtered*n_conn_kernel); 

% Adam solver
beta1 = 0.9;
beta2 = 0.999;
lr = lr * sqrt(B_train/16); % large batch payoff
e = 1e-8;

m = zeros(size(par0));
v = zeros(size(par0));
par = par0;

% GPU acceleration
raster = gpuArray(raster);
predjs_PS = gpuArray(predjs_PS);
predjs_conn = gpuArray(predjs_conn);
raster_test = gpuArray(raster_test);
predjs_PS_test = gpuArray(predjs_PS_test);
predjs_conn_test = gpuArray(predjs_conn_test);
par = gpuArray(par);
facts = factorial(raster);
logfacts = log(facts);
    
beta1_t = 1;
beta2_t = 1;
fprintf("Ready\n");
for epoch=1:max_epoch
    if (mod(epoch, 100)==0 && log_level==1)||log_level==2
        [loss, grad, err] = minuslogL_grad_hess_fun(par,B_train,N_filtered, ...
            n_PS_kernel,n_conn_kernel,raster,predjs_PS,predjs_conn,logfacts,reg); 
        loss_test = minuslogL_grad_hess_fun(par,B_test,N_filtered, ...
            n_PS_kernel,n_conn_kernel,raster_test,predjs_PS_test,predjs_conn_test,logfacts,reg);
        normalized_loss = loss.minuslogL / (B_train * N_filtered);
        normalized_loss_test = loss_test.minuslogL / (B_test * N_filtered);
    else
        [loss, grad] = minuslogL_grad_hess_fun(par,B_train,N_filtered, ...
            n_PS_kernel,n_conn_kernel,raster,predjs_PS,predjs_conn,logfacts,reg); 
        loss_test = NaN;
        normalized_loss = loss.minuslogL / (B_train * N_filtered);
        normalized_loss_test = NaN;
    end
    
    % Adam
    beta1_t = beta1_t * beta1;
    beta2_t = beta2_t * beta2;
    m = beta1*m + (1-beta1)*grad.total;
    v = beta2*v + (1-beta2)*(grad.total).^2;
    mh = m/(1-beta1_t);
    vh = v/(1-beta2_t);
    % grad update
    par = par - lr*mh./(sqrt(vh+e));
    
    if log_level==2
        fprintf("Epoch %d/%d, normalized train L=%f, normalized test L=%f, reg=%f, sparsity=%f\n", epoch, max_epoch, ...
            normalized_loss, normalized_loss_test, loss.reg, sum(abs(par)>1e-2,"all")/numel(par));
        grad_norm_minuslogL = sum(grad.minuslogL.^2, "all");
        grad_norm_reg = sum(grad.reg.^2, "all");
        grad_norm_total = sum(grad.total.^2, "all");
        cos_angle = sum(grad.minuslogL.*grad.reg, "all")/sqrt(grad_norm_minuslogL*grad_norm_reg);
        fprintf("gradient: minuslogL=%f, reg=%f, total=%f, cos_angle=%f\n", ...
            grad_norm_minuslogL, grad_norm_reg, grad_norm_total, cos_angle);
    end
    
    % save model
    if mod(epoch, 100)==0
        if log_level==1
            fprintf("Epoch %d/%d, normalized train L=%f, normalized test L=%f, reg=%f, sparsity=%f\n", epoch, max_epoch, ...
                normalized_loss, normalized_loss_test, loss.reg, sum(abs(par)>1e-2,"all")/numel(par));
            grad_norm_minuslogL = sum(grad.minuslogL.^2, "all");
            grad_norm_reg = sum(grad.reg.^2, "all");
            grad_norm_total = sum(grad.total.^2, "all");
            cos_angle = sum(grad.minuslogL.*grad.reg, "all")/sqrt(grad_norm_minuslogL*grad_norm_reg);
            fprintf("gradient: minuslogL=%f, reg=%f, total=%f, cos_angle=%f\n", ...
                grad_norm_minuslogL, grad_norm_reg, grad_norm_total, cos_angle);
        end

        train_loss.minuslogL = gather(loss.minuslogL);
        train_loss.reg = gather(loss.reg);
        train_loss.total = gather(loss.total);
        test_loss.minuslogL = gather(loss_test.minuslogL);
        test_loss.reg = gather(loss_test.reg);
        test_loss.total = gather(loss_test.total);
        model_par_filtered = gather(par);
        model_err_filtered.minuslogL = gather(err.minuslogL);
        model_err_filtered.total = gather(err.total);

        % model_err_selected = model_err_filtered.minuslogL;
        % model_err_selected = model_err_filtered.total;

        % reconstruct model_par and model_err 
        model_par = zeros(N, 1 + n_PS_kernel + N*n_conn_kernel)*NaN;
        model_par(raster_filter, 1:(1 + n_PS_kernel)) = model_par_filtered(:, 1:(1 + n_PS_kernel));

        model_err.minuslogL = zeros(N, 1 + n_PS_kernel + N*n_conn_kernel)*NaN;
        model_err.total = zeros(N, 1 + n_PS_kernel + N*n_conn_kernel)*NaN;

        model_err.minuslogL(raster_filter, 1:(1 + n_PS_kernel)) = model_err_filtered.minuslogL(:, 1:(1 + n_PS_kernel));
        model_err.total(raster_filter, 1:(1 + n_PS_kernel)) = model_err_filtered.total(:, 1:(1 + n_PS_kernel));

        for i=1:n_conn_kernel
            model_par(raster_filter, [false(1 + n_PS_kernel + (i-1)*N, 1); raster_filter; false(N*(n_conn_kernel-i), 1)]) = ...
                model_par_filtered(:, 1 + n_PS_kernel + (i-1)*N_filtered + (1:N_filtered));
            model_err.minuslogL(raster_filter, [false(1 + n_PS_kernel + (i-1)*N, 1); raster_filter; false(N*(n_conn_kernel-i), 1)]) = ...
                model_err_filtered.minuslogL(:, 1 + n_PS_kernel + (i-1)*N_filtered + (1:N_filtered));
            model_err.total(raster_filter, [false(1 + n_PS_kernel + (i-1)*N, 1); raster_filter; false(N*(n_conn_kernel-i), 1)]) = ...
                model_err_filtered.total(:, 1 + n_PS_kernel + (i-1)*N_filtered + (1:N_filtered));
            % model_err.total_plus(raster_filter, [false(1 + n_PS_kernel + (i-1)*N, 1); raster_filter; false(N*(n_conn_kernel-i), 1)]) = ...
            %     model_err_filtered.total_plus(:, 1 + n_PS_kernel + (i-1)*N_filtered + (1:N_filtered));
            % model_err.total_minus(raster_filter, [false(1 + n_PS_kernel + (i-1)*N, 1); raster_filter; false(N*(n_conn_kernel-i), 1)]) = ...
            %     model_err_filtered.total_minus(:, 1 + n_PS_kernel + (i-1)*N_filtered + (1:N_filtered));
        end

        folder_name = fullfile(root, 'Data', 'Working', 'GLM_models');
        save_name = sprintf('GLM_%s_s%d_shuffle%d_%s_%s_epoch%d_fold%d.mat', dataset_name, session, shuffle_id, kernel_name, ...
            reg.name, epoch, test_fold);
        check_path(folder_name);
        model_path = fullfile(folder_name, save_name);
        fprintf('Saving to: %s\n', model_path);
        save(model_path, 'model_par', 'train_loss', 'test_loss', 'model_err', 'N', "reg", "kernel", ...
           "raster_filter", "N_filtered");
        fprintf('saved\n');
    end
end