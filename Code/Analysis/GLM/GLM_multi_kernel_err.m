function GLM_multi_kernel_err(dataset_name, session, kernel_name, shuffle_seed, max_epoch, reg, log_level, lr)
%% GLM inference.
%%%% required input : (from convolution)
%%%% dataset: "../GLM_data/[dataset_name]/
%%%%   GLMdata_[dataset_name]_[session]_[kernel_name][_shuffle].mat"
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


foldername=['../GLM_data/', dataset_name];
data_filename=[foldername, '/GLMdata_', dataset_name, '_', ...
    int2str(session), '_', kernel_name, '_', int2str(shuffle_seed), '.mat'];
load(data_filename, "raster", "predjs_conn", "predjs_PS", ...
    "conn_kernels", "PS_kernels", ...
    "n_conn_kernel", "n_PS_kernel", "kernel_len", "N", "B"); 

% filter: ignore 0 firing rate neurons in inference
raster_filter = sum(raster, 2)>0;
% N_original = N;
% raster_original = raster;
N_filtered = sum(raster_filter);
raster = raster(raster_filter, :);

% temp solution, remove this.
if isempty(predjs_PS)
    predjs_PS = zeros([N, B, 0]);
end

predjs_PS = predjs_PS(raster_filter, :, :);
predjs_conn = predjs_conn(raster_filter, :, :);

% Initial condition (can be better):
% h ~ (:, 1)
% P_ik ~ (:, k+1)
% J_ijk ~ (:, (N*(k-1) + n_PS_kernel + 2):(N*k + n_PS_kernel + 1))
par0=zeros(N_filtered, 1 + n_PS_kernel + N_filtered*n_conn_kernel); 

% Adam solver
beta1 = 0.9;
beta2 = 0.999;
lr = lr * sqrt(B/16); % large batch payoff
e = 1e-8;

m = zeros(size(par0));
v = zeros(size(par0));
par = par0;

% GPU acceleration
raster = gpuArray(raster);
predjs_PS = gpuArray(predjs_PS);
predjs_conn = gpuArray(predjs_conn);
par = gpuArray(par);
facts = factorial(raster);
logfacts = log(facts);
    
beta1_t = 1;
beta2_t = 1;
fprintf("Ready\n");
for epoch=1:max_epoch
    if (mod(epoch, 100)==0 && log_level==1)||log_level==2
        [loss, grad, err] = minuslogL_grad_hess_fun(par,B,N_filtered, ...
            n_PS_kernel,n_conn_kernel,raster,predjs_PS,predjs_conn,logfacts,reg); 
    else
        [loss, grad] = minuslogL_grad_hess_fun(par,B,N_filtered, ...
            n_PS_kernel,n_conn_kernel,raster,predjs_PS,predjs_conn,logfacts,reg); 
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
        fprintf("Epoch %d/%d, -logL=%f, reg=%f, sparsity=%f\n", epoch, max_epoch, ...
            loss.minuslogL, loss.reg, sum(abs(par)>1e-2,"all")/numel(par));
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
            fprintf("Epoch %d/%d, -logL=%f, reg=%f, sparsity=%f\n", epoch, max_epoch, ...
                loss.minuslogL, loss.reg, sum(abs(par)>1e-2,"all")/numel(par));
            grad_norm_minuslogL = sum(grad.minuslogL.^2, "all");
            grad_norm_reg = sum(grad.reg.^2, "all");
            grad_norm_total = sum(grad.total.^2, "all");
            cos_angle = sum(grad.minuslogL.*grad.reg, "all")/sqrt(grad_norm_minuslogL*grad_norm_reg);
            fprintf("gradient: minuslogL=%f, reg=%f, total=%f, cos_angle=%f\n", ...
                grad_norm_minuslogL, grad_norm_reg, grad_norm_total, cos_angle);
        end
        foldername = ['../GLM_model/', dataset_name];
        check_path(foldername);

        model_path = [foldername, '/GLM_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_', int2str(shuffle_seed), '_', ...
        reg.name, '_', int2str(epoch), '.mat'];

        model_loss.minuslogL = gather(loss.minuslogL);
        model_loss.reg = gather(loss.reg);
        model_loss.total = gather(loss.total);
        model_par_filtered = gather(par);
        model_err_filtered.minuslogL = gather(err.minuslogL);
        model_err_filtered.total = gather(err.total);

        % calc predicted firing prob
        p = GLM_test(par,B,N_filtered, ...
            n_PS_kernel,n_conn_kernel,raster,predjs_PS,predjs_conn,logfacts,reg);
        p = gather(p);
        pred_prob = zeros(N, B); % reconstruct
        pred_prob(raster_filter, :) = p;

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

        % save(model_path, 'model_par', 'model_loss', 'model_err', 'N', "reg", "kernel_len", ...
        %     "PS_kernels", "conn_kernels", "n_PS_kernel", "n_conn_kernel", "raster_filter", ...
        %     "N_filtered", "pred_prob");
        fprintf('Saving to: %s\n', model_path);
        save(model_path, 'model_par', 'model_loss', 'model_err', 'N', "reg", "kernel_len", ...
            "PS_kernels", "conn_kernels", "n_PS_kernel", "n_conn_kernel", "raster_filter", ...
            "N_filtered");
        fprintf('saved');
    end
end

  


end