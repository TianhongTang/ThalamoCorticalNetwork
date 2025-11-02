function [loss, grad, err]=minuslogL_grad_hess_fun(par,B,N, ...
        n_PS_kernel,n_conn_kernel,raster,predjs_PS,predjs_conn,logfacts,reg)
% raster:(N, B), predj:(N, B, n_kernel)

chunk_size_B = 200000;

hi=par(:, 1); % (N, 1)
Pik=par(:, 2:(n_PS_kernel+1)); % (N, n_PS)
Jijk=reshape(par(:, (n_PS_kernel+2):end), N, N, n_conn_kernel); % (N, N, n_conn)

htot_it = repmat(hi, 1, B); % (N, B)

%%%% tensor notations: 
%%%% <A>=elementwise product, 
%%%% [A]=contracted tensor product

% conn kernels
for k=1:n_conn_kernel
    htot_it = htot_it + Jijk(:, :, k) * predjs_conn(:, :, k); % (N, [N])*([N], B)=(N, B)
end
% post-spike kernels
for k=1:n_PS_kernel
    htot_it = htot_it + diag(Pik(:, k)) * predjs_PS(:, :, k); % (<N>)*(<N>, B)=(N, B)
end

% link function for poisson distribution
lambda_it=exp(htot_it);  % (N, B)
% minuslogL= -sum(raster.*log(lambda_it) - log(1+lambda_it) ,'all');

L_it = -raster.*log(lambda_it)+log(1+lambda_it);            % -log(L_it) for Bernoulli
% L_it = lambda_it - (raster.*log(lambda_it)-logfacts);       % -log(L_it) for Poisson

minuslogL = sum(L_it, "all");
loss.minuslogL = minuslogL;

RegLoss = 0;
if reg.l1>0 
    RegLoss = RegLoss + reg.l1*sum(abs(par(:, 2:end)), "all");
end
if reg.l2>0
    RegLoss = RegLoss + reg.l2*sum(par(:, 2:end).^2, "all");
end
loss.reg = RegLoss;
loss.total = loss.minuslogL + RegLoss;

if nargout > 1     
    % fprintf("calc grad\n");
    % todo: fix the grad calc
    grad_minuslogL=zeros(N, 1 + n_PS_kernel + N*n_conn_kernel);
    
    dLdh_it = lambda_it./(1+lambda_it)-raster;  % d(minuslogL)/d(htot_it) for Bernoulli, (N, B)
    % dLdh_it = lambda_it-raster;                 % d(minuslogL)/d(htot_it) for Poisson, (N, B)

    dLdh_i = sum(dLdh_it, 2);% (N, [B])=(N, 1)
    dLdP_ik = sum(dLdh_it.*predjs_PS(:, :, :), 2);% (<N>, [B])*(<N>, [B], n_PS)=(N, n_PS)
    dLdJ_ijk = tensorprod(dLdh_it, predjs_conn, 2, 2); % (N, [B])*(N, [B], n_conn)=(N, N, n_conn)

    grad_minuslogL(:, 1) = dLdh_i;
    grad_minuslogL(:, 2:(n_PS_kernel+1)) = dLdP_ik;
    grad_minuslogL(:, (n_PS_kernel+2):end) = reshape(dLdJ_ijk, N, N*n_conn_kernel);

    % predj_conn_reshaped = permute(predjs_conn, [1, 3, 2]); % (N, n_conn, B)
    % predj_conn_reshaped = reshape(predj_conn_reshaped, 1, N*n_conn_kernel, B); % (1, N*n_conn, B)
    % predj_conn_reshaped = repmat(predj_conn_reshaped, N, 1, 1); % (N, N*n_conn, B)
    % grad_htot_it = [ones(N, B, 1), predjs_PS, predj_conn_reshaped]; % (N, 1 + n_PS + N*n_conn, B)

    % grad_minuslogL = sum(reshape(dLdh_it, N, 1, B).*grad_htot_it, 3); % (<N>, [B])*(<N>, 1 + n_PS + N*n_conn, [B])=(N, 1 + n_PS + N*n_conn)

    % eliminate self-connections in conn kernels
    for i=1:N
        for k=1:n_conn_kernel
            grad_minuslogL(i, i + n_PS_kernel + (k-1)*N + 1)=0;
        end
    end
    grad.minuslogL = grad_minuslogL;

    % grad_minuslogL(:, 1) = dLdh_i;
    % grad_minuslogL(:, 2:(n_PS_kernel+1)) = dLdP_ik;
    % grad_minuslogL(:, (n_PS_kernel+2):end) = reshape(dLdJ_ijk, N, N*n_conn_kernel);
    
    % regularizations
    % if reg.l1>0 
    %     grad_minuslogL = grad_minuslogL + [zeros(N, 1), reg.l1*sign(par(:, 2:end))];
    % end
    % if reg.l2>0
    %     grad_minuslogL = grad_minuslogL + [zeros(N, 1), reg.l2*2*par(:, 2:end)];
    % end
    % if reg.l1>0 
    %     grad_minuslogL = grad_minuslogL + [zeros(N, 1), (1/(reg.l1*B))*sign(par(:, 2:end))];
    % end
    % if reg.l2>0
    %     grad_minuslogL = grad_minuslogL + [zeros(N, 1), (1/(2*reg.l2*B))*2*par(:, 2:end)];
    % end
    grad_reg = zeros(N, 1 + n_PS_kernel + N*n_conn_kernel);
    if reg.l1>0 
        grad_reg = grad_reg + [zeros(N, 1), reg.l1*sign(par(:, 2:end))];
    end
    if reg.l2>0
        grad_reg = grad_reg + [zeros(N, 1), reg.l2*2*par(:, 2:end)];
    end
    
    % eliminate self-connections in conn kernels
    for i=1:N
        for k=1:n_conn_kernel
            grad_reg(i, i + n_PS_kernel + (k-1)*N + 1)=0;
        end
    end

    grad.reg = grad_reg;
    grad.total = grad_minuslogL + grad_reg;

end

% Hessian matrix error (N, 1 + n_PS_kernel + N*n_conn_kernel)
if nargout > 2 
    % fprintf("calc hess\n");
    err.minuslogL = zeros(N, 1 + n_PS_kernel + N*n_conn_kernel);
    err.total = zeros(N, 1 + n_PS_kernel + N*n_conn_kernel);
    % for i=1:N
    %     predj_PS_reshaped = permute(predjs_PS(i, :, :), [1, 3, 2]); % (1, n_PS, B)
    %     predj_PS_reshaped = reshape(predj_PS_reshaped, 1*n_PS_kernel, B); % (n_PS, B)
    %     predj_conn_reshaped = permute(predjs_conn, [1, 3, 2]); % (N, n_conn, B)
    %     predj_conn_reshaped = reshape(predj_conn_reshaped, N*n_conn_kernel, B); % (N*n_conn, B)
    %     grad_htot_it = [ones(1, B); predj_PS_reshaped; predj_conn_reshaped]; % (1 + n_PS + N*n_conn, B)
    %     hess = tensorprod(grad_htot_it, grad_htot_it, 2, 2); % (1 + n_PS + N*n_conn, 1 + n_PS + N*n_conn)
    %     err(i, :) = sqrt(diag(inv(hess)));
    % end
    for i=1:N
        n_chunk = ceil(B/chunk_size_B);
        hess_minuslogL_total = zeros(1 + n_PS_kernel + (N-1)*n_conn_kernel, 1 + n_PS_kernel + (N-1)*n_conn_kernel);
        hess_total_total = zeros(1 + n_PS_kernel + (N-1)*n_conn_kernel, 1 + n_PS_kernel + (N-1)*n_conn_kernel);

        for chunk_idx=1:n_chunk
            chunk_start = (chunk_idx-1)*chunk_size_B + 1;
            chunk_end = min(chunk_idx*chunk_size_B, B);
            chunk_size = chunk_end - chunk_start + 1;

            predj_PS_reshaped = permute(predjs_PS(i, chunk_start:chunk_end, :), [1, 3, 2]); % (1, n_PS, B)
            predj_PS_reshaped = reshape(predj_PS_reshaped, 1*n_PS_kernel, chunk_size); % (n_PS, B)
            predj_conn_reshaped = permute(predjs_conn(:, chunk_start:chunk_end, :), [1, 3, 2]); % (N, n_conn, B)
            predj_conn_reshaped = reshape(predj_conn_reshaped, N*n_conn_kernel, chunk_size); % (N*n_conn, B)
            % remove self-connections
            for k=1:n_conn_kernel
                predj_conn_reshaped(i + (k-1)*(N-1), :) = []; 
            end
        
            grad_htot_it = [ones(1, chunk_size); predj_PS_reshaped; predj_conn_reshaped]; % (1 + n_PS + (N-1)*n_conn, B)
            lambda_factor_it = lambda_it(i, chunk_start:chunk_end)./((1+lambda_it(i, chunk_start:chunk_end)).^2); % (1, B)
            % lambda_factor_it_reshaped = repmat(lambda_factor_it, 1 + n_PS_kernel + (N-1)*n_conn_kernel, 1); % (1 + n_PS + (N-1)*n_conn, B)
            % hess = tensorprod(grad_htot_it, grad_htot_it.*lambda_factor_it, 2, 2); % (1 + n_PS + (N-1)*n_conn, 1 + n_PS + (N-1)*n_conn)
            hess = (grad_htot_it.*lambda_factor_it) * grad_htot_it.';
            hess_minuslogL_total = hess_minuslogL_total + hess;

            % add hessian of the regularization
            hess = hess + diag([zeros(1, 1+n_PS_kernel), 2*reg.l2*ones(1, (N-1)*n_conn_kernel)]); % (1 + n_PS + (N-1)*n_conn, 1 + n_PS + (N-1)*n_conn)
            hess_total_total = hess_total_total + hess;
        end
        
        err_i = sqrt(diag(inv(hess_minuslogL_total))); % (1 + n_PS + (N-1)*n_conn)
        err_i_total = sqrt(diag(inv(hess_total_total))); % (1 + n_PS + (N-1)*n_conn)
        
        for k = 1:n_conn_kernel
            % put back zero to self-connections
            err_i = [err_i(1:(1+n_PS_kernel+(k-1)*N+i-1)); 0; err_i((1+n_PS_kernel+(k-1)*N+i):end)]; % (1 + n_PS + N*n_conn)
            err_i_total = [err_i_total(1:(1+n_PS_kernel+(k-1)*N+i-1)); 0; err_i_total((1+n_PS_kernel+(k-1)*N+i):end)]; % (1 + n_PS + N*n_conn)
        end
        err.minuslogL(i, :) = err_i;
        err.total(i, :) = err_i_total;
    end
end

