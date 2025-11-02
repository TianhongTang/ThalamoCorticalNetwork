function plot_GLM_sorted(dataset_name, session, kernel_name, epoch, reg, shuffle_size, sorting, idx)
model_path_ori = ['../GLM_model/', dataset_name, '/GLM_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_0_', ...
        reg.name, '_', int2str(epoch), '.mat'];
load(model_path_ori, "model_par", "model_err", "PS_kernels", "conn_kernels", "n_PS_kernel", "n_conn_kernel", "kernel_len", "N");

data_path = ['../GLM_data/', dataset_name, '/GLMdata_', dataset_name, '_', ...
        int2str(session),'_', kernel_name, '_0.mat'];
load(data_path, "raster");
firing_rate = mean(raster, 2);

if session<100
    session_border = session;
else
    session_border = floor(session/100);
end

load(['../GLM_data/', dataset_name,'/borders_', dataset_name, '_', ...
        int2str(session_border),'.mat'], "borders");

% for acc-thalamus-vlpfc dataset
A_T_border = borders(1);
T_P_border = borders(2);
    
par_ori = model_par;

par_sfl = zeros([size(par_ori), shuffle_size]);
for i=1:shuffle_size
    model_path_sfl = ['../GLM_model/', dataset_name, '/GLM_', dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_', int2str(i), '_', ...
        reg.name, '_', int2str(epoch), '.mat'];
    par_sfl(:, :, i) = load(model_path_sfl).model_par;
end

% significant = zeros(size(par_ori));

% for i=1:N
%     for j=1:(1+n_PS_kernel+N*n_conn_kernel)
%         % one-sample t-test
%         h = ttest(reshape(par_sfl(i, j, :)-par_ori(i, j), [1, shuffle_size]), ...
%             0, "Alpha", 0.001);
%         significant(i, j) = h;
%         % permutation test (doesn't work!)
%         % p = permutationTest(par_ori(i, j), reshape(par_sfl(i, j, :), [1, shuffle_size]), 0, 'exact', 1);
%         % significant(i, j) = p<0.05;

%     end
% end

significant = abs(par_ori) > model_err.total*2;

% for k = 1:3
%     sig_mat = significant(:, (N*(k-1) + n_PS_kernel + 2):(N*k + n_PS_kernel + 1));
%     J_mat = par_ori(:, (N*(k-1) + n_PS_kernel + 2):(N*k + n_PS_kernel + 1));
%     fprintf("dataset: %s, session: %d, kernel: %s\n", dataset_name, session, int2str(k));
%     within_pos = sum(sig_mat(1:30, 1:30) & J_mat(1:30, 1:30) > 0, "all") + ...
%         sum(sig_mat(31:60, 31:60) & J_mat(31:60, 31:60) > 0, "all");
%     within_neg = sum(sig_mat(1:30, 1:30) & J_mat(1:30, 1:30) < 0, "all") + ...
%         sum(sig_mat(31:60, 31:60) & J_mat(31:60, 31:60) < 0, "all");
%     across_pos = sum(sig_mat(1:30, 31:60) & J_mat(1:30, 31:60) > 0, "all") + ...
%         sum(sig_mat(31:60, 1:30) & J_mat(31:60, 1:30) > 0, "all");
%     across_neg = sum(sig_mat(1:30, 31:60) & J_mat(1:30, 31:60) < 0, "all") + ...
%         sum(sig_mat(31:60, 1:30) & J_mat(31:60, 1:30) < 0, "all");
%     fprintf("within_pos: %d, within_neg: %d, across_pos: %d, across_neg: %d\n", ...
%         within_pos, within_neg, across_pos, across_neg);
% end

par_sig = par_ori;
% par_sig = (par_sig - mean(par_sfl, 3)).*significant;
par_sig = par_sig.*significant;
% save("par_sig.mat", "par_sig");
% par_sig(~significant) = nan;

% clim for J_ijk plottings
clim_ori = max(abs(par_ori(:, (n_PS_kernel+2):end)), [], "all");
clim_sfl = max(abs(par_sfl(:, (n_PS_kernel+2):end, :)), [], "all");
clim_all = max(clim_ori, clim_sfl);
clim_all = 2;

%% sorting
sorting_ranges = [1, A_T_border-0.5;A_T_border+0.5,T_P_border-0.5;T_P_border+0.5, N];
criterion = 1;
sort_idx = zeros(1, N);

criterion_mat = par_sig(:, (N*(criterion-1) + n_PS_kernel + 2):(N*criterion + n_PS_kernel + 1));
criterion_mat(isnan(criterion_mat))=0;
% get sorting index
for i=1:3
    sort_s = sorting_ranges(i, 1);
    sort_e = sorting_ranges(i, 2);

    if sorting == "count"
        s = criterion_mat~=0;
        inward = sum(s, 2).';
        outward = sum(s, 1);
        v = inward+outward;
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e),'descend');
    elseif sorting == "abs"
        s = abs(criterion_mat);
        inward = sum(s, 2).';
        outward = sum(s, 1);
        v = inward+outward;
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e), 'descend');
    elseif sorting == "square"
        s = criterion_mat.^2;
        inward = sum(s, 2).';
        outward = sum(s, 1);
        v = inward+outward;
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e), 'descend');
    elseif sorting == "max"
        s = abs(criterion_mat);
        inward = max(s, [], 2).';
        outward = max(s, [], 1);
        v = max(inward,outward);
        [~, sort_idx(sort_s:sort_e)] = sort(v(sort_s:sort_e), 'descend');
    elseif sorting == "idx"
        [~, sort_idx(sort_s:sort_e)] = sort(idx(sort_s:sort_e), 'ascend');
    end
    sort_idx(sort_s:sort_e) = sort_idx(sort_s:sort_e) + sort_s-1;
end


%% sorting
sort_s = sorting_ranges(i, 1);
sort_e = sorting_ranges(i, 2);
par_ori(:, :) = par_ori(sort_idx, :);
par_sig(:, :) = par_sig(sort_idx, :);
par_sfl(:, :, :) = par_sfl(sort_idx, :, :);
for j=1:n_conn_kernel
    shift = N*(j-1) + n_PS_kernel + 1;
    par_ori(:, shift+1:shift+N) = par_ori(:, shift+sort_idx);
    par_sig(:, shift+1:shift+N) = par_sig(:, shift+sort_idx);
    par_sfl(:, shift+1:shift+N, :) = par_sfl(:, shift+sort_idx, :);
end
firing_rate = firing_rate(sort_idx);

%% save sorting index
sort_path = ['../GLM_data/', dataset_name, '/sortidx_', dataset_name, '_', ...
        int2str(session),'_', kernel_name, '.mat'];
save(sort_path, "sort_idx");


%% plotting
% very large figure show all kernels/models
fig = figure("Visible", "off");
set(fig, 'PaperPosition', [0, 0, (shuffle_size+3)*6+2, (1+n_conn_kernel+n_PS_kernel)*8+2]);
tiles = tiledlayout(1+n_conn_kernel+n_PS_kernel, shuffle_size+3);
cmap=brewermap(256,'*RdBu');

% y: kernel, ori, significant, shuffle1, ..., shuffleN
% x: h, conn_kernels, PS_kernels
for plot_x=1:1+n_conn_kernel+n_PS_kernel
    for plot_y = 1:shuffle_size+3

        t = nexttile((plot_x-1)*(shuffle_size+3) + plot_y);

        % kernel plot
        if plot_y==1
            if plot_x==1
                % h, no plot
                set(gca,'XColor', 'none','YColor','none');
                set(gca, 'color', 'none');
                title('Kernels');

            elseif plot_x<2+n_conn_kernel
                % conn_kernel
                plot(0:(kernel_len-1), conn_kernels{plot_x-1});
                xlabel(['Connection kernel ', int2str(plot_x-1)]);
                axis square;

            else
                % PS_kernel
                plot(0:(kernel_len-1), PS_kernels{plot_x-n_conn_kernel-1});
                xlabel(['Post-spike kernel ', int2str(plot_x-n_conn_kernel-1)]);
                axis square;

            end

        % original data plot
        elseif plot_y==2
            if plot_x==1
                % h
                plot(1:N, par_ori(:, 1));
                axis square;
                title('original');

            elseif plot_x<2+n_conn_kernel
                % conn_kernel
                imagesc(par_ori(:, (N*(plot_x-2) + n_PS_kernel + 2):(N*(plot_x-1) + n_PS_kernel + 1)));
                clim([-clim_all, clim_all]);
                colormap(cmap);
                % colorbar;
                axis square;
                
                % brain area borders
                hold on
                if A_T_border == T_P_border
                    line([0,N+0.5], [A_T_border,A_T_border], 'Color', 'k');
                    line([A_T_border,A_T_border], [0,N+0.5], 'Color', 'k');
                else
                    line([0,N+0.5], [A_T_border,A_T_border], 'Color', 'r');
                    line([0,N+0.5], [T_P_border,T_P_border], 'Color', 'b');
                    line([A_T_border,A_T_border], [0,N+0.5], 'Color', 'r');
                    line([T_P_border,T_P_border], [0,N+0.5], 'Color', 'b');
                end
                hold off

                xlabel('from idx');
                ylabel('to idx');

            else
                % PS_kernel
                plot(1:N, par_ori(:, plot_x-n_conn_kernel));
                axis square;
                
                xlabel('cell idx');
                ylabel('P');
                
            end
        
        % significant data plot
        elseif plot_y==3
            if plot_x==1
                plot(firing_rate);
                axis square;
                title('Firing rate');

            elseif plot_x<2+n_conn_kernel
                % conn_kernel
                im = imagesc(par_sig(:, (N*(plot_x-2) + n_PS_kernel + 2):(N*(plot_x-1) + n_PS_kernel + 1)));
                set(im, 'AlphaData', ~isnan(par_sig(:, (N*(plot_x-2) + n_PS_kernel + 2):(N*(plot_x-1) + n_PS_kernel + 1))))
                colormap(cmap);
                clim([-clim_all, clim_all]);
                % colorbar;
                axis square;

                % brain area borders
                A_T_border = borders(1);
                T_P_border = borders(2);
                hold on
                if A_T_border == T_P_border
                    line([0,N+0.5], [A_T_border,A_T_border], 'Color', 'k');
                    line([A_T_border,A_T_border], [0,N+0.5], 'Color', 'k');
                else
                    line([0,N+0.5], [A_T_border,A_T_border], 'Color', 'r');
                    line([0,N+0.5], [T_P_border,T_P_border], 'Color', 'b');
                    line([A_T_border,A_T_border], [0,N+0.5], 'Color', 'r');
                    line([T_P_border,T_P_border], [0,N+0.5], 'Color', 'b');
                end
                hold off

                xlabel('from idx');
                ylabel('to idx');

            else
                % PS_kernel
                plot(1:N, par_sig(:, plot_x-n_conn_kernel));
                axis square;

                xlabel('cell idx');
                ylabel('P');
                
            end

        % shuffled data plot
        else
            if plot_x==1
                % h
                plot(1:N, par_sfl(:, 1, plot_y-3));
                title(['Shuffled ', int2str(plot_y-3)]);
                axis square;

            elseif plot_x<2+n_conn_kernel
                % conn_kernel
                imagesc(par_sfl(:, (N*(plot_x-2) + n_PS_kernel + 2):(N*(plot_x-1) + n_PS_kernel + 1), plot_y-3));
                colormap(cmap);
                clim([-clim_all, clim_all]);
                % colorbar;
                axis square;

                % brain area borders
                A_T_border = borders(1);
                T_P_border = borders(2);
                hold on
                if A_T_border == T_P_border
                    line([0,N+0.5], [A_T_border,A_T_border], 'Color', 'k');
                    line([A_T_border,A_T_border], [0,N+0.5], 'Color', 'k');
                else
                    line([0,N+0.5], [A_T_border,A_T_border], 'Color', 'r');
                    line([0,N+0.5], [T_P_border,T_P_border], 'Color', 'b');
                    line([A_T_border,A_T_border], [0,N+0.5], 'Color', 'r');
                    line([T_P_border,T_P_border], [0,N+0.5], 'Color', 'b');
                end
                hold off

                xlabel('from idx');
                ylabel('to idx');

            else
                % PS_kernel
                plot(1:N, par_sfl(:, plot_x-n_conn_kernel, plot_y-3));
                axis square;

                xlabel('cell idx');
                ylabel('P');
            end
                
        end
    end
end

fig_path = ['../figures/GLM/', dataset_name];
check_path(fig_path);
fig_file = [fig_path, '/GLMparameters_' dataset_name, '_', ...
        int2str(session), '_', kernel_name, '_', ...
        reg.name, '_', int2str(epoch), '_sorted.png'];
% exportgraphics(fig, fig_file);
print(fig, fig_file,'-dpng', '-r100');

end