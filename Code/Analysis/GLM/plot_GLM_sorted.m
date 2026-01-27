function plot_GLM_sorted(dataset_name, border_name, session, kernel_name, epoch, reg, shuffle_size, sorting, idx)
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
model_folder = fullfile(root, 'Data', 'Working', 'GLM_models');
model_file = sprintf('GLM_%s_s%d_shuffle%d_%s_%s_epoch%d_fold%d.mat', dataset_name, session, 0, kernel_name, ...
                        reg.name, epoch, 0);
model_path = fullfile(model_folder, model_file);

data_folder = fullfile(root, 'Data', 'Working', 'GLM_data');
data_file = sprintf('GLMdata_%s_%d_%d_%s.mat', dataset_name, session, 0, kernel_name);
data_path = fullfile(data_folder, data_file);

border_folder = fullfile(root, 'Data', 'Working', 'border');
border_file = sprintf('borders_%s_%d.mat', border_name, session);
border_path = fullfile(border_folder, border_file);
                
load(model_path, "model_par", "model_err", "kernel", "N");
load(data_path, "folds", "fold_num");
PS_kernels = kernel.PS_kernels;
conn_kernels = kernel.conn_kernels;
n_PS_kernel = kernel.n_PS_kernel;
n_conn_kernel = kernel.n_conn_kernel;
kernel_len = kernel.kernel_len;

spike_total = zeros(N, 1);
t_total = 0;
for fold_idx = 1:fold_num
    spike_total = spike_total + sum(folds{fold_idx}.raster, 2);
    t_total = t_total + size(folds{fold_idx}.raster, 2);
end
firing_rate = 1000 * spike_total/t_total;

% if session<100
%     session_border = session;
% else
%     session_border = floor(session/100);
% end

load(border_path, "borders");
if length(borders)==2
    borders = [borders(1), borders(2), borders(2)];
end

% for new dataset
A_V_border = borders(2)-0.5;
V_T_border = borders(3)-0.5;
    
par_ori = model_par;

par_sfl = zeros([size(par_ori), shuffle_size]);
for i=1:shuffle_size
    model_path_sfl = fullfile(model_folder, sprintf('GLM_%s_s%d_shuffle%d_%s_%s_epoch%d_fold%d.mat', dataset_name, session, i, kernel_name, reg.name, epoch, 0));
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
sorting_ranges = [1, A_V_border-0.5;A_V_border+0.5,V_T_border-0.5;V_T_border+0.5, N];
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
sort_folder = fullfile(root, 'Data', 'Working', "sort_idx");
check_path(sort_folder);
sort_file = sprintf('sortidx_%s_%d_%s.mat', dataset_name, session, kernel_name);
sort_path = fullfile(sort_folder, sort_file);
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
                if A_V_border == V_T_border
                    line([0,N+0.5], [A_V_border,A_V_border], 'Color', 'k');
                    line([A_V_border,A_V_border], [0,N+0.5], 'Color', 'k');
                else
                    line([0,N+0.5], [A_V_border,A_V_border], 'Color', 'r');
                    line([0,N+0.5], [V_T_border,V_T_border], 'Color', 'b');
                    line([A_V_border,A_V_border], [0,N+0.5], 'Color', 'r');
                    line([V_T_border,V_T_border], [0,N+0.5], 'Color', 'b');
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
                hold on
                if A_V_border == V_T_border
                    line([0,N+0.5], [A_V_border,A_V_border], 'Color', 'k');
                    line([A_V_border,A_V_border], [0,N+0.5], 'Color', 'k');
                else
                    line([0,N+0.5], [A_V_border,A_V_border], 'Color', 'r');
                    line([0,N+0.5], [V_T_border,V_T_border], 'Color', 'b');
                    line([A_V_border,A_V_border], [0,N+0.5], 'Color', 'r');
                    line([V_T_border,V_T_border], [0,N+0.5], 'Color', 'b');
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
                hold on
                if A_V_border == V_T_border
                    line([0,N+0.5], [A_V_border,A_V_border], 'Color', 'k');
                    line([A_V_border,A_V_border], [0,N+0.5], 'Color', 'k');
                else
                    line([0,N+0.5], [A_V_border,A_V_border], 'Color', 'r');
                    line([0,N+0.5], [V_T_border,V_T_border], 'Color', 'b');
                    line([A_V_border,A_V_border], [0,N+0.5], 'Color', 'r');
                    line([V_T_border,V_T_border], [0,N+0.5], 'Color', 'b');
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

fig_folder = fullfile(root, 'Figures', 'GLM_models', dataset_name);
check_path(fig_folder);
fig_file = sprintf('GLMpar_%s_%d_%s_%s_epoch%d_sorted.png', dataset_name, session, kernel_name, reg.name, epoch);
fig_path = fullfile(fig_folder, fig_file);
% exportgraphics(fig, fig_file);
print(fig, fig_path,'-dpng', '-r100');

end