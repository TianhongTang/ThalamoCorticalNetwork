%% J_bars_kernels.m - Plot significant J percentage across kernels.
% Three kernels are combined into one 2x3 figure.
% Columns are Kernel 1/2/3. Rows are state or pre/post.
% Only All animals is used; individual Slayer/Emperor columns are removed.

clear;

%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end

addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Main parameters
resting_dur_threshold = 15;
mt = load_meta(root, 'table');

kernel_indices = 1:3;
injection_types = {'Muscimol', 'Saline'};
states = {'RestOpen', 'RestClose'};
state_labels = {'Eyes Open', 'Eyes Closed'};
posneg_labels = {'Positive J', 'Negative J'};
prepost_labels = {'Pre', 'Post'};
alignment = 'Longest';

animal_filter = strcmp(mt.GLM.animal_name, 'Slayer') | strcmp(mt.GLM.animal_name, 'Emperor');

%% Fig1: Pre vs Post. Fig2: Open vs Closed.
for injection_idx = 1:numel(injection_types)
    injection = injection_types{injection_idx};

    %% Figure 1: Pre vs Post, columns = kernels, rows = states.
    f1 = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
    t1 = tiledlayout(f1, 2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for kernel_i = 1:numel(kernel_indices)
        kernel_idx = kernel_indices(kernel_i);

        for state_idx = 1:numel(states)
            state = states{state_idx};
            state_label = state_labels{state_idx};

            J_counts = zeros(2, 2); % pos/neg x pre/post
            total_counts = zeros(2, 2); % total count for pos/neg in pre/post

            for prepost_idx = 1:numel(prepost_labels)
                prepost_str = prepost_labels{prepost_idx};

                condition_filter = animal_filter & ...
                                    strcmp(mt.GLM.state, state) & ...
                                    strcmp(mt.GLM.kernel_name, "DeltaPure") & ...
                                    strcmp(mt.GLM.align, alignment) & ...
                                    strcmp(mt.GLM.area, "Cortex") & ...
                                    strcmp(mt.GLM.injection, injection) & ...
                                    (mt.GLM.epoch == 3000) & ...
                                    (mt.GLM.fold_idx == 0) & ...
                                    (mt.GLM.shuffle_idx == 0) & ...
                                    cellfun(@(x) ~isempty(x) && x == resting_dur_threshold, mt.GLM.resting_dur_threshold) & ...
                                    strcmp(mt.GLM.prepost, prepost_str);

                condition_meta = mt.GLM(condition_filter, :);
                if isempty(condition_meta)
                    warning('No data found for condition: All animals, %s, %s %s, Kernel %d', ...
                        state, prepost_str, injection, kernel_idx);
                    continue;
                else
                    fprintf('Plotting condition: All animals, %s, %s %s, Kernel %d with %d sessions.\n', ...
                        state, prepost_str, injection, kernel_idx, height(condition_meta));
                end

                [pos_count, neg_count, total_count] = count_condition_J(root, condition_meta, kernel_idx);
                J_counts(:, prepost_idx) = [pos_count; neg_count];
                total_counts(:, prepost_idx) = [total_count; total_count];
            end

            [ratios, low_errs, high_errs, sig_labels] = compute_bar_stats(J_counts, total_counts);

            tile_idx = (state_idx - 1) * numel(kernel_indices) + kernel_i;
            tile = nexttile(t1, tile_idx);

            hold(tile, "on");
            bar(tile, (1:2)-0.15, ratios(:, 1), 'FaceColor', 'b', 'BarWidth', 0.3);
            bar(tile, (1:2)+0.15, ratios(:, 2), 'FaceColor', 'r', 'BarWidth', 0.3);
            errorbar(tile, (1:2)-0.15, ratios(:, 1), low_errs(:, 1), high_errs(:, 1), 'k', 'LineStyle', 'none');
            errorbar(tile, (1:2)+0.15, ratios(:, 2), low_errs(:, 2), high_errs(:, 2), 'k', 'LineStyle', 'none');
            text(tile, 1, max(ratios(:, 1))+5, sig_labels{1}, 'HorizontalAlignment', 'center', 'FontSize', 14);
            text(tile, 2, max(ratios(:, 2))+5, sig_labels{2}, 'HorizontalAlignment', 'center', 'FontSize', 14);
            hold(tile, "off");

            xticks(tile, 1:2);
            xticklabels(tile, posneg_labels);
            ylabel(tile, 'Significant J %');
            legend(tile, prepost_labels, 'Location', 'Best');
            title(tile, sprintf('Kernel %d, %s', kernel_idx, state_label));
            ylim(tile, [0, 20]);

            disp(J_counts);
            disp(total_counts);
        end
    end

    sgtitle(f1, sprintf("%s, Pre vs Post", injection), 'FontSize', 14);
    export_figure(root, f1, sprintf('J_bars_prepost_kernels_%s', injection));

    %% Figure 2: Open vs Closed, columns = kernels, rows = Pre/Post.
    f2 = figure('Position', [100, 100, 1200, 800], 'Visible', 'off');
    t2 = tiledlayout(f2, 2, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

    for kernel_i = 1:numel(kernel_indices)
        kernel_idx = kernel_indices(kernel_i);

        for prepost_idx = 1:numel(prepost_labels)
            prepost_str = prepost_labels{prepost_idx};

            J_counts = zeros(2, 2); % pos/neg x open/closed
            total_counts = zeros(2, 2); % total count for pos/neg in open/closed

            for state_idx = 1:numel(states)
                state = states{state_idx};

                condition_filter = animal_filter & ...
                                    strcmp(mt.GLM.state, state) & ...
                                    strcmp(mt.GLM.kernel_name, "DeltaPure") & ...
                                    strcmp(mt.GLM.align, alignment) & ...
                                    strcmp(mt.GLM.area, "Cortex") & ...
                                    strcmp(mt.GLM.injection, injection) & ...
                                    (mt.GLM.epoch == 3000) & ...
                                    (mt.GLM.fold_idx == 0) & ...
                                    (mt.GLM.shuffle_idx == 0) & ...
                                    cellfun(@(x) ~isempty(x) && x == resting_dur_threshold, mt.GLM.resting_dur_threshold) & ...
                                    strcmp(mt.GLM.prepost, prepost_str);

                condition_meta = mt.GLM(condition_filter, :);
                if isempty(condition_meta)
                    warning('No data found for condition: All animals, %s, %s %s, Kernel %d', ...
                        state, prepost_str, injection, kernel_idx);
                    continue;
                else
                    fprintf('Plotting condition: All animals, %s, %s %s, Kernel %d with %d sessions.\n', ...
                        state, prepost_str, injection, kernel_idx, height(condition_meta));
                end

                [pos_count, neg_count, total_count] = count_condition_J(root, condition_meta, kernel_idx);
                J_counts(:, state_idx) = [pos_count; neg_count];
                total_counts(:, state_idx) = [total_count; total_count];
            end

            [ratios, low_errs, high_errs, sig_labels] = compute_bar_stats(J_counts, total_counts);

            tile_idx = (prepost_idx - 1) * numel(kernel_indices) + kernel_i;
            tile = nexttile(t2, tile_idx);

            hold(tile, "on");
            bar(tile, (1:2)-0.15, ratios(:, 1), 'FaceColor', 'b', 'BarWidth', 0.3);
            bar(tile, (1:2)+0.15, ratios(:, 2), 'FaceColor', 'r', 'BarWidth', 0.3);
            errorbar(tile, (1:2)-0.15, ratios(:, 1), low_errs(:, 1), high_errs(:, 1), 'k', 'LineStyle', 'none');
            errorbar(tile, (1:2)+0.15, ratios(:, 2), low_errs(:, 2), high_errs(:, 2), 'k', 'LineStyle', 'none');
            text(tile, 1, max(ratios(:, 1))+5, sig_labels{1}, 'HorizontalAlignment', 'center', 'FontSize', 14);
            text(tile, 2, max(ratios(:, 2))+5, sig_labels{2}, 'HorizontalAlignment', 'center', 'FontSize', 14);
            hold(tile, "off");

            xticks(tile, 1:2);
            xticklabels(tile, posneg_labels);
            ylabel(tile, 'Significant J %');
            legend(tile, state_labels, 'Location', 'Best');
            title(tile, sprintf('Kernel %d, %s', kernel_idx, prepost_str));
            ylim(tile, [0, 20]);

            disp(J_counts);
            disp(total_counts);
        end
    end

    sgtitle(f2, sprintf("%s, Open vs Closed", injection), 'FontSize', 14);
    export_figure(root, f2, sprintf('J_bars_openclosed_kernels_%s', injection));
end

function [pos_count, neg_count, total_count] = count_condition_J(root, condition_meta, kernel_idx)
    pos_count = 0;
    neg_count = 0;
    total_count = 0;

    for session_idx = 1:height(condition_meta)
        meta = condition_meta(session_idx, :);
        meta = table2struct(meta);

        data_folder = fullfile(root, 'Data', 'Working', 'GLM');
        file_name = meta.file_name;
        file_path = fullfile(data_folder, file_name);
        if ~isfile(file_path)
            warning('File not found: %s. Skipping.', file_path);
            continue;
        end
        session_data = load(file_path, 'meta', 'data');
        model_par = session_data.data.model_par;
        model_err = session_data.data.model_err.total;

        border_folder = fullfile(root, 'Data', 'Working', 'border');
        border_file_name = generate_filename('border', meta);
        border_file_path = fullfile(border_folder, border_file_name);
        if ~isfile(border_file_path)
            warning('Border file not found: %s. Skipping.', border_file_path);
            continue;
        end
        border_data = load(border_file_path, 'data', 'meta');
        borders = border_data.data.borders;
        N = border_data.meta.N;
        assert(numel(borders) == 2);
        borders = [borders, N+1];

        selected_areas = {[1, 2], [2, 1]}; % between ACC and VLPFC

        [pos, neg, total] = J_count(model_par, model_err, kernel_idx, borders, selected_areas);
        pos_count = pos_count + pos;
        neg_count = neg_count + neg;
        total_count = total_count + total;
    end
end

function [ratios, low_errs, high_errs, sig_labels] = compute_bar_stats(J_counts, total_counts)
    ratios_raw = J_counts ./ total_counts;
    CIs = zeros(2, 2, 2); % pos/neg x condition x low/high
    sig_labels = cell(1, 2);

    for posneg_idx = 1:2
        count = J_counts(posneg_idx, :);
        total = total_counts(posneg_idx, :);

        for condition_idx = 1:2
            [CI_low, CI_high] = wilsonCI(count(condition_idx), total(condition_idx), 0.05);
            CIs(posneg_idx, condition_idx, :) = [CI_low, CI_high];
        end

        p_val = twoProportionPValue(count(1), total(1), count(2), total(2));
        if p_val < 0.001
            sig_labels{posneg_idx} = '***';
        elseif p_val < 0.01
            sig_labels{posneg_idx} = '**';
        elseif p_val < 0.05
            sig_labels{posneg_idx} = '*';
        else
            sig_labels{posneg_idx} = 'N.S.';
        end
    end

    ratios = ratios_raw * 100;
    low_errs = ratios - squeeze(CIs(:, :, 1)) * 100;
    high_errs = squeeze(CIs(:, :, 2)) * 100 - ratios;
end

function export_figure(root, fig, output_stub)
    save_folder = fullfile(root, 'Figures', 'Paper');
    check_path(save_folder);

    figWidth  = 12.0;
    figHeight = 8.0;
    resolution = 300;

    set(fig, 'Units', 'inches');
    fig.Position(3:4) = [figWidth, figHeight];

    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperSize', [figWidth, figHeight]);
    set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);
    set(fig, 'Color', 'w');

    preview_filename = fullfile(save_folder, [output_stub, '_preview.jpg']);
    exportgraphics(fig, preview_filename, ...
        'ContentType', 'image', ...
        'BackgroundColor', 'white', ...
        'Resolution', resolution);

    pdf_filename = fullfile(save_folder, [output_stub, '.pdf']);
    exportgraphics(fig, pdf_filename, ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white', ...
        'Resolution', resolution);

    close(fig);
end

function [p_low, p_high] = wilsonCI(M, N, alpha)
    if nargin < 3
        alpha = 0.05;
    end

    if N == 0
        p_low = NaN;
        p_high = NaN;
        return;
    end

    p = M / N;
    z = norminv(1 - alpha/2);

    denominator = 1 + (z^2)/N;
    center = p + (z^2)/(2*N);
    radius = z * sqrt((p*(1-p)/N) + (z^2)/(4*N^2));

    p_low = (center - radius) / denominator;
    p_high = (center + radius) / denominator;
end

function p_val = twoProportionPValue(success1, trials1, success2, trials2)
    if trials1 == 0 || trials2 == 0
        p_val = NaN;
        return;
    end

    p1 = success1 / trials1;
    p2 = success2 / trials2;
    pooled = (success1 + success2) / (trials1 + trials2);
    denom = sqrt(pooled * (1 - pooled) * (1 / trials1 + 1 / trials2));
    if denom == 0
        p_val = 1;
        return;
    end

    z = (p1 - p2) / denom;
    abs_z = abs(z);
    Phi = 0.5 * erfc(-abs_z / sqrt(2));
    p_val = 2 * (1 - Phi);
end
