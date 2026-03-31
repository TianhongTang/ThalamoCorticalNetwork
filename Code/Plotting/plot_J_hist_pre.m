%% plot_hist_new.m - Plot histogram of J values for different kernels and distances.

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
mt = load_meta(root, 'table'); % metadata table

f = figure('Position', [100, 100, 1200, 600], 'Visible', 'off');
t = tiledlayout(2, 4, 'TileSpacing', 'Compact', 'Padding', 'Compact');

animals = {'Slayer', 'Emperor'};
injection_types = {'Muscimol', 'Saline'};
states = {'RestOpen', 'RestClose'};
state_labels = {'Eyes Open', 'Eyes Closed'};
posneg_labels = {'Positive J', 'Negative J'};

for animal_idx = 1:2
    animal_name = animals{animal_idx};
    for injection_idx = 1:2
        injection_type = injection_types{injection_idx};
        for posneg_idx = 1:2
            tile_idx = (posneg_idx-1)*4 + (animal_idx-1)*2 + injection_idx;

            open_count  = 0; open_total  = 0;
            close_count = 0; close_total = 0;

            for state_idx = 1:2
                state = states{state_idx};
                state_label = state_labels{state_idx};

                % filter metadata for current condition 
                % TODO: better condition filtering by a function
                condition_filter = strcmp(mt.GLM.animal_name, animal_name) & ...
                                   strcmp(mt.GLM.injection, injection_type) & ...
                                   strcmp(mt.GLM.state, state) & ...
                                   strcmp(mt.GLM.kernel_name, "DeltaPure") & ...
                                   strcmp(mt.GLM.align, "Last") & ...
                                   strcmp(mt.GLM.area, "Full") & ...
                                   (mt.GLM.epoch == 3000) & ...
                                   (mt.GLM.fold_idx == 0) & ...
                                   (mt.GLM.shuffle_idx == 0) & ...
                                   strcmp(mt.GLM.prepost, "Pre");
                condition_meta = mt.GLM(condition_filter, :);
                if isempty(condition_meta)
                    warning('No data found for condition: %s, %s, %s', animal_name, injection_type, state);
                    continue;
                else
                    fprintf('Plotting condition: %s, %s, %s, with %d sessions.\n', animal_name, injection_type, state, height(condition_meta));
                end

                for session_idx = 1:height(condition_meta)
                    meta = condition_meta(session_idx, :);
                    meta = table2struct(meta);

                    % load model parameters and errors
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

                    % load area borders
                    border_folder = fullfile(root, 'Data', 'Working', 'border');
                    border_file_name = generate_filename('border', meta);
                    border_file_path = fullfile(border_folder, border_file_name);
                    if ~isfile(border_file_path)
                        warning('Border file not found: %s. Skipping.', border_file_path);
                        continue;
                    end
                    border_data = load(border_file_path, 'data');
                    borders = border_data.data.borders;

                    selected_areas = {[1, 2], [2, 1]}; % between ACC and VLPFC

                    % count significant J values
                    [pos, neg, total] = J_count(model_par, model_err, 1, borders, selected_areas);
                    counts = [pos, neg];
                    count = counts(posneg_idx);

                    if state_idx == 1
                        open_count = open_count + count;
                        open_total = open_total + total;
                    else
                        close_count = close_count + count;
                        close_total = close_total + total;
                    end

                end
            end
            % error bar and statistical test
            [CI_open_low, CI_open_high] = wilsonCI(open_count, open_total, 0.05);
            [CI_close_low, CI_close_high] = wilsonCI(close_count, close_total, 0.05);
            p_val = twoProportionPValue(open_count, open_total, close_count, close_total);
            if p_val < 0.001
                sig_label = '***';
            elseif p_val < 0.01
                sig_label = '**';
            elseif p_val < 0.05
                sig_label = '*';
            else
                sig_label = 'N.S.';
            end

            % plot histogram
            nexttile(tile_idx);
            ratios = [open_count/open_total, close_count/close_total]*100;
            low_errs = [ratios(1) - CI_open_low*100, ratios(2) - CI_close_low*100];
            high_errs = [CI_open_high*100 - ratios(1), CI_close_high*100 - ratios(2)];

            hold on;
            bar(1:2, ratios);
            errorbar(1:2, ratios, low_errs, high_errs, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
            text(1.5, 25, sig_label, 'HorizontalAlignment', 'center', 'FontSize', 20);
            hold off;

            xticks(1:2);
            xticklabels(state_labels);
            ylabel('Percentage of Significant J (%)');
            title(sprintf('%s, %s, %s', animal_name, injection_type, posneg_labels{posneg_idx}));
            ylim([0, 30]);
        end
    end
end

save_folder = fullfile(root, 'Figures', 'J_hist');
check_path(save_folder);
save_path = fullfile(save_folder, 'J_hist_summary.png');
saveas(f, save_path);


function [p_low, p_high] = wilsonCI(M, N, alpha)
    % Wilson score interval
    % M = successes, N = trials
    % alpha = significance level (e.g., 0.05 for 95% CI)

    if nargin < 3
        alpha = 0.05;
    end

    p = M / N;
    z = norminv(1 - alpha/2); % Z for CI

    denominator = 1 + (z^2)/N;
    center = p + (z^2)/(2*N);
    radius = z * sqrt( (p*(1-p)/N) + (z^2)/(4*N^2) );

    p_low = (center - radius) / denominator;
    p_high = (center + radius) / denominator;
end

function p_val = twoProportionPValue(success1, trials1, success2, trials2)
    % Two-proportion z-test (two-tailed) between pre and post proportions

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
