%% plot_post_pre_effect_from_extracted_counts.m
% Hard-coded counts extracted from the provided bar plots.
% Plot pooled Post - Pre effect for Saline vs Muscimol.
%
% Effect = 100 * (Post significant J ratio - Pre significant J ratio)
% Error bars = normal approximation SE of the within-injection difference
% Stars = z-test for difference-in-differences:
%   (Muscimol_Post - Muscimol_Pre) - (Saline_Post - Saline_Pre)

clear;

%% Data
kernels = 1:3;
kernel_labels = arrayfun(@(k) sprintf('Kernel %d', k), kernels, 'UniformOutput', false);
injection_labels = {'Saline', 'Muscimol'};
state_labels = {'Eyes Open', 'Eyes Closed'};
sign_labels = {'Positive J', 'Negative J'};
prepost_labels = {'Pre', 'Post'};

% Dimensions: state x sign x prepost x kernel
% state: 1=Open, 2=Close
% sign:  1=Positive, 2=Negative
% prepost: 1=Pre, 2=Post
% kernel: 1,2,3

counts = struct();
totals = struct();

% Saline counts extracted from the second provided figure.
totals.Saline = 7442 * ones(2, 2, 2, 3);
counts.Saline = nan(2, 2, 2, 3);

% Saline, Positive J
counts.Saline(1, 1, 1, :) = [526, 443, 339]; % Open Pre
counts.Saline(1, 1, 2, :) = [632, 519, 398]; % Open Post
counts.Saline(2, 1, 1, :) = [876, 783, 489]; % Close Pre
counts.Saline(2, 1, 2, :) = [857, 753, 526]; % Close Post

% Saline, Negative J
counts.Saline(1, 2, 1, :) = [324, 307, 219]; % Open Pre
counts.Saline(1, 2, 2, :) = [396, 337, 226]; % Open Post
counts.Saline(2, 2, 1, :) = [308, 272, 291]; % Close Pre
counts.Saline(2, 2, 2, :) = [418, 407, 337]; % Close Post

% Muscimol counts extracted from the previous Muscimol figure.
totals.Muscimol = 10834 * ones(2, 2, 2, 3);
counts.Muscimol = nan(2, 2, 2, 3);

% Muscimol, Positive J
counts.Muscimol(1, 1, 1, :) = [619, 488, 286];  % Open Pre
counts.Muscimol(1, 1, 2, :) = [900, 745, 500];  % Open Post
counts.Muscimol(2, 1, 1, :) = [1055, 848, 556]; % Close Pre
counts.Muscimol(2, 1, 2, :) = [1154, 959, 685]; % Close Post

% Muscimol, Negative J
counts.Muscimol(1, 2, 1, :) = [331, 278, 130]; % Open Pre
counts.Muscimol(1, 2, 2, :) = [496, 398, 252]; % Open Post
counts.Muscimol(2, 2, 1, :) = [437, 341, 238]; % Close Pre
counts.Muscimol(2, 2, 2, :) = [685, 606, 454]; % Close Post

%% Compute effects and pooled z-tests
n_state = numel(state_labels);
n_sign = numel(sign_labels);
n_kernel = numel(kernels);
n_injection = numel(injection_labels);

effect_pct = nan(n_state, n_sign, n_kernel, n_injection);
se_effect_pct = nan(n_state, n_sign, n_kernel, n_injection);
p_diff = nan(n_state, n_sign, n_kernel);
effect_diff_pct = nan(n_state, n_sign, n_kernel);

for state_idx = 1:n_state
    for sign_idx = 1:n_sign
        for kernel_idx = 1:n_kernel
            delta = nan(1, n_injection);
            se_delta = nan(1, n_injection);

            for inj_idx = 1:n_injection
                inj = injection_labels{inj_idx};
                pre_count = counts.(inj)(state_idx, sign_idx, 1, kernel_idx);
                post_count = counts.(inj)(state_idx, sign_idx, 2, kernel_idx);
                pre_total = totals.(inj)(state_idx, sign_idx, 1, kernel_idx);
                post_total = totals.(inj)(state_idx, sign_idx, 2, kernel_idx);

                p_pre = pre_count / pre_total;
                p_post = post_count / post_total;

                delta(inj_idx) = p_post - p_pre;
                se_delta(inj_idx) = sqrt(p_post * (1 - p_post) / post_total + ...
                                         p_pre  * (1 - p_pre)  / pre_total);

                effect_pct(state_idx, sign_idx, kernel_idx, inj_idx) = 100 * delta(inj_idx);
                se_effect_pct(state_idx, sign_idx, kernel_idx, inj_idx) = 100 * se_delta(inj_idx);
            end

            % Compare Muscimol effect vs Saline effect.
            effect_diff = delta(2) - delta(1);
            se_diff = sqrt(se_delta(1)^2 + se_delta(2)^2);
            if se_diff > 0
                z = effect_diff / se_diff;
                p_diff(state_idx, sign_idx, kernel_idx) = erfc(abs(z) / sqrt(2));
            end
            effect_diff_pct(state_idx, sign_idx, kernel_idx) = 100 * effect_diff;
        end
    end
end

%% Plot, matching the first figure layout
bar_width = 0.30;
bar_offset = 0.18;
colors = struct();
colors.Saline = [0, 0, 1];
colors.Muscimol = [1, 0, 0];

f = figure('Color', 'w', 'Position', [100, 100, 1300, 900]);
t = tiledlayout(f, 2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for state_idx = 1:n_state
    for sign_idx = 1:n_sign
        tile_idx = (state_idx - 1) * 2 + sign_idx;
        ax = nexttile(t, tile_idx);
        hold(ax, 'on');

        x = 1:n_kernel;
        saline_y = squeeze(effect_pct(state_idx, sign_idx, :, 1));
        musc_y   = squeeze(effect_pct(state_idx, sign_idx, :, 2));
        saline_e = squeeze(se_effect_pct(state_idx, sign_idx, :, 1));
        musc_e   = squeeze(se_effect_pct(state_idx, sign_idx, :, 2));

        bar(ax, x - bar_offset, saline_y, bar_width, 'FaceColor', colors.Saline, 'DisplayName', 'Saline');
        bar(ax, x + bar_offset, musc_y,   bar_width, 'FaceColor', colors.Muscimol, 'DisplayName', 'Muscimol');

        errorbar(ax, x - bar_offset, saline_y, saline_e, 'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
        errorbar(ax, x + bar_offset, musc_y,   musc_e,   'k', 'LineStyle', 'none', 'HandleVisibility', 'off');

        yline(ax, 0, 'k-', 'HandleVisibility', 'off');

        all_y = [saline_y(:) - saline_e(:); saline_y(:) + saline_e(:); ...
                 musc_y(:) - musc_e(:); musc_y(:) + musc_e(:); 0];
        y_min = min(all_y, [], 'omitnan');
        y_max = max(all_y, [], 'omitnan');
        y_pad = max(0.35, 0.12 * (y_max - y_min + eps));
        ylim(ax, [y_min - y_pad, y_max + 2.0 * y_pad]);

        for kernel_idx = 1:n_kernel
            y_text = max([saline_y(kernel_idx) + saline_e(kernel_idx), ...
                          musc_y(kernel_idx) + musc_e(kernel_idx)], [], 'omitnan') + 0.45 * y_pad;
            text(ax, x(kernel_idx), y_text, p_to_stars(p_diff(state_idx, sign_idx, kernel_idx)), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 12);
        end

        xticks(ax, x);
        xticklabels(ax, kernel_labels);
        ylabel(ax, 'Post - Pre significant J (percentage points)');
        title(ax, sprintf('%s, %s', state_labels{state_idx}, sign_labels{sign_idx}), 'FontWeight', 'bold');
        legend(ax, 'Location', 'best');
        box(ax, 'off');
        hold(ax, 'off');
    end
end

sgtitle(f, 'Pooled Post - Pre effect: Saline vs Muscimol');

%% Print summary table in command window
summary_rows = {};
for state_idx = 1:n_state
    for sign_idx = 1:n_sign
        for kernel_idx = 1:n_kernel
            summary_rows(end+1, :) = { ... %#ok<SAGROW>
                state_labels{state_idx}, ...
                sign_labels{sign_idx}, ...
                kernel_idx, ...
                effect_pct(state_idx, sign_idx, kernel_idx, 1), ...
                effect_pct(state_idx, sign_idx, kernel_idx, 2), ...
                effect_diff_pct(state_idx, sign_idx, kernel_idx), ...
                p_diff(state_idx, sign_idx, kernel_idx), ...
                p_to_stars(p_diff(state_idx, sign_idx, kernel_idx))};
        end
    end
end
summary_table = cell2table(summary_rows, 'VariableNames', ...
    {'State', 'Sign', 'Kernel', 'Saline_PostMinusPre_pctpt', 'Muscimol_PostMinusPre_pctpt', ...
     'MuscimolMinusSaline_pctpt', 'P_value', 'Stars'});
disp(summary_table);

%% Optional export
% save_folder = fullfile(pwd, 'Figures');
% if ~isfolder(save_folder), mkdir(save_folder); end
% exportgraphics(f, fullfile(save_folder, 'post_pre_effect_extracted_counts_preview.jpg'), ...
%     'ContentType', 'image', 'Resolution', 300, 'BackgroundColor', 'white');
% exportgraphics(f, fullfile(save_folder, 'post_pre_effect_extracted_counts.pdf'), ...
%     'ContentType', 'vector', 'BackgroundColor', 'white');

function s = p_to_stars(p)
    if isnan(p)
        s = 'n/a';
    elseif p < 0.001
        s = '***';
    elseif p < 0.01
        s = '**';
    elseif p < 0.05
        s = '*';
    else
        s = 'N.S.';
    end
end
