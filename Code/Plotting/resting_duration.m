%% resting_duration.m - Histogram showing resting state duration distribution
% 

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
threshold = 15;

% load metadata
metadata_folder = fullfile(root, 'Data', 'Working', 'Meta');  
metadata_path = fullfile(metadata_folder, 'PDS_dataset_info.mat');  
load(metadata_path, 'dataset_num', 'dataset_names', 'session_nums', 'cortex_files', 'thalamus_files', 'eyeID_files');

merge_types = {'Full', 'Cortex'};
prepost_types = {'Pre', 'Post'};
% states = {'RestOpen', 'RestClose', 'Task'};
states = {'RestOpen', 'RestClose'};
state_num = length(states);
align = 'AlignLast';
area_type = 'Cortex';
prepost = 'Post';
mode_str = sprintf('_separate_%d', threshold);
% mode_str = '_overlay';

durations = cell(dataset_num, state_num);
session_count = zeros(dataset_num, state_num);
session_durations = cell(dataset_num, state_num);

for dataset_idx = [1,2,4,5,7,8] % (Slayer Muscimol, Slayer Saline, Zeppelin Muscimol, Zeppelin Saline, Emperor Muscimol, Emperor Saline)
    dataset_name = dataset_names{dataset_idx};
    session_num = session_nums(dataset_idx);
    % if dataset_idx == 7
    %     session_num = 2;
    % end

    for state_idx = 1:state_num
        state = states{state_idx};
        for session_idx = 1:session_num
            % load aligned data
            full_name = [dataset_name, prepost, state, area_type, align];
            data_folder = fullfile(root, 'Data', 'Working', 'raster');
            check_path(data_folder);
            file_name = sprintf('raster_%s_%d.mat', full_name, session_idx);
            data_path = fullfile(data_folder, file_name);
            load(data_path, 'trial_len');
            trial_len = trial_len(trial_len >= (threshold*1000)); % only include trials longer than threshold seconds
            durations{dataset_idx, state_idx} = [durations{dataset_idx, state_idx}, trial_len];
            session_durations{dataset_idx, state_idx} = [session_durations{dataset_idx, state_idx}, sum(trial_len)];
        end
        fprintf('Dataset: %s, State: %s\n', dataset_name, state);
        fprintf('Session Num: %d, trial num: %d, min duration: %d, max duration: %d, mean duration: %.2f\n', session_num, length(durations{dataset_idx, state_idx}), min(durations{dataset_idx, state_idx}), max(durations{dataset_idx, state_idx}), mean(durations{dataset_idx, state_idx}));
        fprintf('All Durations: %s\n', mat2str(durations{dataset_idx, state_idx}));
        fprintf('Session Durations: %s\n', mat2str(session_durations{dataset_idx, state_idx}));
        session_count(dataset_idx, state_idx) = session_num;
    end
end
durations{3, 1} = [durations{1, 1}, durations{2, 1}];
durations{3, 2} = [durations{1, 2}, durations{2, 2}];
session_count(3, :) = session_count(1, :) + session_count(2, :);
durations{6, 1} = [durations{4, 1}, durations{5, 1}]; 
durations{6, 2} = [durations{4, 2}, durations{5, 2}];
session_count(6, :) = session_count(4, :) + session_count(5, :);
durations{9, 1} = [durations{7, 1}, durations{8, 1}];
durations{9, 2} = [durations{7, 2}, durations{8, 2}];
session_count(9, :) = session_count(7, :) + session_count(8, :);

session_durations{3, 1} = [session_durations{1, 1}, session_durations{2, 1}];
session_durations{3, 2} = [session_durations{1, 2}, session_durations{2, 2}];
session_durations{6, 1} = [session_durations{4, 1}, session_durations{5, 1}];
session_durations{6, 2} = [session_durations{4, 2}, session_durations{5, 2}];
session_durations{9, 1} = [session_durations{7, 1}, session_durations{8, 1}];
session_durations{9, 2} = [session_durations{7, 2}, session_durations{8, 2}];

% f = figure('Position', [100, 100, 1200, 600], 'Visible', 'off');
% t = tiledlayout(3, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
% dataset_labels = {'Slayer Muscimol', 'Slayer Saline', 'Slayer', 'Zeppelin Muscimol', 'Zeppelin Saline', 'Zeppelin', 'Emperor Muscimol', 'Emperor Saline', 'Emperor'};


% for dataset_idx = 1:dataset_num
%     nexttile(dataset_idx);
%     % two colors for two states, overlayed histogram with some transparency
%     durations_1 = durations{dataset_idx, 1}/1000;
%     durations_2 = durations{dataset_idx, 2}/1000;

%     % combine all >150s durations into one bin and label it as >150s
%     durations_1(durations_1 > 150) = 172.5;
%     durations_2(durations_2 > 150) = 172.5;

%     % dummy histograms to show legend
%     histogram([], hist_borders, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', states{1}); 
%     hold on;   
%     histogram([], hist_borders, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', states{2});

%     % actual histograms with no edge color and some transparency
%     histogram(durations_1, hist_borders, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
%     histogram(durations_2, hist_borders, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
%     xline(150, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
%     hold off;
%     xlabel('Duration (s)');
%     ylabel('Count');
%     title_str = sprintf('%d sessions, Open: %d trials, Close: %d trials', session_count(dataset_idx, 1), numel(durations{dataset_idx, 1}), numel(durations{dataset_idx, 2}));
%     title({dataset_labels{dataset_idx}, title_str}, 'FontSize', 10);
%     legend('Location', 'northeast');
%     xlim([0, 190]);
%     ylim([0, 50]);
%     xticks([0, 25, 50, 75, 100, 125, 150, 172.5]);
%     xticklabels({'0', '25', '50', '75', '100', '125', '150', '(>150)'});

% end
controls = {'Muscimol', 'Saline', 'All'};

for control_idx = 1:3
    control = controls{control_idx};
    f = figure('Position', [100, 100, 1200, 600], 'Visible', 'off');
    t = tiledlayout(3, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    
    for dataset_idx = 1:3
        nexttile(dataset_idx*2-1);
        % nexttile(1);
        hist_borders = 0:5:190;
        % two colors for two states, overlayed histogram with some transparency
        durations_Slayer = durations{control_idx, 2}/1000;
        durations_Zeppelin = durations{control_idx+3, 2}/1000;
        durations_Emperor = durations{control_idx+6, 2}/1000;

        % combine all >150s durations into one bin and label it as >150s
        durations_Slayer(durations_Slayer > 150) = 172.5;
        durations_Zeppelin(durations_Zeppelin > 150) = 172.5;
        durations_Emperor(durations_Emperor > 150) = 172.5;

        % dummy histograms to show legend
        hold on;   

        % actual histograms with no edge color and some transparency
        switch dataset_idx
            case 1
                histogram([], hist_borders, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'Slayer');
                histogram(durations_Slayer, hist_borders, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                title('Slayer, by Epoch', 'FontSize', 12);
            case 2
                histogram([], hist_borders, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'Zeppelin');
                histogram(durations_Zeppelin, hist_borders, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                title('Zeppelin, by Epoch', 'FontSize', 12);
            case 3
                histogram([], hist_borders, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'DisplayName', 'Emperor');
                histogram(durations_Emperor, hist_borders, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                title('Emperor, by Epoch', 'FontSize', 12);
        end
        % title('All Monkeys, by Epoch', 'FontSize', 12);
        xline(150, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        xline(threshold, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
        hold off;
        xlabel('Duration (s)');
        ylabel('Count');
        % title_str = sprintf('%d sessions, Open: %d trials, Close: %d trials', session_count(dataset_idx, 1), numel(durations{dataset_idx, 1}), numel(durations{dataset_idx, 2}));
        % title({dataset_labels{dataset_idx}, title_str}, 'FontSize', 10);
        legend('Location', 'northeast');
        xlim([0, 190]);
        ylim([0, 50]);
        xticks([0, 25, 50, 75, 100, 125, 150, 172.5]);
        xticklabels({'0', '25', '50', '75', '100', '125', '150', '(>150)'});

        %% by session
        nexttile(dataset_idx*2);
        % nexttile(2);
        hist_borders = 0:25:650;
        % two colors for two states, overlayed histogram with some transparency
        durations_Slayer = session_durations{control_idx, 2}/1000;
        durations_Zeppelin = session_durations{control_idx+3, 2}/1000;
        durations_Emperor = session_durations{control_idx+6, 2}/1000;

        % combine all >150s durations into one bin and label it as >150s
        durations_Slayer(durations_Slayer > 500) = 562.5;
        durations_Zeppelin(durations_Zeppelin > 500) = 562.5;
        durations_Emperor(durations_Emperor > 500) = 562.5;

        % dummy histograms to show legend
        hold on;   

        % actual histograms with no edge color and some transparency
        switch dataset_idx
            case 1
                histogram([], hist_borders, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'Slayer');
                histogram(durations_Slayer, hist_borders, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                title('Slayer, by Session', 'FontSize', 12);
            case 2
                histogram([], hist_borders, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'Zeppelin');
                histogram(durations_Zeppelin, hist_borders, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                title('Zeppelin, by Session', 'FontSize', 12);
            case 3
                histogram([], hist_borders, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'DisplayName', 'Emperor');
                histogram(durations_Emperor, hist_borders, 'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                title('Emperor, by Session', 'FontSize', 12);
        end
        % title('All Monkeys, by Session', 'FontSize', 12);
        xline(500, 'k-', 'LineWidth', 1, 'HandleVisibility', 'off');
        xline(threshold, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
        hold off;
        xlabel('Duration (s)');
        ylabel('Count');
        % title_str = sprintf('%d sessions, Open: %d trials, Close: %d trials', session_count(dataset_idx, 1), numel(durations{dataset_idx, 1}), numel(durations{dataset_idx, 2}));
        % title({dataset_labels{dataset_idx}, title_str}, 'FontSize', 10);
        legend('Location', 'northeast');
        xlim([0, 650]);
        ylim([0, 5]);
        xticks([0, 100, 200, 300, 400, 500, 562.5]);
        xticklabels({'0', '100', '200', '300', '400', '500', '(>500)'});
    end

    % global title
    sgtitle(sprintf('%s sessions, threshold = %d s', control, threshold));

    save_folder = fullfile(root, 'Figures', 'RestingDuration');
    check_path(save_folder);
    save_path = fullfile(save_folder, sprintf('resting_duration_histogram_%s%s.png', control, mode_str));
    saveas(f, save_path);
end


