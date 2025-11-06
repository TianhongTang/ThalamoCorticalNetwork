%% shuffle spike trains and keep the firing rate.
% shuffle_type: "None", "Within trial", "Across trial"

%% Get root folder
code_depth = 4;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
% include code folder and utils
% addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Main
reg_levels = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5];
% reg_names = {'Mix', 'L1', 'L2'};
reg_names = {'L2'};
regs = cell(length(reg_levels) * length(reg_names), 1);
for i = 1:length(reg_levels)
    reg_level = reg_levels(i);
    for j = 1:length(reg_names)
        reg_name = reg_names{j};
        reg = struct();
        switch reg_name
            case 'Mix'
                reg.l1 = reg_level;
                reg.l2 = reg_level;
            case 'L1'
                reg.l1 = reg_level;
                reg.l2 = 0;
            case 'L2'
                reg.l1 = 0;
                reg.l2 = reg_level;
        end
        reg.name = sprintf('%s=%d', reg_name, reg_level*100);
        regs{(i-1)*length(reg_names) + j} = reg;
    end
end
controls = {'Muscimol'};
session_idxs_all = {1:10}; % session indices for each control
area_types = {'Cortex'};
prepost_types = {'Pre', 'Post'};
states = {'RestOpen', 'RestClose', 'Task'};
alignments = {'Last'};
kernel = 'DeltaPure';

tasks = {};
found_count = 0;
total_count = 0;
% check valid files
for control_idx = 1:length(controls)
    control = controls{control_idx};
    sessions = session_idxs_all{control_idx};
    for area_idx = 1:length(area_types)
        area_type = area_types{area_idx};
        for prepost_idx = 1:length(prepost_types)
            prepost = prepost_types{prepost_idx};
            if strcmp(area_type, 'Full') && strcmp(prepost, 'Post')
                % All thalamus data in post sessions are not available
                continue;
            end
            for align_idx = 1:length(alignments)
                alignment = alignments{align_idx};
                for state_idx = 1:length(states)
                    state = states{state_idx};

                    loss_fig = figure('Position', [100, 100, 800, 1200]);
                    % 5 row 2 columns on loss figure
                    t = tiledlayout(5, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

                    for session_idx = sessions
                        reg_train_loss_curve = zeros(1, length(regs));
                        reg_test_loss_curve = zeros(1, length(regs));
                        best_loss_idx = zeros(1, length(regs));

                        all_train_loss_curve = zeros(25, length(regs));
                        all_test_loss_curve = zeros(25, length(regs));

                        for reg_idx = 1:length(regs)
                            train_loss_curve = zeros(25, 3);
                            test_loss_curve = zeros(25, 3);
                            for epoch = 100:100:2500
                                for fold_idx = 1:3
                                    reg = regs{reg_idx};

                                    task = struct();
                                    task.control = control;
                                    task.area_type = area_type;
                                    task.prepost = prepost;
                                    task.state = state;
                                    task.alignment = alignment;
                                    task.name = [control, prepost, state, area_type, 'Align', alignment];
                                    task.session_idx = session_idx;
                                    task.kernel = kernel;
                                    task.reg = reg;
                                    task.shuffle_size = 0; % number of shuffles
                                    task.epoch = epoch;
                                    tasks{end+1} = task;

                                    folder_name = fullfile(root, 'Data', 'Working', 'GLM_models');
                                    file_name = sprintf('GLM_%s_s%d_shuffle%d_%s_%s_epoch%d_fold%d.mat', ...
                                        task.name, task.session_idx, 0, task.kernel, task.reg.name, task.epoch, fold_idx);
                                    file_path = fullfile(folder_name, file_name);
                                    total_count = total_count + 1;
                                    if isfile(file_path)
                                        found_count = found_count + 1;
                                    else
                                        fprintf('Missing file: %s\n', file_path);
                                    end

                                    load(file_path, 'train_loss', 'test_loss');
                                    train_loss_curve(epoch/100, fold_idx) = train_loss.minuslogL;
                                    test_loss_curve(epoch/100, fold_idx) = test_loss.minuslogL;
                                end
                            end
                            train_loss_total = sum(train_loss_curve, 2) / 2;
                            test_loss_total = sum(test_loss_curve, 2);

                            all_train_loss_curve(:, reg_idx) = train_loss_total;
                            all_test_loss_curve(:, reg_idx) = test_loss_total;

                            best_test_loss = min(test_loss_total);
                            best_epoch_idx = find(test_loss_total == best_test_loss, 1, 'first');

                            reg_train_loss_curve(reg_idx) = train_loss_total(best_epoch_idx);
                            reg_test_loss_curve(reg_idx) = best_test_loss;
                            best_loss_idx(reg_idx) = best_epoch_idx * 100;
                        end
                        % save results
                        folder_name = fullfile(root, 'Data', 'Working', 'reg_analysis');
                        check_path(folder_name);
                        save_name = sprintf('reg_analysis_%s_s%d_%s_%s_%s.mat', ...
                            [control, prepost, state, area_type, 'Align', alignment], task.session_idx, task.kernel, 'NoShuffle', state);
                        % save(fullfile(folder_name, save_name), ...
                        %     'all_train_loss_curve', all_train_loss_curve, ...
                        %     'all_test_loss_curve', all_test_loss_curve, ...
                        %     'reg_train_loss_curve', reg_train_loss_curve, ...
                        %     'reg_test_loss_curve', reg_test_loss_curve, ...
                        %     'best_loss_idx', best_loss_idx);

                        % plot results: train loss total and test loss total
                        figure(loss_fig);
                        nexttile;
                        hold on;
                        % plot(1:length(regs), reg_train_loss_curve', 'b');
                        % plot(1:length(regs), reg_test_loss_curve', 'r');
                        plot(1:length(regs), (reg_test_loss_curve' - reg_train_loss_curve'), 'k');
                        hold off;
                        % legend('Best Train Loss', 'Best Test Loss');
                        legend('Test - Train');
                        title(sprintf('Session %d, %s', session_idx, [control, prepost, state, area_type, 'Align', alignment]));
                        xticks(1:length(regs));
                        xticklabels(cellfun(@(r) r.name, regs, 'UniformOutput', false));
                        xlabel('Regularization');
                        ylabel('Loss (minus log-likelihood)');
                    end
                    % save loss figure
                    folder_name = fullfile(root, 'Figures', 'reg_analysis');
                    check_path(folder_name);
                    save_name = sprintf('reg_analysis_diff_%s.png', ...
                        [control, prepost, state, area_type, 'Align', alignment]);
                    saveas(loss_fig, fullfile(folder_name, save_name));
                    close(loss_fig);
                end
            end
        end
    end
end

