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

f = figure('Position', [100, 100, 1600, 800], 'Visible', 'off');
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

            % plot histogram
            nexttile(tile_idx);
            ratios = [open_count/open_total, close_count/close_total]*100;
            bar(1:2, ratios);
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





