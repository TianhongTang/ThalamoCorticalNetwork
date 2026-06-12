%% Figure 3 network-only, all sessions across columns, loop kernels 1-3
% Each column is one selected session.
% Each row is one state:
%   1. Pre  RestOpen
%   2. Pre  RestClose
%   3. Post RestOpen
%   4. Post RestClose
% A separate figure is exported for each kernel_idx = 1:3.

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

%% State definitions
preposts = {'Pre', 'Pre', 'Post', 'Post'};
states   = {'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};
row_labels = {'Pre | RestOpen', 'Pre | RestClose', 'Post | RestOpen', 'Post | RestClose'};
n_state = numel(states);

%% Parameters
kernel_list = 1:3;
network_err_multi = 2;

max_sessions_to_include = inf; % set smaller while debugging.
figure_visible = 'off';
show_row_ylabel_only_first_col = true;

%% Load and filter metadata
mt = load_meta(root, 'table');
mt = mt.GLM;

% -------------------------------------------------------------------------
% EDIT THIS FILTER FOR THE FINAL SESSION SET.
% Examples:
% selected_rows = strcmp(string(mt.animal_name), "Slayer") & strcmp(string(mt.injection), "Muscimol");
% selected_rows = mt.session_idx >= 1 & strcmp(string(mt.injection), "Muscimol");
% -------------------------------------------------------------------------
selected_rows = default_metadata_filter(mt);

selected_mt = mt(selected_rows, :);
meta_array = table2struct(selected_mt);
if isfinite(max_sessions_to_include)
    meta_array = meta_array(1:min(numel(meta_array), max_sessions_to_include));
end

if isempty(meta_array)
    error('No metadata rows selected.');
end

n_session = numel(meta_array);

fprintf('Selected %d sessions for plotting.\n', n_session);
for session_i = 1:n_session
    fprintf('  Column %d: %s\n', session_i, make_column_label(meta_array(session_i)));
end

%% Main loop: one figure per kernel
for kernel_idx = kernel_list
    fprintf('\n================ KERNEL %d ================\n', kernel_idx);

    f = figure('Color', 'w', 'Visible', figure_visible);
    tiledlayout(n_state, n_session, "TileSpacing", "Compact", "Padding", "Compact");

    for session_i = 1:n_session
        meta = meta_array(session_i);
        column_label = make_column_label(meta);
        fprintf('\n===== Loading %s | kernel %d (%d/%d) =====\n', ...
            column_label, kernel_idx, session_i, n_session);

        try
            state_data = struct([]);
            for state_i = 1:n_state
                sd = load_state_connectivity(root, meta, preposts{state_i}, states{state_i}, kernel_idx);
                if state_i == 1
                    state_data = sd;
                else
                    state_data = [state_data; sd]; %#ok<AGROW>
                end
            end

            for state_i = 1:n_state
                tile_idx = (state_i - 1) * n_session + session_i;
                ax = nexttile(tile_idx);

                call_plot_network(ax, ...
                    state_data(state_i).J12, state_data(state_i).J21, ...
                    state_data(state_i).err12, state_data(state_i).err21, ...
                    network_err_multi, [], []);

                if state_i == 1
                    title(ax, column_label, 'Interpreter', 'none', 'FontWeight', 'bold');
                end

                if show_row_ylabel_only_first_col && session_i == 1
                    ylabel(ax, row_labels{state_i}, 'Interpreter', 'none', 'FontWeight', 'bold');
                end
            end

        catch ME
            warning('Failed to load %s | kernel %d: %s', column_label, kernel_idx, ME.message);

            for state_i = 1:n_state
                tile_idx = (state_i - 1) * n_session + session_i;
                ax = nexttile(tile_idx);
                axis(ax, 'off');

                if state_i == 1
                    title(ax, column_label, 'Interpreter', 'none', 'FontWeight', 'bold');
                end

                if show_row_ylabel_only_first_col && session_i == 1
                    ylabel(ax, row_labels{state_i}, 'Interpreter', 'none', 'FontWeight', 'bold');
                end

                text(ax, 0.5, 0.5, sprintf('Load failed\n%s', ME.message), ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Interpreter', 'none', ...
                    'Color', 'r', ...
                    'FontSize', 8);
            end
        end
    end

    sgtitle(sprintf('Network plots | Kernel %d | columns = sessions', kernel_idx), ...
        'Interpreter', 'none');

    %% Export
    fig = gcf;
    save_folder = fullfile(root, 'Figures', 'Paper');
    check_path(save_folder);

    figWidth  = max(4.0 * n_session, 8.0); % inches
    figHeight = 16.0; % inches
    resolution = 300;

    set(fig, 'Units', 'inches');
    fig.Position(3:4) = [figWidth, figHeight];

    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperSize', [figWidth, figHeight]);
    set(fig, 'PaperPosition', [0, 0, figWidth, figHeight]);
    set(fig, 'Color', 'w');

    output_stub = sprintf('Figure3_network_all_sessions_k%d', kernel_idx);

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

%% functions

function selected_rows = default_metadata_filter(mt)
    base_filter = strcmp(mt.kernel_name, "DeltaPure") & ...
                  strcmp(mt.align, 'Last') & ...
                  strcmp(mt.area, "Cortex") & ...
                  strcmp(mt.injection, 'Muscimol') & ...
                  (mt.epoch == 3000) & ...
                  (mt.fold_idx == 0) & ...
                  (mt.shuffle_idx == 0) & ...
                  cellfun(@(x) ~isempty(x) && x == 15, mt.resting_dur_threshold);

    anchor_filter = base_filter & ...
                    strcmp(mt.prepost, 'Pre') & ...
                    strcmp(mt.state, 'RestOpen');

    selected_rows = false(height(mt), 1);
    anchor_idx = find(anchor_filter);

    required_preposts = {'Pre', 'Pre', 'Post', 'Post'};
    required_states   = {'RestOpen', 'RestClose', 'RestOpen', 'RestClose'};

    match_fields = {'animal_name', 'injection', 'align', 'session_idx', ...
                    'resting_dur_threshold', 'area', 'kernel_name', ...
                    'reg_name', 'epoch', 'fold_idx', 'shuffle_idx'};

    for k = 1:numel(anchor_idx)
        idx = anchor_idx(k);

        same_session = base_filter;

        for f = 1:numel(match_fields)
            field = match_fields{f};
            if ~ismember(field, mt.Properties.VariableNames)
                continue;
            end

            anchor_value = mt.(field)(idx);

            if iscell(mt.(field))
                anchor_value = anchor_value{1};
                same_session = same_session & cellfun(@(x) isequal(x, anchor_value), mt.(field));
            else
                same_session = same_session & arrayfun(@(x) isequal(x, anchor_value), mt.(field));
            end
        end

        has_all_required = true;
        for q = 1:numel(required_preposts)
            has_this = any(same_session & ...
                           strcmp(mt.prepost, required_preposts{q}) & ...
                           strcmp(mt.state, required_states{q}));

            if ~has_this
                has_all_required = false;
                anchor_name = mt.file_name{idx};
                warning('Session %s is missing required Pre/Post x Open/Close combination: %s %s. Skipping this session.', ...
                    anchor_name, required_preposts{q}, required_states{q});
                break;
            end
        end

        if has_all_required
            selected_rows(idx) = true;
        end
    end

    fprintf('Default metadata filter selected %d/%d rows:\n', sum(selected_rows), sum(base_filter));
    for k = 1:height(mt)
        if selected_rows(k)
            fprintf('  %s\n', mt.file_name{k});
        end
    end

    if ~any(selected_rows)
        warning('Default metadata filter selected no complete Pre/Post x Open/Close sets. Falling back to Pre RestOpen anchors.');
        selected_rows = anchor_filter;
    end
end

function label = make_column_label(meta)
    animal_str = char(string(meta.animal_name));
    session_str = char(string(meta.session_idx));
    label = sprintf('%s | session %s', animal_str, session_str);
end

function state_struct = load_state_connectivity(root, meta, prepost, state, kernel_idx)
    state_struct = struct();
    state_struct.prepost = prepost;
    state_struct.state = state;

    meta.prepost = prepost;
    meta.state = state;

    meta.filename = generate_filename('raster', meta);
    raster_data = load(fullfile(root, 'Data', 'Working', 'raster', meta.filename));
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);
    fprintf('Trial_len: %d, N: %d\n', raster_data.meta.trial_len, raster_data.meta.N);
    fprintf('Trial_num: %d\n', raster_data.meta.trial_num);

    cell_area = raster_data.data.cell_area;
    filter1 = ismember(cell_area, {'ACC'});
    filter2 = ismember(cell_area, {'VLPFC'});

    meta.filename = generate_filename('GLM', meta);
    GLM_data = load(fullfile(root, 'Data', 'Working', 'GLM', meta.filename));
    fprintf('Loaded GLM data for %s %s %s\n', meta.prepost, meta.state, meta.area);

    N = GLM_data.meta.N;
    J = GLM_data.data.model_par(:, ((2 + N * (kernel_idx - 1)) : (1 + N * kernel_idx)));
    err = GLM_data.data.model_err.total(:, ((2 + N * (kernel_idx - 1)) : (1 + N * kernel_idx)));

    state_struct.cell_area = cell_area;
    state_struct.filter1 = filter1;
    state_struct.filter2 = filter2;
    state_struct.J = J;
    state_struct.err = err;
    state_struct.J12 = J(filter1, filter2);
    state_struct.J21 = J(filter2, filter1);
    state_struct.err12 = err(filter1, filter2);
    state_struct.err21 = err(filter2, filter1);
end

function call_plot_network(ax, J12, J21, err12, err21, err_multi, highlight_i, highlight_j)
    if ~isempty(highlight_i) && ~isempty(highlight_j)
        try
            plot_network(ax, J12, J21, err12, err21, err_multi, highlight_i, highlight_j);
            return;
        catch ME_highlight
            warning('plot_network with highlight failed: %s. Retrying without highlight.', ME_highlight.message);
        end
    end

    try
        plot_network(ax, J12, J21, err12, err21, err_multi);
    catch ME_no_highlight
        if ~isempty(J12)
            warning('plot_network without highlight failed: %s. Retrying with fallback highlight (1,1).', ME_no_highlight.message);
            plot_network(ax, J12, J21, err12, err21, err_multi, 1, 1);
        else
            rethrow(ME_no_highlight);
        end
    end
end
