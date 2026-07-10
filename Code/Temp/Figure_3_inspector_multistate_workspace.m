%% Figure_3_inspector_multistate_workspace.m
% Load corresponding connection vectors for multiple selected states.
% No plotting. The output remains in the MATLAB workspace for downstream analysis.
%
% Main outputs:
%   multi_state_data.J        : n_connection x n_state J matrix
%   multi_state_data.err      : n_connection x n_state error matrix
%   multi_state_data.cat      : n_connection x n_state category matrix (-1, 0, +1)
%   multi_state_data.sig      : n_connection x n_state significant mask
%   connection_table          : metadata for each row/connection
%   state_vectors             : per-state struct views of J/err/cat/sig
%
% Row correspondence is enforced. For every output matrix, row r refers to
% the same session, direction, from-neuron, and to-neuron across all states.

clear;

%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
if ~isfolder(fullfile(root, 'Data')) && isfolder(fullfile(pwd, 'Data'))
    root = pwd;
end
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));

%% Data selection parameters
% Add, remove, or reorder selected states here. Output columns follow this order.
selected_states = [ ...
    make_state_spec('Pre',  'RestOpen',  2, 'Pre_RestOpen_K2'), ...
    make_state_spec('Pre',  'RestClose', 2, 'Pre_RestClose_K2'), ...
    make_state_spec('Post', 'RestOpen',  2, 'Post_RestOpen_K2'), ...
    make_state_spec('Post', 'RestClose', 2, 'Post_RestClose_K2') ...
];

% Connection blocks to extract. The default reproduces the original ACC/VLPFC comparison.
area1_names = {'ACC'};
area2_names = {'VLPFC'};
area1_label = 'ACC';
area2_label = 'VLPFC';

% Metadata filter.
filter_opts = struct();
filter_opts.kernel_name = "DeltaPure";
filter_opts.align = 'Last';
filter_opts.area = "Cortex";
filter_opts.injection = 'Muscimol';
filter_opts.epoch = 3000;
filter_opts.fold_idx = 0;
filter_opts.shuffle_idx = 0;
filter_opts.resting_dur_threshold = 15;

err_multi = 1;
skip_failed_sessions = false;
max_sessions_to_include = inf;

% Optional disk output. The script variables are available in workspace either way.
save_workspace_mat = false;
workspace_mat_name = 'Figure3_multistate_workspace.mat';

%% Load metadata
mt = load_meta(root, 'table');
mt = mt.GLM;

selected_rows = default_metadata_filter(mt, filter_opts, selected_states);
selected_mt = mt(selected_rows, :);
meta_array = table2struct(selected_mt);
if isfinite(max_sessions_to_include)
    meta_array = meta_array(1:min(numel(meta_array), max_sessions_to_include));
end

if isempty(meta_array)
    error('No metadata rows selected.');
end

%% Load states and build aligned connection matrix
multi_state_data = empty_multistate_data(selected_states);
loaded_session_data = struct([]);
valid_session_count = 0;
failed_session_count = 0;

for session_i = 1:numel(meta_array)
    meta = meta_array(session_i);
    session_label = make_session_label(meta);
    fprintf('\n===== Loading %s (%d/%d) =====\n', session_label, session_i, numel(meta_array));

    try
        loaded_states = load_required_state_specs_for_session(root, meta, selected_states, ...
            area1_names, area2_names, area1_label, area2_label);

        session_data = make_multistate_vectors(loaded_states, selected_states, err_multi, meta, session_label);

        valid_session_count = valid_session_count + 1;
        loaded_session_data(valid_session_count).meta = meta; %#ok<SAGROW>
        loaded_session_data(valid_session_count).session_label = session_label; %#ok<SAGROW>
        loaded_session_data(valid_session_count).states = loaded_states; %#ok<SAGROW>
        loaded_session_data(valid_session_count).n_connection = size(session_data.J, 1); %#ok<SAGROW>

        multi_state_data = append_multistate_data(multi_state_data, session_data);
    catch ME
        failed_session_count = failed_session_count + 1;
        if skip_failed_sessions
            warning('Skipping %s because loading/processing failed: %s', session_label, ME.message);
            continue;
        else
            rethrow(ME);
        end
    end
end

if valid_session_count == 0 || isempty(multi_state_data.J)
    error('No valid multi-state connection data were loaded.');
end

multi_state_data.connection_table = make_connection_table(multi_state_data.connection_meta);
multi_state_data.sig = multi_state_data.cat ~= 0;

% Convenience workspace variables.
J = multi_state_data.J;
J_err = multi_state_data.err;
J_cat = multi_state_data.cat;
J_sig = multi_state_data.sig;
state_labels = multi_state_data.state_labels;
state_keys = multi_state_data.state_keys;
connection_table = multi_state_data.connection_table;
state_vectors = make_state_vectors(multi_state_data);

fprintf('\nLoaded sessions: %d. Failed sessions: %d.\n', valid_session_count, failed_session_count);
fprintf('Created aligned matrices: %d connections x %d states.\n', size(J, 1), size(J, 2));
fprintf('Workspace variables: multi_state_data, J, J_err, J_cat, J_sig, state_labels, state_keys, connection_table, state_vectors.\n');

if save_workspace_mat
    save_folder = fullfile(root, 'Data', 'Working', 'workspace');
    check_path(save_folder);
    save_path = fullfile(save_folder, workspace_mat_name);
    save(save_path, 'multi_state_data', 'J', 'J_err', 'J_cat', 'J_sig', ...
        'state_labels', 'state_keys', 'connection_table', 'state_vectors', '-v7.3');
    fprintf('Saved workspace data to:\n%s\n', save_path);
end

% Explicitly assign outputs to base workspace. This is redundant when run as a script,
% but useful if this file is invoked from another function using run().
assignin('base', 'multi_state_data', multi_state_data);
assignin('base', 'J', J);
assignin('base', 'J_err', J_err);
assignin('base', 'J_cat', J_cat);
assignin('base', 'J_sig', J_sig);
assignin('base', 'state_labels', state_labels);
assignin('base', 'state_keys', state_keys);
assignin('base', 'connection_table', connection_table);
assignin('base', 'state_vectors', state_vectors);
assignin('base', 'loaded_session_data', loaded_session_data);

%% Local functions
function spec = make_state_spec(prepost, state, kernel_idx, label)
    if nargin < 4 || isempty(label)
        label = sprintf('%s_%s_K%d', prepost, state, kernel_idx);
    end
    spec = struct();
    spec.prepost = char(string(prepost));
    spec.state = char(string(state));
    spec.kernel_idx = kernel_idx;
    spec.label = char(string(label));
    spec.key = state_key(prepost, state, kernel_idx);
end

function data = empty_multistate_data(selected_states)
    n_state = numel(selected_states);
    data = struct();
    data.state_specs = selected_states;
    data.state_labels = string({selected_states.label}).';
    data.state_keys = string({selected_states.key}).';
    data.J = zeros(0, n_state);
    data.err = zeros(0, n_state);
    data.cat = zeros(0, n_state);
    data.sig = false(0, n_state);
    data.connection_meta = empty_connection_meta();
    data.connection_table = table();
end

function meta = empty_connection_meta()
    meta = struct();
    meta.session_label = strings(0, 1);
    meta.animal_name = strings(0, 1);
    meta.injection = strings(0, 1);
    meta.session_idx = [];
    meta.connection_direction = strings(0, 1);
    meta.from_area = strings(0, 1);
    meta.to_area = strings(0, 1);
    meta.from_cell_area_idx = [];
    meta.to_cell_area_idx = [];
    meta.from_cell_global_idx = [];
    meta.to_cell_global_idx = [];
    meta.from_cell_idx = [];
    meta.to_cell_idx = [];
end

function out = append_multistate_data(out, incoming)
    out.J = [out.J; incoming.J];
    out.err = [out.err; incoming.err];
    out.cat = [out.cat; incoming.cat];
    out.sig = out.cat ~= 0;
    out.connection_meta = append_connection_meta(out.connection_meta, incoming.connection_meta);
end

function out = append_connection_meta(out, incoming)
    fields = fieldnames(out);
    for k = 1:numel(fields)
        field_name = fields{k};
        if isfield(incoming, field_name)
            out.(field_name) = [out.(field_name); incoming.(field_name)];
        end
    end
end

function selected_rows = default_metadata_filter(mt, filter_opts, selected_states)
    base_filter = strcmp(mt.kernel_name, filter_opts.kernel_name) & ...
                  strcmp(mt.align, filter_opts.align) & ...
                  strcmp(mt.area, filter_opts.area) & ...
                  strcmp(mt.injection, filter_opts.injection) & ...
                  (mt.epoch == filter_opts.epoch) & ...
                  (mt.fold_idx == filter_opts.fold_idx) & ...
                  (mt.shuffle_idx == filter_opts.shuffle_idx) & ...
                  cellfun(@(x) ~isempty(x) && x == filter_opts.resting_dur_threshold, mt.resting_dur_threshold);

    anchor = selected_states(1);
    anchor_filter = base_filter & ...
                    strcmp(mt.prepost, anchor.prepost) & ...
                    strcmp(mt.state, anchor.state);

    selected_rows = false(height(mt), 1);
    anchor_idx = find(anchor_filter);

    required_pairs = unique_state_pairs(selected_states);

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

            if iscell(mt.(field))
                anchor_value = mt.(field){idx};
                same_session = same_session & cellfun(@(x) isequal(x, anchor_value), mt.(field));
            else
                anchor_value = mt.(field)(idx);
                same_session = same_session & arrayfun(@(x) isequal(x, anchor_value), mt.(field));
            end
        end

        has_all_required = true;
        for q = 1:numel(required_pairs)
            pair = required_pairs(q);
            has_this = any(same_session & ...
                           strcmp(mt.prepost, pair.prepost) & ...
                           strcmp(mt.state, pair.state));

            if ~has_this
                has_all_required = false;
                break;
            end
        end

        if has_all_required
            selected_rows(idx) = true;
        end
    end

    fprintf('Metadata filter selected %d/%d anchor rows with all requested states.\n', sum(selected_rows), sum(base_filter));

    if ~any(selected_rows)
        warning('No complete requested state sets found. Falling back to anchor rows only. Missing states may fail during loading.');
        selected_rows = anchor_filter;
    end
end

function pairs = unique_state_pairs(selected_states)
    keys = strings(numel(selected_states), 1);
    for i = 1:numel(selected_states)
        keys(i) = sprintf('%s|%s', selected_states(i).prepost, selected_states(i).state);
    end
    [unique_keys, ia] = unique(keys, 'stable'); %#ok<ASGLU>
    pairs = struct([]);
    for i = 1:numel(ia)
        pairs(i).prepost = selected_states(ia(i)).prepost; %#ok<AGROW>
        pairs(i).state = selected_states(ia(i)).state; %#ok<AGROW>
    end
end

function label = make_session_label(meta)
    animal_str = char(string(meta.animal_name));
    injection_str = char(string(meta.injection));
    label = sprintf('%s %s session %s', animal_str, injection_str, char(string(meta.session_idx)));
end

function loaded_states = load_required_state_specs_for_session(root, meta, selected_states, area1_names, area2_names, area1_label, area2_label)
    loaded_states = struct();
    for s = 1:numel(selected_states)
        spec = selected_states(s);
        key = spec.key;
        if isfield(loaded_states, key)
            continue;
        end
        loaded_states.(key) = load_state_connectivity(root, meta, spec.prepost, spec.state, spec.kernel_idx, ...
            area1_names, area2_names, area1_label, area2_label);
    end
end

function key = state_key(prepost, state, kernel_idx)
    key = sprintf('k%d_%s_%s', kernel_idx, char(string(prepost)), char(string(state)));
    key = matlab.lang.makeValidName(key);
end

function state_struct = load_state_connectivity(root, meta, prepost, state, kernel_idx, area1_names, area2_names, area1_label, area2_label)
    state_struct = struct();
    state_struct.prepost = prepost;
    state_struct.state = state;
    state_struct.kernel_idx = kernel_idx;
    state_struct.area1_names = area1_names;
    state_struct.area2_names = area2_names;
    state_struct.area1_label = area1_label;
    state_struct.area2_label = area2_label;

    meta.prepost = prepost;
    meta.state = state;

    meta.filename = generate_filename('raster', meta);
    raster_path = fullfile(root, 'Data', 'Working', 'raster', meta.filename);
    if ~isfile(raster_path)
        error('Raster file not found: %s', raster_path);
    end
    raster_data = load(raster_path);
    fprintf('Loaded raster data for %s %s %s\n', meta.prepost, meta.state, meta.area);

    cell_area = raster_data.data.cell_area;
    filter1 = ismember(cell_area, area1_names);
    filter2 = ismember(cell_area, area2_names);
    if ~any(filter1)
        error('No cells found for area1: %s', join_area_names(area1_names));
    end
    if ~any(filter2)
        error('No cells found for area2: %s', join_area_names(area2_names));
    end

    meta.filename = generate_filename('GLM', meta);
    glm_path = fullfile(root, 'Data', 'Working', 'GLM', meta.filename);
    if ~isfile(glm_path)
        error('GLM file not found: %s', glm_path);
    end
    GLM_data = load(glm_path);
    fprintf('Loaded GLM data for %s %s %s K%d\n', meta.prepost, meta.state, meta.area, kernel_idx);

    N = GLM_data.meta.N;
    J = GLM_data.data.model_par(:, ((2 + N * (kernel_idx - 1)) : (1 + N * kernel_idx)));
    err = GLM_data.data.model_err.total(:, ((2 + N * (kernel_idx - 1)) : (1 + N * kernel_idx)));

    idx1 = find(filter1);
    idx2 = find(filter2);

    state_struct.cell_area = cell_area;
    state_struct.filter1 = filter1;
    state_struct.filter2 = filter2;
    state_struct.idx1 = idx1;
    state_struct.idx2 = idx2;
    state_struct.J = J;
    state_struct.err = err;
    state_struct.J12 = J(filter1, filter2);
    state_struct.J21 = J(filter2, filter1);
    state_struct.err12 = err(filter1, filter2);
    state_struct.err21 = err(filter2, filter1);
end

function data = make_multistate_vectors(loaded_states, selected_states, err_multi, meta, session_label)
    n_state = numel(selected_states);

    first_state = loaded_states.(selected_states(1).key);
    [J0, err0, connection_meta] = vectorize_state_connectivity(first_state, meta, session_label);
    n_connection = numel(J0);

    J = nan(n_connection, n_state);
    err = nan(n_connection, n_state);
    J(:, 1) = J0;
    err(:, 1) = err0;

    for s = 2:n_state
        state_s = loaded_states.(selected_states(s).key);
        validate_matching_filters(first_state, state_s, selected_states(1).label, selected_states(s).label);
        [Js, errs, connection_meta_s] = vectorize_state_connectivity(state_s, meta, session_label);
        validate_connection_identity(connection_meta, connection_meta_s, selected_states(1).label, selected_states(s).label);
        J(:, s) = Js;
        err(:, s) = errs;
    end

    valid = all(isfinite(J), 2) & all(isfinite(err), 2);
    if ~all(valid)
        fprintf('Dropping %d/%d connections with NaN/Inf in at least one selected state for %s.\n', ...
            sum(~valid), n_connection, session_label);
    end

    J = J(valid, :);
    err = err(valid, :);
    connection_meta = subset_connection_meta(connection_meta, valid);

    cat = zeros(size(J));
    cat(J >  err_multi * err) = 1;
    cat(J < -err_multi * err) = -1;

    data = struct();
    data.J = J;
    data.err = err;
    data.cat = cat;
    data.sig = cat ~= 0;
    data.connection_meta = connection_meta;
end

function [J_vec, err_vec, meta_out] = vectorize_state_connectivity(state_struct, meta, session_label)
    J12 = state_struct.J12(:);
    J21 = state_struct.J21(:);
    err12 = state_struct.err12(:);
    err21 = state_struct.err21(:);

    % J(i,j) represents the connection from neuron j to neuron i.
    [to12_global, from12_global, to12_area, from12_area] = make_connection_index_vectors(state_struct.idx1, state_struct.idx2);
    [to21_global, from21_global, to21_area, from21_area] = make_connection_index_vectors(state_struct.idx2, state_struct.idx1);

    J_vec = [J12; J21];
    err_vec = [err12; err21];

    from_cell_global_idx = [from12_global; from21_global];
    to_cell_global_idx = [to12_global; to21_global];
    from_cell_area_idx = [from12_area; from21_area];
    to_cell_area_idx = [to12_area; to21_area];

    n12 = numel(J12);
    n21 = numel(J21);
    from_area = [repmat(string(state_struct.area2_label), n12, 1); repmat(string(state_struct.area1_label), n21, 1)];
    to_area = [repmat(string(state_struct.area1_label), n12, 1); repmat(string(state_struct.area2_label), n21, 1)];
    dir12 = sprintf('%s_to_%s', sanitize_token(state_struct.area2_label), sanitize_token(state_struct.area1_label));
    dir21 = sprintf('%s_to_%s', sanitize_token(state_struct.area1_label), sanitize_token(state_struct.area2_label));
    connection_direction = [repmat(string(dir12), n12, 1); repmat(string(dir21), n21, 1)];

    n = numel(J_vec);
    meta_out = empty_connection_meta();
    meta_out.session_label = repmat(string(session_label), n, 1);
    meta_out.animal_name = repmat(string(meta.animal_name), n, 1);
    meta_out.injection = repmat(string(meta.injection), n, 1);
    meta_out.session_idx = repmat(double(meta.session_idx), n, 1);
    meta_out.connection_direction = connection_direction;
    meta_out.from_area = from_area;
    meta_out.to_area = to_area;
    meta_out.from_cell_area_idx = from_cell_area_idx;
    meta_out.to_cell_area_idx = to_cell_area_idx;
    meta_out.from_cell_global_idx = from_cell_global_idx;
    meta_out.to_cell_global_idx = to_cell_global_idx;
    meta_out.from_cell_idx = from_cell_global_idx;
    meta_out.to_cell_idx = to_cell_global_idx;
end

function [row_global_vec, col_global_vec, row_area_vec, col_area_vec] = make_connection_index_vectors(row_idx, col_idx)
    [row_area_grid, col_area_grid] = ndgrid((1:numel(row_idx)).', (1:numel(col_idx)).');

    row_area_vec = row_area_grid(:);
    col_area_vec = col_area_grid(:);

    row_global_vec = row_idx(row_area_vec);
    col_global_vec = col_idx(col_area_vec);

    row_global_vec = row_global_vec(:);
    col_global_vec = col_global_vec(:);
end

function meta_out = subset_connection_meta(meta_in, mask)
    meta_out = struct();
    fields = fieldnames(meta_in);
    for k = 1:numel(fields)
        field = fields{k};
        values = meta_in.(field);
        if numel(values) == numel(mask)
            meta_out.(field) = values(mask);
        else
            meta_out.(field) = values;
        end
    end
end

function tbl = make_connection_table(connection_meta)
    tbl = table(connection_meta.session_label, ...
                connection_meta.animal_name, ...
                connection_meta.injection, ...
                connection_meta.session_idx, ...
                connection_meta.connection_direction, ...
                connection_meta.from_area, ...
                connection_meta.to_area, ...
                connection_meta.from_cell_area_idx, ...
                connection_meta.to_cell_area_idx, ...
                connection_meta.from_cell_global_idx, ...
                connection_meta.to_cell_global_idx, ...
                'VariableNames', {'session_label', 'animal_name', 'injection', 'session_idx', ...
                                  'connection_direction', 'from_area', 'to_area', ...
                                  'from_cell_area_idx', 'to_cell_area_idx', ...
                                  'from_cell_global_idx', 'to_cell_global_idx'});
end

function state_vectors = make_state_vectors(data)
    state_vectors = struct();
    for s = 1:numel(data.state_labels)
        field_name = matlab.lang.makeValidName(char(data.state_labels(s)));
        state_vectors.(field_name) = struct();
        state_vectors.(field_name).label = data.state_labels(s);
        state_vectors.(field_name).key = data.state_keys(s);
        state_vectors.(field_name).J = data.J(:, s);
        state_vectors.(field_name).err = data.err(:, s);
        state_vectors.(field_name).cat = data.cat(:, s);
        state_vectors.(field_name).sig = data.cat(:, s) ~= 0;
        state_vectors.(field_name).connection_table = data.connection_table;
    end
end

function validate_matching_filters(state_a, state_b, label_a, label_b)
    if numel(state_a.filter1) ~= numel(state_b.filter1) || any(state_a.filter1(:) ~= state_b.filter1(:)) || ...
       numel(state_a.filter2) ~= numel(state_b.filter2) || any(state_a.filter2(:) ~= state_b.filter2(:))
        error('Area filters do not match between %s and %s.', label_a, label_b);
    end

    if ~isequal(size(state_a.J12), size(state_b.J12)) || ~isequal(size(state_a.J21), size(state_b.J21))
        error('Connectivity block sizes do not match between %s and %s.', label_a, label_b);
    end
end

function validate_connection_identity(meta_a, meta_b, label_a, label_b)
    fields_to_check = {'session_label', 'connection_direction', 'from_area', 'to_area', ...
                       'from_cell_area_idx', 'to_cell_area_idx', ...
                       'from_cell_global_idx', 'to_cell_global_idx'};
    for f = 1:numel(fields_to_check)
        field = fields_to_check{f};
        a = meta_a.(field);
        b = meta_b.(field);
        if ~isequal(a, b)
            error('Connection identity field "%s" does not match between %s and %s.', field, label_a, label_b);
        end
    end
end

function s = sanitize_token(name)
    s = regexprep(char(string(name)), '[^A-Za-z0-9]+', '_');
    s = regexprep(s, '^_+|_+$', '');
end

function joined = join_area_names(area_names)
    joined = strjoin(cellstr(string(area_names)), ', ');
end
