%% simulation_hybrid_network_v1.m
% Hybrid GLM network simulation.
%
% Two cortical-network modes are supported in one file:
%
%   config.network_source = 'load'
%       Load cortical baseline and cortical connections from an inferred GLM,
%       normally a Post-Muscimol model.
%
%   config.network_source = 'generate'
%       Load only cell counts/area labels, then generate all cortical
%       baselines and cortical connections from configurable random rules.
%
% Thalamocortical connections are controlled independently:
%
%   config.tc_source = 'load_pre'
%       Attempt to copy thalamocortical connections from a Pre GLM model.
%       This requires the Pre model to contain identifiable thalamic cells.
%
%   config.tc_source = 'generate'
%       Generate thalamocortical connections from configurable random rules.
%
%   config.tc_source = 'none'
%       Do not add thalamocortical connections.
%
% IMPORTANT:
%   Project-specific file formats are isolated in adapter functions marked
%   "PROJECT ADAPTER". Review those functions before a production run.
%
% Output:
%   Data/Working/simulation_v2/<simulation_tag>/
%
% The saved MAT file contains:
%   config      - complete simulation configuration
%   sim_meta    - simulation metadata and area labels
%   network     - h, J, masks, and provenance
%   results     - raster and firing-rate outputs for each simulated condition
%
% This first version saves a self-contained simulation file. It does not
% automatically register the result as a project raster dataset because the
% exact current raster-registration schema is not known.

clear;
clc;

run_tic = tic;
progress_log('SCRIPT', 'Started.');

%% Project root
code_depth = 4;
script_path = mfilename('fullpath');
root = script_path;
for depth_i = 1:code_depth
    root = fileparts(root);
end
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));
progress_log('SCRIPT', 'Project root: %s', root);

%% User configuration
config = default_simulation_config();

% -------------------------------------------------------------------------
% REFERENCE SESSION
% Fill these fields so generate_filename() can locate the current project
% raster and GLM files.
% -------------------------------------------------------------------------
config.reference.animal_name = 'Slayer';
config.reference.injection = 'Muscimol';
config.reference.session_idx = 1;
config.reference.state = 'RestOpen';
config.reference.area = 'Cortex';
config.reference.align = 'Last';
config.reference.resting_dur_threshold = 15;

% GLM hyperparameters used when loading inferred connections.
config.reference.kernel_name = 'DeltaPure';
config.reference.kernel_num = 3;
config.reference.reg_name = 'L2=0_2';
config.reference.epoch = 3000;
config.reference.fold_idx = 0;
config.reference.shuffle_idx = 0;

% -------------------------------------------------------------------------
% MAIN MODE SWITCHES
% -------------------------------------------------------------------------
% 'load'     = load cortical h and J from the reference Post GLM.
% 'generate' = copy cell counts only and generate cortical h and J.
config.network_source = 'load';

% 'load_pre' = copy thalamocortical J from a Pre GLM if available.
% 'generate' = generate thalamocortical J.
% 'none'     = no thalamocortical J.
config.tc_source = 'generate';

% 'data'   = derive ACC/VLPFC counts and ordering from the reference raster.
% 'manual' = use config.manual_cell_counts.
config.cell_count_source = 'data';

% -------------------------------------------------------------------------
% CELL COUNTS AND AREA LABELS
% -------------------------------------------------------------------------
config.manual_cell_counts.ACC = 40;
config.manual_cell_counts.VLPFC = 40;
config.thalamus_cells_per_area = 10;

% Labels used to identify cells in project raster files.
config.area_labels.ACC = {'ACC'};
config.area_labels.VLPFC = {'VLPFC'};
config.area_labels.thalamus_ACC = {'Thal-ACC', 'Thalamus-ACC'};
config.area_labels.thalamus_VLPFC = {'Thal-VLPFC', 'Thalamus-VLPFC'};

% -------------------------------------------------------------------------
% SIMULATION SIZE
% Start with small values while validating file adapters and model scaling.
% Production values can be increased to n_trials=30 and n_time_bins=30000.
% -------------------------------------------------------------------------
config.simulation.n_trials = 5;
config.simulation.n_time_bins = 5000;
config.simulation.parallel_trials = false;
config.simulation.base_seed = 137;
config.simulation.save_cortex_only = true;

% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
config.output.simulation_tag = sprintf( ...
    '%s_s%d_%s_%s_%s', ...
    config.reference.animal_name, ...
    config.reference.session_idx, ...
    config.reference.state, ...
    config.network_source, ...
    config.tc_source);
config.output.overwrite = false;
config.output.save_rasters = true;
config.output.save_network = true;

%% Validate configuration and run
validate_simulation_config(config);
progress_log('CONFIG', ...
    'network_source=%s, tc_source=%s, kernel=%s, reg=%s, K=%d.', ...
    config.network_source, config.tc_source, ...
    config.reference.kernel_name, config.reference.reg_name, ...
    config.reference.kernel_num);

[sim_meta, network, kernels] = build_simulation_network(root, config);
results = simulate_all_conditions(config, sim_meta, network, kernels);
save_simulation_output(root, config, sim_meta, network, results);

progress_log('SUMMARY', 'Finished in %.1f s.', toc(run_tic));

%% =========================================================================
%  CONFIGURATION
%  =========================================================================
function config = default_simulation_config()
    config = struct();

    config.reference = struct();
    config.reference.animal_name = '';
    config.reference.injection = 'Muscimol';
    config.reference.session_idx = 1;
    config.reference.state = 'RestOpen';
    config.reference.area = 'Cortex';
    config.reference.align = 'Last';
    config.reference.resting_dur_threshold = 15;
    config.reference.kernel_name = 'DeltaPure';
    config.reference.kernel_num = 3;
    config.reference.reg_name = 'L2=0_2';
    config.reference.epoch = 3000;
    config.reference.fold_idx = 0;
    config.reference.shuffle_idx = 0;

    config.network_source = 'load';
    config.tc_source = 'generate';
    config.cell_count_source = 'data';

    config.manual_cell_counts = struct('ACC', 40, 'VLPFC', 40);
    config.thalamus_cells_per_area = 10;

    config.area_labels = struct();
    config.area_labels.ACC = {'ACC'};
    config.area_labels.VLPFC = {'VLPFC'};
    config.area_labels.thalamus_ACC = {'Thal-ACC', 'Thalamus-ACC'};
    config.area_labels.thalamus_VLPFC = {'Thal-VLPFC', 'Thalamus-VLPFC'};

    % GLM column layout.
    %
    % Current analysis scripts generally use:
    %   baseline = model_par(:, 1)
    %   kernel k connection block =
    %       model_par(:, 2 + N*(k-1) : 1 + N*k)
    %
    % If the stored model contains post-spike-history columns before the
    % connection blocks, change connection_start_col accordingly.
    config.glm_layout.baseline_col = 1;
    config.glm_layout.connection_start_col = 2;

    % Kernel-basis file adapter.
    config.kernel_basis.source = 'file'; % 'file' or 'manual'
    config.kernel_basis.manual = {};
    config.kernel_basis.manual_dt = 1;

    % Generated cortical network.
    %
    % Matrices use target-area rows and source-area columns:
    %   [ACC<-ACC,     ACC<-VLPFC
    %    VLPFC<-ACC,   VLPFC<-VLPFC]
    config.generate.cortex.connection_probability = [0.08, 0.04; 0.04, 0.08];
    config.generate.cortex.positive_probability = [0.75, 0.65; 0.65, 0.75];
    config.generate.cortex.weight_abs_mean = [1.5, 0.8; 0.8, 1.5];
    config.generate.cortex.weight_cv = 0.20;
    config.generate.cortex.kernel_scale = [1.0, 0.6, 0.3];
    config.generate.cortex.baseline_mean = [-4.0, -4.0];
    config.generate.cortex.baseline_sd = [0.4, 0.4];

    % Generated thalamic cells and thalamocortical connections.
    config.generate.thalamus.baseline_mean = -4.5;
    config.generate.thalamus.baseline_sd = 0.3;
    config.generate.tc.connection_probability = 0.10;
    config.generate.tc.weight_mean = 3.0;
    config.generate.tc.weight_cv = 0.20;
    config.generate.tc.kernel_scale = [1.0, 0.6, 0.3];

    % Optional generated connections not used by default.
    config.generate.allow_corticothalamic = false;
    config.generate.allow_thalamus_recurrent = false;

    % Simulated conditions.
    conditions = struct([]);
    conditions(1).name = 'FittedAsync';
    conditions(1).enable_tc = true;
    conditions(1).sync_thalamus = false;
    conditions(1).sync_period = 100;
    conditions(1).sync_on_fraction = 1/8;
    conditions(1).sync_on_offset = 1.5;
    conditions(1).sync_off_offset = -0.6;

    conditions(2) = conditions(1);
    conditions(2).name = 'FittedSync';
    conditions(2).sync_thalamus = true;

    conditions(3) = conditions(1);
    conditions(3).name = 'FittedNoInput';
    conditions(3).enable_tc = false;
    conditions(3).sync_thalamus = false;

    config.conditions = conditions;

    config.simulation.n_trials = 5;
    config.simulation.n_time_bins = 5000;
    config.simulation.parallel_trials = false;
    config.simulation.base_seed = 137;
    config.simulation.save_cortex_only = true;

    config.output.simulation_tag = 'simulation';
    config.output.overwrite = false;
    config.output.save_rasters = true;
    config.output.save_network = true;
end

function validate_simulation_config(config)
    valid_network_sources = {'load', 'generate'};
    valid_tc_sources = {'load_pre', 'generate', 'none'};
    valid_count_sources = {'data', 'manual'};

    if ~ismember(config.network_source, valid_network_sources)
        error('Unknown config.network_source: %s', config.network_source);
    end
    if ~ismember(config.tc_source, valid_tc_sources)
        error('Unknown config.tc_source: %s', config.tc_source);
    end
    if ~ismember(config.cell_count_source, valid_count_sources)
        error('Unknown config.cell_count_source: %s', config.cell_count_source);
    end
    if config.reference.kernel_num < 1 || ...
            config.reference.kernel_num ~= round(config.reference.kernel_num)
        error('config.reference.kernel_num must be a positive integer.');
    end
    if config.thalamus_cells_per_area < 0 || ...
            config.thalamus_cells_per_area ~= round(config.thalamus_cells_per_area)
        error('config.thalamus_cells_per_area must be a non-negative integer.');
    end
    if config.simulation.n_trials < 1 || ...
            config.simulation.n_trials ~= round(config.simulation.n_trials)
        error('config.simulation.n_trials must be a positive integer.');
    end
    if config.simulation.n_time_bins < 1 || ...
            config.simulation.n_time_bins ~= round(config.simulation.n_time_bins)
        error('config.simulation.n_time_bins must be a positive integer.');
    end
end

%% =========================================================================
%  NETWORK CONSTRUCTION
%  =========================================================================
function [sim_meta, network, kernels] = build_simulation_network(root, config)
    progress_log('NETWORK', 'Building cell layout.');

    cell_info = build_cell_layout(root, config);
    kernels = load_connection_kernels(root, config, cell_info.N_total);

    K = config.reference.kernel_num;
    N = cell_info.N_total;

    network = struct();
    network.h = zeros(N, 1);
    network.J = zeros(N, N, K);
    network.provenance = struct();
    network.provenance.network_source = config.network_source;
    network.provenance.tc_source = config.tc_source;

    % Cortical baseline and cortical connections.
    switch config.network_source
        case 'load'
            progress_log('NETWORK', 'Loading cortical h and J from Post GLM.');
            post_meta = make_reference_meta(config, 'Post');
            loaded_post = load_reference_glm_network(root, post_meta, config);
            [network.h, network.J] = insert_loaded_cortex_network( ...
                network.h, network.J, loaded_post, cell_info, config);

        case 'generate'
            progress_log('NETWORK', 'Generating cortical h and J from random rules.');
            [network.h, network.J] = generate_cortical_network( ...
                network.h, network.J, cell_info, config);
    end

    % Thalamic baselines are generated in both modes because the loaded Post
    % cortical model normally has no artificial thalamic cells.
    network.h = generate_thalamic_baseline(network.h, cell_info, config);

    % Thalamocortical connections.
    switch config.tc_source
        case 'load_pre'
            progress_log('NETWORK', ...
                'Attempting to load thalamocortical J from Pre GLM.');
            pre_meta = make_reference_meta(config, 'Pre');
            loaded_pre = load_reference_glm_network(root, pre_meta, config);
            network.J = insert_loaded_tc_network( ...
                network.J, loaded_pre, cell_info, config);

        case 'generate'
            progress_log('NETWORK', 'Generating thalamocortical J.');
            network.J = generate_tc_network(network.J, cell_info, config);

        case 'none'
            progress_log('NETWORK', 'No thalamocortical connections added.');
    end

    % Remove self-connections for every kernel.
    for k = 1:K
        Jk = network.J(:, :, k);
        Jk(1:N+1:end) = 0;
        network.J(:, :, k) = Jk;
    end

    network.cell_area = cell_info.cell_area;
    network.cell_id = cell_info.cell_id;
    network.masks = cell_info.masks;
    network.N = N;
    network.K = K;

    sim_meta = struct();
    sim_meta.reference = config.reference;
    sim_meta.cell_area = cell_info.cell_area;
    sim_meta.cell_id = cell_info.cell_id;
    sim_meta.N_total = cell_info.N_total;
    sim_meta.N_cortex = cell_info.N_cortex;
    sim_meta.N_thalamus = cell_info.N_thalamus;
    sim_meta.kernel_name = config.reference.kernel_name;
    sim_meta.kernel_num = K;
    sim_meta.reg_name = config.reference.reg_name;
    sim_meta.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    progress_log('NETWORK', ...
        'Network ready: N=%d, cortex=%d, thalamus=%d, K=%d.', ...
        sim_meta.N_total, sim_meta.N_cortex, sim_meta.N_thalamus, K);
end

function cell_info = build_cell_layout(root, config)
    switch config.cell_count_source
        case 'data'
            post_meta = make_reference_meta(config, 'Post');
            raster_info = load_reference_raster_cell_info(root, post_meta);

            is_acc = ismember(string(raster_info.cell_area), ...
                string(config.area_labels.ACC));
            is_vlpfc = ismember(string(raster_info.cell_area), ...
                string(config.area_labels.VLPFC));

            if ~any(is_acc)
                error('No ACC cells found in the reference raster.');
            end
            if ~any(is_vlpfc)
                error('No VLPFC cells found in the reference raster.');
            end

            % Reorder cortical cells into ACC then VLPFC. Loaded GLM matrices
            % are reordered consistently in insert_loaded_cortex_network().
            source_cortex_indices = [find(is_acc); find(is_vlpfc)];
            n_acc = sum(is_acc);
            n_vlpfc = sum(is_vlpfc);
            source_cell_id = raster_info.cell_id(source_cortex_indices);

        case 'manual'
            n_acc = config.manual_cell_counts.ACC;
            n_vlpfc = config.manual_cell_counts.VLPFC;
            source_cortex_indices = (1:(n_acc + n_vlpfc)).';
            source_cell_id = compose('manual_cortex_%d', source_cortex_indices);
    end

    n_thal_acc = config.thalamus_cells_per_area;
    n_thal_vlpfc = config.thalamus_cells_per_area;

    n_cortex = n_acc + n_vlpfc;
    n_thalamus = n_thal_acc + n_thal_vlpfc;
    n_total = n_cortex + n_thalamus;

    idx_acc = (1:n_acc).';
    idx_vlpfc = (n_acc + (1:n_vlpfc)).';
    idx_thal_acc = (n_cortex + (1:n_thal_acc)).';
    idx_thal_vlpfc = (n_cortex + n_thal_acc + (1:n_thal_vlpfc)).';

    cell_area = strings(n_total, 1);
    cell_area(idx_acc) = "ACC";
    cell_area(idx_vlpfc) = "VLPFC";
    cell_area(idx_thal_acc) = "Thal-ACC";
    cell_area(idx_thal_vlpfc) = "Thal-VLPFC";

    cell_id = strings(n_total, 1);
    cell_id(1:n_cortex) = string(source_cell_id(:));
    cell_id(idx_thal_acc) = compose('sim_thal_acc_%d', (1:n_thal_acc).');
    cell_id(idx_thal_vlpfc) = compose('sim_thal_vlpfc_%d', (1:n_thal_vlpfc).');

    masks = struct();
    masks.ACC = false(n_total, 1);
    masks.VLPFC = false(n_total, 1);
    masks.thalamus_ACC = false(n_total, 1);
    masks.thalamus_VLPFC = false(n_total, 1);
    masks.ACC(idx_acc) = true;
    masks.VLPFC(idx_vlpfc) = true;
    masks.thalamus_ACC(idx_thal_acc) = true;
    masks.thalamus_VLPFC(idx_thal_vlpfc) = true;
    masks.cortex = masks.ACC | masks.VLPFC;
    masks.thalamus = masks.thalamus_ACC | masks.thalamus_VLPFC;

    cell_info = struct();
    cell_info.N_total = n_total;
    cell_info.N_cortex = n_cortex;
    cell_info.N_thalamus = n_thalamus;
    cell_info.n_acc = n_acc;
    cell_info.n_vlpfc = n_vlpfc;
    cell_info.n_thal_acc = n_thal_acc;
    cell_info.n_thal_vlpfc = n_thal_vlpfc;
    cell_info.cell_area = cell_area;
    cell_info.cell_id = cell_id;
    cell_info.masks = masks;
    cell_info.source_cortex_indices = source_cortex_indices(:);

    progress_log('CELLS', ...
        'ACC=%d, VLPFC=%d, Thal-ACC=%d, Thal-VLPFC=%d.', ...
        n_acc, n_vlpfc, n_thal_acc, n_thal_vlpfc);
end

function meta = make_reference_meta(config, prepost)
    meta = config.reference;
    meta.prepost = prepost;
end

%% =========================================================================
%  PROJECT ADAPTERS: REVIEW THESE FUNCTIONS
%  =========================================================================
function raster_info = load_reference_raster_cell_info(root, meta)
    % PROJECT ADAPTER 1:
    % Assumed current raster format:
    %   raster_file.meta.N
    %   raster_file.data.cell_area
    % Optional:
    %   raster_file.data.cell_id / cell_ids / unit_id / unit_ids
    %
    % If the current project stores these fields elsewhere, modify only this
    % function.

    meta.filename = generate_filename('raster', meta);
    raster_path = fullfile(root, 'Data', 'Working', 'raster', meta.filename);
    progress_log('LOAD-RASTER', 'Loading cell information: %s', raster_path);

    if ~isfile(raster_path)
        error('Reference raster file not found: %s', raster_path);
    end

    loaded = load(raster_path);
    if ~isfield(loaded, 'data') || ~isfield(loaded.data, 'cell_area')
        error(['Raster adapter expected data.cell_area but it was not found: %s. ' ...
               'Modify load_reference_raster_cell_info().'], raster_path);
    end

    cell_area = string(loaded.data.cell_area(:));
    n_cells = numel(cell_area);

    cell_id = strings(n_cells, 1);
    id_fields = {'cell_id', 'cell_ids', 'unit_id', 'unit_ids', ...
        'cell_name', 'cell_names'};
    found_id = false;
    for i = 1:numel(id_fields)
        field = id_fields{i};
        if isfield(loaded.data, field) && ...
                numel(loaded.data.(field)) == n_cells
            cell_id = string(loaded.data.(field)(:));
            found_id = true;
            break;
        end
    end
    if ~found_id
        cell_id = compose('raster_cell_%d', (1:n_cells).');
        progress_log('LOAD-RASTER', ...
            'No explicit cell IDs found; generated positional IDs.');
    end

    raster_info = struct();
    raster_info.cell_area = cell_area;
    raster_info.cell_id = cell_id;
    raster_info.N = n_cells;
end

function loaded_network = load_reference_glm_network(root, meta, config)
    % PROJECT ADAPTER 2:
    % Assumed current GLM format:
    %   glm_file.meta.N
    %   glm_file.data.model_par
    %
    % Assumed parameter layout:
    %   h = model_par(:, baseline_col)
    %   J(:,:,k) =
    %       model_par(:, connection_start_col + N*(k-1) :
    %                       connection_start_col + N*k - 1)
    %
    % This adapter does not use model_err because simulation needs the
    % inferred parameter values, not their significance classification.
    %
    % If the stored GLM contains post-spike-history parameters before the
    % connection blocks, update config.glm_layout.connection_start_col.

    meta.filename = generate_filename('GLM', meta);
    glm_path = fullfile(root, 'Data', 'Working', 'GLM', meta.filename);
    progress_log('LOAD-GLM', 'Loading inferred network: %s', glm_path);

    if ~isfile(glm_path)
        error('Reference GLM file not found: %s', glm_path);
    end

    loaded = load(glm_path);
    if ~isfield(loaded, 'meta') || ~isfield(loaded.meta, 'N') || ...
            ~isfield(loaded, 'data') || ...
            ~isfield(loaded.data, 'model_par')
        error(['GLM adapter expected meta.N and data.model_par: %s. ' ...
               'Modify load_reference_glm_network().'], glm_path);
    end

    N = double(loaded.meta.N);
    model_par = loaded.data.model_par;
    K = config.reference.kernel_num;

    if size(model_par, 1) ~= N
        error('model_par rows (%d) do not match meta.N (%d): %s', ...
            size(model_par, 1), N, glm_path);
    end

    baseline_col = config.glm_layout.baseline_col;
    connection_start_col = config.glm_layout.connection_start_col;
    required_last_col = connection_start_col + N * K - 1;
    if baseline_col > size(model_par, 2) || ...
            required_last_col > size(model_par, 2)
        error(['Configured GLM layout exceeds model_par columns. ' ...
               'baseline_col=%d, connection_start_col=%d, N=%d, K=%d, ' ...
               'available columns=%d. Review config.glm_layout.'], ...
            baseline_col, connection_start_col, N, K, size(model_par, 2));
    end

    h = double(model_par(:, baseline_col));
    h(~isfinite(h)) = -20;

    J = zeros(N, N, K);
    for k = 1:K
        col_start = connection_start_col + N * (k - 1);
        col_end = col_start + N - 1;
        Jk = double(model_par(:, col_start:col_end));
        Jk(~isfinite(Jk)) = 0;
        J(:, :, k) = Jk;
    end

    % Load corresponding raster labels so the GLM can be reordered.
    raster_info = load_reference_raster_cell_info(root, meta);
    if raster_info.N ~= N
        error(['Reference GLM N=%d but corresponding raster has %d cells. ' ...
               'The current adapter assumes one raster label per GLM cell. ' ...
               'Modify load_reference_glm_network() if the GLM is a subset.'], ...
            N, raster_info.N);
    end

    loaded_network = struct();
    loaded_network.N = N;
    loaded_network.h = h;
    loaded_network.J = J;
    loaded_network.cell_area = raster_info.cell_area;
    loaded_network.cell_id = raster_info.cell_id;
    loaded_network.source_path = glm_path;
end

function [h_out, J_out] = insert_loaded_cortex_network( ...
    h_out, J_out, loaded_post, cell_info, config)

    is_acc = ismember(string(loaded_post.cell_area), ...
        string(config.area_labels.ACC));
    is_vlpfc = ismember(string(loaded_post.cell_area), ...
        string(config.area_labels.VLPFC));
    source_idx = [find(is_acc); find(is_vlpfc)];

    if numel(source_idx) ~= cell_info.N_cortex
        error(['Loaded cortical cell count (%d) does not match simulation ' ...
               'cortical count (%d). Check raster area labels and GLM scope.'], ...
            numel(source_idx), cell_info.N_cortex);
    end

    target_idx = find(cell_info.masks.cortex);
    h_out(target_idx) = loaded_post.h(source_idx);

    for k = 1:size(J_out, 3)
        J_out(target_idx, target_idx, k) = ...
            loaded_post.J(source_idx, source_idx, k);
    end

    progress_log('NETWORK', ...
        'Inserted loaded cortical network from %s.', ...
        loaded_post.source_path);
end

function J_out = insert_loaded_tc_network( ...
    J_out, loaded_pre, cell_info, config)
    % PROJECT ADAPTER 3:
    % This implementation works only if the Pre GLM itself contains
    % identifiable thalamic cells and uses the same target-row/source-column
    % convention as the cortical GLM.
    %
    % If the current Pre model stores thalamocortical parameters in a
    % separate object or predictor block, replace this function.

    areas = string(loaded_pre.cell_area);
    src_thal_acc = find(ismember(areas, ...
        string(config.area_labels.thalamus_ACC)));
    src_thal_vlpfc = find(ismember(areas, ...
        string(config.area_labels.thalamus_VLPFC)));
    src_acc = find(ismember(areas, string(config.area_labels.ACC)));
    src_vlpfc = find(ismember(areas, string(config.area_labels.VLPFC)));

    if isempty(src_thal_acc) || isempty(src_thal_vlpfc)
        error(['tc_source=''load_pre'' was requested, but the Pre GLM/raster ' ...
               'does not expose identifiable Thal-ACC and Thal-VLPFC cells. ' ...
               'Modify insert_loaded_tc_network() for the actual storage ' ...
               'format, or use tc_source=''generate''.']);
    end

    dst_acc = find(cell_info.masks.ACC);
    dst_vlpfc = find(cell_info.masks.VLPFC);
    dst_thal_acc = find(cell_info.masks.thalamus_ACC);
    dst_thal_vlpfc = find(cell_info.masks.thalamus_VLPFC);

    if numel(src_acc) ~= numel(dst_acc) || ...
            numel(src_vlpfc) ~= numel(dst_vlpfc)
        error('Pre GLM cortex counts do not match the simulation cortex counts.');
    end
    if numel(src_thal_acc) ~= numel(dst_thal_acc) || ...
            numel(src_thal_vlpfc) ~= numel(dst_thal_vlpfc)
        error(['Pre GLM thalamic counts do not match configured artificial ' ...
               'thalamic counts. Either change thalamus_cells_per_area or ' ...
               'implement a resampling rule in insert_loaded_tc_network().']);
    end

    for k = 1:size(J_out, 3)
        J_out(dst_acc, dst_thal_acc, k) = ...
            loaded_pre.J(src_acc, src_thal_acc, k);
        J_out(dst_vlpfc, dst_thal_vlpfc, k) = ...
            loaded_pre.J(src_vlpfc, src_thal_vlpfc, k);
    end

    progress_log('NETWORK', ...
        'Inserted loaded Pre thalamocortical connections from %s.', ...
        loaded_pre.source_path);
end

function kernels = load_connection_kernels(root, config, N)
    % PROJECT ADAPTER 4:
    % Assumed kernel file:
    %   Data/Working/kernel/kernel_<kernel_name>.mat
    % with:
    %   conn_kernels - cell array, one temporal kernel per connection block
    %
    % Each conn_kernels{k} may be:
    %   1 x L, L x 1, or N x L.
    %
    % If the current kernel file uses a different field or orientation,
    % modify this function.

    K = config.reference.kernel_num;
    switch config.kernel_basis.source
        case 'file'
            kernel_path = fullfile(root, 'Data', 'Working', 'kernel', ...
                sprintf('kernel_%s.mat', config.reference.kernel_name));
            progress_log('LOAD-KERNEL', 'Loading: %s', kernel_path);
            if ~isfile(kernel_path)
                error('Kernel file not found: %s', kernel_path);
            end
            loaded = load(kernel_path);
            if ~isfield(loaded, 'conn_kernels')
                error(['Kernel adapter expected conn_kernels: %s. ' ...
                       'Modify load_connection_kernels().'], kernel_path);
            end
            raw_kernels = loaded.conn_kernels;

        case 'manual'
            raw_kernels = config.kernel_basis.manual;

        otherwise
            error('Unknown kernel_basis.source: %s', ...
                config.kernel_basis.source);
    end

    if numel(raw_kernels) < K
        error('Only %d temporal kernels available, but kernel_num=%d.', ...
            numel(raw_kernels), K);
    end

    kernels = struct();
    kernels.values = cell(1, K);
    kernels.lengths = zeros(1, K);
    for k = 1:K
        kernel_value = double(raw_kernels{k});
        if isvector(kernel_value)
            kernel_value = reshape(kernel_value, 1, []);
        elseif size(kernel_value, 1) ~= N
            error(['Kernel %d must be a vector or an N x L matrix. ' ...
                   'Received size %s for N=%d.'], ...
                k, mat2str(size(kernel_value)), N);
        end
        if any(~isfinite(kernel_value(:)))
            error('Kernel %d contains nonfinite values.', k);
        end
        kernels.values{k} = kernel_value(:, end:-1:1);
        kernels.lengths(k) = size(kernel_value, 2);
    end
    kernels.max_length = max(kernels.lengths);
end

%% =========================================================================
%  GENERATED NETWORK RULES
%  =========================================================================
function [h, J] = generate_cortical_network(h, J, cell_info, config)
    rng(config.simulation.base_seed, 'twister');

    acc_idx = find(cell_info.masks.ACC);
    vlpfc_idx = find(cell_info.masks.VLPFC);
    area_indices = {acc_idx, vlpfc_idx};

    baseline_mean = config.generate.cortex.baseline_mean;
    baseline_sd = config.generate.cortex.baseline_sd;
    h(acc_idx) = baseline_mean(1) + baseline_sd(1) * randn(numel(acc_idx), 1);
    h(vlpfc_idx) = baseline_mean(2) + baseline_sd(2) * randn(numel(vlpfc_idx), 1);

    P = config.generate.cortex.connection_probability;
    Ppos = config.generate.cortex.positive_probability;
    Wmean = config.generate.cortex.weight_abs_mean;
    weight_cv = config.generate.cortex.weight_cv;
    kernel_scale = expand_kernel_scale( ...
        config.generate.cortex.kernel_scale, size(J, 3), ...
        'generate.cortex.kernel_scale');

    for target_area = 1:2
        target_idx = area_indices{target_area};
        for source_area = 1:2
            source_idx = area_indices{source_area};

            connection_mask = rand(numel(target_idx), numel(source_idx)) < ...
                P(target_area, source_area);
            positive_mask = rand(numel(target_idx), numel(source_idx)) < ...
                Ppos(target_area, source_area);
            signs = ones(size(connection_mask));
            signs(~positive_mask) = -1;

            base_abs = Wmean(target_area, source_area) .* ...
                max(1 + weight_cv * randn(size(connection_mask)), 0.05);
            base_weights = connection_mask .* signs .* base_abs;

            for k = 1:size(J, 3)
                J(target_idx, source_idx, k) = ...
                    kernel_scale(k) .* base_weights;
            end
        end
    end

    progress_log('GENERATE', 'Generated cortical baselines and connections.');
end

function h = generate_thalamic_baseline(h, cell_info, config)
    rng(config.simulation.base_seed + 1, 'twister');
    thal_idx = find(cell_info.masks.thalamus);
    mean_h = config.generate.thalamus.baseline_mean;
    sd_h = config.generate.thalamus.baseline_sd;
    h(thal_idx) = mean_h + sd_h * randn(numel(thal_idx), 1);
end

function J = generate_tc_network(J, cell_info, config)
    rng(config.simulation.base_seed + 2, 'twister');

    kernel_scale = expand_kernel_scale( ...
        config.generate.tc.kernel_scale, size(J, 3), ...
        'generate.tc.kernel_scale');

    J = generate_one_tc_projection( ...
        J, find(cell_info.masks.ACC), ...
        find(cell_info.masks.thalamus_ACC), ...
        config, kernel_scale);
    J = generate_one_tc_projection( ...
        J, find(cell_info.masks.VLPFC), ...
        find(cell_info.masks.thalamus_VLPFC), ...
        config, kernel_scale);
end

function J = generate_one_tc_projection( ...
    J, target_idx, source_idx, config, kernel_scale)

    connection_mask = rand(numel(target_idx), numel(source_idx)) < ...
        config.generate.tc.connection_probability;
    base_weights = config.generate.tc.weight_mean .* ...
        max(1 + config.generate.tc.weight_cv .* ...
        randn(size(connection_mask)), 0.05);
    base_weights = connection_mask .* base_weights;

    for k = 1:size(J, 3)
        J(target_idx, source_idx, k) = ...
            kernel_scale(k) .* base_weights;
    end
end

function scale = expand_kernel_scale(value, K, field_name)
    value = double(value(:).');
    if isscalar(value)
        scale = repmat(value, 1, K);
    elseif numel(value) >= K
        scale = value(1:K);
    else
        error('%s has %d values but K=%d.', field_name, numel(value), K);
    end
end

%% =========================================================================
%  SIMULATION
%  =========================================================================
function results = simulate_all_conditions(config, sim_meta, network, kernels)
    n_conditions = numel(config.conditions);
    results = struct([]);

    for condition_i = 1:n_conditions
        condition = config.conditions(condition_i);
        progress_log('SIMULATION', ...
            '[%d/%d] START condition=%s.', ...
            condition_i, n_conditions, condition.name);

        condition_network = network;
        if ~condition.enable_tc
            thalamus_sources = condition_network.masks.thalamus;
            condition_network.J(:, thalamus_sources, :) = 0;
        end

        condition_result = simulate_condition( ...
            config, sim_meta, condition_network, kernels, ...
            condition, condition_i);
        results(condition_i) = condition_result;

        progress_log('SIMULATION', ...
            '[%d/%d] DONE condition=%s.', ...
            condition_i, n_conditions, condition.name);
    end
end

function result = simulate_condition( ...
    config, sim_meta, network, kernels, condition, condition_i)

    n_trials = config.simulation.n_trials;
    n_time_bins = config.simulation.n_time_bins;
    N = network.N;

    trial_rasters = cell(1, n_trials);
    firing_rates = zeros(N, n_trials);

    if config.simulation.parallel_trials
        progress_log('SIMULATION', ...
            'Using parfor for %d trials.', n_trials);
        parfor trial_i = 1:n_trials
            trial_seed = config.simulation.base_seed + ...
                100000 * condition_i + trial_i;
            raster = simulate_one_trial( ...
                network, kernels, condition, n_time_bins, trial_seed);
            trial_rasters{trial_i} = raster;
            firing_rates(:, trial_i) = mean(raster, 2);
        end
    else
        for trial_i = 1:n_trials
            trial_tic = tic;
            trial_seed = config.simulation.base_seed + ...
                100000 * condition_i + trial_i;
            trial_rasters{trial_i} = simulate_one_trial( ...
                network, kernels, condition, n_time_bins, trial_seed);
            firing_rates(:, trial_i) = ...
                mean(trial_rasters{trial_i}, 2);
            progress_log('TRIAL', ...
                '%s trial %d/%d done in %.1f s.', ...
                condition.name, trial_i, n_trials, toc(trial_tic));
        end
    end

    saved_cell_mask = true(N, 1);
    if config.simulation.save_cortex_only
        saved_cell_mask = network.masks.cortex;
        for trial_i = 1:n_trials
            trial_rasters{trial_i} = ...
                trial_rasters{trial_i}(saved_cell_mask, :);
        end
        firing_rates = firing_rates(saved_cell_mask, :);
    end

    result = struct();
    result.condition = condition;
    result.N = sum(saved_cell_mask);
    result.cell_area = sim_meta.cell_area(saved_cell_mask);
    result.cell_id = sim_meta.cell_id(saved_cell_mask);
    result.trial_num = n_trials;
    result.trial_len = repmat(n_time_bins, 1, n_trials);
    result.firing_rates = firing_rates;

    if config.output.save_rasters
        result.rasters = trial_rasters;
    else
        result.rasters = {};
    end
end

function raster = simulate_one_trial( ...
    network, kernels, condition, n_time_bins, seed)

    rng(seed, 'twister');

    N = network.N;
    K = network.K;
    raster = false(N, n_time_bins);

    for t = 1:n_time_bins
        h_total = network.h;

        if condition.sync_thalamus
            phase = mod(t - 1, condition.sync_period);
            if phase < condition.sync_period * condition.sync_on_fraction
                h_total(network.masks.thalamus) = ...
                    h_total(network.masks.thalamus) + ...
                    condition.sync_on_offset;
            else
                h_total(network.masks.thalamus) = ...
                    h_total(network.masks.thalamus) + ...
                    condition.sync_off_offset;
            end
        end

        for k = 1:K
            L = kernels.lengths(k);
            if t > L
                history = double(raster(:, (t - L):(t - 1)));
                kernel_value = kernels.values{k};
                convolved = sum(history .* kernel_value, 2);
                h_total = h_total + network.J(:, :, k) * convolved;
            end
        end

        % Numerically stable logistic probability.
        h_total = max(min(h_total, 30), -30);
        spike_probability = 1 ./ (1 + exp(-h_total));
        raster(:, t) = rand(N, 1) < spike_probability;
    end
end

%% =========================================================================
%  OUTPUT
%  =========================================================================
function save_simulation_output(root, config, sim_meta, network, results)
    output_folder = fullfile(root, 'Data', 'Working', ...
        'simulation_v2', sanitize_for_path(config.output.simulation_tag));
    check_path(output_folder);

    output_path = fullfile(output_folder, ...
        sprintf('simulation_%s.mat', ...
        sanitize_for_path(config.output.simulation_tag)));

    if isfile(output_path) && ~config.output.overwrite
        error(['Simulation output already exists: %s. ' ...
               'Set config.output.overwrite=true to replace it.'], ...
            output_path);
    end

    if ~config.output.save_network
        network_to_save = rmfield(network, {'h', 'J'});
    else
        network_to_save = network;
    end

    network = network_to_save; %#ok<NASGU>
    progress_log('SAVE', 'Writing: %s', output_path);
    save(output_path, 'config', 'sim_meta', 'network', 'results', '-v7.3');

    % PROJECT ADAPTER 5:
    % To feed simulated rasters directly into the current GLM training
    % pipeline, add a conversion/registration block here. The exact required
    % current schema for raster metadata, data fields, and metadata indexing
    % must be supplied before this can be implemented safely.
end

%% =========================================================================
%  UTILITIES
%  =========================================================================
function progress_log(stage, format_string, varargin)
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    message = sprintf(format_string, varargin{:});
    fprintf('[%s][%s] %s\n', timestamp, stage, message);
end

function out = sanitize_for_path(value)
    out = char(string(value));
    out = regexprep(out, '[^A-Za-z0-9._=-]+', '_');
    out = regexprep(out, '^_+|_+$', '');
    if isempty(out)
        out = 'simulation';
    end
end
