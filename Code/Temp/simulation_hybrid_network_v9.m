%% simulation_hybrid_network_v9.m
% Hybrid GLM network simulation using the project meta/data file format.
%
% Model-size source:
%   config.model_size_source = 'match_data'
%       Match cell counts and dt to project raster data. Temporal kernels
%       are loaded from the externally selected project kernel file.
%
%   config.model_size_source = 'fixed'
%       Use fixed cell counts and fixed dt. No reference raster or GLM file
%       is loaded. Temporal kernels are still loaded from the externally
%       selected project kernel file. Cortical and TC connections are
%       generated from the configured random rules.
%
% Cortical network source:
%   config.network_source = 'load'
%       Load cortical h, post-spike weights P, and cortical J from a fitted
%       Post GLM model.
%
%   config.network_source = 'generate'
%       Copy cortical cell counts/labels only, then generate h, P, and J.
%
% Thalamocortical source:
%   config.tc_source = 'load_pre'
%       Load thalamocortical J from a fitted Pre GLM model.
%
%   config.tc_source = 'generate'
%       Generate thalamocortical J.
%
%   config.tc_source = 'none'
%       Do not add thalamocortical J.
%
% Project-compatible output:
%   For every simulated condition, two raster/border pairs are written:
%       area='Full'   - ACC, VLPFC, and Thalamus
%       area='Cortex' - ACC and VLPFC only
%
%   Data/Working/raster/raster_##.mat
%       Contains structs named meta and data.
%
%   Data/Working/border/border_##.mat
%       Contains structs named meta and data.
%
%   Data/Working/simulation_v2/<simulation_tag>/simulation_manifest_*.mat
%       Contains meta and data with configuration, network, kernel
%       definition, provenance, and paths of generated raster/border files.
%
% The input/output adapters follow README.md. Simulated thalamic neurons
% use cell_area='Thalamus'. The optional data.cell_type field records
% generated projection classes ('Thal-ACC' or 'Thal-VLPFC') when needed.

clear;
clc;

run_tic = tic;
progress_log('SCRIPT', 'Started.');

%% Project root
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for depth_i = 1:code_depth
    root = fileparts(root);
end
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));
progress_log('SCRIPT', 'Project root: %s', root);

%% User configuration
% This is the only section intended for routine parameter editing.
% Each user-facing config field is assigned exactly once.
config = struct();

% -------------------------------------------------------------------------
% MODEL SIZE / DATA SOURCE
% -------------------------------------------------------------------------
% 'match_data' = obtain cell counts/dt from reference raster data.
% 'fixed'      = use fixed cell counts/dt without reference raster/GLM data.
config.model_size_source = 'fixed';

% 'load' or 'generate'. Fixed mode requires 'generate'.
config.network_source = 'generate';

% 'load_pre', 'generate', or 'none'. Fixed mode cannot use 'load_pre'.
config.tc_source = 'generate';

% 'data' or 'manual'. Fixed mode requires 'manual'.
config.cell_count_source = 'manual';

% -------------------------------------------------------------------------
% REFERENCE SESSION
% Used only by match_data mode or when loading fitted GLM parameters.
% -------------------------------------------------------------------------
config.reference.animal_name = 'Slayer';
config.reference.injection = 'Muscimol';
config.reference.session_idx = 1;
config.reference.state = 'RestOpen';
config.reference.align = 'Last';
config.reference.resting_dur_threshold = 15;
config.reference.cortex_area = 'Cortex';
config.reference.tc_area = 'Full';
config.reference.reg_name = 'L2=0_2';
config.reference.epoch = 3000;
config.reference.fold_idx = 0;
config.reference.shuffle_idx = 0;

% -------------------------------------------------------------------------
% EXTERNAL TEMPORAL KERNEL
% Name, expected number of connection kernels, and generated J scales.
% -------------------------------------------------------------------------
config.kernel.name = 'LongExp320';
config.kernel.n_conn_kernel = 1;
config.kernel.cortex_scale = 1.0;
config.kernel.tc_scale = 1.0;

% -------------------------------------------------------------------------
% FIXED MODEL SIZE
% Used only when model_size_source='fixed'.
% -------------------------------------------------------------------------
config.fixed_model.cell_counts.ACC = 40;
config.fixed_model.cell_counts.VLPFC = 40;
config.fixed_model.cell_counts.thalamus_ACC = 10;
config.fixed_model.cell_counts.thalamus_VLPFC = 10;
config.fixed_model.dt = 0.001;
config.fixed_model.reg_name = 'Generated';

% -------------------------------------------------------------------------
% MANUAL COUNTS FOR match_data + cell_count_source='manual'
% -------------------------------------------------------------------------
config.manual_cell_counts.ACC = 40;
config.manual_cell_counts.VLPFC = 40;
config.thalamus_cells_per_area = 10;

% -------------------------------------------------------------------------
% AREA LABELS
% Canonical output cell_area values are ACC, VLPFC, and Thalamus.
% -------------------------------------------------------------------------
config.area_labels.ACC = {'ACC'};
config.area_labels.VLPFC = {'VLPFC'};
config.area_labels.thalamus_ACC = {'Thal-ACC', 'Thalamus-ACC'};
config.area_labels.thalamus_VLPFC = {'Thal-VLPFC', 'Thalamus-VLPFC'};
config.area_labels.thalamus_generic = {'Thalamus', 'Thal'};

% -------------------------------------------------------------------------
% GENERATED CORTICAL NETWORK PARAMETERS
% Matrix convention: target area in rows, source area in columns.
% -------------------------------------------------------------------------
config.generate.cortex.connection_probability = [0.08, 0.04; 0.04, 0.08];
config.generate.cortex.positive_probability = [0.75, 0.65; 0.65, 0.75];
config.generate.cortex.weight_abs_mean = [1.5, 0.8; 0.8, 1.5];
config.generate.cortex.weight_cv = 0.20;
config.generate.cortex.baseline_mean = [-4.0, -4.0];
config.generate.cortex.baseline_sd = [0.4, 0.4];
config.generate.cortex.post_spike_weight_mean = -2.0;
config.generate.cortex.post_spike_weight_sd = 0.2;

% -------------------------------------------------------------------------
% GENERATED THALAMIC INTRINSIC PARAMETERS
% -------------------------------------------------------------------------
config.generate.thalamus.baseline_mean = -4.5;
config.generate.thalamus.baseline_sd = 0.3;
config.generate.thalamus.post_spike_weight_mean = -2.0;
config.generate.thalamus.post_spike_weight_sd = 0.2;

% -------------------------------------------------------------------------
% GENERATED THALAMOCORTICAL CONNECTION PARAMETERS
% -------------------------------------------------------------------------
config.generate.tc.connection_probability = 0.10;
config.generate.tc.weight_mean = 3.0;
config.generate.tc.weight_cv = 0.20;

% -------------------------------------------------------------------------
% CONDITIONS
% All conditions in one simulated session share one session_idx and must
% therefore have unique output_state values.
% -------------------------------------------------------------------------
condition_template = struct( ...
    'name', '', ...
    'output_state', '', ...
    'enable_tc', true, ...
    'sync_thalamus', false, ...
    'sync_period', 100, ...
    'sync_on_fraction', 1/8, ...
    'sync_on_offset', 1.5, ...
    'sync_off_offset', -0.6);
config.conditions = repmat(condition_template, 1, 3);

config.conditions(1).name = 'FittedAsync';
config.conditions(1).output_state = 'FittedAsync';

config.conditions(2).name = 'FittedSync';
config.conditions(2).output_state = 'FittedSync';
config.conditions(2).sync_thalamus = true;

config.conditions(3).name = 'FittedNoInput';
config.conditions(3).output_state = 'FittedNoInput';
config.conditions(3).enable_tc = false;

% -------------------------------------------------------------------------
% SIMULATION SIZE AND RANDOMNESS
% These values are not redefined anywhere else in the script.
% -------------------------------------------------------------------------
config.simulation.n_sessions = 5;
config.simulation.n_trials = 10;
config.simulation.n_time_bins = 30000;
config.simulation.dt = [];          % [] = use input dt or default_dt.
config.simulation.default_dt = 0.001;
config.simulation.parallel_trials = false;
config.simulation.base_seed = 137;
config.simulation.session_seed_stride = 1000000;

% -------------------------------------------------------------------------
% OUTPUT
% Nonempty values are preserved exactly. Empty identity fields are derived
% once by finalize_user_config() and are not overwritten later.
% -------------------------------------------------------------------------
config.output.simulation_tag = '';  % '' = derive automatically.
config.output.animal_name = '';     % '' = SimFixed or Sim<reference animal>.
config.output.injection = 'No injection';
config.output.prepost = 'Post';
config.output.align = 'None';
config.output.base_session_idx = [];          % [] = derive by mode.
config.output.resting_dur_threshold = [];     % [] = derive by mode.
config.output.overwrite = true;
config.output.save_project_raster = true;
config.output.save_project_border = true;
config.output.save_manifest = true;
config.output.save_network_in_manifest = true;

% Resolve only explicitly empty derived output fields. No simulation-size,
% kernel, network, or cell-count parameter is modified here.
config = finalize_user_config(config);

%% Validate and run
validate_simulation_config(config);
progress_log('CONFIG', ...
    ['model_size_source=%s, network_source=%s, tc_source=%s, ' ...
     'kernel=%s, reg=%s, animal=%s, first_session=%d, ' ...
     'n_sessions=%d, n_trials=%d, n_time_bins=%d.'], ...
    config.model_size_source, config.network_source, config.tc_source, ...
    config.kernel.name, ...
    resolve_effective_reg_name(config), ...
    config.output.animal_name, ...
    config.output.base_session_idx, ...
    config.simulation.n_sessions, ...
    config.simulation.n_trials, ...
    config.simulation.n_time_bins);

output_manifest = run_simulation_sessions(root, config); %#ok<NASGU>

progress_log('SUMMARY', 'Finished in %.1f s.', toc(run_tic));

%% =========================================================================
%  CONFIGURATION FINALIZATION
%  =========================================================================
function config = finalize_user_config(config)
    % Incompatible combinations raise errors instead of being silently
    % overwritten. The top configuration remains the source of truth.
    switch config.model_size_source
        case 'fixed'
            if ~strcmp(config.network_source, 'generate')
                error(['fixed mode requires network_source=''generate''; ' ...
                       'the script no longer overrides this setting.']);
            end
            if ~strcmp(config.cell_count_source, 'manual')
                error(['fixed mode requires cell_count_source=''manual''; ' ...
                       'the script no longer overrides this setting.']);
            end
            if strcmp(config.tc_source, 'load_pre')
                error(['fixed mode cannot use tc_source=''load_pre''; ' ...
                       'choose ''generate'' or ''none''.']);
            end

        case 'match_data'
            % No mode-derived changes are required.

        otherwise
            error('Unknown config.model_size_source: %s', ...
                config.model_size_source);
    end

    % Fill only empty output identity fields. Explicit user values are kept.
    if strlength(string(config.output.animal_name)) == 0
        if strcmp(config.model_size_source, 'fixed')
            config.output.animal_name = 'SimFixed';
        else
            config.output.animal_name = sprintf( ...
                'Sim%s', config.reference.animal_name);
        end
    end

    if isempty(config.output.base_session_idx)
        if strcmp(config.model_size_source, 'fixed')
            config.output.base_session_idx = 1;
        else
            config.output.base_session_idx = ...
                config.reference.session_idx;
        end
    end

    if isempty(config.output.resting_dur_threshold)
        if strcmp(config.model_size_source, 'fixed')
            config.output.resting_dur_threshold = 0;
        else
            config.output.resting_dur_threshold = ...
                config.reference.resting_dur_threshold;
        end
    end

    if strlength(string(config.output.simulation_tag)) == 0
        if strcmp(config.model_size_source, 'fixed')
            config.output.simulation_tag = sprintf( ...
                ['fixed_ACC%d_VLPFC%d_ThalACC%d_ThalVLPFC%d_' ...
                 'kernel_%s_%s_nSess%d'], ...
                config.fixed_model.cell_counts.ACC, ...
                config.fixed_model.cell_counts.VLPFC, ...
                config.fixed_model.cell_counts.thalamus_ACC, ...
                config.fixed_model.cell_counts.thalamus_VLPFC, ...
                config.kernel.name, ...
                config.tc_source, ...
                config.simulation.n_sessions);
        else
            config.output.simulation_tag = sprintf( ...
                '%s_s%d_%s_%s_%s_nSess%d', ...
                config.reference.animal_name, ...
                config.reference.session_idx, ...
                config.reference.state, ...
                config.network_source, ...
                config.tc_source, ...
                config.simulation.n_sessions);
        end
    end
end

function validate_simulation_config(config)
    if ~ismember(config.model_size_source, {'match_data', 'fixed'})
        error('Unknown config.model_size_source: %s', ...
            config.model_size_source);
    end
    if ~ismember(config.network_source, {'load', 'generate'})
        error('Unknown config.network_source: %s', config.network_source);
    end
    if ~ismember(config.tc_source, {'load_pre', 'generate', 'none'})
        error('Unknown config.tc_source: %s', config.tc_source);
    end
    if ~ismember(config.cell_count_source, {'data', 'manual'})
        error('Unknown config.cell_count_source: %s', config.cell_count_source);
    end
    if strlength(string(config.kernel.name)) == 0
        error(['config.kernel.name must identify an external kernel file ' ...
               'in every model-size mode.']);
    end
    if ~isscalar(config.kernel.n_conn_kernel) || ...
            config.kernel.n_conn_kernel < 1 || ...
            config.kernel.n_conn_kernel ~= round(config.kernel.n_conn_kernel)
        error('config.kernel.n_conn_kernel must be a positive integer.');
    end
    validate_manual_kernel_scale( ...
        config.kernel.cortex_scale, config.kernel.n_conn_kernel, ...
        'config.kernel.cortex_scale');
    validate_manual_kernel_scale( ...
        config.kernel.tc_scale, config.kernel.n_conn_kernel, ...
        'config.kernel.tc_scale');
    if strcmp(config.model_size_source, 'fixed')
        fixed_counts = [ ...
            config.fixed_model.cell_counts.ACC, ...
            config.fixed_model.cell_counts.VLPFC, ...
            config.fixed_model.cell_counts.thalamus_ACC, ...
            config.fixed_model.cell_counts.thalamus_VLPFC];
        if any(fixed_counts < 0) || any(fixed_counts ~= round(fixed_counts)) || ...
                fixed_counts(1) < 1 || fixed_counts(2) < 1
            error(['Fixed model cell counts must be non-negative integers, ' ...
                   'with at least one ACC and one VLPFC cell.']);
        end
        if ~isscalar(config.fixed_model.dt) || ...
                ~isfinite(config.fixed_model.dt) || config.fixed_model.dt <= 0
            error('config.fixed_model.dt must be a positive scalar in seconds.');
        end
        if ~strcmp(config.network_source, 'generate') || ...
                strcmp(config.tc_source, 'load_pre')
            error(['Fixed mode must use generated cortical connections and ' ...
                   'cannot use tc_source=''load_pre''.']);
        end
    end
    if config.thalamus_cells_per_area < 0 || ...
            config.thalamus_cells_per_area ~= round(config.thalamus_cells_per_area)
        error('config.thalamus_cells_per_area must be a non-negative integer.');
    end
    if config.simulation.n_sessions < 1 || ...
            config.simulation.n_sessions ~= round(config.simulation.n_sessions)
        error('config.simulation.n_sessions must be a positive integer.');
    end
    if config.simulation.n_trials < 1 || ...
            config.simulation.n_trials ~= round(config.simulation.n_trials)
        error('config.simulation.n_trials must be a positive integer.');
    end
    if config.simulation.n_time_bins < 1 || ...
            config.simulation.n_time_bins ~= round(config.simulation.n_time_bins)
        error('config.simulation.n_time_bins must be a positive integer.');
    end
    if ~isempty(config.simulation.dt) && ...
            (~isscalar(config.simulation.dt) || config.simulation.dt <= 0)
        error('config.simulation.dt must be empty or a positive scalar in seconds.');
    end
    if ~isscalar(config.simulation.default_dt) || ...
            ~isfinite(config.simulation.default_dt) || ...
            config.simulation.default_dt <= 0
        error('config.simulation.default_dt must be a positive scalar in seconds.');
    end
    if ~isscalar(config.simulation.session_seed_stride) || ...
            ~isfinite(config.simulation.session_seed_stride) || ...
            config.simulation.session_seed_stride < 1 || ...
            config.simulation.session_seed_stride ~= ...
                round(config.simulation.session_seed_stride)
        error(['config.simulation.session_seed_stride must be a positive ' ...
               'integer.']);
    end
    if config.output.base_session_idx < 0 || ...
            config.output.base_session_idx ~= round(config.output.base_session_idx)
        error('config.output.base_session_idx must be a non-negative integer.');
    end
    if isempty(config.conditions)
        error('config.conditions must contain at least one condition.');
    end
    output_states = string({config.conditions.output_state});
    if any(strlength(output_states) == 0) || ...
            numel(unique(output_states)) ~= numel(output_states)
        error(['Every condition must have a nonempty, unique output_state ' ...
               'because conditions share one session_idx.']);
    end
end

%% =========================================================================
%  MULTI-SESSION DRIVER
%  =========================================================================
function output_manifest = run_simulation_sessions(root, config)
    n_sessions = config.simulation.n_sessions;
    n_conditions = numel(config.conditions);
    output_areas = {'Full', 'Cortex'};
    n_areas = numel(output_areas);

    output_manifest = struct();
    output_manifest.session_indices = ...
        config.output.base_session_idx + (0:(n_sessions - 1));
    output_manifest.area_labels = string(output_areas);
    output_manifest.condition_states = ...
        string({config.conditions.output_state});
    output_manifest.raster_paths = ...
        strings(n_sessions, n_conditions, n_areas);
    output_manifest.border_paths = strings(n_sessions, n_areas);
    output_manifest.session_seeds = zeros(1, n_sessions);
    output_manifest.manifest_path = "";

    session_records = cell(1, n_sessions);
    run_kernels = [];

    for session_i = 1:n_sessions
        output_session_idx = ...
            config.output.base_session_idx + session_i - 1;
        session_seed = config.simulation.base_seed + ...
            (session_i - 1) * config.simulation.session_seed_stride;

        session_config = config;
        session_config.simulation.base_seed = session_seed;

        progress_log('SESSION', ...
            '[%d/%d] START session_idx=%d, seed=%d.', ...
            session_i, n_sessions, output_session_idx, session_seed);

        [sim_meta, network, kernels, source_info] = ...
            build_simulation_network(root, session_config);
        if session_i == 1
            run_kernels = kernels;
        end
        sim_meta.simulation_session_number = session_i;
        sim_meta.output_session_idx = output_session_idx;
        sim_meta.session_seed = session_seed;

        results = simulate_all_conditions( ...
            session_config, sim_meta, network, kernels);
        session_output = save_simulation_session_output( ...
            root, session_config, sim_meta, results, ...
            output_session_idx);

        output_manifest.raster_paths(session_i, :, :) = ...
            reshape(session_output.raster_paths, ...
            1, n_conditions, n_areas);
        output_manifest.border_paths(session_i, :) = ...
            session_output.border_paths;
        output_manifest.session_seeds(session_i) = session_seed;

        session_record = struct();
        session_record.session_number = session_i;
        session_record.session_idx = output_session_idx;
        session_record.session_seed = session_seed;
        session_record.sim_meta = sim_meta;
        session_record.source_info = source_info;
        if config.output.save_network_in_manifest
            session_record.network = network;
        else
            session_record.network = rmfield(network, {'h', 'P', 'J'});
        end
        session_records{session_i} = session_record;

        clear results network
        progress_log('SESSION', ...
            '[%d/%d] DONE session_idx=%d.', ...
            session_i, n_sessions, output_session_idx);
    end

    if config.output.save_manifest
        output_manifest.manifest_path = save_simulation_run_manifest( ...
            root, config, run_kernels, session_records, output_manifest);
    end
end

%% =========================================================================
%  NETWORK CONSTRUCTION
%  =========================================================================
function [sim_meta, network, kernels, source_info] = ...
    build_simulation_network(root, config)

    is_fixed = strcmp(config.model_size_source, 'fixed');
    kernels = load_kernel_definition(root, config);
    K_conn = kernels.n_conn_kernel;
    K_ps = kernels.n_PS_kernel;

    if config.kernel.n_conn_kernel ~= K_conn
        error(['Configured config.kernel.n_conn_kernel=%d, but imported ' ...
               'kernel metadata reports n_conn_kernel=%d.'], ...
            config.kernel.n_conn_kernel, K_conn);
    end

    post_raster_info = [];
    full_count_source = [];
    loaded_post = [];
    loaded_pre = [];

    if is_fixed
        progress_log('NETWORK', ...
            ['Fixed mode: skipping reference raster and GLM files; ' ...
             'using external kernel file %s.'], kernels.source_path);
        cell_info = build_cell_layout(config, [], [], []);
    else
        post_raster_meta = make_reference_meta(config, 'Post', ...
            config.reference.cortex_area);
        post_raster_info = load_reference_raster_cell_info( ...
            root, post_raster_meta, config.simulation.default_dt);

        if strcmp(config.network_source, 'load')
            progress_log('NETWORK', 'Loading Post cortical GLM.');
            loaded_post = load_reference_glm_network( ...
                root, post_raster_meta, kernels, config.simulation.default_dt);
            cortical_source = loaded_post;
        else
            cortical_source = post_raster_info;
        end

        % If cell counts are matched to data, obtain the thalamic count from
        % a Pre Full source. When TC weights are loaded, the filtered Pre GLM
        % itself defines the simulated thalamic population. Otherwise only
        % the Pre Full raster is loaded and generic Thalamus cells are split
        % approximately in half for generated TC projections.
        if strcmp(config.cell_count_source, 'data')
            pre_full_meta = make_reference_meta(config, 'Pre', ...
                config.reference.tc_area);
            if strcmp(config.tc_source, 'load_pre')
                progress_log('NETWORK', ...
                    'Loading Pre Full GLM for TC weights and cell counts.');
                loaded_pre = load_reference_glm_network( ...
                    root, pre_full_meta, kernels, config.simulation.default_dt);
                full_count_source = loaded_pre;
            else
                progress_log('NETWORK', ...
                    'Loading Pre Full raster for thalamic cell count.');
                full_count_source = load_reference_raster_cell_info( ...
                    root, pre_full_meta, config.simulation.default_dt);
            end
        elseif strcmp(config.tc_source, 'load_pre')
            pre_full_meta = make_reference_meta(config, 'Pre', ...
                config.reference.tc_area);
            progress_log('NETWORK', ...
                'Loading Pre Full GLM for TC weights.');
            loaded_pre = load_reference_glm_network( ...
                root, pre_full_meta, kernels, ...
                config.simulation.default_dt);
        end

        cell_info = build_cell_layout( ...
            config, cortical_source, post_raster_info, full_count_source);
    end

    N = cell_info.N_total;
    network = struct();
    network.h = zeros(N, 1);
    network.P = zeros(N, K_ps);
    network.J = zeros(N, N, K_conn);
    network.provenance = struct();
    network.provenance.model_size_source = config.model_size_source;
    network.provenance.network_source = config.network_source;
    network.provenance.tc_source = config.tc_source;

    switch config.network_source
        case 'load'
            [network.h, network.P, network.J] = ...
                insert_loaded_cortex_network( ...
                network.h, network.P, network.J, ...
                loaded_post, cell_info);

        case 'generate'
            [network.h, network.P, network.J] = ...
                generate_cortical_network( ...
                network.h, network.P, network.J, ...
                cell_info, config);
    end

    [network.h, network.P] = generate_thalamic_intrinsic_parameters( ...
        network.h, network.P, cell_info, config);

    switch config.tc_source
        case 'load_pre'
            if isempty(loaded_pre)
                pre_meta = make_reference_meta(config, 'Pre', ...
                    config.reference.tc_area);
                progress_log('NETWORK', ...
                    'Loading Pre Full GLM for TC connections.');
                loaded_pre = load_reference_glm_network( ...
                    root, pre_meta, kernels, config.simulation.default_dt);
            end
            network.J = insert_loaded_tc_network( ...
                network.J, loaded_pre, cell_info, config);

        case 'generate'
            network.J = generate_tc_network(network.J, cell_info, config);

        case 'none'
            progress_log('NETWORK', 'No thalamocortical connections added.');
    end

    % Remove self-connections.
    for k = 1:K_conn
        Jk = network.J(:, :, k);
        Jk(1:N+1:end) = 0;
        network.J(:, :, k) = Jk;
    end

    network.cell_area = cell_info.cell_area;
    network.cell_type = cell_info.cell_type;
    network.cell_id = cell_info.cell_id;
    network.channel = cell_info.channel;
    network.masks = cell_info.masks;
    network.N = N;
    network.n_conn_kernel = K_conn;
    network.n_PS_kernel = K_ps;

    if is_fixed
        dt = config.fixed_model.dt;
    else
        dt = config.simulation.dt;
        if isempty(dt)
            dt = get_optional_struct_field( ...
                post_raster_info.meta, 'dt', config.simulation.default_dt);
        end
    end
    if isempty(dt) || ~isscalar(dt) || ~isfinite(dt) || dt <= 0
        error(['Simulation dt could not be resolved. Set config.simulation.dt ' ...
               'or config.fixed_model.dt to a positive value in seconds.']);
    end

    sim_meta = struct();
    sim_meta.model_size_source = config.model_size_source;
    if is_fixed
        sim_meta.reference = struct();
    else
        sim_meta.reference = config.reference;
    end
    sim_meta.cell_area = cell_info.cell_area;
    sim_meta.cell_type = cell_info.cell_type;
    sim_meta.cell_id = cell_info.cell_id;
    sim_meta.channel = cell_info.channel;
    sim_meta.N_total = cell_info.N_total;
    sim_meta.N_cortex = cell_info.N_cortex;
    sim_meta.N_thalamus = cell_info.N_thalamus;
    sim_meta.dt = dt;
    sim_meta.kernel_name = config.kernel.name;
    sim_meta.n_conn_kernel = K_conn;
    sim_meta.n_PS_kernel = K_ps;
    sim_meta.reg_name = resolve_effective_reg_name(config);
    sim_meta.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    source_info = struct();
    source_info.uses_external_reference = ~is_fixed;
    source_info.uses_external_data_reference = ~is_fixed;
    source_info.uses_external_kernel = true;
    source_info.post_raster_path = get_optional_struct_field( ...
        post_raster_info, 'source_path', '');
    source_info.full_count_source_path = get_optional_struct_field( ...
        full_count_source, 'source_path', '');
    source_info.post_glm_path = get_optional_struct_field( ...
        loaded_post, 'source_path', '');
    source_info.pre_glm_path = get_optional_struct_field( ...
        loaded_pre, 'source_path', '');
    source_info.kernel_path = kernels.source_path;

    progress_log('NETWORK', ...
        ['Network ready: mode=%s, N=%d, cortex=%d, thalamus=%d, ' ...
         'connection kernels=%d, PS kernels=%d, dt=%g s.'], ...
        config.model_size_source, sim_meta.N_total, ...
        sim_meta.N_cortex, sim_meta.N_thalamus, ...
        K_conn, K_ps, dt);
end
function cell_info = build_cell_layout( ...
    config, cortical_source, raster_source, thalamus_source)

    is_fixed = strcmp(config.model_size_source, 'fixed');
    use_loaded_tc = strcmp(config.tc_source, 'load_pre');
    source_thalamus_id = strings(0, 1);

    if is_fixed
        n_acc = config.fixed_model.cell_counts.ACC;
        n_vlpfc = config.fixed_model.cell_counts.VLPFC;
        n_thal_acc = config.fixed_model.cell_counts.thalamus_ACC;
        n_thal_vlpfc = config.fixed_model.cell_counts.thalamus_VLPFC;
        n_thal_generic = 0;
        source_cell_id = [ ...
            compose('fixed_ACC_%d', (1:n_acc).'); ...
            compose('fixed_VLPFC_%d', (1:n_vlpfc).')];
        reference_raster_meta = struct();
    else
        if strcmp(config.cell_count_source, 'manual')
            n_acc = config.manual_cell_counts.ACC;
            n_vlpfc = config.manual_cell_counts.VLPFC;
            source_cell_id = compose('manual_cortex_%d', ...
                (1:(n_acc + n_vlpfc)).');
            total_thalamus = 2 * config.thalamus_cells_per_area;
        else
            areas = string(cortical_source.cell_area(:));
            is_acc = ismember(areas, string(config.area_labels.ACC));
            is_vlpfc = ismember(areas, string(config.area_labels.VLPFC));
            source_idx = [find(is_acc); find(is_vlpfc)];

            if isempty(source_idx) || ~any(is_acc) || ~any(is_vlpfc)
                error(['Could not identify both ACC and VLPFC cells. ' ...
                       'Review config.area_labels.']);
            end

            n_acc = sum(is_acc);
            n_vlpfc = sum(is_vlpfc);
            source_cell_id = string(cortical_source.cell_id(source_idx));

            if isempty(thalamus_source) || ...
                    ~isfield(thalamus_source, 'cell_area')
                error(['cell_count_source=''data'' requires a Pre Full raster ' ...
                       'or GLM source containing Thalamus cells.']);
            end
            thal_areas = string(thalamus_source.cell_area(:));
            is_thalamus = ismember(thal_areas, ...
                string(config.area_labels.thalamus_generic)) | ...
                ismember(thal_areas, string(config.area_labels.thalamus_ACC)) | ...
                ismember(thal_areas, string(config.area_labels.thalamus_VLPFC));
            total_thalamus = sum(is_thalamus);
            source_thalamus_id = string( ...
                thalamus_source.cell_id(is_thalamus));
            if total_thalamus == 0
                error(['No Thalamus cells were found in the Pre Full count ' ...
                       'source. Review config.area_labels.']);
            end
        end

        if use_loaded_tc
            % Loaded TC weights already encode projection targets. Keep one
            % generic Thalamus population and do not invent a split.
            n_thal_acc = 0;
            n_thal_vlpfc = 0;
            n_thal_generic = total_thalamus;
        else
            % Generated TC weights need projection classes. Split generic
            % Thalamus cells approximately in half; the extra odd cell is
            % assigned to VLPFC.
            n_thal_acc = floor(total_thalamus / 2);
            n_thal_vlpfc = total_thalamus - n_thal_acc;
            n_thal_generic = 0;
            if mod(total_thalamus, 2) ~= 0
                progress_log('CELLS', ...
                    ['Odd Thalamus count=%d: assigning %d to Thal-ACC ' ...
                     'and %d to Thal-VLPFC.'], ...
                    total_thalamus, n_thal_acc, n_thal_vlpfc);
            end
        end
        reference_raster_meta = raster_source.meta;
    end

    n_cortex = n_acc + n_vlpfc;
    n_thalamus = n_thal_acc + n_thal_vlpfc + n_thal_generic;
    n_total = n_cortex + n_thalamus;

    idx_acc = (1:n_acc).';
    idx_vlpfc = (n_acc + (1:n_vlpfc)).';
    idx_thal_acc = (n_cortex + (1:n_thal_acc)).';
    idx_thal_vlpfc = (n_cortex + n_thal_acc + ...
        (1:n_thal_vlpfc)).';
    idx_thal_generic = (n_cortex + n_thal_acc + n_thal_vlpfc + ...
        (1:n_thal_generic)).';

    cell_area = strings(n_total, 1);
    cell_area(idx_acc) = "ACC";
    cell_area(idx_vlpfc) = "VLPFC";
    cell_area(n_cortex + 1:end) = "Thalamus";

    cell_type = strings(n_total, 1);
    cell_type(idx_acc) = "ACC";
    cell_type(idx_vlpfc) = "VLPFC";
    cell_type(idx_thal_acc) = "Thal-ACC";
    cell_type(idx_thal_vlpfc) = "Thal-VLPFC";
    cell_type(idx_thal_generic) = "Thalamus";

    cell_id = strings(n_total, 1);
    cell_id(1:n_cortex) = source_cell_id(:);
    cell_id(idx_thal_acc) = compose('sim_thalamus_%d', ...
        (1:n_thal_acc).');
    cell_id(idx_thal_vlpfc) = compose('sim_thalamus_%d', ...
        n_thal_acc + (1:n_thal_vlpfc).');
    if n_thal_generic > 0 && ...
            numel(source_thalamus_id) == n_thal_generic
        cell_id(idx_thal_generic) = source_thalamus_id(:);
    else
        cell_id(idx_thal_generic) = compose('sim_thalamus_%d', ...
            n_thal_acc + n_thal_vlpfc + (1:n_thal_generic).');
    end

    % Simulated data are already stored in the intended order, so channel is
    % a simple one-based ordering index rather than a physical probe channel.
    channel = (1:n_total).';

    masks = struct();
    masks.ACC = false(n_total, 1);
    masks.VLPFC = false(n_total, 1);
    masks.thalamus_ACC = false(n_total, 1);
    masks.thalamus_VLPFC = false(n_total, 1);
    masks.thalamus_generic = false(n_total, 1);
    masks.ACC(idx_acc) = true;
    masks.VLPFC(idx_vlpfc) = true;
    masks.thalamus_ACC(idx_thal_acc) = true;
    masks.thalamus_VLPFC(idx_thal_vlpfc) = true;
    masks.thalamus_generic(idx_thal_generic) = true;
    masks.cortex = masks.ACC | masks.VLPFC;
    masks.thalamus = masks.thalamus_ACC | ...
        masks.thalamus_VLPFC | masks.thalamus_generic;

    cell_info = struct();
    cell_info.N_total = n_total;
    cell_info.N_cortex = n_cortex;
    cell_info.N_thalamus = n_thalamus;
    cell_info.n_acc = n_acc;
    cell_info.n_vlpfc = n_vlpfc;
    cell_info.n_thal_acc = n_thal_acc;
    cell_info.n_thal_vlpfc = n_thal_vlpfc;
    cell_info.n_thal_generic = n_thal_generic;
    cell_info.cell_area = cell_area;
    cell_info.cell_type = cell_type;
    cell_info.cell_id = cell_id;
    cell_info.channel = channel;
    cell_info.masks = masks;
    cell_info.reference_raster_meta = reference_raster_meta;

    progress_log('CELLS', ...
        ['mode=%s, ACC=%d, VLPFC=%d, Thalamus=%d ' ...
         '(Thal-ACC=%d, Thal-VLPFC=%d, generic=%d).'], ...
        config.model_size_source, n_acc, n_vlpfc, n_thalamus, ...
        n_thal_acc, n_thal_vlpfc, n_thal_generic);
end
function meta = make_reference_meta(config, prepost, area)
    meta = struct();
    meta.animal_name = config.reference.animal_name;
    meta.injection = config.reference.injection;
    meta.prepost = prepost;
    meta.state = config.reference.state;
    meta.area = area;
    meta.align = config.reference.align;
    meta.session_idx = config.reference.session_idx;
    meta.resting_dur_threshold = config.reference.resting_dur_threshold;
    meta.kernel_name = config.kernel.name;
    meta.reg_name = config.reference.reg_name;
    meta.epoch = config.reference.epoch;
    meta.fold_idx = config.reference.fold_idx;
    meta.shuffle_idx = config.reference.shuffle_idx;
end

%% =========================================================================
%  README-BASED INPUT ADAPTERS
%  =========================================================================
function raster_info = load_reference_raster_cell_info(root, meta, default_dt)
    meta.filename = generate_filename('raster', meta);
    raster_path = fullfile(root, 'Data', 'Working', 'raster', meta.filename);
    progress_log('LOAD-RASTER', 'Loading: %s', raster_path);

    if ~isfile(raster_path)
        error('Reference raster file not found: %s', raster_path);
    end

    loaded = load(raster_path, 'meta', 'data');
    if ~isfield(loaded, 'meta') || ~isfield(loaded, 'data')
        error('Raster file must contain structs meta and data: %s', raster_path);
    end
    if ~isfield(loaded.meta, 'N') || ...
            ~isfield(loaded.data, 'cell_area') || ...
            ~isfield(loaded.data, 'cell_id')
        error(['Raster file is missing required fields ' ...
               'meta.N/data.cell_area/data.cell_id: %s'], raster_path);
    end

    if ~isfield(loaded.meta, 'dt') || ...
            ~isscalar(loaded.meta.dt) || ...
            ~isfinite(loaded.meta.dt) || loaded.meta.dt <= 0
        loaded.meta.dt = default_dt;
        progress_log('LOAD-RASTER', ...
            'meta.dt missing/invalid; using default dt=%g s.', default_dt);
    end

    N = double(loaded.meta.N);
    cell_area = string(loaded.data.cell_area(:));
    cell_id = string(loaded.data.cell_id(:));
    if numel(cell_area) ~= N || numel(cell_id) ~= N
        error('Raster cell metadata lengths do not match meta.N: %s', raster_path);
    end

    if isfield(loaded.data, 'cell_type') && ...
            numel(loaded.data.cell_type) == N
        cell_type = string(loaded.data.cell_type(:));
    else
        cell_type = cell_area;
    end

    if isfield(loaded.data, 'channel') && numel(loaded.data.channel) == N
        channel = double(loaded.data.channel(:));
    else
        channel = (1:N).';
        progress_log('LOAD-RASTER', ...
            'data.channel missing or wrong length; using sequential channels.');
    end

    raster_info = struct();
    raster_info.meta = loaded.meta;
    raster_info.cell_area = cell_area;
    raster_info.cell_type = cell_type;
    raster_info.cell_id = cell_id;
    raster_info.channel = channel;
    raster_info.N = N;
    raster_info.source_path = raster_path;
end
function kernels = load_kernel_definition(root, config)
    % Temporal kernels are always loaded from the project kernel file.
    % model_size_source controls only cell-size/dt and raster/GLM references.
    kernel_meta = struct();
    kernel_meta.kernel_name = config.kernel.name;
    kernel_filename = generate_filename('kernel', kernel_meta);
    kernel_path = fullfile(root, 'Data', 'Working', 'kernel', kernel_filename);
    progress_log('LOAD-KERNEL', 'Loading: %s', kernel_path);

    if ~isfile(kernel_path)
        error('Kernel file not found: %s', kernel_path);
    end

    loaded = load(kernel_path, 'meta', 'data');
    if ~isfield(loaded, 'meta') || ~isfield(loaded, 'data')
        error('Kernel file must contain structs meta and data: %s', kernel_path);
    end

    required_meta = {'kernel_name', 'n_conn_kernel', 'n_PS_kernel', 'kernel_len'};
    required_data = {'conn_kernels', 'PS_kernels'};
    assert_struct_fields(loaded.meta, required_meta, 'kernel meta');
    assert_struct_fields(loaded.data, required_data, 'kernel data');

    imported_kernel_name = string(loaded.meta.kernel_name);
    configured_kernel_name = string(config.kernel.name);
    if ~isscalar(imported_kernel_name) || ...
            imported_kernel_name ~= configured_kernel_name
        error(['Imported kernel metadata name (%s) does not match ' ...
               'config.kernel.name (%s).'], ...
            char(imported_kernel_name), char(configured_kernel_name));
    end

    n_conn = double(loaded.meta.n_conn_kernel);
    n_ps = double(loaded.meta.n_PS_kernel);
    if n_conn ~= config.kernel.n_conn_kernel
        error(['Imported kernel %s contains %d connection kernels, but ' ...
               'config.kernel.n_conn_kernel=%d.'], ...
            config.kernel.name, n_conn, config.kernel.n_conn_kernel);
    end
    if numel(loaded.data.conn_kernels) ~= n_conn
        error('conn_kernels count does not match meta.n_conn_kernel.');
    end
    if numel(loaded.data.PS_kernels) ~= n_ps
        error(['PS_kernels count does not match meta.n_PS_kernel. ' ...
               'README table may contain a shape typo; the code follows n_PS_kernel.']);
    end

    kernels = struct();
    kernels.meta = loaded.meta;
    kernels.conn_values = normalize_temporal_kernel_cells( ...
        loaded.data.conn_kernels, n_conn, 'connection');
    kernels.PS_values = normalize_temporal_kernel_cells( ...
        loaded.data.PS_kernels, n_ps, 'post-spike');
    kernels.conn_lengths = cellfun(@numel, kernels.conn_values);
    kernels.PS_lengths = cellfun(@numel, kernels.PS_values);
    kernels.n_conn_kernel = n_conn;
    kernels.n_PS_kernel = n_ps;
    kernels.kernel_len = double(loaded.meta.kernel_len);
    kernels.source_path = kernel_path;

    progress_log('LOAD-KERNEL', ...
        ['Validated kernel=%s, n_conn=%d, cortex_scale=%s, ' ...
         'tc_scale=%s.'], ...
        config.kernel.name, n_conn, ...
        mat2str(config.kernel.cortex_scale), ...
        mat2str(config.kernel.tc_scale));
end


function values = normalize_temporal_kernel_cells(raw_values, expected_n, label)
    values = cell(1, expected_n);
    for k = 1:expected_n
        value = double(raw_values{k});
        if ~isvector(value)
            error('%s kernel %d must be one-dimensional.', label, k);
        end
        value = reshape(value, 1, []);
        if any(~isfinite(value))
            error('%s kernel %d contains nonfinite values.', label, k);
        end
        values{k} = value(:, end:-1:1);
    end
end

function loaded_network = load_reference_glm_network( ...
    root, meta, kernels, default_dt)
    meta.filename = generate_filename('GLM', meta);
    glm_path = fullfile(root, 'Data', 'Working', 'GLM', meta.filename);
    progress_log('LOAD-GLM', 'Loading: %s', glm_path);

    if ~isfile(glm_path)
        error('Reference GLM file not found: %s', glm_path);
    end

    loaded = load(glm_path, 'meta', 'data');
    if ~isfield(loaded, 'meta') || ~isfield(loaded, 'data')
        error('GLM file must contain structs meta and data: %s', glm_path);
    end
    assert_struct_fields(loaded.data, {'model_par'}, 'GLM data');

    model_par = double(loaded.data.model_par);
    n_target = size(model_par, 1);
    n_ps = kernels.n_PS_kernel;
    n_conn = kernels.n_conn_kernel;

    % README layout:
    %   col 1                         = h
    %   cols 2 : 1+n_PS_kernel        = P
    %   remaining N*n_conn_kernel     = J
    connection_start_col = 2 + n_ps;
    connection_column_count = size(model_par, 2) - connection_start_col + 1;
    if connection_column_count <= 0 || ...
            mod(connection_column_count, n_conn) ~= 0
        error(['model_par columns are incompatible with README layout: ' ...
               'size=%s, n_PS_kernel=%d, n_conn_kernel=%d.'], ...
            mat2str(size(model_par)), n_ps, n_conn);
    end
    n_source = connection_column_count / n_conn;
    if n_source ~= n_target
        error(['The current simulator requires a square GLM connection matrix, ' ...
               'but model_par implies %d targets and %d sources. ' ...
               'A project-specific filtered-neuron mapping is required.'], ...
            n_target, n_source);
    end

    h = model_par(:, 1);
    h(~isfinite(h)) = -20;

    P = zeros(n_target, n_ps);
    if n_ps > 0
        P = model_par(:, 2:(1 + n_ps));
        P(~isfinite(P)) = 0;
    end

    J = zeros(n_target, n_source, n_conn);
    for k = 1:n_conn
        col_start = connection_start_col + n_source * (k - 1);
        col_end = col_start + n_source - 1;
        Jk = model_par(:, col_start:col_end);
        Jk(~isfinite(Jk)) = 0;
        J(:, :, k) = Jk;
    end

    raster_info = load_reference_raster_cell_info(root, meta, default_dt);
    raster_indices = resolve_glm_to_raster_indices( ...
        loaded.meta, loaded.data, raster_info.N, n_target);

    loaded_network = struct();
    loaded_network.meta = loaded.meta;
    loaded_network.N = n_target;
    loaded_network.h = h;
    loaded_network.P = P;
    loaded_network.J = J;
    loaded_network.cell_area = raster_info.cell_area(raster_indices);
    loaded_network.cell_type = raster_info.cell_type(raster_indices);
    loaded_network.cell_id = raster_info.cell_id(raster_indices);
    loaded_network.channel = raster_info.channel(raster_indices);
    loaded_network.raster_indices = raster_indices;
    loaded_network.source_path = glm_path;
end

function raster_indices = resolve_glm_to_raster_indices( ...
    glm_meta, glm_data, raster_N, model_N)

    if raster_N == model_N
        raster_indices = (1:raster_N).';
        return;
    end

    if isfield(glm_data, 'filter') && numel(glm_data.filter) == raster_N
        filter_mask = logical(glm_data.filter(:));
        if sum(filter_mask) == model_N
            raster_indices = find(filter_mask);
            return;
        end
    end

    if isfield(glm_meta, 'N_filtered') && glm_meta.N_filtered == model_N
        error(['GLM meta.N_filtered matches model rows, but data.filter does ' ...
               'not provide a usable mapping to the raster. Supply the exact ' ...
               'filtered-neuron index mapping in resolve_glm_to_raster_indices().']);
    end

    error(['Cannot map a %d-neuron GLM to a %d-neuron raster. ' ...
           'Review data.filter and meta.N_filtered.'], model_N, raster_N);
end

%% =========================================================================
%  NETWORK INSERTION AND GENERATION
%  =========================================================================
function [h_out, P_out, J_out] = insert_loaded_cortex_network( ...
    h_out, P_out, J_out, loaded_post, cell_info)

    source_area = string(loaded_post.cell_area(:));
    source_idx = [find(source_area == "ACC"); find(source_area == "VLPFC")];
    target_idx = find(cell_info.masks.cortex);

    if numel(source_idx) ~= numel(target_idx)
        error(['Loaded cortical cell count (%d) does not match the simulation ' ...
               'cortical count (%d).'], numel(source_idx), numel(target_idx));
    end

    h_out(target_idx) = loaded_post.h(source_idx);
    P_out(target_idx, :) = loaded_post.P(source_idx, :);
    for k = 1:size(J_out, 3)
        J_out(target_idx, target_idx, k) = ...
            loaded_post.J(source_idx, source_idx, k);
    end

    progress_log('NETWORK', 'Inserted loaded cortical h, P, and J.');
end

function J_out = insert_loaded_tc_network( ...
    J_out, loaded_pre, cell_info, config)

    % Loaded TC weights already contain their target-specific structure.
    % Therefore all source neurons with cell_area='Thalamus' are treated as
    % one population; no Thal-ACC/Thal-VLPFC classification is invented.
    areas = string(loaded_pre.cell_area(:));
    all_source_ids = string(loaded_pre.cell_id(:));
    source_cortex_candidates = find( ...
        ismember(areas, string(config.area_labels.ACC)) | ...
        ismember(areas, string(config.area_labels.VLPFC)));
    source_thalamus_candidates = find( ...
        ismember(areas, string(config.area_labels.thalamus_generic)) | ...
        ismember(areas, string(config.area_labels.thalamus_ACC)) | ...
        ismember(areas, string(config.area_labels.thalamus_VLPFC)));

    if isempty(source_cortex_candidates)
        error('Pre Full GLM contains no identifiable cortical target cells.');
    end
    if isempty(source_thalamus_candidates)
        error(['tc_source=''load_pre'' requested, but no Thalamus cells ' ...
               'were found in the Pre Full GLM.']);
    end

    dst_cortex = find(cell_info.masks.cortex);
    dst_thalamus = find(cell_info.masks.thalamus);

    % Prefer explicit cell-ID alignment. This protects against different
    % neuron ordering between the Post cortical model and the Pre Full model.
    dst_cortex_ids = string(cell_info.cell_id(dst_cortex));
    [cortex_matched, src_cortex] = ismember( ...
        dst_cortex_ids, all_source_ids);
    src_cortex = src_cortex(:);
    cortex_ids_are_manual = all( ...
        startsWith(dst_cortex_ids, "manual_") | ...
        startsWith(dst_cortex_ids, "fixed_"));
    if ~all(cortex_matched) || ...
            any(~ismember(src_cortex, source_cortex_candidates))
        if cortex_ids_are_manual
            progress_log('NETWORK', ...
                ['Manual cortical identities cannot be matched by cell_id; ' ...
                 'using ACC/VLPFC area order.']);
            src_acc = find(ismember(areas, string(config.area_labels.ACC)));
            src_vlpfc = find(ismember(areas, string(config.area_labels.VLPFC)));
            src_cortex = [src_acc; src_vlpfc];
        else
            error(['Pre Full and simulated cortical cell_id values do not ' ...
                   'match. A safe TC target mapping cannot be constructed.']);
        end
    end

    dst_thalamus_ids = string(cell_info.cell_id(dst_thalamus));
    [thalamus_matched, src_thalamus] = ismember( ...
        dst_thalamus_ids, all_source_ids);
    src_thalamus = src_thalamus(:);
    thalamus_ids_are_simulated = all(startsWith( ...
        dst_thalamus_ids, "sim_thalamus_"));
    if ~all(thalamus_matched) || ...
            any(~ismember(src_thalamus, source_thalamus_candidates))
        if thalamus_ids_are_simulated
            progress_log('NETWORK', ...
                ['Simulated Thalamus identities cannot be matched by cell_id; ' ...
                 'using source-file Thalamus order.']);
            src_thalamus = source_thalamus_candidates;
        else
            error(['Pre Full and simulated Thalamus cell_id values do not ' ...
                   'match. A safe TC source mapping cannot be constructed.']);
        end
    end

    assert_matching_count(src_cortex, dst_cortex, 'cortical');
    assert_matching_count(src_thalamus, dst_thalamus, 'Thalamus');

    for k = 1:size(J_out, 3)
        J_out(dst_cortex, dst_thalamus, k) = ...
            loaded_pre.J(src_cortex, src_thalamus, k);
    end

    progress_log('NETWORK', ...
        'Inserted loaded Pre thalamocortical J without projection splitting.');
end
function assert_matching_count(source_idx, target_idx, label)
    if numel(source_idx) ~= numel(target_idx)
        error(['Loaded %s count (%d) does not match configured simulation ' ...
               'count (%d). A resampling or cell-matching rule is required.'], ...
            label, numel(source_idx), numel(target_idx));
    end
end

function [h, P, J] = generate_cortical_network( ...
    h, P, J, cell_info, config)

    rng(config.simulation.base_seed, 'twister');
    area_indices = {find(cell_info.masks.ACC), ...
                    find(cell_info.masks.VLPFC)};

    baseline_mean = config.generate.cortex.baseline_mean;
    baseline_sd = config.generate.cortex.baseline_sd;
    for area_i = 1:2
        idx = area_indices{area_i};
        h(idx) = baseline_mean(area_i) + ...
            baseline_sd(area_i) * randn(numel(idx), 1);
    end

    if size(P, 2) > 0
        P(cell_info.masks.cortex, :) = ...
            config.generate.cortex.post_spike_weight_mean + ...
            config.generate.cortex.post_spike_weight_sd .* ...
            randn(sum(cell_info.masks.cortex), size(P, 2));
    end

    kernel_scale = config.kernel.cortex_scale(:).';

    for target_area = 1:2
        target_idx = area_indices{target_area};
        for source_area = 1:2
            source_idx = area_indices{source_area};

            connection_mask = rand(numel(target_idx), numel(source_idx)) < ...
                config.generate.cortex.connection_probability( ...
                target_area, source_area);
            positive_mask = rand(numel(target_idx), numel(source_idx)) < ...
                config.generate.cortex.positive_probability( ...
                target_area, source_area);
            signs = ones(size(connection_mask));
            signs(~positive_mask) = -1;

            base_abs = config.generate.cortex.weight_abs_mean( ...
                target_area, source_area) .* ...
                max(1 + config.generate.cortex.weight_cv .* ...
                randn(size(connection_mask)), 0.05);
            base_weights = connection_mask .* signs .* base_abs;

            for k = 1:size(J, 3)
                J(target_idx, source_idx, k) = ...
                    kernel_scale(k) .* base_weights;
            end
        end
    end

    progress_log('GENERATE', 'Generated cortical h, P, and J.');
end

function [h, P] = generate_thalamic_intrinsic_parameters( ...
    h, P, cell_info, config)

    rng(config.simulation.base_seed + 1, 'twister');
    thal_idx = find(cell_info.masks.thalamus);
    h(thal_idx) = config.generate.thalamus.baseline_mean + ...
        config.generate.thalamus.baseline_sd .* randn(numel(thal_idx), 1);

    if size(P, 2) > 0
        P(thal_idx, :) = config.generate.thalamus.post_spike_weight_mean + ...
            config.generate.thalamus.post_spike_weight_sd .* ...
            randn(numel(thal_idx), size(P, 2));
    end
end

function J = generate_tc_network(J, cell_info, config)
    rng(config.simulation.base_seed + 2, 'twister');
    kernel_scale = config.kernel.tc_scale(:).';

    thal_acc_idx = find(cell_info.masks.thalamus_ACC);
    thal_vlpfc_idx = find(cell_info.masks.thalamus_VLPFC);
    if isempty(thal_acc_idx) || isempty(thal_vlpfc_idx)
        error(['Generated TC connections require Thalamus cells to be split ' ...
               'into cell_type Thal-ACC and Thal-VLPFC.']);
    end

    J = generate_one_tc_projection(J, ...
        find(cell_info.masks.ACC), thal_acc_idx, ...
        config, kernel_scale);
    J = generate_one_tc_projection(J, ...
        find(cell_info.masks.VLPFC), thal_vlpfc_idx, ...
        config, kernel_scale);

    progress_log('GENERATE', 'Generated thalamocortical J.');
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
        J(target_idx, source_idx, k) = kernel_scale(k) .* base_weights;
    end
end


%% =========================================================================
%  SIMULATION
%  =========================================================================
function results = simulate_all_conditions(config, sim_meta, network, kernels)
    n_conditions = numel(config.conditions);
    if n_conditions == 0
        error('config.conditions must contain at least one condition.');
    end

    for condition_i = 1:n_conditions
        condition = config.conditions(condition_i);
        progress_log('SIMULATION', '[%d/%d] START %s.', ...
            condition_i, n_conditions, condition.name);

        condition_network = network;
        if ~condition.enable_tc
            condition_network.J(:, condition_network.masks.thalamus, :) = 0;
        end

        condition_result = simulate_condition( ...
            config, sim_meta, condition_network, kernels, ...
            condition, condition_i);

        % MATLAB requires all elements of a structure array to have the same
        % field set. Assign the first populated result to the whole variable
        % so it establishes the structure-array schema. Indexed assignments
        % are safe after that first initialization.
        if condition_i == 1
            results = condition_result;
        else
            results(condition_i) = condition_result;
        end

        progress_log('SIMULATION', '[%d/%d] DONE %s.', ...
            condition_i, n_conditions, condition.name);
    end
end

function result = simulate_condition( ...
    config, sim_meta, network, kernels, condition, condition_i)

    n_trials = config.simulation.n_trials;
    n_time_bins = config.simulation.n_time_bins;
    N = network.N;

    trial_rasters = cell(1, n_trials);
    trial_firing_rates = cell(1, n_trials);

    if config.simulation.parallel_trials
        progress_log('SIMULATION', 'Using parfor for %d trials.', n_trials);
        parfor trial_i = 1:n_trials
            trial_seed = config.simulation.base_seed + ...
                100000 * condition_i + trial_i;
            raster = simulate_one_trial( ...
                network, kernels, condition, n_time_bins, trial_seed);
            trial_rasters{trial_i} = raster;
            trial_firing_rates{trial_i} = ...
                sum(double(raster), 2) ./ (n_time_bins * sim_meta.dt);
        end
    else
        for trial_i = 1:n_trials
            trial_tic = tic;
            trial_seed = config.simulation.base_seed + ...
                100000 * condition_i + trial_i;
            raster = simulate_one_trial( ...
                network, kernels, condition, n_time_bins, trial_seed);
            trial_rasters{trial_i} = raster;
            trial_firing_rates{trial_i} = ...
                sum(double(raster), 2) ./ (n_time_bins * sim_meta.dt);
            progress_log('TRIAL', '%s trial %d/%d done in %.1f s.', ...
                condition.name, trial_i, n_trials, toc(trial_tic));
        end
    end

    % Keep one complete simulation result. Full and Cortex-only project
    % raster views are derived later without rerunning the stochastic model.
    for trial_i = 1:n_trials
        trial_rasters{trial_i} = uint8(trial_rasters{trial_i});
        trial_firing_rates{trial_i} = trial_firing_rates{trial_i}.';
    end

    result = struct();
    result.condition = condition;
    result.N = N;
    result.dt = sim_meta.dt;
    result.cell_area = sim_meta.cell_area.';
    result.cell_type = sim_meta.cell_type.';
    result.cell_id = sim_meta.cell_id.';
    result.channel = sim_meta.channel.';
    result.trial_num = n_trials;
    result.trial_len = repmat(n_time_bins, 1, n_trials);
    result.rasters = trial_rasters;
    result.firing_rates = trial_firing_rates;
    result.cortex_mask = network.masks.cortex;
end

function raster = simulate_one_trial( ...
    network, kernels, condition, n_time_bins, seed)

    rng(seed, 'twister');
    N = network.N;
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

        % Connection-history terms.
        for k = 1:kernels.n_conn_kernel
            L = kernels.conn_lengths(k);
            if t > L
                history = double(raster(:, (t - L):(t - 1)));
                convolved = sum(history .* kernels.conn_values{k}, 2);
                h_total = h_total + network.J(:, :, k) * convolved;
            end
        end

        % Post-spike-history terms from README model_par layout.
        for k = 1:kernels.n_PS_kernel
            L = kernels.PS_lengths(k);
            if t > L
                history = double(raster(:, (t - L):(t - 1)));
                convolved = sum(history .* kernels.PS_values{k}, 2);
                h_total = h_total + network.P(:, k) .* convolved;
            end
        end

        h_total = max(min(h_total, 30), -30);
        spike_probability = 1 ./ (1 + exp(-h_total));
        raster(:, t) = rand(N, 1) < spike_probability;
    end
end

%% =========================================================================
%  README-BASED OUTPUT
%  =========================================================================
function session_output = save_simulation_session_output( ...
    root, config, sim_meta, results, output_session_idx)

    raster_folder = fullfile(root, 'Data', 'Working', 'raster');
    border_folder = fullfile(root, 'Data', 'Working', 'border');
    check_path(raster_folder);
    check_path(border_folder);

    output_areas = {'Full', 'Cortex'};
    n_conditions = numel(results);
    n_areas = numel(output_areas);

    session_output = struct();
    session_output.raster_paths = strings(n_conditions, n_areas);
    session_output.border_paths = strings(1, n_areas);
    session_output.area_labels = string(output_areas);
    session_output.session_idx = output_session_idx;

    for condition_i = 1:n_conditions
        full_result = results(condition_i);

        for area_i = 1:n_areas
            area_label = output_areas{area_i};
            result_view = build_result_area_view(full_result, area_label);

            raster_meta = build_output_raster_meta( ...
                config, sim_meta, result_view, ...
                output_session_idx, area_label);
            raster_data = build_output_raster_data(result_view);

            if config.output.save_project_raster
                raster_meta.file_name = generate_filename( ...
                    'raster', raster_meta);
                raster_path = fullfile( ...
                    raster_folder, raster_meta.file_name);
                assert_output_write_allowed( ...
                    raster_path, config.output.overwrite);

                meta = raster_meta; %#ok<NASGU>
                data = raster_data; %#ok<NASGU>
                progress_log('SAVE-RASTER', ...
                    ['session=%d, state=%s, area=%s: %s'], ...
                    output_session_idx, raster_meta.state, ...
                    area_label, raster_path);
                save(raster_path, 'meta', 'data', '-v7.3');
                session_output.raster_paths(condition_i, area_i) = ...
                    string(raster_path);
            end
        end
    end

    % Border metadata does not include state. All conditions in one
    % simulated session share the same cell layout, so save one border per
    % area and reuse it for every state.
    border_source_result = results(1);
    for area_i = 1:n_areas
        area_label = output_areas{area_i};
        result_view = build_result_area_view( ...
            border_source_result, area_label);
        raster_meta = build_output_raster_meta( ...
            config, sim_meta, result_view, ...
            output_session_idx, area_label);
        [border_meta, border_data] = build_output_border( ...
            raster_meta, result_view);

        if config.output.save_project_border
            border_meta.file_name = generate_filename( ...
                'border', border_meta);
            border_path = fullfile( ...
                border_folder, border_meta.file_name);
            assert_output_write_allowed( ...
                border_path, config.output.overwrite);

            meta = border_meta; %#ok<NASGU>
            data = border_data; %#ok<NASGU>
            progress_log('SAVE-BORDER', ...
                'session=%d, area=%s: %s', ...
                output_session_idx, area_label, border_path);
            save(border_path, 'meta', 'data');
            session_output.border_paths(area_i) = string(border_path);
        end
    end
end

function manifest_path = save_simulation_run_manifest( ...
    root, config, kernels, session_records, output_manifest)

    manifest_folder = fullfile(root, 'Data', 'Working', ...
        'simulation_v2', sanitize_for_path(config.output.simulation_tag));
    check_path(manifest_folder);

    manifest_path = fullfile(manifest_folder, sprintf( ...
        'simulation_manifest_%s.mat', ...
        sanitize_for_path(config.output.simulation_tag)));
    assert_output_write_allowed( ...
        manifest_path, config.output.overwrite);

    meta = struct();
    meta.simulation_tag = config.output.simulation_tag;
    meta.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    meta.model_size_source = config.model_size_source;
    meta.network_source = config.network_source;
    meta.tc_source = config.tc_source;
    meta.kernel_name = config.kernel.name;
    meta.reg_name = resolve_effective_reg_name(config);
    meta.n_sessions = config.simulation.n_sessions;
    meta.first_session_idx = config.output.base_session_idx;
    meta.last_session_idx = config.output.base_session_idx + ...
        config.simulation.n_sessions - 1;

    output_manifest.manifest_path = string(manifest_path);

    data = struct();
    data.config = config;
    data.kernels = kernels;
    data.session_records = session_records;
    data.output_manifest = output_manifest;

    progress_log('SAVE-MANIFEST', 'Writing: %s', manifest_path);
    save(manifest_path, 'meta', 'data', '-v7.3');
    manifest_path = string(manifest_path);
end


function result_view = build_result_area_view(full_result, area_label)
    switch area_label
        case 'Full'
            cell_mask = true(full_result.N, 1);
        case 'Cortex'
            cell_mask = logical(full_result.cortex_mask(:));
        otherwise
            error('Unknown output area view: %s', area_label);
    end

    result_view = full_result;
    result_view.N = sum(cell_mask);
    result_view.cell_area = full_result.cell_area(cell_mask);
    result_view.cell_type = full_result.cell_type(cell_mask);
    result_view.cell_id = full_result.cell_id(cell_mask);
    result_view.channel = 1:result_view.N;

    for trial_i = 1:full_result.trial_num
        result_view.rasters{trial_i} = ...
            full_result.rasters{trial_i}(cell_mask, :);
        result_view.firing_rates{trial_i} = ...
            full_result.firing_rates{trial_i}(cell_mask);
    end

    result_view.area_label = area_label;
    result_view.cell_mask = cell_mask;
end

function meta = build_output_raster_meta( ...
    config, sim_meta, result, output_session_idx, area_label)
    condition = result.condition;

    meta = struct();
    meta.animal_name = config.output.animal_name;
    meta.injection = config.output.injection;
    meta.prepost = config.output.prepost;
    meta.state = condition.output_state;
    meta.area = area_label;
    meta.align = config.output.align;
    meta.session_idx = output_session_idx;
    meta.date = datestr(now, 'mmddyyyy');
    meta.N = result.N;
    meta.dt = result.dt;
    meta.trial_num = result.trial_num;
    meta.trial_len = result.trial_len;
    meta.max_len = max(result.trial_len);
    meta.min_len = min(result.trial_len);
    meta.total_len = sum(result.trial_len);
    meta.resting_dur_threshold = ...
        config.output.resting_dur_threshold;

    % Extra provenance fields are allowed but are not required identifiers.
    meta.simulation_tag = config.output.simulation_tag;
    meta.model_size_source = config.model_size_source;
    meta.simulation_session_number = ...
        sim_meta.simulation_session_number;
    meta.session_seed = sim_meta.session_seed;
    if strcmp(config.model_size_source, 'fixed')
        meta.source_animal_name = '';
        meta.source_session_idx = NaN;
    else
        meta.source_animal_name = config.reference.animal_name;
        meta.source_session_idx = config.reference.session_idx;
    end
    meta.network_source = config.network_source;
    meta.tc_source = config.tc_source;
    meta.kernel_name = sim_meta.kernel_name;
    meta.reg_name = sim_meta.reg_name;
end

function data = build_output_raster_data(result)
    data = struct();
    data.rasters = result.rasters;
    data.spikes = raster_to_spike_times( ...
        result.rasters, result.dt);
    data.trial_len = reshape(result.trial_len, 1, []);
    data.cell_id = reshape(string(result.cell_id), 1, []);
    data.cell_area = reshape(string(result.cell_area), 1, []);
    data.cell_type = reshape(string(result.cell_type), 1, []);
    data.channel = reshape(double(result.channel), 1, []);
    data.cuetype = repmat({[]}, 1, result.trial_num);
    data.firing_rates = result.firing_rates;
end
function spikes = raster_to_spike_times(rasters, dt)
    % Output format:
    %   spikes{cell_i, trial_i} = spike_count x 1 double vector in seconds.
    %
    % A spike in raster bin b is placed at the center of the bin:
    %   time = (b - 0.5) * dt.
    n_trials = numel(rasters);
    if n_trials == 0
        spikes = cell(0, 0);
        return;
    end

    n_cells = size(rasters{1}, 1);
    spikes = cell(n_cells, n_trials);

    for trial_i = 1:n_trials
        raster = rasters{trial_i};
        if size(raster, 1) ~= n_cells
            error('All trial rasters must contain the same number of cells.');
        end

        for cell_i = 1:n_cells
            spike_bins = find(raster(cell_i, :) ~= 0);
            spike_times = (double(spike_bins(:)) - 0.5) .* dt;
            spikes{cell_i, trial_i} = spike_times;
        end
    end
end

function [meta, data] = build_output_border(raster_meta, result)
    meta = struct();
    meta.animal_name = raster_meta.animal_name;
    meta.injection = raster_meta.injection;
    meta.prepost = raster_meta.prepost;
    meta.area = raster_meta.area;
    meta.align = raster_meta.align;
    meta.session_idx = raster_meta.session_idx;
    meta.N = raster_meta.N;

    areas = string(result.cell_area(:));
    run_starts = [1; find(areas(2:end) ~= areas(1:end-1)) + 1];
    borders = [run_starts.', result.N + 1];
    meta.area_num = numel(run_starts);

    data = struct();
    data.borders = borders;
end


function validate_manual_kernel_scale(scale, expected_n, field_name)
    scale = double(scale(:).');
    if numel(scale) ~= expected_n
        error('%s must contain exactly %d values; received %d.', ...
            field_name, expected_n, numel(scale));
    end
    if any(~isfinite(scale))
        error('%s contains nonfinite values.', field_name);
    end
end


function assert_output_write_allowed(path, overwrite)
    if isfile(path) && ~overwrite
        error(['Output already exists: %s. Set ' ...
               'config.output.overwrite=true to replace it.'], path);
    end
end

%% =========================================================================
%  UTILITIES
%  =========================================================================
function assert_struct_fields(value, fields, label)
    missing = fields(~isfield(value, fields));
    if ~isempty(missing)
        error('%s is missing fields: %s', label, strjoin(missing, ', '));
    end
end

function name = resolve_effective_reg_name(config)
    if strcmp(config.model_size_source, 'fixed')
        name = config.fixed_model.reg_name;
    else
        name = config.reference.reg_name;
    end
end

function value = get_optional_struct_field(source, field, fallback)
    if isstruct(source) && isfield(source, field)
        value = source.(field);
    else
        value = fallback;
    end
end

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
