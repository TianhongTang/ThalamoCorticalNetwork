%% simulation_hybrid_network_v2.m
% Hybrid GLM network simulation using the project meta/data file format.
%
% Model-size source:
%   config.model_size_source = 'match_data'
%       Match cell counts, dt, and temporal kernels to project data.
%
%   config.model_size_source = 'fixed'
%       Use only fixed settings defined in config.fixed_model. No raster,
%       GLM, or kernel data files are loaded. Cortical and TC networks are
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
% The input/output adapters follow README.md. Remaining project-specific
% assumptions are marked "MANUAL CONFIRMATION".

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
% MODEL SIZE / EXTERNAL REFERENCE SWITCH
% -------------------------------------------------------------------------
% 'match_data' = use reference raster/kernel size and optionally GLM values.
% 'fixed'      = use config.fixed_model only; load no raster/GLM/kernel data.
config.model_size_source = 'match_data';

% -------------------------------------------------------------------------
% REFERENCE SESSION
% -------------------------------------------------------------------------
config.reference.animal_name = 'Slayer';
config.reference.injection = 'Muscimol';
config.reference.session_idx = 1;
config.reference.state = 'RestOpen';
config.reference.align = 'Last';
config.reference.resting_dur_threshold = 15;

% The cortical Post model is normally fitted with area='Cortex'.
config.reference.cortex_area = 'Cortex';

% A Pre model containing thalamic cells is normally expected to use
% area='Full'. Change this if the project uses another identifier.
config.reference.tc_area = 'Full';

config.reference.kernel_name = 'DeltaPure';
config.reference.kernel_num = []; % [] = read n_conn_kernel from kernel metadata.
config.reference.reg_name = 'L2=0_2';
config.reference.epoch = 3000;
config.reference.fold_idx = 0;
config.reference.shuffle_idx = 0;

% -------------------------------------------------------------------------
% FIXED MODEL PRESET
% Used only when config.model_size_source = 'fixed'. Time constants and
% kernel length are expressed in simulation bins.
% -------------------------------------------------------------------------
config.fixed_model.cell_counts.ACC = 40;
config.fixed_model.cell_counts.VLPFC = 40;
config.fixed_model.cell_counts.thalamus_ACC = 10;
config.fixed_model.cell_counts.thalamus_VLPFC = 10;
config.fixed_model.dt = 0.001;
config.fixed_model.kernel_name = 'FixedExpTriple';
config.fixed_model.reg_name = 'Generated';
config.fixed_model.kernel_len = 100;
config.fixed_model.conn_time_constants_bins = [5, 20, 80];
config.fixed_model.PS_time_constants_bins = 10;
config.fixed_model.normalize_temporal_kernels = true;

% -------------------------------------------------------------------------
% MAIN MODE SWITCHES
% -------------------------------------------------------------------------
config.network_source = 'load';       % 'load' or 'generate'
config.tc_source = 'generate';        % 'load_pre', 'generate', or 'none'
config.cell_count_source = 'data';    % 'data' or 'manual'

% -------------------------------------------------------------------------
% CELL COUNTS AND LABELS
% -------------------------------------------------------------------------
config.manual_cell_counts.ACC = 40;
config.manual_cell_counts.VLPFC = 40;
config.thalamus_cells_per_area = 10;

config.area_labels.ACC = {'ACC'};
config.area_labels.VLPFC = {'VLPFC'};
config.area_labels.thalamus_ACC = {'Thal-ACC', 'Thalamus-ACC'};
config.area_labels.thalamus_VLPFC = {'Thal-VLPFC', 'Thalamus-VLPFC'};
config.area_labels.thalamus_generic = {'Thalamus', 'Thal'};

% -------------------------------------------------------------------------
% SIMULATION SIZE
% Start small while validating network scaling and adapters.
% -------------------------------------------------------------------------
config.simulation.n_trials = 5;
config.simulation.n_time_bins = 5000;
config.simulation.dt = []; % [] = copy meta.dt from reference raster.
config.simulation.parallel_trials = false;
config.simulation.base_seed = 137;
config.simulation.save_cortex_only = true;

% Apply fixed-mode constraints before resolving output identity. In fixed
% mode this forces generated networks and prevents external data loading.
config = apply_model_size_mode(config);

% -------------------------------------------------------------------------
% OUTPUT IDENTITY
% Each condition gets a different session_idx so border file names cannot
% collide even if border naming does not include state.
% -------------------------------------------------------------------------
if strcmp(config.model_size_source, 'fixed')
    config.output.simulation_tag = sprintf( ...
        'fixed_ACC%d_VLPFC%d_ThalACC%d_ThalVLPFC%d_%s', ...
        config.fixed_model.cell_counts.ACC, ...
        config.fixed_model.cell_counts.VLPFC, ...
        config.fixed_model.cell_counts.thalamus_ACC, ...
        config.fixed_model.cell_counts.thalamus_VLPFC, ...
        config.tc_source);
    config.output.animal_name = 'SimFixed';
    config.output.base_session_idx = 1;
    config.output.resting_dur_threshold = 0;
else
    config.output.simulation_tag = sprintf('%s_s%d_%s_%s_%s', ...
        config.reference.animal_name, ...
        config.reference.session_idx, ...
        config.reference.state, ...
        config.network_source, ...
        config.tc_source);
    config.output.animal_name = sprintf('Sim%s', ...
        config.reference.animal_name);
    config.output.base_session_idx = config.reference.session_idx;
    config.output.resting_dur_threshold = ...
        config.reference.resting_dur_threshold;
end
config.output.injection = 'No injection';
config.output.prepost = 'Post';
config.output.align = 'None';
config.output.session_index_stride = 10;
config.output.overwrite = false;
config.output.save_project_raster = true;
config.output.save_project_border = true;
config.output.save_manifest = true;
config.output.save_network_in_manifest = true;
config.output.update_metadata = false;

%% Validate and run
validate_simulation_config(config);
progress_log('CONFIG', ...
    ['model_size_source=%s, network_source=%s, tc_source=%s, ' ...
     'kernel=%s, reg=%s.'], ...
    config.model_size_source, config.network_source, config.tc_source, ...
    resolve_effective_kernel_name(config), ...
    resolve_effective_reg_name(config));

[sim_meta, network, kernels, source_info] = ...
    build_simulation_network(root, config);
results = simulate_all_conditions(config, sim_meta, network, kernels);
output_manifest = save_simulation_output( ...
    root, config, sim_meta, network, kernels, source_info, results); %#ok<NASGU>

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
    config.reference.align = 'Last';
    config.reference.resting_dur_threshold = 15;
    config.reference.cortex_area = 'Cortex';
    config.reference.tc_area = 'Full';
    config.reference.kernel_name = 'DeltaPure';
    config.reference.kernel_num = [];
    config.reference.reg_name = 'L2=0_2';
    config.reference.epoch = 3000;
    config.reference.fold_idx = 0;
    config.reference.shuffle_idx = 0;

    config.model_size_source = 'match_data';
    config.fixed_model = struct();
    config.fixed_model.cell_counts = struct( ...
        'ACC', 40, 'VLPFC', 40, ...
        'thalamus_ACC', 10, 'thalamus_VLPFC', 10);
    config.fixed_model.dt = 0.001;
    config.fixed_model.kernel_name = 'FixedExpTriple';
    config.fixed_model.reg_name = 'Generated';
    config.fixed_model.kernel_len = 100;
    config.fixed_model.conn_time_constants_bins = [5, 20, 80];
    config.fixed_model.PS_time_constants_bins = 10;
    config.fixed_model.normalize_temporal_kernels = true;

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
    config.area_labels.thalamus_generic = {'Thalamus', 'Thal'};

    % Generated cortical network. Matrices use target rows and source cols:
    %   [ACC<-ACC,     ACC<-VLPFC
    %    VLPFC<-ACC,   VLPFC<-VLPFC]
    config.generate.cortex.connection_probability = [0.08, 0.04; 0.04, 0.08];
    config.generate.cortex.positive_probability = [0.75, 0.65; 0.65, 0.75];
    config.generate.cortex.weight_abs_mean = [1.5, 0.8; 0.8, 1.5];
    config.generate.cortex.weight_cv = 0.20;
    config.generate.cortex.kernel_scale = [1.0, 0.6, 0.3];
    config.generate.cortex.baseline_mean = [-4.0, -4.0];
    config.generate.cortex.baseline_sd = [0.4, 0.4];
    config.generate.cortex.post_spike_weight_mean = -2.0;
    config.generate.cortex.post_spike_weight_sd = 0.2;

    config.generate.thalamus.baseline_mean = -4.5;
    config.generate.thalamus.baseline_sd = 0.3;
    config.generate.thalamus.post_spike_weight_mean = -2.0;
    config.generate.thalamus.post_spike_weight_sd = 0.2;

    config.generate.tc.connection_probability = 0.10;
    config.generate.tc.weight_mean = 3.0;
    config.generate.tc.weight_cv = 0.20;
    config.generate.tc.kernel_scale = [1.0, 0.6, 0.3];

    % Simulated conditions.
    conditions = struct([]);
    conditions(1).name = 'FittedAsync';
    conditions(1).output_state = 'FittedAsync';
    conditions(1).enable_tc = true;
    conditions(1).sync_thalamus = false;
    conditions(1).sync_period = 100;
    conditions(1).sync_on_fraction = 1/8;
    conditions(1).sync_on_offset = 1.5;
    conditions(1).sync_off_offset = -0.6;

    conditions(2) = conditions(1);
    conditions(2).name = 'FittedSync';
    conditions(2).output_state = 'FittedSync';
    conditions(2).sync_thalamus = true;

    conditions(3) = conditions(1);
    conditions(3).name = 'FittedNoInput';
    conditions(3).output_state = 'FittedNoInput';
    conditions(3).enable_tc = false;
    conditions(3).sync_thalamus = false;

    config.conditions = conditions;

    config.simulation.n_trials = 5;
    config.simulation.n_time_bins = 5000;
    config.simulation.dt = [];
    config.simulation.parallel_trials = false;
    config.simulation.base_seed = 137;
    config.simulation.save_cortex_only = true;

    config.output = struct();
    config.output.simulation_tag = 'simulation';
    config.output.animal_name = 'Simulation';
    config.output.injection = 'No injection';
    config.output.prepost = 'Post';
    config.output.align = 'None';
    config.output.base_session_idx = 1;
    config.output.resting_dur_threshold = 0;
    config.output.session_index_stride = 10;
    config.output.overwrite = false;
    config.output.save_project_raster = true;
    config.output.save_project_border = true;
    config.output.save_manifest = true;
    config.output.save_network_in_manifest = true;
    config.output.update_metadata = false;
end

function config = apply_model_size_mode(config)
    switch config.model_size_source
        case 'match_data'
            % Keep user-selected data-matching options unchanged.

        case 'fixed'
            % Fixed mode must be self-contained. A single switch is enough
            % to guarantee that no raster, GLM, or kernel data are loaded.
            if ~strcmp(config.network_source, 'generate')
                progress_log('CONFIG', ...
                    'Fixed mode overrides network_source: %s -> generate.', ...
                    config.network_source);
            end
            config.network_source = 'generate';

            if ~strcmp(config.cell_count_source, 'manual')
                progress_log('CONFIG', ...
                    'Fixed mode ignores data cell counts and uses fixed presets.');
            end
            config.cell_count_source = 'manual';

            if strcmp(config.tc_source, 'load_pre')
                progress_log('CONFIG', ...
                    'Fixed mode overrides tc_source: load_pre -> generate.');
                config.tc_source = 'generate';
            end
            config.simulation.dt = config.fixed_model.dt;

        otherwise
            error('Unknown config.model_size_source: %s', ...
                config.model_size_source);
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
        if config.fixed_model.kernel_len < 1 || ...
                config.fixed_model.kernel_len ~= round(config.fixed_model.kernel_len)
            error('config.fixed_model.kernel_len must be a positive integer.');
        end
        if isempty(config.fixed_model.conn_time_constants_bins) || ...
                any(config.fixed_model.conn_time_constants_bins <= 0) || ...
                isempty(config.fixed_model.PS_time_constants_bins) || ...
                any(config.fixed_model.PS_time_constants_bins <= 0)
            error('Fixed temporal-kernel time constants must be positive.');
        end
        if ~strcmp(config.network_source, 'generate') || ...
                strcmp(config.tc_source, 'load_pre')
            error(['Fixed mode must use generated cortical connections and ' ...
                   'cannot use tc_source=''load_pre''.']);
        end
    end
    if ~isempty(config.reference.kernel_num) && ...
            (config.reference.kernel_num < 1 || ...
             config.reference.kernel_num ~= round(config.reference.kernel_num))
        error('config.reference.kernel_num must be empty or a positive integer.');
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
    if ~isempty(config.simulation.dt) && ...
            (~isscalar(config.simulation.dt) || config.simulation.dt <= 0)
        error('config.simulation.dt must be empty or a positive scalar in seconds.');
    end
    if config.output.base_session_idx < 0 || ...
            config.output.base_session_idx ~= round(config.output.base_session_idx)
        error('config.output.base_session_idx must be a non-negative integer.');
    end
    if config.output.session_index_stride < numel(config.conditions)
        error(['config.output.session_index_stride must be at least the ' ...
               'number of simulated conditions.']);
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

    if ~is_fixed && ~isempty(config.reference.kernel_num) && ...
            config.reference.kernel_num ~= K_conn
        error(['Configured kernel_num=%d but kernel metadata reports ' ...
               'n_conn_kernel=%d.'], ...
            config.reference.kernel_num, K_conn);
    end

    post_raster_info = [];
    loaded_post = [];
    loaded_pre = [];

    if is_fixed
        progress_log('NETWORK', ...
            'Fixed mode: skipping all raster, GLM, and kernel data references.');
        cell_info = build_cell_layout(config, [], []);
    else
        post_raster_meta = make_reference_meta(config, 'Post', ...
            config.reference.cortex_area);
        post_raster_info = load_reference_raster_cell_info( ...
            root, post_raster_meta);

        if strcmp(config.network_source, 'load')
            progress_log('NETWORK', 'Loading Post cortical GLM.');
            loaded_post = load_reference_glm_network( ...
                root, post_raster_meta, kernels);
            cortical_source = loaded_post;
        else
            cortical_source = post_raster_info;
        end
        cell_info = build_cell_layout( ...
            config, cortical_source, post_raster_info);
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
            pre_meta = make_reference_meta(config, 'Pre', ...
                config.reference.tc_area);
            progress_log('NETWORK', 'Loading Pre Full GLM for TC connections.');
            loaded_pre = load_reference_glm_network(root, pre_meta, kernels);
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
            dt = post_raster_info.meta.dt;
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
    sim_meta.cell_id = cell_info.cell_id;
    sim_meta.channel = cell_info.channel;
    sim_meta.N_total = cell_info.N_total;
    sim_meta.N_cortex = cell_info.N_cortex;
    sim_meta.N_thalamus = cell_info.N_thalamus;
    sim_meta.dt = dt;
    sim_meta.kernel_name = resolve_effective_kernel_name(config);
    sim_meta.n_conn_kernel = K_conn;
    sim_meta.n_PS_kernel = K_ps;
    sim_meta.reg_name = resolve_effective_reg_name(config);
    sim_meta.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');

    source_info = struct();
    source_info.uses_external_reference = ~is_fixed;
    source_info.post_raster_path = get_optional_struct_field( ...
        post_raster_info, 'source_path', '');
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

function cell_info = build_cell_layout(config, cortical_source, raster_source)
    is_fixed = strcmp(config.model_size_source, 'fixed');

    if is_fixed
        n_acc = config.fixed_model.cell_counts.ACC;
        n_vlpfc = config.fixed_model.cell_counts.VLPFC;
        n_thal_acc = config.fixed_model.cell_counts.thalamus_ACC;
        n_thal_vlpfc = config.fixed_model.cell_counts.thalamus_VLPFC;
        source_cell_id = [ ...
            compose('fixed_ACC_%d', (1:n_acc).'); ...
            compose('fixed_VLPFC_%d', (1:n_vlpfc).')];
        source_channel = nan(n_acc + n_vlpfc, 1);
        reference_raster_meta = struct();
    else
        if strcmp(config.cell_count_source, 'manual')
            n_acc = config.manual_cell_counts.ACC;
            n_vlpfc = config.manual_cell_counts.VLPFC;
            source_cell_id = compose('manual_cortex_%d', ...
                (1:(n_acc + n_vlpfc)).');
            source_channel = nan(n_acc + n_vlpfc, 1);
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
            source_channel = double(cortical_source.channel(source_idx));
        end
        n_thal_acc = config.thalamus_cells_per_area;
        n_thal_vlpfc = config.thalamus_cells_per_area;
        reference_raster_meta = raster_source.meta;
    end

    n_cortex = n_acc + n_vlpfc;
    n_thalamus = n_thal_acc + n_thal_vlpfc;
    n_total = n_cortex + n_thalamus;

    idx_acc = (1:n_acc).';
    idx_vlpfc = (n_acc + (1:n_vlpfc)).';
    idx_thal_acc = (n_cortex + (1:n_thal_acc)).';
    idx_thal_vlpfc = (n_cortex + n_thal_acc + ...
        (1:n_thal_vlpfc)).';

    cell_area = strings(n_total, 1);
    cell_area(idx_acc) = "ACC";
    cell_area(idx_vlpfc) = "VLPFC";
    cell_area(idx_thal_acc) = "Thal-ACC";
    cell_area(idx_thal_vlpfc) = "Thal-VLPFC";

    cell_id = strings(n_total, 1);
    cell_id(1:n_cortex) = source_cell_id(:);
    cell_id(idx_thal_acc) = compose('sim_thal_acc_%d', ...
        (1:n_thal_acc).');
    cell_id(idx_thal_vlpfc) = compose('sim_thal_vlpfc_%d', ...
        (1:n_thal_vlpfc).');

    channel = nan(n_total, 1);
    channel(1:n_cortex) = source_channel(:);

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
    cell_info.channel = channel;
    cell_info.masks = masks;
    cell_info.reference_raster_meta = reference_raster_meta;

    progress_log('CELLS', ...
        'mode=%s, ACC=%d, VLPFC=%d, Thal-ACC=%d, Thal-VLPFC=%d.', ...
        config.model_size_source, n_acc, n_vlpfc, ...
        n_thal_acc, n_thal_vlpfc);
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
    meta.kernel_name = config.reference.kernel_name;
    meta.reg_name = config.reference.reg_name;
    meta.epoch = config.reference.epoch;
    meta.fold_idx = config.reference.fold_idx;
    meta.shuffle_idx = config.reference.shuffle_idx;
end

%% =========================================================================
%  README-BASED INPUT ADAPTERS
%  =========================================================================
function raster_info = load_reference_raster_cell_info(root, meta)
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
    if ~isfield(loaded.meta, 'N') || ~isfield(loaded.meta, 'dt') || ...
            ~isfield(loaded.data, 'cell_area') || ...
            ~isfield(loaded.data, 'cell_id')
        error(['Raster file is missing README-required fields ' ...
               'meta.N/meta.dt/data.cell_area/data.cell_id: %s'], raster_path);
    end

    N = double(loaded.meta.N);
    cell_area = string(loaded.data.cell_area(:));
    cell_id = string(loaded.data.cell_id(:));
    if numel(cell_area) ~= N || numel(cell_id) ~= N
        error('Raster cell metadata lengths do not match meta.N: %s', raster_path);
    end

    if isfield(loaded.data, 'channel') && numel(loaded.data.channel) == N
        channel = double(loaded.data.channel(:));
    else
        channel = nan(N, 1);
        progress_log('LOAD-RASTER', ...
            'data.channel missing or wrong length; using NaN channels.');
    end

    raster_info = struct();
    raster_info.meta = loaded.meta;
    raster_info.cell_area = cell_area;
    raster_info.cell_id = cell_id;
    raster_info.channel = channel;
    raster_info.N = N;
    raster_info.source_path = raster_path;
end

function kernels = load_kernel_definition(root, config)
    if strcmp(config.model_size_source, 'fixed')
        kernels = build_fixed_kernel_definition(config);
        return;
    end

    kernel_meta = struct();
    kernel_meta.kernel_name = config.reference.kernel_name;
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

    n_conn = double(loaded.meta.n_conn_kernel);
    n_ps = double(loaded.meta.n_PS_kernel);
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
end

function kernels = build_fixed_kernel_definition(config)
    fixed = config.fixed_model;
    kernel_len = fixed.kernel_len;
    t = 0:(kernel_len - 1);

    conn_taus = double(fixed.conn_time_constants_bins(:).');
    ps_taus = double(fixed.PS_time_constants_bins(:).');
    conn_values = cell(1, numel(conn_taus));
    PS_values = cell(1, numel(ps_taus));

    for k = 1:numel(conn_taus)
        value = exp(-t / conn_taus(k));
        if fixed.normalize_temporal_kernels
            value = value / sum(value);
        end
        conn_values{k} = value(:, end:-1:1);
    end
    for k = 1:numel(ps_taus)
        value = exp(-t / ps_taus(k));
        if fixed.normalize_temporal_kernels
            value = value / sum(value);
        end
        PS_values{k} = value(:, end:-1:1);
    end

    kernels = struct();
    kernels.meta = struct();
    kernels.meta.kernel_name = fixed.kernel_name;
    kernels.meta.n_conn_kernel = numel(conn_values);
    kernels.meta.n_PS_kernel = numel(PS_values);
    kernels.meta.kernel_len = kernel_len;
    kernels.conn_values = conn_values;
    kernels.PS_values = PS_values;
    kernels.conn_lengths = cellfun(@numel, conn_values);
    kernels.PS_lengths = cellfun(@numel, PS_values);
    kernels.n_conn_kernel = numel(conn_values);
    kernels.n_PS_kernel = numel(PS_values);
    kernels.kernel_len = kernel_len;
    kernels.source_path = '';

    progress_log('KERNEL', ...
        'Built fixed kernels in memory: connection=%d, PS=%d, length=%d.', ...
        kernels.n_conn_kernel, kernels.n_PS_kernel, kernel_len);
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

function loaded_network = load_reference_glm_network(root, meta, kernels)
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

    raster_info = load_reference_raster_cell_info(root, meta);
    raster_indices = resolve_glm_to_raster_indices( ...
        loaded.meta, loaded.data, raster_info.N, n_target);

    loaded_network = struct();
    loaded_network.meta = loaded.meta;
    loaded_network.N = n_target;
    loaded_network.h = h;
    loaded_network.P = P;
    loaded_network.J = J;
    loaded_network.cell_area = raster_info.cell_area(raster_indices);
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
    % MANUAL CONFIRMATION:
    % This function assumes the Pre Full GLM contains separately labelled
    % Thal-ACC and Thal-VLPFC cells. README only guarantees the generic area
    % label "Thalamus"; it does not define how thalamic cells are assigned to
    % their cortical projection target.

    areas = string(loaded_pre.cell_area(:));
    src_acc = find(ismember(areas, string(config.area_labels.ACC)));
    src_vlpfc = find(ismember(areas, string(config.area_labels.VLPFC)));
    src_thal_acc = find(ismember(areas, ...
        string(config.area_labels.thalamus_ACC)));
    src_thal_vlpfc = find(ismember(areas, ...
        string(config.area_labels.thalamus_VLPFC)));

    if isempty(src_thal_acc) || isempty(src_thal_vlpfc)
        generic_thalamus = find(ismember(areas, ...
            string(config.area_labels.thalamus_generic)));
        if ~isempty(generic_thalamus)
            error(['Pre GLM contains generic Thalamus cells but README does ' ...
                   'not specify how to split them into Thal-ACC and ' ...
                   'Thal-VLPFC. Implement this mapping in ' ...
                   'insert_loaded_tc_network().']);
        end
        error(['tc_source=''load_pre'' requested, but no identifiable ' ...
               'Thal-ACC/Thal-VLPFC cells were found.']);
    end

    dst_acc = find(cell_info.masks.ACC);
    dst_vlpfc = find(cell_info.masks.VLPFC);
    dst_thal_acc = find(cell_info.masks.thalamus_ACC);
    dst_thal_vlpfc = find(cell_info.masks.thalamus_VLPFC);

    assert_matching_count(src_acc, dst_acc, 'ACC');
    assert_matching_count(src_vlpfc, dst_vlpfc, 'VLPFC');
    assert_matching_count(src_thal_acc, dst_thal_acc, 'Thal-ACC');
    assert_matching_count(src_thal_vlpfc, dst_thal_vlpfc, 'Thal-VLPFC');

    for k = 1:size(J_out, 3)
        J_out(dst_acc, dst_thal_acc, k) = ...
            loaded_pre.J(src_acc, src_thal_acc, k);
        J_out(dst_vlpfc, dst_thal_vlpfc, k) = ...
            loaded_pre.J(src_vlpfc, src_thal_vlpfc, k);
    end

    progress_log('NETWORK', 'Inserted loaded Pre thalamocortical J.');
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

    kernel_scale = expand_kernel_scale( ...
        config.generate.cortex.kernel_scale, size(J, 3), ...
        'generate.cortex.kernel_scale');

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
    kernel_scale = expand_kernel_scale( ...
        config.generate.tc.kernel_scale, size(J, 3), ...
        'generate.tc.kernel_scale');

    J = generate_one_tc_projection(J, ...
        find(cell_info.masks.ACC), ...
        find(cell_info.masks.thalamus_ACC), ...
        config, kernel_scale);
    J = generate_one_tc_projection(J, ...
        find(cell_info.masks.VLPFC), ...
        find(cell_info.masks.thalamus_VLPFC), ...
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

function scale = expand_kernel_scale(value, K, field_name)
    value = double(value(:).');
    if isscalar(value)
        scale = repmat(value, 1, K);
    elseif numel(value) >= K
        scale = value(1:K);
    else
        error('%s has %d values but %d are required.', ...
            field_name, numel(value), K);
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
        progress_log('SIMULATION', '[%d/%d] START %s.', ...
            condition_i, n_conditions, condition.name);

        condition_network = network;
        if ~condition.enable_tc
            condition_network.J(:, condition_network.masks.thalamus, :) = 0;
        end

        results(condition_i) = simulate_condition( ...
            config, sim_meta, condition_network, kernels, ...
            condition, condition_i);

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

    saved_cell_mask = true(N, 1);
    if config.simulation.save_cortex_only
        saved_cell_mask = network.masks.cortex;
    end

    for trial_i = 1:n_trials
        trial_rasters{trial_i} = uint8( ...
            trial_rasters{trial_i}(saved_cell_mask, :));
        trial_firing_rates{trial_i} = ...
            trial_firing_rates{trial_i}(saved_cell_mask).';
    end

    result = struct();
    result.condition = condition;
    result.N = sum(saved_cell_mask);
    result.dt = sim_meta.dt;
    result.cell_area = sim_meta.cell_area(saved_cell_mask).';
    result.cell_id = sim_meta.cell_id(saved_cell_mask).';
    result.channel = sim_meta.channel(saved_cell_mask).';
    result.trial_num = n_trials;
    result.trial_len = repmat(n_time_bins, 1, n_trials);
    result.rasters = trial_rasters;
    result.firing_rates = trial_firing_rates;
    result.saved_cell_mask = saved_cell_mask;
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
function output_manifest = save_simulation_output( ...
    root, config, sim_meta, network, kernels, source_info, results)

    raster_folder = fullfile(root, 'Data', 'Working', 'raster');
    border_folder = fullfile(root, 'Data', 'Working', 'border');
    manifest_folder = fullfile(root, 'Data', 'Working', ...
        'simulation_v2', sanitize_for_path(config.output.simulation_tag));
    check_path(raster_folder);
    check_path(border_folder);
    check_path(manifest_folder);

    output_manifest = struct();
    output_manifest.raster_paths = strings(1, numel(results));
    output_manifest.border_paths = strings(1, numel(results));

    for condition_i = 1:numel(results)
        result = results(condition_i);
        raster_meta = build_output_raster_meta( ...
            config, sim_meta, result, condition_i);
        raster_data = build_output_raster_data(result);

        if config.output.save_project_raster
            raster_meta.file_name = generate_filename('raster', raster_meta);
            raster_path = fullfile(raster_folder, raster_meta.file_name);
            assert_output_write_allowed(raster_path, config.output.overwrite);
            meta = raster_meta; %#ok<NASGU>
            data = raster_data; %#ok<NASGU>
            progress_log('SAVE-RASTER', 'Writing: %s', raster_path);
            save(raster_path, 'meta', 'data', '-v7.3');
            output_manifest.raster_paths(condition_i) = string(raster_path);
        end

        if config.output.save_project_border
            [border_meta, border_data] = build_output_border( ...
                raster_meta, result);
            border_meta.file_name = generate_filename('border', border_meta);
            border_path = fullfile(border_folder, border_meta.file_name);
            assert_output_write_allowed(border_path, config.output.overwrite);
            meta = border_meta; %#ok<NASGU>
            data = border_data; %#ok<NASGU>
            progress_log('SAVE-BORDER', 'Writing: %s', border_path);
            save(border_path, 'meta', 'data');
            output_manifest.border_paths(condition_i) = string(border_path);
        end
    end

    if config.output.save_manifest
        manifest_path = fullfile(manifest_folder, sprintf( ...
            'simulation_manifest_%s.mat', ...
            sanitize_for_path(config.output.simulation_tag)));
        assert_output_write_allowed(manifest_path, config.output.overwrite);

        meta = struct();
        meta.simulation_tag = config.output.simulation_tag;
        meta.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        meta.model_size_source = config.model_size_source;
        meta.network_source = config.network_source;
        meta.tc_source = config.tc_source;
        meta.kernel_name = sim_meta.kernel_name;
        meta.reg_name = sim_meta.reg_name;

        data = struct();
        data.config = config;
        data.sim_meta = sim_meta;
        data.kernels = kernels;
        data.source_info = source_info;
        data.output_manifest = output_manifest;
        if config.output.save_network_in_manifest
            data.network = network;
        else
            data.network = rmfield(network, {'h', 'P', 'J'});
        end

        progress_log('SAVE-MANIFEST', 'Writing: %s', manifest_path);
        save(manifest_path, 'meta', 'data', '-v7.3');
        output_manifest.manifest_path = string(manifest_path);
    else
        output_manifest.manifest_path = "";
    end

    if config.output.update_metadata
        register_simulation_metadata(root, output_manifest);
    end
end

function meta = build_output_raster_meta(config, sim_meta, result, condition_i)
    condition = result.condition;
    output_session_idx = config.output.base_session_idx * ...
        config.output.session_index_stride + condition_i;

    meta = struct();
    meta.animal_name = config.output.animal_name;
    meta.injection = config.output.injection;
    meta.prepost = config.output.prepost;
    meta.state = condition.output_state;
    if config.simulation.save_cortex_only
        meta.area = 'Cortex';
    else
        meta.area = 'Full';
    end
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

    % MANUAL CONFIRMATION:
    % README defines data.spikes as one vector per trial but does not define
    % how neuron identity is encoded. Empty cells are written to preserve the
    % field without inventing an incompatible representation.
    data.spikes = repmat({[]}, 1, result.trial_num);

    data.trial_len = reshape(result.trial_len, 1, []);
    data.cell_id = reshape(string(result.cell_id), 1, []);
    data.cell_area = reshape(string(result.cell_area), 1, []);
    data.channel = reshape(double(result.channel), 1, []);
    data.cuetype = repmat({[]}, 1, result.trial_num);
    data.firing_rates = result.firing_rates;
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

function register_simulation_metadata(~, ~)
    % MANUAL CONFIRMATION:
    % README describes metadata.mat conceptually but does not define the
    % write/update API, duplicate handling, locking, or whether metadata
    % stores relative paths or full paths. Do not modify metadata.mat without
    % the project's registration helper.
    error(['config.output.update_metadata=true, but metadata registration ' ...
           'is intentionally not implemented. Supply the current project ' ...
           'metadata registration function.']);
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

function name = resolve_effective_kernel_name(config)
    if strcmp(config.model_size_source, 'fixed')
        name = config.fixed_model.kernel_name;
    else
        name = config.reference.kernel_name;
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
