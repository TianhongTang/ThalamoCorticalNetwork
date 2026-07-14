%% Supplement Figure: Open/Close and Pre/Post compound-pattern transition matrices
% This script pools multiple kernel_name models into configurable kernel setups.
% Each setup produces two 2x7 figures for every injection x align combination.
%
% Export behavior is controlled independently by export_jpg and export_pdf.
% Default: export JPG only.
%
% Metadata behavior:
%   1. On the first run, load the full GLM metadata, retain epoch = 3000 and
%      fixed analysis filters, verify complete Pre/Post x Open/Close sessions,
%      and save a compact session index.
%   2. On later runs, load only the compact index and do not load the full
%      metadata file.
%
% Figure versions for each pooled kernel setup:
%   1. Pre -> Post transition of Open/Close compound patterns.
%   2. Open vs Close comparison of Pre/Post compound patterns.
%
% Output organization:
%   Figures/Paper/Fig3_Supp/<hyperparameter-folder>/
%
% Layout for every saved figure:
%   Row 1 = 9x9 compound-category table including non-significant category.
%   Row 2 = 4x4 compound-category table excluding any non-significant member.
%
% Columns:
%   1. Raw count.
%   2. Normalized by row.
%   3. Normalized by column.
%   4. Transition difference = J_ij - J_ji.
%   5. Normalized difference = (J_ij - J_ji) / (J_ij + J_ji).
%   6. O/E ratio under independence baseline.
%   7. Standardized residual under independence baseline.
%
% Compound category labels use +, o, -.
%   For Open/Close patterns: first symbol = RestOpen, second symbol = RestClose.
%   For Pre/Post patterns:   first symbol = Pre,      second symbol = Post.

clear;
run_tic = tic;
progress_log('SCRIPT', 'Started.');

%% Get root folder
code_depth = 3;
script_path = mfilename('fullpath');
root = script_path;
for i = 1:code_depth
    root = fileparts(root);
end
addpath(fileparts(script_path));
addpath(fullfile(root, 'Code', 'Utils'));
progress_log('SCRIPT', 'Project root: %s', root);

%% Kernel setup and hyperparameter configurations
% Edit this block to change the number of pooled kernel setups or their members.
% Only kernel_name is pooled. reg_name, injection, align, and resting-duration
% threshold remain separate hyperparameter dimensions.
%
% The listed Long* kernels are separate one-kernel GLM models, so kernel block 1
% is loaded from each model file before pooling its connections.
kernel_setups = struct([]);

kernel_setups(1).name = 'fast_kernels';
kernel_setups(1).kernel_names = { ...
    'LongExp5',       'LongExp10',       'LongExp20', ...
    'LongGaussC5',    'LongGaussC10',    'LongGaussC20', ...
    'LongGDeriv5',    'LongGDeriv10',    'LongGDeriv20', ...
    'LongStepB0_5',   'LongStepB5_10',   'LongStepB10_20' ...
};

kernel_setups(2).name = 'slow_kernels';
kernel_setups(2).kernel_names = { ...
    'LongExp40',      'LongExp80',      'LongExp160', ...
    'LongExp320',     'LongExp640',     'LongExp1000', ...
    'LongGaussC40',   'LongGaussC80',   'LongGaussC160', ...
    'LongGaussC320',  'LongGaussC640',  'LongGaussC1000', ...
    'LongGDeriv40',   'LongGDeriv80',   'LongGDeriv160', ...
    'LongGDeriv320',  'LongGDeriv640',  'LongGDeriv1000', ...
    'LongStepB20_40', 'LongStepB40_80', 'LongStepB80_160', ...
    'LongStepB160_320', 'LongStepB320_640', ...
    'LongStepB640_1000', 'LongStepB1000_3000' ...
};

validate_kernel_setups(kernel_setups);
progress_log('CONFIG', 'Generating grouped hyperparameter configurations.');
hyperparam_configs = build_grouped_hyperparam_configs(kernel_setups);
progress_log('CONFIG', 'Generated %d configurations from %d pooled kernel setups.', ...
    numel(hyperparam_configs), numel(kernel_setups));

%% Plot and analysis parameters
err_multi = 1; % threshold for significant J, in multiples of the GLM error estimate.
category_labels_9 = {'++', '+o', '+-', 'o+', 'oo', 'o-', '-+', '-o', '--'};
category_labels_4 = {'++', '+-', '-+', '--'};

n_row = 2;
n_col = 7;
figure_visible = 'off';

skip_failed_sessions = false;
max_sessions_to_include = inf; % Limits kernel-session records; set smaller only for debugging.

% Independent export switches.
export_jpg = true;
export_pdf = false;
export_resolution = 300; % dpi.

% Metadata index controls.
metadata_index_version = 2;
force_rebuild_metadata_index = false;
metadata_index_filename = fullfile(root, 'Data', 'Working', 'metadata_index', ...
    'GLM_epoch3000_complete_session_index.mat');

base_params = struct();
base_params.err_multi = err_multi;
base_params.category_labels_9 = category_labels_9;
base_params.category_labels_4 = category_labels_4;
base_params.n_row = n_row;
base_params.n_col = n_col;
base_params.figure_visible = figure_visible;
base_params.skip_failed_sessions = skip_failed_sessions;
base_params.max_sessions_to_include = max_sessions_to_include;
base_params.export_jpg = export_jpg;
base_params.export_pdf = export_pdf;
base_params.export_resolution = export_resolution;

progress_log('CONFIG', 'Export settings: JPG=%d, PDF=%d, resolution=%d dpi.', ...
    export_jpg, export_pdf, export_resolution);

%% Load or build compact metadata index
progress_log('INDEX', 'Index file: %s', metadata_index_filename);
[session_index, metadata_index_info] = load_or_build_metadata_session_index( ...
    root, metadata_index_filename, metadata_index_version, force_rebuild_metadata_index);
progress_log('INDEX', 'Index ready: %d complete sessions; created %s.', ...
    numel(session_index), metadata_index_info.created_at);

fig3_root_folder = fullfile(root, 'Figures', 'Paper', 'Fig3_Supp');
check_path(fig3_root_folder);
run_issue_log = empty_run_issue_log();
run_issue_log_filename = fullfile(fig3_root_folder, 'run_issue_log.txt');

total_rendered_figure_count = 0;
total_valid_session_count = 0;
total_failed_session_count = 0;

%% Run each hyperparameter configuration
for hp_i = 1:numel(hyperparam_configs)
    hp_tic = tic;
    hp_title_for_log = sprintf('hyperparameter set %d/%d', hp_i, numel(hyperparam_configs));

    try
        hp = fill_default_hyperparams(hyperparam_configs(hp_i));
        hp_title_for_log = make_hyperparam_title(hp);
        model_kernel_idx = 1;

        params = base_params;
        params.hyperparams = hp;
        params.hyperparam_title = hp_title_for_log;
        params.output_folder = fullfile(fig3_root_folder, make_hyperparam_folder_name(hp));

        progress_log('HP', '[%d/%d] START: %s', hp_i, numel(hyperparam_configs), params.hyperparam_title);
        progress_log('HP', '[%d/%d][1/4] Selecting sessions from compact index.', hp_i, numel(hyperparam_configs));

        figure_configs = build_pattern_figure_configs();

        %% Select indexed sessions for this hyperparameter set
        meta_array = select_sessions_from_index(session_index, hp);
        report_kernel_setup_coverage(meta_array, hp);
        if isfinite(params.max_sessions_to_include)
            meta_array = meta_array(1:min(numel(meta_array), params.max_sessions_to_include));
        end

        progress_log('HP', ...
            '[%d/%d][1/4] Matched %d kernel-session records. Output: %s', ...
            hp_i, numel(hyperparam_configs), numel(meta_array), params.output_folder);
        if isempty(meta_array)
            message = 'No indexed kernel-session records selected. Skipping this hyperparameter set.';
            warning('%s Hyperparameter set %d: %s', message, hp_i, params.hyperparam_title);
            run_issue_log = append_run_issue(run_issue_log, hp_i, params.hyperparam_title, 'NoMetadata', message);
            write_run_issue_log(run_issue_log_filename, run_issue_log);
            progress_log('HP', '[%d/%d] SKIPPED: no indexed kernel-session records.', hp_i, numel(hyperparam_configs));
            continue;
        end

        %% Initialize pooled storage
        progress_log('HP', '[%d/%d][2/4] Initializing pooled count matrices for %d figures.', ...
            hp_i, numel(hyperparam_configs), numel(figure_configs));
        n_fig = numel(figure_configs);
        pooled = struct([]);
        for fig_i = 1:n_fig
            pooled(fig_i).counts9 = zeros(numel(category_labels_9), numel(category_labels_9));
            pooled(fig_i).counts4 = zeros(numel(category_labels_4), numel(category_labels_4));
            pooled(fig_i).valid_session_count = 0;
        end

        valid_session_count = 0;
        failed_session_count = 0;

        %% Pool data across all selected sessions for this hyperparameter set
        progress_log('HP', '[%d/%d][3/4] Pooling %d kernel-session records.', ...
            hp_i, numel(hyperparam_configs), numel(meta_array));

        for session_i = 1:numel(meta_array)
            session_tic = tic;
            meta = meta_array(session_i);
            session_label = make_session_label(meta);
            progress_log('SESSION', '[HP %d/%d][%d/%d] START: %s', ...
                hp_i, numel(hyperparam_configs), session_i, numel(meta_array), session_label);

            try
                loaded_states = load_all_required_states(root, meta, model_kernel_idx);

                % Accumulate into temporary session storage first. This prevents
                % partial session results from entering pooled counts if a later
                % figure fails and skip_failed_sessions is enabled.
                session_counts9 = cell(n_fig, 1);
                session_counts4 = cell(n_fig, 1);

                for fig_i = 1:n_fig
                    cfg = figure_configs(fig_i);
                    k = model_kernel_idx;
                    progress_log('ANALYSIS', ...
                        '[HP %d/%d][Record %d/%d][Figure %d/%d] %s, kernel_name=%s.', ...
                        hp_i, numel(hyperparam_configs), session_i, numel(meta_array), ...
                        fig_i, n_fig, cfg.analysis_type, char(string(meta.kernel_name)));

                    pre_open_state   = loaded_states.(state_key('Pre',  'RestOpen',  k));
                    pre_close_state  = loaded_states.(state_key('Pre',  'RestClose', k));
                    post_open_state  = loaded_states.(state_key('Post', 'RestOpen',  k));
                    post_close_state = loaded_states.(state_key('Post', 'RestClose', k));

                    switch cfg.analysis_type
                        case 'pre_to_post_oc_pattern'
                            [counts9, counts4] = make_pre_to_post_oc_pattern_transition_counts( ...
                                pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi);

                        case 'open_vs_close_prepost_pattern'
                            [counts9, counts4] = make_open_vs_close_prepost_pattern_counts( ...
                                pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi);

                        otherwise
                            error('Unknown analysis_type: %s', cfg.analysis_type);
                    end

                    session_counts9{fig_i} = counts9;
                    session_counts4{fig_i} = counts4;
                end

                for fig_i = 1:n_fig
                    pooled(fig_i).counts9 = pooled(fig_i).counts9 + session_counts9{fig_i};
                    pooled(fig_i).counts4 = pooled(fig_i).counts4 + session_counts4{fig_i};
                end

                valid_session_count = valid_session_count + 1;
                progress_log('SESSION', '[HP %d/%d][%d/%d] DONE in %.1f s: %s', ...
                    hp_i, numel(hyperparam_configs), session_i, numel(meta_array), ...
                    toc(session_tic), session_label);
            catch ME_session
                failed_session_count = failed_session_count + 1;
                progress_log('SESSION', '[HP %d/%d][%d/%d] FAILED after %.1f s: %s | %s', ...
                    hp_i, numel(hyperparam_configs), session_i, numel(meta_array), ...
                    toc(session_tic), session_label, ME_session.message);
                if params.skip_failed_sessions
                    warning('Skipping %s because loading/processing failed: %s', session_label, ME_session.message);
                    continue;
                else
                    rethrow(ME_session);
                end
            end
        end

        total_valid_session_count = total_valid_session_count + valid_session_count;
        total_failed_session_count = total_failed_session_count + failed_session_count;

        if valid_session_count == 0
            message = sprintf('No valid kernel-session records were loaded. Failed sessions: %d. Skipping rendering for this hyperparameter set.', failed_session_count);
            warning('%s Hyperparameter set %d: %s', message, hp_i, params.hyperparam_title);
            run_issue_log = append_run_issue(run_issue_log, hp_i, params.hyperparam_title, 'NoValidSessions', message);
            write_run_issue_log(run_issue_log_filename, run_issue_log);
            progress_log('HP', '[%d/%d] SKIPPED: no valid kernel-session records.', hp_i, numel(hyperparam_configs));
            continue;
        end

        progress_log('HP', '[%d/%d][3/4] Pooling complete. Valid=%d, failed=%d.', ...
            hp_i, numel(hyperparam_configs), valid_session_count, failed_session_count);

        for fig_i = 1:n_fig
            pooled(fig_i).valid_session_count = valid_session_count;
        end

        %% Render and save figures for this hyperparameter set
        progress_log('HP', '[%d/%d][4/4] Rendering %d figures.', ...
            hp_i, numel(hyperparam_configs), n_fig);
        for fig_i = 1:n_fig
            progress_log('RENDER', '[HP %d/%d][Figure %d/%d] START: %s', ...
                hp_i, numel(hyperparam_configs), fig_i, n_fig, figure_configs(fig_i).output_stub);
            render_pattern_transition_figure(root, figure_configs(fig_i), pooled(fig_i), params);
            total_rendered_figure_count = total_rendered_figure_count + 1;
            progress_log('RENDER', '[HP %d/%d][Figure %d/%d] DONE.', ...
                hp_i, numel(hyperparam_configs), fig_i, n_fig);
        end

        progress_log('HP', '[%d/%d] DONE in %.1f s. Valid kernel-session records=%d, failed=%d, figures=%d.', ...
            hp_i, numel(hyperparam_configs), toc(hp_tic), valid_session_count, failed_session_count, n_fig);

    catch ME_hp
        message = sprintf('%s: %s', ME_hp.identifier, ME_hp.message);
        warning('Skipping hyperparameter set %d/%d after error: %s\n%s', ...
            hp_i, numel(hyperparam_configs), hp_title_for_log, message);
        run_issue_log = append_run_issue(run_issue_log, hp_i, hp_title_for_log, 'Error', message);
        write_run_issue_log(run_issue_log_filename, run_issue_log);
        progress_log('HP', '[%d/%d] FAILED after %.1f s: %s', ...
            hp_i, numel(hyperparam_configs), toc(hp_tic), message);
        continue;
    end
end

if total_rendered_figure_count == 0
    warning('No figures were rendered. Check hyperparameter values, metadata coverage, and run_issue_log.txt.');
else
    progress_log('SUMMARY', 'Rendered %d figures across all hyperparameter sets.', total_rendered_figure_count);
end

if ~isempty(run_issue_log)
    write_run_issue_log(run_issue_log_filename, run_issue_log);
    progress_log('SUMMARY', 'Run issue log written to: %s', run_issue_log_filename);
end

progress_log('SUMMARY', 'Finished in %.1f s. Total valid kernel-session records=%d, failed records=%d.', ...
    toc(run_tic), total_valid_session_count, total_failed_session_count);


function hyperparam_configs = build_grouped_hyperparam_configs(kernel_setups)
    reg_name = "L2=0_2";
    resting_dur_threshold = 15;
    injections = {"Muscimol", "Saline"};
    aligns = {"Last", "Longest"};

    hyperparam_configs = struct([]);
    cfg_i = 0;
    for setup_i = 1:numel(kernel_setups)
        for inj_i = 1:numel(injections)
            for align_i = 1:numel(aligns)
                cfg_i = cfg_i + 1;
                hyperparam_configs(cfg_i).kernel_setup_name = ...
                    kernel_setups(setup_i).name; %#ok<AGROW>
                hyperparam_configs(cfg_i).kernel_names = ...
                    kernel_setups(setup_i).kernel_names;
                hyperparam_configs(cfg_i).reg_name = reg_name;
                hyperparam_configs(cfg_i).resting_dur_threshold = ...
                    resting_dur_threshold;
                hyperparam_configs(cfg_i).injection = injections{inj_i};
                hyperparam_configs(cfg_i).align = aligns{align_i};
            end
        end
    end

    progress_log('CONFIG', ...
        'Generated %d configurations: %d setups x %d injections x %d aligns.', ...
        numel(hyperparam_configs), numel(kernel_setups), ...
        numel(injections), numel(aligns));
end

function validate_kernel_setups(kernel_setups)
    if isempty(kernel_setups) || ~isstruct(kernel_setups)
        error('kernel_setups must be a nonempty struct array.');
    end

    setup_names = strings(numel(kernel_setups), 1);
    all_kernel_names = strings(0, 1);
    for setup_i = 1:numel(kernel_setups)
        if ~isfield(kernel_setups, 'name') || ...
                ~isfield(kernel_setups, 'kernel_names')
            error('Every kernel setup requires name and kernel_names fields.');
        end

        setup_name = string(kernel_setups(setup_i).name);
        members = string(kernel_setups(setup_i).kernel_names(:));
        if ~isscalar(setup_name) || ismissing(setup_name) || strlength(setup_name) == 0
            error('Kernel setup %d has an invalid name.', setup_i);
        end
        if isempty(members) || any(ismissing(members)) || any(strlength(members) == 0)
            error('Kernel setup %s has an empty or invalid member.', char(setup_name));
        end
        if numel(unique(members)) ~= numel(members)
            error('Kernel setup %s contains duplicate kernel names.', char(setup_name));
        end

        setup_names(setup_i) = setup_name;
        all_kernel_names = [all_kernel_names; members]; %#ok<AGROW>
        progress_log('CONFIG', 'Kernel setup %s: %d members: %s', ...
            char(setup_name), numel(members), strjoin(cellstr(members), ', '));
    end

    if numel(unique(setup_names)) ~= numel(setup_names)
        error('kernel_setups contains duplicate setup names.');
    end
    if numel(unique(all_kernel_names)) ~= numel(all_kernel_names)
        warning(['One or more kernel names occur in multiple setups. Their connections ' ...
            'will be pooled independently into every setup containing them.']);
    end
end

function [session_index, index_info] = load_or_build_metadata_session_index( ...
    root, index_filename, expected_version, force_rebuild)

    index_folder = fileparts(index_filename);
    check_path(index_folder);

    if isfile(index_filename) && ~force_rebuild
        load_tic = tic;
        progress_log('INDEX', 'Existing index detected. Loading compact index only.');
        loaded = load(index_filename, 'session_index', 'index_info');

        if ~isfield(loaded, 'session_index') || ~isfield(loaded, 'index_info')
            error(['Metadata index is missing session_index or index_info. ' ...
                   'Delete the index file or set force_rebuild_metadata_index=true.']);
        end

        needs_resave = false;
        if istable(loaded.session_index)
            progress_log('INDEX', 'Migrating legacy compact table index to struct format.');
            session_index = table2struct(loaded.session_index);
            needs_resave = true;
        elseif isstruct(loaded.session_index)
            session_index = loaded.session_index;
        else
            error('Metadata index session_index must be a struct array or a legacy table.');
        end

        index_info = loaded.index_info;
        if ~isfield(index_info, 'version') || index_info.version ~= expected_version
            if isfield(index_info, 'version') && index_info.version > expected_version
                error('Metadata index version %d is newer than supported version %d.', ...
                    index_info.version, expected_version);
            end
            progress_log('INDEX', 'Updating compact index version metadata to %d.', expected_version);
            index_info.version = expected_version;
            index_info.format = 'struct';
            index_info.migrated_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
            needs_resave = true;
        end
        if ~isfield(index_info, 'created_at')
            index_info.created_at = 'unknown';
        end

        if needs_resave
            progress_log('INDEX', 'Saving migrated struct index: %s', index_filename);
            save(index_filename, 'session_index', 'index_info');
        end

        progress_log('INDEX', 'Loaded %d indexed sessions in %.1f s without loading full metadata.', ...
            numel(session_index), toc(load_tic));
        return;
    end

    if force_rebuild && isfile(index_filename)
        progress_log('INDEX', 'Forced rebuild requested. Existing index will be replaced.');
    else
        progress_log('INDEX', 'No index found. Full metadata will be loaded once to build it.');
    end

    build_tic = tic;
    [session_index, index_info] = build_metadata_session_index(root, expected_version);

    progress_log('INDEX', 'Saving compact struct index: %s', index_filename);
    save(index_filename, 'session_index', 'index_info');
    progress_log('INDEX', 'Index saved; total build time %.1f s.', toc(build_tic));
end

function [session_index, index_info] = build_metadata_session_index(root, index_version)
    target_epoch = 3000;
    extraction_chunk_size = 100000;

    progress_log('INDEX-BUILD', '[1/7] Loading full metadata as struct. This occurs only when building the index.');
    load_tic = tic;
    loaded_meta = load_meta(root, 'struct');
    if ~isfield(loaded_meta, 'GLM') || ~isstruct(loaded_meta.GLM)
        error('load_meta(root, ''struct'') did not return a GLM struct array.');
    end
    mt = loaded_meta.GLM;
    clear loaded_meta;
    source_row_count = numel(mt);
    progress_log('INDEX-BUILD', '[1/7] Full metadata loaded: %d struct records in %.1f s.', ...
        source_row_count, toc(load_tic));

    required_fields = {'epoch', 'area', 'fold_idx', 'shuffle_idx', ...
        'prepost', 'state', 'animal_name', 'injection', 'align', ...
        'session_idx', 'resting_dur_threshold', 'kernel_name', 'reg_name'};
    missing_fields = required_fields(~isfield(mt, required_fields));
    if ~isempty(missing_fields)
        error('Metadata is missing required fields: %s', strjoin(missing_fields, ', '));
    end

    progress_log('INDEX-BUILD', '[2/7] Extracting epoch in chunks of %d records.', extraction_chunk_size);
    epoch_tic = tic;
    epoch_values = struct_field_numeric(mt, 1:source_row_count, 'epoch', extraction_chunk_size);
    candidate_idx = find(epoch_values == target_epoch);
    epoch_candidate_count = numel(candidate_idx);
    clear epoch_values;
    progress_log('INDEX-BUILD', '[2/7] epoch=%d candidates: %d/%d records in %.1f s.', ...
        target_epoch, epoch_candidate_count, source_row_count, toc(epoch_tic));

    if isempty(candidate_idx)
        error('No metadata records found for epoch %d.', target_epoch);
    end

    progress_log('INDEX-BUILD', '[3/7] Extracting fixed-filter and grouping fields for epoch candidates.');
    field_tic = tic;
    area = struct_field_string(mt, candidate_idx, 'area', extraction_chunk_size);
    fold_idx = struct_field_numeric(mt, candidate_idx, 'fold_idx', extraction_chunk_size);
    shuffle_idx = struct_field_numeric(mt, candidate_idx, 'shuffle_idx', extraction_chunk_size);
    prepost = struct_field_string(mt, candidate_idx, 'prepost', extraction_chunk_size);
    state = struct_field_string(mt, candidate_idx, 'state', extraction_chunk_size);

    animal_name = struct_field_string(mt, candidate_idx, 'animal_name', extraction_chunk_size);
    injection = struct_field_string(mt, candidate_idx, 'injection', extraction_chunk_size);
    align_name = struct_field_string(mt, candidate_idx, 'align', extraction_chunk_size);
    session_idx = struct_field_string(mt, candidate_idx, 'session_idx', extraction_chunk_size);
    resting_dur = struct_field_numeric(mt, candidate_idx, 'resting_dur_threshold', extraction_chunk_size);
    kernel_name = struct_field_string(mt, candidate_idx, 'kernel_name', extraction_chunk_size);
    reg_name = struct_field_string(mt, candidate_idx, 'reg_name', extraction_chunk_size);
    progress_log('INDEX-BUILD', '[3/7] Candidate fields extracted in %.1f s.', toc(field_tic));

    progress_log('INDEX-BUILD', '[4/7] Applying fixed filters without copying the full struct array.');
    keep = area == "Cortex" & ...
           fold_idx == 0 & ...
           shuffle_idx == 0 & ...
           ismember(prepost, ["Pre", "Post"]) & ...
           ismember(state, ["RestOpen", "RestClose"]);

    filtered_row_count = sum(keep);
    progress_log('INDEX-BUILD', ...
        '[4/7] Fixed filters retained %d/%d epoch candidates.', ...
        filtered_row_count, numel(candidate_idx));
    if filtered_row_count == 0
        error('No metadata records remain after fixed index filters.');
    end

    candidate_idx = candidate_idx(keep);
    prepost = prepost(keep);
    state = state(keep);
    animal_name = animal_name(keep);
    injection = injection(keep);
    align_name = align_name(keep);
    session_idx = session_idx(keep);
    resting_dur = resting_dur(keep);
    kernel_name = kernel_name(keep);
    reg_name = reg_name(keep);
    clear area fold_idx shuffle_idx keep;

    progress_log('INDEX-BUILD', '[5/7] Encoding conditions and removing missing grouping values.');
    condition_code = zeros(filtered_row_count, 1, 'uint8');
    condition_code(prepost == "Pre"  & state == "RestOpen")  = 1;
    condition_code(prepost == "Pre"  & state == "RestClose") = 2;
    condition_code(prepost == "Post" & state == "RestOpen")  = 3;
    condition_code(prepost == "Post" & state == "RestClose") = 4;

    valid_group_values = condition_code > 0 & ...
        ~ismissing(animal_name) & ~ismissing(injection) & ...
        ~ismissing(align_name) & ~ismissing(session_idx) & ...
        isfinite(resting_dur) & ~ismissing(kernel_name) & ~ismissing(reg_name);

    if ~all(valid_group_values)
        progress_log('INDEX-BUILD', 'Discarding %d records with invalid conditions or missing grouping values.', ...
            sum(~valid_group_values));
        candidate_idx = candidate_idx(valid_group_values);
        condition_code = condition_code(valid_group_values);
        animal_name = animal_name(valid_group_values);
        injection = injection(valid_group_values);
        align_name = align_name(valid_group_values);
        session_idx = session_idx(valid_group_values);
        resting_dur = resting_dur(valid_group_values);
        kernel_name = kernel_name(valid_group_values);
        reg_name = reg_name(valid_group_values);
    end

    if isempty(candidate_idx)
        error('No valid metadata records remain after grouping-value checks.');
    end

    progress_log('INDEX-BUILD', '[6/7] Grouping records by hyperparameters and session identity.');
    group_tic = tic;
    group_id = findgroups(animal_name, injection, align_name, session_idx, ...
        resting_dur, kernel_name, reg_name);
    if any(~isfinite(group_id))
        error('findgroups returned invalid group IDs after missing-value filtering.');
    end
    n_groups = max(group_id);
    progress_log('INDEX-BUILD', '[6/7] Created %d candidate session groups in %.1f s.', ...
        n_groups, toc(group_tic));

    condition_counts = accumarray( ...
        [double(group_id), double(condition_code)], 1, [n_groups, 4], @sum, 0);
    complete_group = all(condition_counts == 1, 2);
    missing_group = any(condition_counts == 0, 2);
    duplicate_group = any(condition_counts > 1, 2);

    progress_log('INDEX-BUILD', ...
        '[6/7] Complete=%d; missing-condition=%d; duplicate-condition=%d groups.', ...
        sum(complete_group), sum(missing_group), sum(duplicate_group));
    if ~any(complete_group)
        error('No complete Pre/Post x RestOpen/RestClose session groups were found.');
    end

    anchor_mask = condition_code == 1 & complete_group(group_id);
    anchor_source_idx = candidate_idx(anchor_mask);
    if numel(anchor_source_idx) ~= sum(complete_group)
        error('Internal index construction error: anchor count does not match complete group count.');
    end

    progress_log('INDEX-BUILD', '[7/7] Copying only %d complete anchor records into the compact struct index.', ...
        numel(anchor_source_idx));
    session_index = mt(anchor_source_idx);
    session_index = sort_session_index_struct(session_index);
    clear mt;

    progress_log('INDEX-BUILD', '[7/7] Compact struct index constructed with %d session records.', ...
        numel(session_index));

    index_info = struct();
    index_info.version = index_version;
    index_info.format = 'struct';
    index_info.created_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    index_info.target_epoch = target_epoch;
    index_info.fixed_area = 'Cortex';
    index_info.fixed_fold_idx = 0;
    index_info.fixed_shuffle_idx = 0;
    index_info.source_row_count = source_row_count;
    index_info.epoch_candidate_count = epoch_candidate_count;
    index_info.filtered_row_count = filtered_row_count;
    index_info.candidate_group_count = n_groups;
    index_info.complete_session_count = numel(session_index);
    index_info.missing_group_count = sum(missing_group);
    index_info.duplicate_group_count = sum(duplicate_group);
end

function selected = select_sessions_from_index(session_index, hp)
    n_sessions = numel(session_index);
    if n_sessions == 0
        selected = session_index;
        return;
    end

    all_idx = 1:n_sessions;
    kernel_name = struct_field_string(session_index, all_idx, ...
        'kernel_name', n_sessions);
    align_name = struct_field_string(session_index, all_idx, ...
        'align', n_sessions);
    injection = struct_field_string(session_index, all_idx, ...
        'injection', n_sessions);
    resting_dur = struct_field_numeric(session_index, all_idx, ...
        'resting_dur_threshold', n_sessions);

    requested_kernels = string(hp.kernel_names(:));
    mask = ismember(kernel_name, requested_kernels) & ...
           align_name == string(hp.align) & ...
           injection == string(hp.injection) & ...
           resting_dur == double(hp.resting_dur_threshold);

    if ~is_empty_filter_value(hp.reg_name)
        reg_name = struct_field_string(session_index, all_idx, ...
            'reg_name', n_sessions);
        mask = mask & reg_name == string(hp.reg_name);
    end

    selected = session_index(mask);
end

function report_kernel_setup_coverage(meta_array, hp)
    requested = string(hp.kernel_names(:));
    if isempty(meta_array)
        progress_log('CONFIG', 'Setup %s matched no kernel-session records.', ...
            char(string(hp.kernel_setup_name)));
        return;
    end

    matched = struct_field_string(meta_array, 1:numel(meta_array), ...
        'kernel_name', numel(meta_array));
    progress_log('CONFIG', 'Setup %s selected %d kernel-session records.', ...
        char(string(hp.kernel_setup_name)), numel(meta_array));
    for kernel_i = 1:numel(requested)
        n_records = sum(matched == requested(kernel_i));
        progress_log('CONFIG', '  kernel=%s: %d records.', ...
            char(requested(kernel_i)), n_records);
    end

    missing = requested(~ismember(requested, unique(matched)));
    if ~isempty(missing)
        warning('Setup %s has no matching records for: %s', ...
            char(string(hp.kernel_setup_name)), strjoin(cellstr(missing), ', '));
    end
end

function values = struct_field_numeric(records, indices, field, chunk_size)
    if nargin < 4 || isempty(chunk_size)
        chunk_size = 100000;
    end
    if ~isstruct(records) || ~isfield(records, field)
        error('Metadata struct is missing required field: %s', field);
    end

    indices = indices(:);
    n = numel(indices);
    values = nan(n, 1);

    for first_i = 1:chunk_size:n
        last_i = min(first_i + chunk_size - 1, n);
        source_idx = indices(first_i:last_i);
        target_idx = first_i:last_i;

        used_fast_path = false;
        try
            raw = [records(source_idx).(field)];
            if numel(raw) == numel(source_idx) && (isnumeric(raw) || islogical(raw))
                values(target_idx) = double(raw(:));
                used_fast_path = true;
            end
        catch
            used_fast_path = false;
        end

        if ~used_fast_path
            for j = 1:numel(source_idx)
                value = unwrap_cell_scalar(records(source_idx(j)).(field));
                if isempty(value)
                    continue;
                elseif (isnumeric(value) || islogical(value)) && isscalar(value)
                    values(target_idx(j)) = double(value);
                elseif (ischar(value) || isstring(value)) && isscalar(string(value))
                    parsed = str2double(string(value));
                    if ~isnan(parsed)
                        values(target_idx(j)) = parsed;
                    else
                        error('Field %s contains a nonnumeric value at struct record %d.', ...
                            field, source_idx(j));
                    end
                else
                    error('Field %s is not scalar numeric at struct record %d.', ...
                        field, source_idx(j));
                end
            end
        end
    end
end

function values = struct_field_string(records, indices, field, chunk_size)
    if nargin < 4 || isempty(chunk_size)
        chunk_size = 100000;
    end
    if ~isstruct(records) || ~isfield(records, field)
        error('Metadata struct is missing required field: %s', field);
    end

    indices = indices(:);
    n = numel(indices);
    values = strings(n, 1);

    for first_i = 1:chunk_size:n
        last_i = min(first_i + chunk_size - 1, n);
        source_idx = indices(first_i:last_i);
        target_idx = first_i:last_i;

        used_fast_path = false;
        try
            raw = {records(source_idx).(field)};
            converted = string(raw(:));
            if numel(converted) == numel(source_idx)
                values(target_idx) = converted;
                used_fast_path = true;
            end
        catch
            used_fast_path = false;
        end

        if ~used_fast_path
            for j = 1:numel(source_idx)
                value = unwrap_cell_scalar(records(source_idx(j)).(field));
                if isempty(value)
                    values(target_idx(j)) = missing;
                else
                    converted = string(value);
                    if ~isscalar(converted)
                        error('Field %s is not scalar-valued at struct record %d.', ...
                            field, source_idx(j));
                    end
                    values(target_idx(j)) = converted;
                end
            end
        end
    end
end

function session_index = sort_session_index_struct(session_index)
    n = numel(session_index);
    if n <= 1
        return;
    end

    idx = 1:n;
    separator = string(char(31));
    kernel_name = struct_field_string(session_index, idx, 'kernel_name', n);
    reg_name = struct_field_string(session_index, idx, 'reg_name', n);
    resting_dur = struct_field_numeric(session_index, idx, 'resting_dur_threshold', n);
    injection = struct_field_string(session_index, idx, 'injection', n);
    align_name = struct_field_string(session_index, idx, 'align', n);
    animal_name = struct_field_string(session_index, idx, 'animal_name', n);
    session_idx = struct_field_string(session_index, idx, 'session_idx', n);

    sort_key = kernel_name + separator + reg_name + separator + ...
        compose('%020.10g', resting_dur) + separator + injection + separator + ...
        align_name + separator + animal_name + separator + session_idx;
    [~, order] = sort(sort_key);
    session_index = session_index(order);
end

function progress_log(stage, format_string, varargin)
    timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    message = sprintf(format_string, varargin{:});
    fprintf('[%s][%s] %s\n', timestamp, stage, message);
end


function issue_log = empty_run_issue_log()
    issue_log = struct('hp_index', {}, 'hyperparam_title', {}, 'issue_type', {}, 'message', {});
end

function issue_log = append_run_issue(issue_log, hp_index, hyperparam_title, issue_type, message)
    entry_i = numel(issue_log) + 1;
    issue_log(entry_i).hp_index = hp_index;
    issue_log(entry_i).hyperparam_title = char(string(hyperparam_title));
    issue_log(entry_i).issue_type = char(string(issue_type));
    issue_log(entry_i).message = char(string(message));
end

function write_run_issue_log(filename, issue_log)
    fid = fopen(filename, 'w');
    if fid < 0
        warning('Could not open run issue log for writing: %s', filename);
        return;
    end
    cleanup_obj = onCleanup(@() fclose(fid));

    fprintf(fid, 'hp_index\tissue_type\thyperparam_title\tmessage\n');
    for i = 1:numel(issue_log)
        fprintf(fid, '%d\t%s\t%s\t%s\n', ...
            issue_log(i).hp_index, ...
            escape_log_field(issue_log(i).issue_type), ...
            escape_log_field(issue_log(i).hyperparam_title), ...
            escape_log_field(issue_log(i).message));
    end
end

function out = escape_log_field(value)
    out = char(string(value));
    out = strrep(out, sprintf('\t'), ' ');
    out = strrep(out, sprintf('\n'), ' ');
    out = strrep(out, sprintf('\r'), ' ');
end


function figure_configs = build_pattern_figure_configs()
    figure_configs = struct([]);

    figure_configs(1).analysis_type = 'pre_to_post_oc_pattern';
    figure_configs(1).output_stub = ...
        'SuppFig_OC_pattern_pre_to_post_pooled_kernels';
    figure_configs(1).figure_title = ...
        'Pooled kernels: Pre to Post transition of Open/Close patterns';
    figure_configs(1).x_label = 'Pre Open/Close pattern';
    figure_configs(1).y_label = 'Post Open/Close pattern';
    figure_configs(1).row1_title = '9x9 Open/Close pattern transition';
    figure_configs(1).row2_title = ...
        '4x4 Open/Close pattern transition, significant only';

    figure_configs(2).analysis_type = 'open_vs_close_prepost_pattern';
    figure_configs(2).output_stub = ...
        'SuppFig_PrePost_pattern_open_vs_close_pooled_kernels';
    figure_configs(2).figure_title = ...
        'Pooled kernels: Open vs Close comparison of Pre/Post patterns';
    figure_configs(2).x_label = 'Open Pre/Post pattern';
    figure_configs(2).y_label = 'Close Pre/Post pattern';
    figure_configs(2).row1_title = '9x9 Pre/Post pattern comparison';
    figure_configs(2).row2_title = ...
        '4x4 Pre/Post pattern comparison, significant only';
end

function loaded_states = load_all_required_states(root, meta, kernel_indices)
    loaded_states = struct();
    preposts = {'Pre', 'Post'};
    states = {'RestOpen', 'RestClose'};
    n_loads = numel(kernel_indices) * numel(preposts) * numel(states);
    load_i = 0;
    progress_log('LOAD', 'Loading %d state/kernel combinations for %s.', ...
        n_loads, make_session_label(meta));
    for k_i = 1:numel(kernel_indices)
        kernel_idx = kernel_indices(k_i);
        for pp_i = 1:numel(preposts)
            for st_i = 1:numel(states)
                load_i = load_i + 1;
                pp = preposts{pp_i};
                st = states{st_i};
                key = state_key(pp, st, kernel_idx);
                progress_log('LOAD', '[%d/%d] kernel=%d, prepost=%s, state=%s.', ...
                    load_i, n_loads, kernel_idx, pp, st);
                loaded_states.(key) = load_state_connectivity(root, meta, pp, st, kernel_idx);
            end
        end
    end
end

function key = state_key(prepost, state, kernel_idx)
    key = sprintf('k%d_%s_%s', kernel_idx, prepost, state);
end

function [counts9, counts4] = make_pre_to_post_oc_pattern_transition_counts( ...
    pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi)

    [pre_open_cat, pre_close_cat, post_open_cat, post_close_cat] = make_all_condition_categories( ...
        pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi);

    % 9x9: Pre Open/Close compound pattern -> Post Open/Close compound pattern.
    pre_code9  = compound_category_code_9(pre_open_cat,  pre_close_cat);
    post_code9 = compound_category_code_9(post_open_cat, post_close_cat);
    counts9 = make_transition_counts(pre_code9, post_code9, 9);

    % 4x4: same transition, but only connections with no non-significant member
    % in either the Pre or Post Open/Close pair are retained.
    pre_code4  = compound_category_code_4(pre_open_cat,  pre_close_cat);
    post_code4 = compound_category_code_4(post_open_cat, post_close_cat);
    counts4 = make_transition_counts(pre_code4, post_code4, 4);
end

function [counts9, counts4] = make_open_vs_close_prepost_pattern_counts( ...
    pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi)

    [pre_open_cat, pre_close_cat, post_open_cat, post_close_cat] = make_all_condition_categories( ...
        pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi);

    % 9x9: Open Pre/Post compound pattern -> Close Pre/Post compound pattern.
    open_code9  = compound_category_code_9(pre_open_cat,  post_open_cat);
    close_code9 = compound_category_code_9(pre_close_cat, post_close_cat);
    counts9 = make_transition_counts(open_code9, close_code9, 9);

    % 4x4: same comparison, but only connections with no non-significant member
    % in either the Open or Close Pre/Post pair are retained.
    open_code4  = compound_category_code_4(pre_open_cat,  post_open_cat);
    close_code4 = compound_category_code_4(pre_close_cat, post_close_cat);
    counts4 = make_transition_counts(open_code4, close_code4, 4);
end

function [pre_open_cat, pre_close_cat, post_open_cat, post_close_cat] = make_all_condition_categories( ...
    pre_open_state, pre_close_state, post_open_state, post_close_state, err_multi)

    validate_matching_filters(pre_open_state, pre_close_state);
    validate_matching_filters(pre_open_state, post_open_state);
    validate_matching_filters(pre_open_state, post_close_state);

    pre_open_cat12   = classify_connections(pre_open_state.J12(:),   pre_open_state.err12(:),   err_multi);
    pre_close_cat12  = classify_connections(pre_close_state.J12(:),  pre_close_state.err12(:),  err_multi);
    post_open_cat12  = classify_connections(post_open_state.J12(:),  post_open_state.err12(:),  err_multi);
    post_close_cat12 = classify_connections(post_close_state.J12(:), post_close_state.err12(:), err_multi);

    pre_open_cat21   = classify_connections(pre_open_state.J21(:),   pre_open_state.err21(:),   err_multi);
    pre_close_cat21  = classify_connections(pre_close_state.J21(:),  pre_close_state.err21(:),  err_multi);
    post_open_cat21  = classify_connections(post_open_state.J21(:),  post_open_state.err21(:),  err_multi);
    post_close_cat21 = classify_connections(post_close_state.J21(:), post_close_state.err21(:), err_multi);

    pre_open_cat   = [pre_open_cat12;   pre_open_cat21];
    pre_close_cat  = [pre_close_cat12;  pre_close_cat21];
    post_open_cat  = [post_open_cat12;  post_open_cat21];
    post_close_cat = [post_close_cat12; post_close_cat21];

    valid = isfinite(pre_open_cat) & isfinite(pre_close_cat) & ...
            isfinite(post_open_cat) & isfinite(post_close_cat);

    pre_open_cat   = pre_open_cat(valid);
    pre_close_cat  = pre_close_cat(valid);
    post_open_cat  = post_open_cat(valid);
    post_close_cat = post_close_cat(valid);
end

function code = compound_category_code_9(first_cat, second_cat)
    % Category order:
    %   ++, +o, +-, o+, oo, o-, -+, -o, --
    % first_cat and second_cat must use: +1 = positive, 0 = non-significant, -1 = negative.
    code = nan(size(first_cat));

    code(first_cat ==  1 & second_cat ==  1) = 1;
    code(first_cat ==  1 & second_cat ==  0) = 2;
    code(first_cat ==  1 & second_cat == -1) = 3;

    code(first_cat ==  0 & second_cat ==  1) = 4;
    code(first_cat ==  0 & second_cat ==  0) = 5;
    code(first_cat ==  0 & second_cat == -1) = 6;

    code(first_cat == -1 & second_cat ==  1) = 7;
    code(first_cat == -1 & second_cat ==  0) = 8;
    code(first_cat == -1 & second_cat == -1) = 9;
end

function code = compound_category_code_4(first_cat, second_cat)
    % Category order:
    %   ++, +-, -+, --
    % Non-significant members are returned as NaN and excluded from counts.
    code = nan(size(first_cat));

    code(first_cat ==  1 & second_cat ==  1) = 1;
    code(first_cat ==  1 & second_cat == -1) = 2;
    code(first_cat == -1 & second_cat ==  1) = 3;
    code(first_cat == -1 & second_cat == -1) = 4;
end

function counts = make_transition_counts(x_code, y_code, n_cat)
    valid = isfinite(x_code) & isfinite(y_code);
    x_code = x_code(valid);
    y_code = y_code(valid);

    counts = zeros(n_cat, n_cat);
    for i_cat = 1:n_cat
        for j_cat = 1:n_cat
            counts(i_cat, j_cat) = sum(x_code == i_cat & y_code == j_cat);
        end
    end
end

function render_pattern_transition_figure(root, cfg, pooled_one, params) %#ok<INUSD>
    render_tic = tic;
    f = figure('Color', 'w', 'Visible', params.figure_visible);
    figure_cleanup = onCleanup(@() close_figure_if_valid(f)); %#ok<NASGU>

    tiledlayout(params.n_row, params.n_col, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    tile_idx = @(row, col) rowcol_to_panel_index(row, col, params.n_col);

    plot_modes = {'raw', 'row', 'column', 'transition_difference', 'normalized_difference', 'oe_ratio', 'standardized_residual'};
    column_titles = {'raw count', 'normalized by row', 'normalized by column', ...
                     'transition difference', 'normalized difference', 'O/E ratio', 'standardized residual'};
    colorbar_labels = {'count', 'row fraction', 'column fraction', ...
                       'count difference', 'normalized difference', 'O/E ratio', 'std residual'};
    value_modes = {'integer', 'fraction', 'fraction', 'signed', 'signed', 'ratio', 'signed'};

    progress_log('RENDER', 'Building 14 panels for %s.', cfg.output_stub);
    for col_i = 1:numel(plot_modes)
        mode = plot_modes{col_i};

        counts = pooled_one.counts9;
        [agreement, kappa, n_valid] = summarize_counts(counts);
        mat = transform_transition_matrix(counts, mode);
        ax = nexttile(tile_idx(1, col_i));
        plot_transition_table_generic(ax, mat, params.category_labels_9, agreement, kappa, n_valid, ...
            cfg.row1_title, cfg.x_label, cfg.y_label, column_titles{col_i}, colorbar_labels{col_i}, value_modes{col_i});
        add_panel_label(ax, 1, col_i, params.n_col);

        counts = pooled_one.counts4;
        [agreement, kappa, n_valid] = summarize_counts(counts);
        mat = transform_transition_matrix(counts, mode);
        ax = nexttile(tile_idx(2, col_i));
        plot_transition_table_generic(ax, mat, params.category_labels_4, agreement, kappa, n_valid, ...
            cfg.row2_title, cfg.x_label, cfg.y_label, column_titles{col_i}, colorbar_labels{col_i}, value_modes{col_i});
        add_panel_label(ax, 2, col_i, params.n_col);
    end

    sgtitle(sprintf('Supplement Fig: %s; %s; valid kernel-session records = %d', ...
        cfg.figure_title, params.hyperparam_title, pooled_one.valid_session_count), ...
        'Interpreter', 'none');

    save_folder = params.output_folder;
    check_path(save_folder);

    figWidth = 28.0;  % inches.
    figHeight = 8.5;  % inches.

    set(f, 'Units', 'inches');
    f.Position(3:4) = [figWidth, figHeight];
    set(f, 'PaperUnits', 'inches');
    set(f, 'PaperSize', [figWidth, figHeight]);
    set(f, 'PaperPosition', [0, 0, figWidth, figHeight]);
    set(f, 'Color', 'w');

    if params.export_jpg
        jpg_filename = fullfile(save_folder, [cfg.output_stub, '_supp_preview.jpg']);
        progress_log('EXPORT', 'Writing JPG: %s', jpg_filename);
        exportgraphics(f, jpg_filename, 'ContentType', 'image', ...
            'BackgroundColor', 'white', 'Resolution', params.export_resolution);
    else
        progress_log('EXPORT', 'JPG export disabled for %s.', cfg.output_stub);
    end

    if params.export_pdf
        pdf_filename = fullfile(save_folder, [cfg.output_stub, '_supp.pdf']);
        progress_log('EXPORT', 'Writing PDF: %s', pdf_filename);
        exportgraphics(f, pdf_filename, 'ContentType', 'vector', ...
            'BackgroundColor', 'white');
    else
        progress_log('EXPORT', 'PDF export disabled for %s.', cfg.output_stub);
    end

    progress_log('RENDER', 'Completed %s in %.1f s.', cfg.output_stub, toc(render_tic));
end

function close_figure_if_valid(f)
    if ~isempty(f) && isgraphics(f)
        close(f);
    end
end

function [agreement, kappa, n_valid] = summarize_counts(counts)
    agreement = compute_raw_agreement(counts);
    kappa = compute_cohen_kappa(counts);
    n_valid = sum(counts(:));
end

function mat_out = transform_transition_matrix(counts, mode)
    counts = double(counts);
    total_n = sum(counts(:));

    switch lower(mode)
        case 'raw'
            mat_out = counts;

        case 'row'
            denom = sum(counts, 2);
            mat_out = zeros(size(counts));
            valid = denom > 0;
            mat_out(valid, :) = bsxfun(@rdivide, counts(valid, :), denom(valid));

        case 'column'
            denom = sum(counts, 1);
            mat_out = zeros(size(counts));
            valid = denom > 0;
            mat_out(:, valid) = bsxfun(@rdivide, counts(:, valid), denom(valid));

        case 'transition_difference'
            mat_out = counts - counts.';

        case 'normalized_difference'
            pair_sum = counts + counts.';
            mat_out = nan(size(counts));
            valid = pair_sum > 0;
            diff_mat = counts - counts.';
            mat_out(valid) = diff_mat(valid) ./ pair_sum(valid);

        case 'oe_ratio'
            expected_prob = expected_probability_under_independence(counts);
            observed_prob = zeros(size(counts));
            if total_n > 0
                observed_prob = counts / total_n;
            end
            mat_out = nan(size(counts));
            valid = expected_prob > 0;
            mat_out(valid) = observed_prob(valid) ./ expected_prob(valid);

        case 'standardized_residual'
            expected_counts = expected_counts_under_independence(counts);
            mat_out = nan(size(counts));
            valid = expected_counts > 0;
            mat_out(valid) = (counts(valid) - expected_counts(valid)) ./ sqrt(expected_counts(valid));

        otherwise
            error('Unknown transform mode: %s', mode);
    end
end

function expected_prob = expected_probability_under_independence(counts)
    counts = double(counts);
    total_n = sum(counts(:));
    if total_n <= 0
        expected_prob = zeros(size(counts));
        return;
    end
    x_marginal_prob = sum(counts, 2) / total_n;
    y_marginal_prob = sum(counts, 1) / total_n;
    expected_prob = x_marginal_prob * y_marginal_prob;
end

function expected_counts = expected_counts_under_independence(counts)
    counts = double(counts);
    total_n = sum(counts(:));
    if total_n <= 0
        expected_counts = zeros(size(counts));
        return;
    end
    expected_counts = total_n * expected_probability_under_independence(counts);
end

function plot_transition_table_generic(ax, mat, category_labels, agreement, kappa, n_valid, title_text, xlabel_text, ylabel_text, normalization_label, colorbar_label, value_mode)
    plot_mat = mat.';
    imagesc(ax, plot_mat);
    set(ax, 'YDir', 'normal');
    axis(ax, 'square');

    finite_vals = plot_mat(isfinite(plot_mat));
    if isempty(finite_vals)
        caxis(ax, [0, 1]);
    elseif strcmp(value_mode, 'signed')
        max_abs = max(abs(finite_vals));
        if max_abs > 0
            caxis(ax, [-max_abs, max_abs]);
        end
    end

    n_cat = numel(category_labels);
    xticks(ax, 1:n_cat);
    yticks(ax, 1:n_cat);
    xticklabels(ax, category_labels);
    yticklabels(ax, category_labels);
    xtickangle(ax, 30);

    xlabel(ax, xlabel_text);
    ylabel(ax, ylabel_text);

    cb = colorbar(ax);
    ylabel(cb, colorbar_label);

    max_val = max(finite_vals);
    min_val = min(finite_vals);
    for y_idx = 1:n_cat
        for x_idx = 1:n_cat
            this_val = plot_mat(y_idx, x_idx);
            if ~isfinite(this_val)
                value_str = 'NaN';
                text_color = 'k';
            else
                this_ratio = (this_val - min_val) / (max_val - min_val + eps);
                if this_ratio < 0.2
                    text_color = 'w';
                else
                    text_color = 'k';
                end

                switch value_mode
                    case 'integer'
                        value_str = sprintf('%d', round(this_val));
                    case 'fraction'
                        value_str = sprintf('%.2f', this_val);
                    case 'ratio'
                        value_str = sprintf('%.2f', this_val);
                    case 'signed'
                        value_str = sprintf('%.2f', this_val);
                    otherwise
                        value_str = sprintf('%.2f', this_val);
                end
            end

            text(ax, x_idx, y_idx, value_str, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold', 'Color', text_color, 'FontSize', choose_cell_font_size(n_cat));
        end
    end

    if isnan(agreement)
        agreement_str = 'NaN';
    else
        agreement_str = sprintf('%.3f', agreement);
    end
    if isnan(kappa)
        kappa_str = 'NaN';
    else
        kappa_str = sprintf('%.3f', kappa);
    end

    title(ax, sprintf('%s\n%s\nAgreement = %s, kappa = %s, n = %d', ...
        title_text, normalization_label, agreement_str, kappa_str, n_valid), 'Interpreter', 'none');
end

function font_size = choose_cell_font_size(n_cat)
    if n_cat <= 4
        font_size = 8;
    else
        font_size = 6;
    end
end

function add_panel_label(ax, row, col, n_col)
    panel_index = rowcol_to_panel_index(row, col, n_col);
    panel_label = panel_index_to_letters(panel_index);
    text(ax, -0.10, 1.08, panel_label, 'Units', 'normalized', ...
        'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', 'Interpreter', 'none', 'Clipping', 'off');
end

function panel_index = rowcol_to_panel_index(row, col, n_col)
    panel_index = (row - 1) * n_col + col;
end

function label = panel_index_to_letters(panel_index)
    if panel_index < 1 || panel_index ~= floor(panel_index)
        error('panel_index must be a positive integer.');
    end
    label = '';
    n = panel_index;
    while n > 0
        rem0 = mod(n - 1, 26);
        label = [char(double('A') + rem0), label]; %#ok<AGROW>
        n = floor((n - 1) / 26);
    end
end

function hp = fill_default_hyperparams(hp)
    default_hp = struct();
    default_hp.kernel_setup_name = "default_setup";
    default_hp.kernel_names = {'DeltaPure'};
    default_hp.reg_name = [];
    default_hp.resting_dur_threshold = 15;
    default_hp.injection = "Muscimol";
    default_hp.align = "Last";

    fields = fieldnames(default_hp);
    for f = 1:numel(fields)
        field = fields{f};
        if ~isfield(hp, field)
            hp.(field) = default_hp.(field);
        end
    end

    setup_name = string(hp.kernel_setup_name);
    kernel_names = string(hp.kernel_names(:));
    if ~isscalar(setup_name) || ismissing(setup_name) || strlength(setup_name) == 0
        error('kernel_setup_name must be a nonempty scalar string.');
    end
    if isempty(kernel_names) || any(ismissing(kernel_names)) || ...
            any(strlength(kernel_names) == 0)
        error('kernel_names must contain at least one valid kernel name.');
    end
    hp.kernel_setup_name = char(setup_name);
    hp.kernel_names = cellstr(kernel_names).';
end

function value = unwrap_cell_scalar(value)
    while iscell(value) && isscalar(value)
        value = value{1};
    end
end

function tf = is_empty_filter_value(value)
    if isempty(value)
        tf = true;
        return;
    end
    if isstring(value) && all(strlength(value) == 0)
        tf = true;
        return;
    end
    if iscell(value) && isempty(value)
        tf = true;
        return;
    end
    tf = false;
end

function title_str = make_hyperparam_title(hp)
    title_str = sprintf(['kernel_setup=%s (%d kernels), reg_name=%s, ' ...
        'resting_dur_threshold=%s, injection=%s, align=%s'], ...
        value_to_display_string(hp.kernel_setup_name, 'all'), ...
        numel(hp.kernel_names), ...
        value_to_display_string(hp.reg_name, 'all'), ...
        value_to_display_string(hp.resting_dur_threshold, 'all'), ...
        value_to_display_string(hp.injection, 'all'), ...
        value_to_display_string(hp.align, 'all'));
end
function folder_name = make_hyperparam_folder_name(hp)
    parts = { ...
        ['kernelSetup_', sanitize_for_path(value_to_display_string( ...
            hp.kernel_setup_name, 'all'))], ...
        ['reg_', sanitize_for_path(value_to_display_string(hp.reg_name, 'all'))], ...
        ['restDur_', sanitize_for_path(value_to_display_string( ...
            hp.resting_dur_threshold, 'all'))], ...
        ['inj_', sanitize_for_path(value_to_display_string(hp.injection, 'all'))], ...
        ['align_', sanitize_for_path(value_to_display_string(hp.align, 'all'))] ...
    };
    folder_name = strjoin(parts, '__');
end

function out = value_to_display_string(value, empty_label)
    if nargin < 2
        empty_label = 'all';
    end

    if is_empty_filter_value(value)
        out = empty_label;
        return;
    end

    value = unwrap_cell_scalar(value);

    if iscell(value)
        parts = cell(size(value));
        for i = 1:numel(value)
            parts{i} = value_to_display_string(value{i}, empty_label);
        end
        out = strjoin(parts(:).', '+');
    elseif isstring(value)
        parts = cell(1, numel(value));
        for i = 1:numel(value)
            parts{i} = char(value(i));
        end
        out = strjoin(parts, '+');
    elseif ischar(value)
        out = value;
    elseif isnumeric(value) || islogical(value)
        if isscalar(value)
            out = num2str(value);
        else
            parts = cell(1, numel(value));
            for i = 1:numel(value)
                parts{i} = num2str(value(i));
            end
            out = strjoin(parts, '+');
        end
    else
        out = char(string(value));
    end
end

function out = sanitize_for_path(str_in)
    out = char(str_in);
    out = regexprep(out, '[^A-Za-z0-9._=-]+', '_');
    out = regexprep(out, '^_+|_+$', '');
    if isempty(out)
        out = 'empty';
    end
end

function label = make_session_label(meta)
    animal_str = char(string(meta.animal_name));
    injection_str = char(string(meta.injection));
    session_str = char(string(meta.session_idx));
    if isfield(meta, 'kernel_name')
        kernel_str = char(string(meta.kernel_name));
    else
        kernel_str = 'unknown_kernel';
    end
    label = sprintf('%s %s session %s, kernel %s', ...
        animal_str, injection_str, session_str, kernel_str);
end

function state_struct = load_state_connectivity(root, meta, prepost, state, kernel_idx)
    state_tic = tic;
    state_struct = struct();
    state_struct.prepost = prepost;
    state_struct.state = state;
    state_struct.kernel_idx = kernel_idx;

    meta.prepost = prepost;
    meta.state = state;

    meta.filename = generate_filename('raster', meta);
    raster_filename = fullfile(root, 'Data', 'Working', 'raster', meta.filename);
    progress_log('LOAD-RASTER', 'Loading: %s', raster_filename);
    raster_data = load(raster_filename);
    progress_log('LOAD-RASTER', 'Loaded %s %s %s; trial_len=%d, N=%d, trials=%d.', ...
        meta.prepost, meta.state, meta.area, raster_data.meta.trial_len, ...
        raster_data.meta.N, raster_data.meta.trial_num);

    cell_area = raster_data.data.cell_area;
    filter1 = ismember(cell_area, {'ACC'});
    filter2 = ismember(cell_area, {'VLPFC'});

    meta.filename = generate_filename('GLM', meta);
    glm_filename = fullfile(root, 'Data', 'Working', 'GLM', meta.filename);
    progress_log('LOAD-GLM', 'Loading: %s', glm_filename);
    GLM_data = load(glm_filename);
    N = GLM_data.meta.N;

    if numel(cell_area) ~= N
        error('Raster cell_area length (%d) does not match GLM N (%d): %s', ...
            numel(cell_area), N, glm_filename);
    end

    col_start = 2 + N * (kernel_idx - 1);
    col_end = 1 + N * kernel_idx;
    if col_end > size(GLM_data.data.model_par, 2)
        error('Kernel %d requires model_par columns %d:%d, but only %d columns exist.', ...
            kernel_idx, col_start, col_end, size(GLM_data.data.model_par, 2));
    end
    if col_end > size(GLM_data.data.model_err.total, 2)
        error('Kernel %d requires model_err.total columns %d:%d, but only %d columns exist.', ...
            kernel_idx, col_start, col_end, size(GLM_data.data.model_err.total, 2));
    end

    J = GLM_data.data.model_par(:, col_start:col_end);
    err = GLM_data.data.model_err.total(:, col_start:col_end);

    if size(J, 1) ~= N || size(err, 1) ~= N || ~isequal(size(J), size(err))
        error('GLM parameter/error matrix dimensions are inconsistent for %s.', glm_filename);
    end

    state_struct.cell_area = cell_area;
    state_struct.filter1 = filter1;
    state_struct.filter2 = filter2;
    state_struct.J12 = J(filter1, filter2);
    state_struct.J21 = J(filter2, filter1);
    state_struct.err12 = err(filter1, filter2);
    state_struct.err21 = err(filter2, filter1);

    progress_log('LOAD-GLM', ...
        'Loaded %s %s kernel=%d; ACC=%d, VLPFC=%d, elapsed=%.1f s.', ...
        meta.prepost, meta.state, kernel_idx, sum(filter1), sum(filter2), toc(state_tic));
end

function cat = classify_connections(J, err, err_multi)
    cat = nan(size(J));
    valid = isfinite(J) & isfinite(err);
    cat(valid) = 0;
    cat(valid & (J >  err_multi * err)) = 1;
    cat(valid & (J < -err_multi * err)) = -1;
end

function agreement = compute_raw_agreement(counts)
    total_n = sum(counts(:));
    if total_n == 0
        agreement = NaN;
    else
        agreement = trace(counts) / total_n;
    end
end

function kappa = compute_cohen_kappa(counts)
    total_n = sum(counts(:));
    if total_n == 0
        kappa = NaN;
        return;
    end

    p_observed = trace(counts) / total_n;
    row_marginal = sum(counts, 2) / total_n;
    col_marginal = sum(counts, 1) / total_n;
    p_expected = row_marginal.' * col_marginal.';

    if abs(1 - p_expected) < eps
        kappa = NaN;
    else
        kappa = (p_observed - p_expected) / (1 - p_expected);
    end
end

function validate_matching_filters(state_a, state_b)
    if numel(state_a.filter1) ~= numel(state_b.filter1) || any(state_a.filter1(:) ~= state_b.filter1(:)) || ...
       numel(state_a.filter2) ~= numel(state_b.filter2) || any(state_a.filter2(:) ~= state_b.filter2(:))
        error('State filters do not match between the two comparison members.');
    end

    if ~isequal(size(state_a.J12), size(state_b.J12)) || ~isequal(size(state_a.J21), size(state_b.J21))
        error('State connectivity matrices do not have matching sizes.');
    end
end
