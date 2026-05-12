%% generate_filename.m - Generate standardized file names from metadata.
function filename = generate_filename(data_type, meta)
switch data_type
    case 'raster'
        required_fields = {'animal_name', 'injection', 'prepost', 'state', 'area', 'align', 'session_idx'};
        check_required_fields(meta, required_fields);

        animal_code = get_code('animal_name', meta.animal_name);
        injection_code = get_code('injection', meta.injection);
        prepost_code = get_code('prepost', meta.prepost);
        state_code = get_code('state', meta.state);
        area_code = get_code('area', meta.area);
        align_code = get_code('align', meta.align);

        if isfield(meta, 'resting_dur_threshold')
            dur_code = sprintf('dur%d', meta.resting_dur_threshold);
        else
            dur_code = '';
        end
        
        filename = sprintf('raster_%s%s%s%s%s%s%d%s.mat', animal_code, injection_code, prepost_code, state_code, area_code, align_code, meta.session_idx, dur_code);
    case 'border'
        required_fields = {'animal_name', 'injection', 'prepost', 'area', 'session_idx'};
        check_required_fields(meta, required_fields);

        animal_code = get_code('animal_name', meta.animal_name);
        injection_code = get_code('injection', meta.injection);
        prepost_code = get_code('prepost', meta.prepost);
        area_code = get_code('area', meta.area);
        % align_code = get_code('align', meta.align);
        special_code = 'None';

        filename = sprintf('border_%s%s%s%s%s%d.mat', animal_code, injection_code, prepost_code, area_code, special_code, meta.session_idx);
    
    case 'sortidx'
        required_fields = {'animal_name', 'injection', 'prepost', 'state', 'area', 'align', 'session_idx', 'criterion'};
        check_required_fields(meta, required_fields);

        animal_code = get_code('animal_name', meta.animal_name);
        injection_code = get_code('injection', meta.injection);
        prepost_code = get_code('prepost', meta.prepost);
        state_code = get_code('state', meta.state);
        area_code = get_code('area', meta.area);
        align_code = get_code('align', meta.align);
        criterion_code = get_code('criterion', meta.criterion);

        filename = sprintf('sortidx_%s%s%s%s%s%s%s%d.mat', animal_code, injection_code, prepost_code, state_code, area_code, align_code, criterion_code, meta.session_idx);
    
    case 'kernel'
        required_fields = {'kernel_name'};
        check_required_fields(meta, required_fields);
        filename = sprintf('kernel_%s.mat', meta.kernel_name);

    case 'shuffled'
        required_fields = {'animal_name', 'injection', 'prepost', 'state', 'area', 'align', 'resting_dur_threshold', 'session_idx', 'shuffle_idx'};
        check_required_fields(meta, required_fields);

        animal_code = get_code('animal_name', meta.animal_name);
        injection_code = get_code('injection', meta.injection);
        prepost_code = get_code('prepost', meta.prepost);
        state_code = get_code('state', meta.state);
        area_code = get_code('area', meta.area);
        align_code = get_code('align', meta.align);
        dur_code = get_code('resting_dur_threshold', meta.resting_dur_threshold);

        filename = sprintf('shuffled_%s%s%s%s%s%s%d%s_%d.mat', animal_code, injection_code, prepost_code, state_code, area_code, align_code, meta.session_idx, dur_code, meta.shuffle_idx);

    case 'crossval'
        required_fields = {'animal_name', 'injection', 'prepost', 'state', 'area', 'align', 'resting_dur_threshold', 'session_idx', 'shuffle_idx'};
        check_required_fields(meta, required_fields);
        
        animal_code = get_code('animal_name', meta.animal_name);
        injection_code = get_code('injection', meta.injection);
        prepost_code = get_code('prepost', meta.prepost);
        state_code = get_code('state', meta.state);
        area_code = get_code('area', meta.area);
        align_code = get_code('align', meta.align);
        dur_code = get_code('resting_dur_threshold', meta.resting_dur_threshold);

        filename = sprintf('crossval_%s%s%s%s%s%s%d%s_%d.mat', animal_code, injection_code, prepost_code, state_code, area_code, align_code, meta.session_idx, dur_code, meta.shuffle_idx);

    case 'GLMdata'
        required_fields = {'animal_name', 'injection', 'prepost', 'state', 'area', 'align', 'resting_dur_threshold', 'session_idx', 'shuffle_idx', 'kernel_name'};
        check_required_fields(meta, required_fields);
        
        animal_code = get_code('animal_name', meta.animal_name);
        injection_code = get_code('injection', meta.injection);
        prepost_code = get_code('prepost', meta.prepost);
        state_code = get_code('state', meta.state);
        area_code = get_code('area', meta.area);
        align_code = get_code('align', meta.align);
        dur_code = get_code('resting_dur_threshold', meta.resting_dur_threshold);

        filename = sprintf('GLMdata_%s%s%s%s%s%s%d%s_%d_%s.mat', animal_code, injection_code, prepost_code, state_code, area_code, align_code, meta.session_idx, dur_code, meta.shuffle_idx, meta.kernel_name);
    case 'GLM'
        required_fields = {'animal_name', 'injection', 'prepost', ...
        'state', 'area', 'align', 'resting_dur_threshold', 'session_idx', 'shuffle_idx', ...
        'kernel_name', 'reg_name', 'epoch', 'fold_idx'};
        check_required_fields(meta, required_fields);

        animal_code = get_code('animal_name', meta.animal_name);
        injection_code = get_code('injection', meta.injection);
        prepost_code = get_code('prepost', meta.prepost);
        state_code = get_code('state', meta.state);
        area_code = get_code('area', meta.area);
        align_code = get_code('align', meta.align);
        dur_code = get_code('resting_dur_threshold', meta.resting_dur_threshold);

        filename = sprintf('GLM_%s%s%s%s%s%s%d%s_%d_%s_%s_epoch%d_fold%d.mat', ...
            animal_code, injection_code, prepost_code, state_code, area_code, align_code, ...
            meta.session_idx, dur_code, meta.shuffle_idx, meta.kernel_name, meta.reg_name, meta.epoch, meta.fold_idx);
    
    case 'tuning'
        required_fields = {'animal_name', 'injection', 'area', 'session_idx'};
        check_required_fields(meta, required_fields);

        animal_code = get_code('animal_name', meta.animal_name);
        injection_code = get_code('injection', meta.injection);
        area_code = get_code('area', meta.area);

        filename = sprintf('tuning_%s%s%s%d.mat', animal_code, injection_code, area_code, meta.session_idx);

    otherwise
        error('Unknown data type: %s', data_type);
end
end % of function

function check_required_fields(meta, required_fields)
    for i = 1:numel(required_fields)
        if ~isfield(meta, required_fields{i})
            error('Missing required field: %s', required_fields{i});
        end
    end
end

function code = get_code(field_name, field_value)
    switch field_name
        case 'animal_name'
            switch field_value
                case 'Slayer'
                    code = 'Slayer';
                case 'Zeppelin'
                    code = 'Zeppelin';
                case 'Emperor'
                    code = 'Emperor';
                otherwise
                    error('Unknown animal name: %s', field_value);
            end
        case 'injection'
            switch field_value
                case 'Muscimol'
                    code = 'Mus';
                case 'Saline'
                    code = 'Sal';
                case 'No injection'
                    code = 'Noinj';
                otherwise
                    error('Unknown injection type: %s', field_value);
            end
        case 'prepost'
            switch field_value
                case 'Pre'
                    code = 'Pre';
                case 'Post'
                    code = 'Post';
                otherwise
                    error('Unknown prepost type: %s', field_value);
            end
        case 'state'
            switch field_value
                case 'Task'
                    code = 'Task';
                case 'RestOpen'
                    code = 'Open';
                case 'RestClose'
                    code = 'Closed';
                otherwise
                    error('Unknown state type: %s', field_value);
            end
        case 'area'
            switch field_value
                case 'ACC'
                    code = 'ACC';
                case 'VLPFC'
                    code = 'PFC';
                case 'Thalamus'
                    code = 'Thal';
                case 'Cortex'
                    code = 'Ctx';
                case 'Full'
                    code = 'Full';
                otherwise
                    error('Unknown area type: %s', field_value);
            end
        case 'align'
            switch field_value
                case 'None'
                    code = 'None';
                case 'First'
                    code = 'First';
                case 'Last'
                    code = 'Last';
                case 'Random'
                    code = 'Rand';
                case 'Longest'
                    code = 'Longest';
                otherwise
                    error('Unknown align type: %s', field_value);
            end
        case 'criterion'
            switch field_value
                case 'channel'
                    code = 'Channel';
                otherwise
                    error('Unknown criterion type: %s', field_value);
            end
        case 'resting_dur_threshold'
            code = sprintf('dur%d', field_value);
        otherwise
            error('Unknown field name: %s', field_name);
    end
end