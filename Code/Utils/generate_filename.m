%% generate_filename.m - Generate standardized file names from metadata.
function filename = generate_filename(data_type, meta)
switch data_type
    case 'raster'
        required_fields = {'animal_name', 'injection', 'prepost', 'state', 'area', 'align', 'session_idx'};
        check_required_fields(meta, required_fields);

        switch meta.animal_name
            case 'Slayer'
                animal_code = 'Slayer';
            case 'Zeppelin'
                animal_code = 'Zeppelin';
            case 'Emperor'
                animal_code = 'Emperor';
            otherwise
                error('Unknown animal name: %s', meta.animal_name);
        end

        switch meta.injection
            case 'Muscimol'
                injection_code = 'Mus';
            case 'Saline'
                injection_code = 'Sal';
            case 'No injection'
                injection_code = 'Noinj';
            otherwise
                error('Unknown injection type: %s', meta.injection);
        end

        switch meta.prepost
            case 'Pre'
                prepost_code = 'Pre';
            case 'Post'
                prepost_code = 'Post';
            otherwise
                error('Unknown prepost type: %s', meta.prepost);
        end

        switch meta.state
            case 'Task'
                state_code = 'Task';
            case 'RestOpen'
                state_code = 'Open';
            case 'RestClose'
                state_code = 'Closed';
            otherwise
                error('Unknown state type: %s', meta.state);
        end

        switch meta.area
            case 'ACC'
                area_code = 'ACC';
            case 'VLPFC'
                area_code = 'PFC';
            case 'Thalamus'
                area_code = 'Thal';
            case 'Cortex'
                area_code = 'Ctx';
            case 'Full'
                area_code = 'Full';
            otherwise
                error('Unknown area type: %s', meta.area);
        end

        switch meta.align
            case 'None'
                align_code = 'None';
            case 'First'
                align_code = 'First';
            case 'Last'
                align_code = 'Last';
            case 'Random'
                align_code = 'Rand';
            case 'Longest'
                align_code = 'Longest';
            otherwise
                error('Unknown align type: %s', meta.align);
        end
        
        filename = sprintf('raster_%s%s%s%s%s%s%d.mat', animal_code, injection_code, prepost_code, state_code, area_code, align_code, meta.session_idx);
    case 'border'
        required_fields = {'animal_name', 'injection', 'prepost', 'area', 'align', 'session_idx'};
        check_required_fields(meta, required_fields);

        switch meta.animal_name
            case 'Slayer'
                animal_code = 'Slayer';
            case 'Zeppelin'
                animal_code = 'Zeppelin';
            case 'Emperor'
                animal_code = 'Emperor';
            otherwise
                error('Unknown animal name: %s', meta.animal_name);
        end

        switch meta.injection
            case 'Muscimol'
                injection_code = 'Mus';
            case 'Saline'
                injection_code = 'Sal';
            case 'No injection'
                injection_code = 'Noinj';
            otherwise
                error('Unknown injection type: %s', meta.injection);
        end

        switch meta.prepost
            case 'Pre'
                prepost_code = 'Pre';
            case 'Post'
                prepost_code = 'Post';
            otherwise
                error('Unknown prepost type: %s', meta.prepost);
        end

        switch meta.area
            case 'ACC'
                area_code = 'ACC';
            case 'VLPFC'
                area_code = 'PFC';
            case 'Thalamus'
                area_code = 'Thal';
            case 'Cortex'
                area_code = 'Ctx';
            case 'Full'
                area_code = 'Full';
            otherwise
                error('Unknown area type: %s', meta.area);
        end

        switch meta.align
            case 'None'
                align_code = 'None';
            case 'First'
                align_code = 'First';
            case 'Last'
                align_code = 'Last';
            case 'Random'
                align_code = 'Rand';
            case 'Longest'
                align_code = 'Longest';
            otherwise
                error('Unknown align type: %s', meta.align);
        end

        filename = sprintf('border_%s%s%s%s%s%d.mat', animal_code, injection_code, prepost_code, area_code, align_code, meta.session_idx);
    
    case 'sortidx'
        required_fields = {'animal_name', 'injection', 'prepost', 'state', 'area', 'align', 'session_idx', 'criterion'};
        check_required_fields(meta, required_fields);

        switch meta.animal_name
            case 'Slayer'
                animal_code = 'Slayer';
            case 'Zeppelin'
                animal_code = 'Zeppelin';
            case 'Emperor'
                animal_code = 'Emperor';
            otherwise
                error('Unknown animal name: %s', meta.animal_name);
        end

        switch meta.injection
            case 'Muscimol'
                injection_code = 'Mus';
            case 'Saline'
                injection_code = 'Sal';
            case 'No injection'
                injection_code = 'Noinj';
            otherwise
                error('Unknown injection type: %s', meta.injection);
        end

        switch meta.prepost
            case 'Pre'
                prepost_code = 'Pre';
            case 'Post'
                prepost_code = 'Post';
            otherwise
                error('Unknown prepost type: %s', meta.prepost);
        end

        switch meta.state
            case 'Task'
                state_code = 'Task';
            case 'RestOpen'
                state_code = 'Open';
            case 'RestClose'
                state_code = 'Closed';
            otherwise
                error('Unknown state type: %s', meta.state);
        end

        switch meta.area
            case 'ACC'
                area_code = 'ACC';
            case 'VLPFC'
                area_code = 'PFC';
            case 'Thalamus'
                area_code = 'Thal';
            case 'Cortex'
                area_code = 'Ctx';
            case 'Full'
                area_code = 'Full';
            otherwise
                error('Unknown area type: %s', meta.area);
        end

        switch meta.align
            case 'None'
                align_code = 'None';
            case 'First'
                align_code = 'First';
            case 'Last'
                align_code = 'Last';
            case 'Random'
                align_code = 'Rand';
            case 'Longest'
                align_code = 'Longest';
            otherwise
                error('Unknown align type: %s', meta.align);
        end

        switch meta.criterion
            case 'channel'
                criterion_code = 'Channel';
            otherwise
                error('Unknown criterion type: %s', meta.criterion);
        end

        filename = sprintf('sortidx_%s%s%s%s%s%s%s%d.mat', animal_code, injection_code, prepost_code, state_code, area_code, align_code, criterion_code, meta.session_idx);
    
    case 'kernel'
        required_fields = {'kernel_name'};
        check_required_fields(meta, required_fields);
        filename = sprintf('kernel_%s.mat', meta.kernel_name);

    case 'GLM_data'
    case 'GLM_model'
    case 'crossval'

    otherwise
        error('Unknown data type: %s', data_type);
end

function check_required_fields(meta, required_fields)
for i = 1:numel(required_fields)
    if ~isfield(meta, required_fields{i})
        error('Missing required field: %s', required_fields{i});
    end
end