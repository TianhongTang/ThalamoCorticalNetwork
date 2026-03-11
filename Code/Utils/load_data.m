function data = load_data(root, data_type, param)

    switch data_type
        case 'raster'
            path = fullfile(root, 'Data', 'Working', 'Spikes', sprintf('%s.mat', param));
            load(path, 'data');
        case ''
            path = fullfile(root, 'Data', 'Working', 'LFP', sprintf('%s.mat', param));
            load(path, 'data');
        otherwise
            error('Unknown data type: %s', data_type);
    end
end