function save_meta(root, meta)
    path = fullfile(root, 'Data', 'Working', 'Meta', 'metadata.mat');
    save(path, 'meta');
end