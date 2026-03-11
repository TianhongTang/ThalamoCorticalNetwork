function meta = load_meta(root)
    path = fullfile(root, 'Data', 'Working', 'Meta', 'metadata.mat');
    load(path, 'meta');
end