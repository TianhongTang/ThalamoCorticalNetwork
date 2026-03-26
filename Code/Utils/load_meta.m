function meta = load_meta(root, mode)
    if ~exist('mode', 'var')
        mode = 'struct';
    end

    path = fullfile(root, 'Data', 'Working', 'Meta', 'metadata.mat');
    load(path, 'metadata');
    meta = metadata;
    if strcmp(mode, 'table')
        for field_name = fieldnames(meta)'
            fname = field_name{1};
            meta.(fname) = struct2table(meta.(fname));
        end
    end
end