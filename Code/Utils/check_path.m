function check_path(path)
    if ~exist(path, 'dir')
        mkdir(path);
    end
end