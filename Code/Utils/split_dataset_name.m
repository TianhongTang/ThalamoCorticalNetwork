%% split_dataset_name.m - Split dataset name (old) into animal name and injection type.
function [animal_name, injection] = split_dataset_name(dataset_name)
    if contains(dataset_name, 'Slayer')
        animal_name = 'Slayer';
    elseif contains(dataset_name, 'Zeppelin')
        animal_name = 'Zeppelin';
    elseif contains(dataset_name, 'Emperor')
        animal_name = 'Emperor';
    end
    if contains(dataset_name, 'Mus')
        injection = 'Muscimol';
    elseif contains(dataset_name, 'Sal')
        injection = 'Saline';
    elseif contains(dataset_name, 'Noinj')
        injection = 'No injection';
    end
end
%% split_dataset_name.m - Split dataset name (old) into animal name and injection type.
function [animal_name, injection] = split_dataset_name(dataset_name)
    if contains(dataset_name, 'Slayer')
        animal_name = 'Slayer';
    elseif contains(dataset_name, 'Zeppelin')
        animal_name = 'Zeppelin';
    elseif contains(dataset_name, 'Emperor')
        animal_name = 'Emperor';
    end
    if contains(dataset_name, 'Mus')
        injection = 'Muscimol';
    elseif contains(dataset_name, 'Sal')
        injection = 'Saline';
    elseif contains(dataset_name, 'Noinj')
        injection = 'No injection';
    end
end