function [pos, neg, total] = J_count(J, err, kernel_idx, borders, selected_areas)
    N = size(J, 1);
    J_mat = J(:, (2+(kernel_idx-1)*N):(1+kernel_idx*N));
    err_mat = err(:, (2+(kernel_idx-1)*N):(1+kernel_idx*N));

    area_ranges = cell(1, numel(borders)-1);
    for i=1:(numel(borders)-1)
        area_ranges{i} = borders(i):(borders(i+1)-1);
    end

    pos = 0; 
    neg = 0;
    total = 0;
    for area_idx = 1:numel(selected_areas)
        selected_area = selected_areas{area_idx};
        area_i = selected_area(1);
        area_j = selected_area(2);

        J_area     = J_mat(area_ranges{area_i}, area_ranges{area_j});
        err_area   = err_mat(area_ranges{area_i}, area_ranges{area_j});
        pos_area   = sum(J_area(:) > err_area(:));
        neg_area   = sum(J_area(:) < -err_area(:));
        total_area = numel(J_area);
        pos        = pos + pos_area;
        neg        = neg + neg_area;
        total      = total + total_area;
    end
end