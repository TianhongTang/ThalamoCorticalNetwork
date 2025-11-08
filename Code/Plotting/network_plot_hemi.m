function network_plot_hemi(ax, mat, err, borders, node_colors, area_names)
    % Plot directed connections between areas arranged on a circle.

    area_num = numel(area_names);
    if numel(borders) ~= area_num + 1
        error('network_plot_hemi:InvalidBorders', 'borders must have length area_num + 1');
    end

    N = size(mat, 1);
    if size(mat, 2) ~= N || size(err, 1) ~= N || size(err, 2) ~= N
        error('network_plot_hemi:InvalidMatrixSize', 'mat and err must be square N-by-N matrices');
    end

    if borders(end) ~= N + 1
        error('network_plot_hemi:InvalidBorders', 'borders(end) must equal N + 1');
    end

    if size(node_colors, 1) < N
        error('network_plot_hemi:InvalidNodeColors', 'node_colors must have at least N rows');
    end

    cla(ax);
    hold(ax, 'on');

    start_angle = pi / 2; % top of the circle
    area_angle = 2 * pi / area_num;

    x = zeros(1, N);
    y = zeros(1, N);
    area_center_angles = zeros(1, area_num);
    area_centers = zeros(area_num, 2);

    for area_idx = 1:area_num
        node_start = borders(area_idx);
        node_end = borders(area_idx + 1) - 1;
        if node_end < node_start
            continue;
        end
    area_nodes = node_start:node_end;
        node_count = numel(area_nodes);

        area_angle_start = start_angle - (area_idx - 1) * area_angle;
        area_angle_end = area_angle_start - area_angle;
        center_angle = area_angle_start - area_angle / 2;
        area_center_angles(area_idx) = center_angle;

        if node_count == 1
            node_angles = center_angle;
        else
            node_angles = linspace(area_angle_start - area_angle / (node_count + 1), ...
                                   area_angle_end + area_angle / (node_count + 1), node_count);
        end

        base_x = cos(node_angles);
        base_y = sin(node_angles);
        offset_vec = 0.5 * [cos(center_angle), sin(center_angle)];
        area_centers(area_idx, :) = offset_vec;

        x(area_nodes) = base_x + offset_vec(1);
        y(area_nodes) = base_y + offset_vec(2);

    end

    for i = 1:N
        for j = 1:N
            if mat(i, j) > err(i, j)
                quiver(ax, x(i), y(i), x(j) - x(i), y(j) - y(i), 0, 'Color', 'r', 'LineStyle', '-');
            elseif mat(i, j) < -err(i, j)
                quiver(ax, x(i), y(i), x(j) - x(i), y(j) - y(i), 0, 'Color', 'b', 'LineStyle', '-');
            end
        end
    end

    for area_idx = 1:area_num
        node_start = borders(area_idx);
        node_end = borders(area_idx + 1) - 1;
        if node_end < node_start
            continue;
        end
        area_nodes = node_start:node_end;
        area_offset = node_start - 1;
        for node_id = area_nodes
            plot(ax, x(node_id), y(node_id), 'o', 'MarkerSize', 3, ...
                'MarkerFaceColor', node_colors(node_id, :), 'MarkerEdgeColor', 'k');
            % text(ax, x(node_id) * 1.1, y(node_id) * 1.1, ...
            %     num2str(node_id - area_offset), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        end
    end

    for area_idx = 1:area_num
        center_angle = area_center_angles(area_idx);
        offset_vec = area_centers(area_idx, :);
        label_point = offset_vec + 1.3 * [cos(center_angle), sin(center_angle)];
        text(ax, label_point(1), label_point(2), area_names{area_idx}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontWeight', 'bold');
    end

    hold(ax, 'off');
    axis(ax, 'equal');
    ax.XLim = [-2, 2];
    ax.YLim = [-2, 2];
    xticks(ax, []);
    yticks(ax, []);
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
end