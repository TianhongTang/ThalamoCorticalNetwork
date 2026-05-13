function ax = plot_network(f, J12, J21, err12, err21, errScale, emphasize1, emphasize2)
% plot_network(f, J12, J21, err12, err21, errScale)
%
% J12: m x n, connections from area 2 to area 1
%      J12(i,j): area 2 neuron j -> area 1 neuron i
%
% J21: n x m, connections from area 1 to area 2
%      J21(j,i): area 1 neuron i -> area 2 neuron j
%
% errScale: scaling factor for the error threshold. Default is 1.
%
% emphasize1: If provided, emphasize the specified neurons in area 1.
% emphasize2: If provided, emphasize the specified neurons in area 2.
%
% Only connections satisfying abs(J) > errScale * err are plotted.
% Positive connections are shown in red.
% Negative connections are shown in blue.
%
% For a unidirectional connection, a dot is placed near the receiving end.
% For a bidirectional connection, dots are placed near both ends.
% If a bidirectional connection has opposite signs in the two directions,
% the two halves of the line are colored separately according to the
% connection received by the nearby area.

    if nargin < 6 || isempty(errScale)
        errScale = 1;
    end

    [m, n] = size(J12);

    if ~isequal(size(J21), [n, m])
        error('J21 must be n x m.');
    end

    if ~isequal(size(err12), [m, n])
        error('err12 must be m x n, same as J12.');
    end

    if ~isequal(size(err21), [n, m])
        error('err21 must be n x m, same as J21.');
    end

    % f can be either a figure handle or an axes handle.
    if isgraphics(f, 'axes')
        ax = f;
    elseif isgraphics(f, 'figure')
        ax = axes('Parent', f);
    else
        error('f must be a figure handle or axes handle.');
    end

    cla(ax);
    hold(ax, 'on');

    % Node layout
    xL = 0;
    xR = 1;

    if m == 1
        yL = 0.5;
    else
        yL = linspace(1, 0, m);
    end

    if n == 1
        yR = 0.5;
    else
        yR = linspace(1, 0, n);
    end

    red  = [0.85, 0.10, 0.10];
    blue = [0.10, 0.25, 0.90];

    lineWidth = 0.5;
    nodeSize = 2;
    dotSize = 2;
    dotFrac = 0;
    nodeShift = 0.02; % Shift nodes slightly away from the center to avoid overlap with connection lines

    % % Vertical guide lines for the two areas
    % plot(ax, [xL xL], [min(yL) max(yL)], '-', ...
    %     'Color', [0.75 0.75 0.75], 'LineWidth', 1);

    % plot(ax, [xR xR], [min(yR) max(yR)], '-', ...
    %     'Color', [0.75 0.75 0.75], 'LineWidth', 1);

    % Draw connections
    for i = 1:m
        for j = 1:n

            % area 2 neuron j -> area 1 neuron i
            sig12 = abs(J12(i,j)) > errScale * err12(i,j);

            % area 1 neuron i -> area 2 neuron j
            sig21 = abs(J21(j,i)) > errScale * err21(j,i);

            if ~sig12 && ~sig21
                continue;
            end

            x1 = xL;
            y1 = yL(i);

            x2 = xR;
            y2 = yR(j);

            if sig12 && sig21
                % Bidirectional connection
                xm = (x1 + x2) / 2;
                ym = (y1 + y2) / 2;

                % Left half is colored by the connection received by area 1.
                cLeft = sign_color(J12(i,j), red, blue);

                % Right half is colored by the connection received by area 2.
                cRight = sign_color(J21(j,i), red, blue);

                plot(ax, [x1 xm], [y1 ym], '-', ...
                    'Color', cLeft, 'LineWidth', lineWidth);

                plot(ax, [xm x2], [ym y2], '-', ...
                    'Color', cRight, 'LineWidth', lineWidth);

                % Dots near both receiving ends
                draw_dot(ax, x1, y1, x2, y2, dotFrac, dotSize, cLeft);
                draw_dot(ax, x1, y1, x2, y2, 1 - dotFrac, dotSize, cRight);

            elseif sig12
                % Only area 2 -> area 1
                c = sign_color(J12(i,j), red, blue);

                plot(ax, [x1 x2], [y1 y2], '-', ...
                    'Color', c, 'LineWidth', lineWidth);

                % Receiver is area 1, on the left side.
                draw_dot(ax, x1, y1, x2, y2, dotFrac, dotSize, c);

            elseif sig21
                % Only area 1 -> area 2
                c = sign_color(J21(j,i), red, blue);

                plot(ax, [x1 x2], [y1 y2], '-', ...
                    'Color', c, 'LineWidth', lineWidth);

                % Receiver is area 2, on the right side.
                draw_dot(ax, x1, y1, x2, y2, 1 - dotFrac, dotSize, c);
            end
        end
    end

    % Draw neuron nodes on top of connection lines
    plot(ax, xL * ones(1,m) - nodeShift, yL, 'ko', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', nodeSize);

    plot(ax, xR * ones(1,n) + nodeShift, yR, 'ko', ...
        'MarkerFaceColor', 'k', ...
        'MarkerSize', nodeSize);

    text(ax, xL, 1.08, 'ACC', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');

    text(ax, xR, 1.08, 'VLPFC', ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold');
    
    % Emphasize specified neurons by drawing a triangle pointing to them
    if nargin >= 7 && ~isempty(emphasize1)
        for idx = emphasize1
            if idx >= 1 && idx <= m
                plot(ax, xL - 0.05, yL(idx), '>', ...
                    'MarkerFaceColor', 'm', ...
                    'MarkerEdgeColor', 'm', ...
                    'MarkerSize', nodeSize + 2);
            end
        end
    end
    if nargin >= 8 && ~isempty(emphasize2)
        for idx = emphasize2
            if idx >= 1 && idx <= n
                plot(ax, xR + 0.05, yR(idx), '<', ...
                    'MarkerFaceColor', 'm', ...
                    'MarkerEdgeColor', 'm', ...
                    'MarkerSize', nodeSize + 2);
            end
        end
    end

    axis(ax, 'equal');
    axis(ax, 'off');
    xlim(ax, [-0.15, 1.15]);
    ylim(ax, [-0.10, 1.15]);

end

function c = sign_color(v, red, blue)
% Return the line color based on the sign of the connection strength.

    if v > 0
        c = red;
    else
        c = blue;
    end

end

function draw_dot(ax, x1, y1, x2, y2, t, dotSize, color)
% Draw a dot at a fractional position along the connection line.

    xd = x1 + t * (x2 - x1);
    yd = y1 + t * (y2 - y1);

    plot(ax, xd, yd, 'o', ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', color, ...
        'MarkerSize', dotSize);

end