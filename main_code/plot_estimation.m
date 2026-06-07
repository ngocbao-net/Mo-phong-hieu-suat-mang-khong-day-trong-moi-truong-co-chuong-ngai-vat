function plot_estimation(L, R, total_shadow, data, mesh, estimation, title_string, buildings)
    plot([L,-L,-L,L,L], [L,L,-L,-L,L], '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
    hold on;
    xlim((1.1)*[-R R]); ylim((1.1)*[-R R]); pbaspect([1 1 1]);

    title(title_string, 'Interpreter', 'none');

    % Border of working disk.
    plot(R*cos(linspace(0,2*pi,1000)), R*sin(linspace(0,2*pi,1000)), 'k-', 'LineWidth', 1.2);

    % Shadow regions.
    plot(total_shadow, 'FaceColor', 'black', 'FaceAlpha', 0.25, 'LineStyle', 'none');

    % Buildings / obstacle edges.
    if exist('buildings','var') && ~isempty(buildings)
        for i = 1:size(buildings,1)
            plot([buildings(i,1),buildings(i,3)],[buildings(i,2),buildings(i,4)],'k-','LineWidth',1.0);
        end
    end

    % Data points.
    if ~isempty(data)
        plot(data(:,1), data(:,2), 'xk', 'MarkerSize', 5, 'LineWidth', 1.0);
    end

    % Estimated mesh labels.
    if ~isempty(mesh) && ~isempty(estimation)
        estimation = logical(estimation(:));
        scatter(mesh(estimation,1), mesh(estimation,2), 28, 'o', ...
            'MarkerEdgeColor', [0 0.4470 0.7410], 'LineWidth', 0.9);
        scatter(mesh(~estimation,1), mesh(~estimation,2), 28, 'o', ...
            'MarkerEdgeColor', [0.8500 0.3250 0.0980], 'LineWidth', 0.9);
    end

    % Base station marker.
    plot(0, 0, 'p', 'MarkerSize', 12, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0);

    % White-out outside the disk.
    in_circle_x = R*cos(linspace(0,2*pi,1000));
    in_circle_y = R*sin(linspace(0,2*pi,1000));
    inner_circle = polyshape(in_circle_x, in_circle_y);
    out_circle_x = 2*R*cos(linspace(0,2*pi,1000));
    out_circle_y = 2*R*sin(linspace(0,2*pi,1000));
    outer_circle = polyshape(out_circle_x, out_circle_y);
    exterior_polygon = subtract(outer_circle, inner_circle);
    plot(exterior_polygon,'FaceColor','white', 'FaceAlpha', 1, 'LineWidth', 1.2);

    xlim((1.1)*[-R R]); ylim((1.1)*[-R R]); pbaspect([1 1 1]);
    box on; grid on;
    hold off;
end
