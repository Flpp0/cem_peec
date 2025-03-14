function [rect_centers, rect_areas, rect_vertices] = generateCircularConductorMesh(radius, n_lines)
    rect_centers = [];
    rect_areas = [];
    rect_vertices = {}; 

    theta = linspace(0, pi, n_lines + 2);

    x_spacing = radius * cos(theta); % Vertical lines
    y_spacing = radius * cos(theta); % Horizontal lines

    % Remove the boundary points
    x_spacing = x_spacing(2:end-1);
    y_spacing = y_spacing(2:end-1);

    [x_grid, y_grid] = meshgrid(x_spacing, y_spacing);

    % Filter points inside the circle
    mask = (x_grid.^2 + y_grid.^2) <= radius^2;
    x_grid = x_grid(mask);
    y_grid = y_grid(mask);

    boundary_x = [];
    boundary_y = [];

    % For each x line (vertical line) find the upper and lower points
    unique_x = unique(x_grid);
    for i = 1:length(unique_x)
        x_val = unique_x(i);
        y_vals = y_grid(x_grid == x_val);
        if ~isempty(y_vals)
            boundary_x = [boundary_x; x_val; x_val];
            boundary_y = [boundary_y; max(y_vals); min(y_vals)];
        end
    end

    % For each y line (horizontal line) find the most left and most right points
    unique_y = unique(y_grid);
    for i = 1:length(unique_y)
        y_val = unique_y(i);
        x_vals = x_grid(y_grid == y_val);
        if ~isempty(x_vals)
            boundary_x = [boundary_x; max(x_vals); min(x_vals)];
            boundary_y = [boundary_y; y_val; y_val];
        end
    end

    boundary_points = unique([boundary_x, boundary_y], 'rows'); % Remove duplicates if present
    boundary_x = boundary_points(:,1);
    boundary_y = boundary_points(:,2);

    % Order the boundary points counterclockwise
    angles = atan2(boundary_y, boundary_x);
    [~, order] = sort(angles);
    boundary_x = boundary_x(order);
    boundary_y = boundary_y(order);

    % Find rectangles over the grid
    for i = 1:length(x_spacing)-1 % -1 because we consider successive points
        for j = 1:length(y_spacing)-1 % -1 because we consider successive points
            x1 = x_spacing(i);
            x2 = x_spacing(i+1);
            y1 = y_spacing(j);
            y2 = y_spacing(j+1);

            % All corners must be inside the circle
            if (x1^2 + y1^2 <= radius^2) && (x1^2 + y2^2 <= radius^2) && ...
               (x2^2 + y1^2 <= radius^2) && (x2^2 + y2^2 <= radius^2)

                rect_x = [x1, x1, x2, x2];
                rect_y = [y1, y2, y2, y1];

                % Compute and store center
                center_x = mean(rect_x);
                center_y = mean(rect_y);
                rect_centers = [rect_centers; center_x, center_y];

                % Compute and store area
                rect_area = polyarea(rect_x, rect_y);
                rect_areas = [rect_areas; rect_area];
                
                % Store vertices in a cell array
                rect_vertices{end+1} = [rect_x(:), rect_y(:)];
            end
        end
    end

    % For the refining mesh step (outer boundary) - similar changes:
    num_outermost_points = length(boundary_x);
    for k = 1:num_outermost_points
        idx1 = k;
        if k < num_outermost_points
            idx2 = k + 1;
        else
            idx2 = 1;
        end

        x1 = boundary_x(idx1);
        y1 = boundary_y(idx1);
        x2 = boundary_x(idx2);
        y2 = boundary_y(idx2);

        mid_x = (x1 + x2) / 2;
        mid_y = (y1 + y2) / 2;
        slope = (y2 - y1) / (x2 - x1);

        if slope ~= 0
            perpendicular_slope = -1 / slope;
        else
            perpendicular_slope = Inf;
        end

        theta_perpendicular = atan(perpendicular_slope);
        new_x = radius * cos(theta_perpendicular);
        new_y = radius * sin(theta_perpendicular);

        if abs(new_x) == 1 || abs(new_y) == 1
            continue
        end

        if dot([new_x, new_y], [mid_x, mid_y]) < 0
            new_x = -new_x;
            new_y = -new_y;
        end

        % Refining the mesh: define new rectangle vertices based on quadrant
        if mid_x > 0 && mid_y > 0 % First quadrant
            new_xA = x2;
            new_yA = new_y;
            new_xB = x2;
            new_yB = y1;
            new_xC = new_x;
            new_yC = y1;
            rect_x = [new_xA, new_xB, new_xC, new_x];
            rect_y = [new_yA, new_yB, new_yC, new_y];
        elseif mid_x < 0 && mid_y > 0 % Second quadrant
            new_xB = new_x;
            new_yB = y2;
            new_xC = x1;
            new_yC = y2;
            new_xD = x1;
            new_yD = new_y;
            rect_x = [new_x, new_xB, new_xC, new_xD];
            rect_y = [new_y, new_yB, new_yC, new_yD];
        elseif mid_x < 0 && mid_y < 0 % Third quadrant
            new_xA = new_x;
            new_yA = y1;
            new_xC = x2;
            new_yC = new_y;
            new_xD = x2;
            new_yD = y1;
            rect_x = [new_xA, new_x, new_xC, new_xD];
            rect_y = [new_yA, new_y, new_yC, new_yD];
        elseif mid_x > 0 && mid_y < 0 % Fourth quadrant
            new_xA = x1;
            new_yA = y2;
            new_xB = x1;
            new_yB = new_y;
            new_xD = new_x;
            new_yD = y2;
            rect_x = [new_xA, new_xB, new_x, new_xD];
            rect_y = [new_yA, new_yB, new_y, new_yD];
        end

        center_x = mean(rect_x);
        center_y = mean(rect_y);
        rect_centers = [rect_centers; center_x, center_y];

        rect_area = polyarea(rect_x, rect_y);
        rect_areas = [rect_areas; rect_area];
        
        % Store the vertices for the refined rectangle
        rect_vertices{end+1} = [rect_x(:), rect_y(:)];
    end
end
