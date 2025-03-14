function [X, Y, Areas, dw, dt] = generateNonUniformMesh_Rect(width, height, numPointsX, numPointsY)
    % Non-uniform mesh using a sin function 

    x = nonUniformGrid(-width/2, width/2, numPointsX);
    y = nonUniformGrid(-height/2, height/2, numPointsY);

    [X, Y] = meshgrid(x, y);

    % Calculate the width (dw) and thickness (dt) of each sub-element
    dw = diff(x); % Widths of sub-elements in the x-direction
    dt = diff(y); % Thicknesses of sub-elements in the y-direction

    % Create a grid to multiply in each element of the grid dx and dy
    [DX, DY] = meshgrid(dw, dt); 
    Areas = DX .* DY;

    function grid = nonUniformGrid(minVal, maxVal, numPoints)
        t = linspace(0, pi, numPoints);
        % (sin(t - pi/2) + 1) / 2 -> [0 , 1]
        % (maxVal - minVal) * (sin(t - pi/2) + 1) / 2 -> [0, width (or height)]
        % minVal + (maxVal - minVal) * (sin(t - pi/2) + 1) / 2 -> [-width/2 (or height/2), width/2 (or height/2)]
        grid = minVal + (maxVal - minVal) * (sin(t - pi/2) + 1) / 2;
    end
end
