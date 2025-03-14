
%% generateUniformMesh function
function [X, Y, Areas] = generateUniformMesh(width, height, nx, ny)
    x = linspace(0, width, nx);
    y = linspace(0, height, ny);
    [X, Y] = meshgrid(x, y);
    dx = diff(x);
    dy = diff(y);
    [DX, DY] = meshgrid(dx, dy);
    Areas = DX .* DY;
end