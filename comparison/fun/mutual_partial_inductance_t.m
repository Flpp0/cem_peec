function L = mutual_partial_inductance_t(i, j, x_center, y_center, mu0, l)
    % Get the center coordinates of elements i and j
    xi = x_center(i);
    yi = y_center(i);
    xj = x_center(j);
    yj = y_center(j);

    % Distance between centers
    dij = sqrt((xi - xj)^2 + (yi - yj)^2);

    % Calculate the mutual inductance
    if dij == 0
        L = 0; % Self inductance is handled separately
    else
        L = mu0 * l / (2 * pi) * (log(l / dij + sqrt((l / dij)^2 + 1)) - sqrt(1 + (dij / l)^2) + dij / l);
    end
end