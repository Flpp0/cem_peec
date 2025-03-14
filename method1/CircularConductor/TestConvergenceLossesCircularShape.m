clear all; close all; clc;

format short e
%% Constants and Parameters
mu0 = 4 * pi * 1e-7;            % Permeability of free space (H/m)
R = 0.1;                        % Reference radius (m)
rho = 1/5.998e7;                % Resistivity of copper (Ohm*m)
f = 50e2;                       % Frequency (Hz)
omega = 2 * pi * f;             % Angular frequency (rad/s)
l = 1;                          % Length of the conductor (m): scale the results

radius = 5e-3;                  % Radius of the circular conductor
numConductors = 3;              % Number of conductors

% Define different discretization levels
n_lines_values = [10:10:60];
num_discretizations = length(n_lines_values);

reference_loss = 1.057375767873015e-3;

% Store the losses for each discretization level
losses = zeros(num_discretizations, 1);
errors = zeros(num_discretizations, 1);

for k = 1:num_discretizations
    n_lines = n_lines_values(k);

    totalCenters = [];
    totalAreas = [];

    for j = 1:numConductors
        %% Generate mesh for the circular conductor
        [rect_centers, rect_areas, ~] = generateCircularConductorMesh(radius, n_lines);
        
        % Shift the positions of the conductors
        xShift = (j - 2) * (5 * radius); % Adjust this factor based on desired spacing
        rect_centers(:, 1) = rect_centers(:, 1) + xShift;
        
        totalCenters = [totalCenters; rect_centers];
        totalAreas = [totalAreas; rect_areas];
    end

    N = length(totalAreas);         % Total number of filaments

    %% Calculate Resistance
    Resistance = diag(rho * l ./ totalAreas);

    %% Calculate Inductance
    xi = totalCenters / R; % Using centers for distance calculation
    Inductance = zeros(N, N);

    % Vectorized distance calculation
    [X1, X2] = meshgrid(xi(:,1), xi(:,1));
    [Y1, Y2] = meshgrid(xi(:,2), xi(:,2));
    distances = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);
    norm_xi = sqrt(xi(:,1).^2 + xi(:,2).^2);

    % Precompute log values for mutual inductance
    norm_xi_sq = norm_xi.^2;
    norms_product = norm_xi * norm_xi';
    G_mutual_log = log(distances) - 0.5 * log(norms_product.^2 - 2 * (xi * xi') + 1);
    G_mutual_log(distances == 0) = 0; % avoid log(0)

    % Compute mutual inductance, excluding self-inductance diagonal
    Inductance(~eye(N)) = -mu0 * l / (2 * pi) * G_mutual_log(~eye(N));

    % Compute self-inductance
    rho_i = sqrt(totalAreas / pi) / R;
    log_rho_i = log(rho_i);
    norms_squared = norm_xi.^2;
    G_self = -mu0 * l / (2 * pi) * (log_rho_i - 0.5 * log((1 - norms_squared).^2 + norms_squared.^2 .* rho_i.^2));
    Inductance(1:size(G_self,1)+1:end) = G_self;

    %% Combined impedance matrix
    Z_Lambda = Resistance + 1i * omega * Inductance;

    %% Connectivity matrix and current setup
    I = [1; 0.8; 1.2];  %[A]
    C = kron(eye(numConductors), ones(N/numConductors,1)); % Assign each filament to its conductor

    %% Compute terminal impedance
    Z_terminal = inv(C' * inv(Z_Lambda) * C);

    %% Compute the current in each filament
    i_lambda = Z_Lambda \ (C * Z_terminal * I);

    %% Compute losses in all the conductors
    P_total = real(I' * Z_terminal * conj(I)) * 0.5;

    % Store the computed losses
    losses(k) = P_total;

    % Compute the absolute error
    errors(k) = abs(P_total - reference_loss);
end

% Plot the losses as a function of the number of lines
figure;
subplot(2, 1, 1);
plot(n_lines_values, losses, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Lines for Discretization');
ylabel('Total Losses (W)');
title('Convergence of Total Losses with Increasing Discretization');
grid on;

% Plot the errors as a function of the number of lines
subplot(2, 1, 2);
plot(n_lines_values, errors, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Number of Lines for Discretization');
ylabel('Absolute Error (W)');
title('Error in Total Losses with Increasing Discretization');
grid on;

%% Plot the Mesh for Different Discretizations for All Conductors

numPlots = numel(n_lines_values);

% For each discretization
for k = 1:numPlots

    n_lines = n_lines_values(k); % level of refinement

    figure;
    hold on;
    
    % For each conductor in the geometry
    for j = 1:numConductors 

        [rect_centers, ~, rect_vertices] = generateCircularConductorMesh(radius, n_lines);
        
        % Shift the conductor at its place
        xShift = (j - 2) * (5 * radius);
        rect_centers(:, 1) = rect_centers(:, 1) + xShift;

        for m = 1:length(rect_vertices)
            rect_vertices{m}(:,1) = rect_vertices{m}(:,1) + xShift;
        end
        
        % Circular conductor
        theta = linspace(0, 2*pi, 200);
        plot(radius * cos(theta) + xShift, radius * sin(theta), 'b', 'LineWidth', 2);
        
        % Elements of the mesh
        for m = 1:length(rect_vertices)
            patch(rect_vertices{m}(:,1), rect_vertices{m}(:,2), 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
        end
        
        % Centers
        scatter(rect_centers(:,1), rect_centers(:,2), 10, 'k', 'filled');

    end

    title(sprintf('Mesh for n-lines = %d', n_lines));
    xlabel('x (m)');
    ylabel('y (m)');
    axis equal;
    grid on;
    hold off;

end