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

% Define different discretization levels
n_lines_values = [10, 20, 30, 40, 50, 60];
num_discretizations = length(n_lines_values);

% Reference value for losses
reference_loss = 3.170205031911107e-4; % (W)

% Store the losses for each discretization level
losses = zeros(num_discretizations, 1);
errors = zeros(num_discretizations, 1);

for k = 1:num_discretizations
    n_lines = n_lines_values(k);

    %% Generate mesh for the circular conductor
    [rect_centers, rect_areas, ~] = generateCircularConductorMesh(radius, n_lines);
    N = length(rect_areas);         % Total number of filaments

    %% Calculate Resistance
    Resistance = diag(rho * l ./ rect_areas);

    %% Calculate Inductance 
    xi = rect_centers / R; % Using centers for distance calculation
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
    rho_i = sqrt(rect_areas / pi) / radius;
    log_rho_i = log(rho_i);
    norms_squared = norm_xi.^2;
    G_self = -mu0 * l / (2 * pi) * (log_rho_i - 0.5 * log((1 - norms_squared).^2 + norms_squared.^2 .* rho_i.^2));
    Inductance(1:size(G_self,1)+1:end) = G_self;

    %% Combined impedance matrix
    Z_Lambda = Resistance + 1i * omega * Inductance;

    %% Connectivity matrix and current setup
    I = 1; % Current through the conductor [A]
    C = ones(N, 1); % All filaments belong to a single conductor

    %% Compute terminal impedance
    Z_terminal = inv(C' * inv(Z_Lambda) * C);

    %% Compute the current in each filament
    i_lambda = Z_Lambda \ (C * Z_terminal * I);

    %% Compute losses in the conductor
    P_total = sum(diag(Resistance) .* abs(i_lambda).^2 * 0.5);

    losses(k) = P_total;
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


%% Plot the Composed Mesh for Each Discretization Level

figure;
numPlots = numel(n_lines_values);
% For 6 plots use this scheme. Adapt if needed
nrows = 3;
ncols = 3;

for k = 1:numPlots
    n_lines_plot = n_lines_values(k);
    [rect_centers, rect_areas, rect_vertices] = generateCircularConductorMesh(radius, n_lines_plot);
    
    subplot(nrows, ncols, k);
    hold on;
    
    % Conductor
    theta = linspace(0, 2*pi, 200);
    plot(radius * cos(theta), radius * sin(theta), 'b', 'LineWidth', 2);
    
    % Rectangles 
    for i = 1:length(rect_vertices)
        patch(rect_vertices{i}(:,1), rect_vertices{i}(:,2), 'r', ...
              'FaceAlpha', 0.3, 'EdgeColor', 'k');
    end
    
    % Centers of each element: each rectangle
    scatter(rect_centers(:,1), rect_centers(:,2), 25, 'k', 'filled');
    
    title(sprintf('n-lines = %d', n_lines_plot));
    xlabel('x (m)');
    ylabel('y (m)');
    axis equal;
    grid on;
    hold off;
end

sgtitle('Circular Conductor Mesh Visualizations for Different Discretizations');

