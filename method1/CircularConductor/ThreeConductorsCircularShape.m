clear all; close all; clc;

format short e
%% Constants and Parameters
mu0 = 4 * pi * 1e-7;            % Permeability of free space (H/m)
R = 0.1;                        % Reference radius (m)
rho = 1/5.998e7;                % Resistivity of copper (Ohm*m)
f = 500;                       % Frequency (Hz)
omega = 2 * pi * f;             % Angular frequency (rad/s)
l = 1;                          % Length of the conductor (m): scale the results

radius = 5e-3;                  % Radius of the circular conductor
numConductors = 3;              % Number of conductors
n_lines = 50;                   % Number of lines for discretization

totalPositions = [];
totalCenters = [];
totalAreas = [];

for k = 1:numConductors
    [rect_centers, rect_areas, ~] = generateCircularConductorMesh(radius, n_lines);
    
    % Shift the positions of the conductors
    xShift = (k - 2) * (5 * radius); % Shift along the x direction: same y
    rect_centers(:, 1) = rect_centers(:, 1) + xShift;
    
    totalCenters = [totalCenters; rect_centers];
    totalAreas = [totalAreas; rect_areas];
end

N = length(totalAreas);  

%% Calculate Resistance
Resistance = diag(rho * l ./ totalAreas);

%% Calculate Inductance
xi = totalCenters / R; 
Inductance = zeros(N, N);

[X1, X2] = meshgrid(xi(:,1), xi(:,1));
[Y1, Y2] = meshgrid(xi(:,2), xi(:,2));

distances = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);
norm_xi = sqrt(xi(:,1).^2 + xi(:,2).^2);

% Log values fot mutual inductances
norm_xi_sq = norm_xi.^2;
norms_product = norm_xi * norm_xi';
G_mutual_log = log(distances) - 0.5 * log(norms_product.^2 - 2 * (xi * xi') + 1);
G_mutual_log(distances == 0) = 0; % Correct -Inf because of log(0)

% Mutual Inductances -> excluding self-inductance along the diagonal
Inductance(~eye(N)) = -mu0 * l / (2 * pi) * G_mutual_log(~eye(N));

% Self Inductances
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

disp('Total Losses in the Conductor (W):');
disp(P_total);

%% Plot Resistance, Inductance, and Impedance matrices
figure;
subplot(1, 3, 1);
imagesc(Resistance);
colorbar;
title('Resistance Matrix');
xlabel('Filament Index');
ylabel('Filament Index');

subplot(1, 3, 2);
imagesc(abs(Inductance));
colorbar;
title('Inductance Matrix (Magnitude)');
xlabel('Filament Index');
ylabel('Filament Index');

subplot(1, 3, 3);
imagesc(abs(Z_Lambda));
colorbar;
title('Impedance Matrix');
xlabel('Filament Index');
ylabel('Filament Index');

%% Compute the current densities
current_density = abs(i_lambda) ./ totalAreas;

%% Plot Current Distribution and Current Density for each conductor
max_current = max(abs(i_lambda));
max_density = max(current_density);

figure;
for k = 1:numConductors
    % Current Distribution (A)
    subplot(2, numConductors, k); % In the first row
    idx = (k-1)*(N/numConductors) + (1:(N/numConductors));
    scatter(totalCenters(idx, 1), totalCenters(idx, 2), 20, abs(i_lambda(idx)), 'filled');
    colorbar;
    caxis([0 max_current]); 
    axis equal tight;
    title(sprintf('Current Distribution in Conductor %d (A)', k));
    xlabel('X Position (m)');
    ylabel('Y Position (m)');

    % Current Density (A/m^2)
    subplot(2, numConductors, numConductors + k); % In the second row
    scatter(totalCenters(idx, 1), totalCenters(idx, 2), 20, current_density(idx), 'filled');
    colorbar;
    caxis([0 max_density]);  
    axis equal tight;
    title(sprintf('Current Density in Conductor %d (A/m^2)', k));
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
end
