clear all; close all; clc;

%% Constants and Parameters
mu0 = 4 * pi * 1e-7;             % Permeability of free space (H/m)
R = 0.1;                         % Reference radius (m)
rho = 1.68e-8;                   % Resistivity of copper (Ohm*m)
f = 50e2;                        % Frequency (Hz)
omega = 2 * pi * f;              % Angular frequency (rad/s)
l = 1;                           % Length of the conductor (m)

reference_losses = 4.808550807587721e-3; % From COMSOL: see ThreeConductors.mph

a = 0.003;                       % Width of the conductor (m)
b = 0.002;                       % Height of the conductor (m)
nx = 50;                         % Number of filaments along x
ny = 50;                         % Number of filaments along y

dx = a / nx;
dy = b / ny;

r_hat = sqrt(a * b / pi);        % Equivalent radius

numConductors = 3;               % Number of conductors
N = nx * ny;                     % Total number of filaments per conductor

%% Matrices
Resistance = rho * l / (dx * dy) * eye(N * numConductors); % Resistance per unit length
Inductance = zeros(N * numConductors, N * numConductors);

%% Filament positions and inductance calculation
positions = zeros(N * numConductors, 2); 

xSpacing = 5 * a; 
for k = 1:numConductors
    [X, Y] = meshgrid(linspace(-a/2 + dx/2, a/2 - dx/2, nx), linspace(-b/2 + dy/2, b/2 - dy/2, ny));
    xShift = (k - 2) * xSpacing; % Shift x positions for each conductor
    positions((k-1)*N + 1:k*N, :) = [X(:) + xShift, Y(:)];
end

% Calculate mutual and self inductances using original formulas
for i = 1:N * numConductors
    xi = positions(i,:) / R;
    for j = 1:N * numConductors
        xj = positions(j,:) / R;
        distance = norm(xi - xj);
        xi_dot_xj = dot(xi, xj);

        if i ~= j
            G_mutual = log(distance) - 0.5 * log(norm(xi)^2 * norm(xj)^2 - 2 * xi_dot_xj + 1);
            Inductance(i, j) = -mu0 * l / (2 * pi) * G_mutual;
        end
    end
    rho_i = r_hat / R;
    G_self = log(rho_i) - 0.5 * log((1 - norm(xi)^2))^2 + (norm(xi)^2 * rho_i)^2;
    Inductance(i, i) = -mu0 * l / (2 * pi) * G_self;
end

%% Combined impedance matrix
Z_Lambda = Resistance + 1i * omega * Inductance;

%% Current for each conductor
I = [1; 0.8; 1.2]; %[A]

%% Connectivity matrix 
C = kron(eye(numConductors), ones(N,1)); 

%% Terminal current and voltage calculations
U = Z_Lambda * C * I; % Terminal voltage

%% Terminal Impedance
Z_terminal = inv(C' * inv(Z_Lambda) * C);

%% Display some results
disp('Terminal Impedance:');
disp(Z_terminal);

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

%% Calculate filament currents and reshape for plotting
i_lambda = Z_Lambda^(-1) * C * Z_terminal * I;
i_lambda_matrix = reshape(i_lambda, ny, nx, numConductors);

%% Calculate current density
current_density = abs(i_lambda_matrix) / (dx * dy);

%% Determine common color limits for consistent plotting
max_current = max(abs(i_lambda_matrix(:)));
max_density = max(abs(current_density(:)));

%% Plot Current Distribution and Current Density for each conductor
figure;
for k = 1:numConductors
    % Current Distribution
    subplot(2, numConductors, k);
    imagesc(abs(i_lambda_matrix(:, :, k)));
    colorbar;
    caxis([0 max_current]); % Consistent color limits
    axis equal tight;
    title(sprintf('Current Distribution in Conductor %d (A)', k));
    xlabel('Filament Index X');
    ylabel('Filament Index Y');
    set(gca, 'Position', get(gca, 'Position') + [0+0.025*(k-1) 0 0 0]); % Adjust subplot spacing

    % Current Density
    subplot(2, numConductors, numConductors + k);
    imagesc(current_density(:, :, k));
    colorbar;
    caxis([0 max_density]); % Consistent color limits
    axis equal tight;
    title(sprintf('Current Density in Conductor %d (A/m^2)', k));
    xlabel('Filament Index X');
    ylabel('Filament Index Y');
    set(gca, 'Position', get(gca, 'Position') + [0+0.025*(k-1) 0 0 0]); % Adjust subplot spacing
end

%% Calculate the losses in all the conductors

P_total = real(I'*Z_terminal*conj(I))*0.5;
disp('Total Losses in the Conductor (W):');
disp(P_total);

%% Scatter plot to visualize the conductors in space

figure()
scatter(positions(:,1),positions(:,2),'ro','SizeData',0.2);
title('Conductors in space');
xlabel('x');
ylabel('y');
xlim([min(positions(:,1))-a/2, max(positions(:,1))+a/2]);
axis equal 

%% Calculate Error in Total Losses and Percentage Error

% Compute absolute error and percentage error
error_loss = abs(P_total - reference_losses);
percentage_error = (error_loss / reference_losses) * 100;

disp('Error in Total Losses (W):');
disp(error_loss);
disp('Percentage Error (%):');
disp(percentage_error);