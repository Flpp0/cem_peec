clear all; close all; clc;

format short e
%% Constants and Parameters
mu0 = 4 * pi * 1e-7;            % Permeability of free space (H/m)
R = 0.1;                        % Reference radius (m)
rho = 1/5.998e7;                % Resistivity of copper (Ohm*m)
f = 50e2;                       % Frequency (Hz)
omega = 2 * pi * f;             % Angular frequency (rad/s)
l = 1;                          % Length of the conductor (m): scale the results

reference_losses = 4.808550807587721e-3; % From COMSOL: see ThreeConductors.mph

a = 0.003;                      % Width of the conductor (m)
b = 0.002;                      % Height of the conductor (m)
nx = 50;                        % Number of points along x
ny = 50;                        % Number of points along y
numConductors = 3;              % Number of conductors
N = (nx - 1) * (ny - 1);        % Total number of filaments per conductor

%% Non - uniform mesh
totalPositions = [];
totalCenters = [];             
totalAreas = [];

for k = 1:numConductors
    [X, Y, Areas] = generateNonUniformMesh(a, b, nx, ny);
    xShift = (k - 2) * (5 * a);  % Shift to the left and to the right the conductors (along x-axis)
    X = X(1:end, 1:end) + xShift; 
    Y = Y(1:end, 1:end);
    
    % Centers of each filament
    CentersX = (X(1:end-1, 1:end-1) + X(2:end, 2:end)) / 2;
    CentersY = (Y(1:end-1, 1:end-1) + Y(2:end, 2:end)) / 2;
    
    totalPositions = [totalPositions; [X(:), Y(:)]];
    totalCenters = [totalCenters; [CentersX(:), CentersY(:)]];
    totalAreas = [totalAreas; Areas(:)];
end

%% Calculate Resistance
Resistance = diag(rho * l ./ totalAreas);

%% Calculate Inductance 
xi = totalCenters / R; % Using centers for distance calculation
Inductance = zeros(N * numConductors, N * numConductors);

distances = sqrt((xi(:,1) - xi(:,1)').^2 + (xi(:,2) - xi(:,2)').^2);

for i = 1:size(xi,1)
    for j = 1:size(xi,1)
        if i ~= j
            G_mutual = log(distances(i,j)) - 0.5 * log(norm(xi(i,:))^2 * norm(xi(j,:))^2 - 2 * dot(xi(i,:), xi(j,:)) + 1);
            Inductance(i, j) = -mu0 * l / (2 * pi) * G_mutual;
        else
            rho_i = sqrt(totalAreas(i) / pi) / R;
            G_self = -mu0 * l / (2 * pi) * (log(rho_i) - 0.5 * log((1 - norm(xi(i,:))^2))^2 + (norm(xi(i,:))^2 * rho_i)^2);
            Inductance(i, i) = G_self;
        end
    end
end

%% Combined impedance matrix
Z_Lambda = Resistance + 1i * omega * Inductance;

%% Connectivity matrix and current setup
I = [1; 0.8; 1.2];  %[A]
C = kron(eye(numConductors), ones(N,1)); % Assign each filament to its conductor 

%% Compute terminal impedance
U = Z_Lambda * C * I; % Not used; only to see the terminal quantities together 

Z_terminal = inv(C' * inv(Z_Lambda) * C);

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

%% Compute the current in each filament !!
i_lambda = Z_Lambda^(-1) * C * Z_terminal * I;

i_lambda_matrix = reshape(i_lambda, ny-1, nx-1, numConductors); % 3D matrix, each layer for conductor  

%% Compute the current densities
current_density = zeros(ny-1, nx-1, numConductors);

for k = 1:numConductors
    area_matrix = reshape(totalAreas((k-1)*N+1:k*N), ny-1, nx-1);  

    current_density(:, :, k) = abs(i_lambda_matrix(:, :, k)) ./ area_matrix;
end

%% Plot Current Distribution and Current Density for each conductor

max_current = max(abs(i_lambda_matrix(:)));
max_density = max(current_density(:));

figure;
for k = 1:numConductors
    % Current Distribution (A)
    subplot(2, numConductors, k); % In the first row
    imagesc(abs(i_lambda_matrix(:, :, k)));
    colorbar;
    caxis([0 max_current]); 
    axis equal tight;
    title(sprintf('Current Distribution in Conductor %d (A)', k));
    xlabel('Filament Index X');
    ylabel('Filament Index Y');

    % Current Density (A/m^2)
    subplot(2, numConductors, numConductors + k); % In the second row 
    imagesc(current_density(:, :, k));
    colorbar;
    caxis([0 max_density]);  
    axis equal tight;
    title(sprintf('Current Density in Conductor %d (A/m^2)', k));
    xlabel('Filament Index X');
    ylabel('Filament Index Y');
end

%% Compute losses in all the conductors
% Two different approaches to compare the results using the terminal
% quantities and the currents in each filament
P_total = real(I' * Z_terminal * conj(I)) * 0.5; 

disp('Total Losses in the Conductor (W):');
disp(P_total);

P_total = 0; 
for i = 1:1:size(i_lambda,1)
    P_total = P_total + Resistance(i,i)*(i_lambda(i)*conj(i_lambda(i)))*0.5;
end

disp('Total Losses in the Conductor (W):');
disp(P_total);

%% Calculate Error in Total Losses and Percentage Error

% Compute absolute error and percentage error
error_loss = abs(P_total - reference_losses);
percentage_error = (error_loss / reference_losses) * 100;

disp('Error in Total Losses (W):');
disp(error_loss);
disp('Percentage Error (%):');
disp(percentage_error);
