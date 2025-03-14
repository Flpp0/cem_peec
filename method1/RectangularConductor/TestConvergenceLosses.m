clear all; close all; clc;

format short e
%% Constants and Parameters
mu0 = 4 * pi * 1e-7;            % Permeability of free space (H/m)
R = 0.1;                        % Reference radius (m)
rho = 1/5.998e7;                % Resistivity of copper (Ohm*m)
f = 50e2;                       % Frequency (Hz)
omega = 2 * pi * f;             % Angular frequency (rad/s)
l = 1;                          % Length of the conductor (m): scale the results

a = 0.003;                      % Width of the conductor (m)
b = 0.002;                      % Height of the conductor (m)
numConductors = 3;              % Number of conductors

reference_losses = 4.808550807587721e-3; % From COMSOL: see ThreeConductors.mph
% As reference from COMSOL. 
% Here we should insert the comsol file from which we calculated the
% losses. 

%% Iteration parameters
startingElements = 10;
numIterations = 5;            
elementIncrements = 5;          
errors_vec = zeros(numIterations, 1);
elements = zeros(numIterations, 1);

%% Vectorized version without parallelization
tic;
for iter = 1:numIterations
    nx = startingElements + iter * elementIncrements; 
    ny = startingElements + iter * elementIncrements;  
    N = (nx - 1) * (ny - 1);             
    
    %% Non-uniform mesh
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
    
    % Precompute log and distance based values
    xi_norms = vecnorm(xi, 2, 2);
    log_distances = log(distances);
    log_norms = log(xi_norms);
    norms_product = xi_norms * xi_norms';

    % Vectorized computation of G_mutual
    G_mutual_matrix = log_distances - 0.5 * log(norms_product.^2 - 2 * (xi * xi') + 1);
    G_mutual_matrix(distances == 0) = 0; % avoid log(0)

    % Inductance matrix off-diagonal elements
    Inductance = -mu0 * l / (2 * pi) * G_mutual_matrix;

    % Self inductance (diagonal elements)
    rho_i = sqrt(totalAreas / pi) / R;
    log_rho_i = log(rho_i);
    xi_norms = vecnorm(xi, 2, 2);
    norms_squared = xi_norms.^2;
    G_self = -mu0 * l / (2 * pi) * (log_rho_i - 0.5 * log((1 - norms_squared).^2 + norms_squared.^2 .* rho_i.^2));
    Inductance(1:size(G_self,1)+1:end) = G_self;

    %% Combined impedance matrix
    Z_Lambda = Resistance + 1i * omega * Inductance;

    %% Connectivity matrix and current setup
    I = [1; 0.8; 1.2];  %[A]
    C = kron(eye(numConductors), ones(N,1)); % Assign each filament to its conductor 

    %% Compute terminal impedance
    Z_terminal = inv(C' * inv(Z_Lambda) * C);

    %% Compute the current in each filament
    i_lambda = inv(Z_Lambda)*C*Z_terminal*I;  % More efficient way to solve

    %% Compute losses in all the conductors
    P_total = sum(diag(Resistance) .* abs(i_lambda).^2 * 0.5);

    %% Calculate error
    error = abs(P_total - reference_losses);
    errors_vec(iter) = error;
    elements(iter) = N * numConductors;

    disp(['Iteration ', num2str(iter), ' (Vectorized): Error = ', num2str(error)]);
end
execution_time_vec = toc;

%% Display execution time
disp(['Total execution time: ', num2str(execution_time_vec), ' seconds']);

%% Plot the error as the number of elements varies
figure;
plot(elements, errors_vec, '-o');
title('Error in Total Losses vs. Number of Elements');
xlabel('Number of Elements');
ylabel('Error in Total Losses (W)');
legend('Vectorized');
grid on;
