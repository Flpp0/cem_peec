clear all; close all; clc;

% Constants
mu0 = 4 * pi * 1e-7;       % Permeability of free space (H/m)
rho = 1 / 5.998e7;         % Resistivity of copper (Ohm*m)
f = 50e2;                  % Frequency (Hz)
omega_f = 2 * pi * f;      % Angular frequency (rad/s)
l = 1;                     % Length of the conductor (m)

reference_losses = 4.808550807587721e-3; % From COMSOL: see ThreeConductors.mph

format short e

% Conductor dimensions
a = 0.003;          % Width (m)
b = 0.002;          % Height (m)
nw = 30;            % Number of points along x
nt = 30;            % Number of points along y
numConductors = 3;  % Number of conductors
N = (nw - 1) * (nt - 1); % Filaments per conductor

% Generate uniform positions, centers, and areas of the filaments
totalPositions = [];
totalCenters = [];
totalAreas = [];

dw = a / (nw - 1);
dt = b / (nt - 1);

for k = 1:numConductors
    [X, Y] = meshgrid(0:dw:a, 0:dt:b);
    xShift = (k - 2) * (5 * a);
    X = X + xShift;
    
    CentersX = (X(1:end-1, 1:end-1) + X(2:end, 2:end)) / 2;
    CentersY = (Y(1:end-1, 1:end-1) + Y(2:end, 2:end)) / 2;
    
    Areas = dw * dt * ones(size(CentersX));
    
    totalPositions = [totalPositions; [X(:), Y(:)]];
    totalCenters = [totalCenters; [CentersX(:), CentersY(:)]];
    totalAreas = [totalAreas; Areas(:)];
end

Resistance = diag(rho * l ./ totalAreas);

Inductance = zeros(N * numConductors, N * numConductors);
for i = 1:N * numConductors
    for j = 1:N * numConductors
        if i == j
            % For self-inductance, compute parameters and call the external L_pii
            u_val = l / dw;
            omega_val = dt / dw;
            Lpi = L_pii(mu0, omega_val, u_val, l);  % Direct call to the existing function
            Inductance(i, j) = Lpi / l;
        else
            x1 = totalCenters(i, 1); y1 = totalCenters(i, 2);
            x2 = totalCenters(j, 1); y2 = totalCenters(j, 2);
            d_ij = sqrt((x1 - x2)^2 + (y1 - y2)^2);
            Inductance(i, j) = (mu0 * l / (2 * pi)) * ( log(l / d_ij + sqrt(1 + (l / d_ij)^2)) - sqrt(1 + (d_ij / l)^2) + d_ij / l );
        end
    end
end

Z_Lambda = Resistance + 1i * omega_f * Inductance;

I = [1; 0.8; 1.2];  % [A]
C = kron(eye(numConductors), ones(N,1));

Z_terminal = inv(C' * inv(Z_Lambda) * C);
disp('Terminal Impedance:');
disp(Z_terminal);

i_lambda = Z_Lambda^(-1) * C * Z_terminal * I;
i_lambda_matrix = reshape(i_lambda, nt-1, nw-1, numConductors);

current_density = zeros(nt-1, nw-1, numConductors);
for k = 1:numConductors
    area_matrix = reshape(totalAreas((k-1)*N+1:k*N), nt-1, nw-1);  
    current_density(:, :, k) = abs(i_lambda_matrix(:, :, k)) ./ area_matrix;
end

max_current = max(abs(i_lambda_matrix(:)));
max_density = max(current_density(:));

figure;
for k = 1:numConductors
    subplot(2, numConductors, k);
    imagesc(abs(i_lambda_matrix(:, :, k)));
    colorbar;
    caxis([0 max_current]); 
    axis equal tight;
    title(sprintf('Current Distribution in Conductor %d (A)', k));
    xlabel('Filament Index X');
    ylabel('Filament Index Y');

    subplot(2, numConductors, numConductors + k);
    imagesc(current_density(:, :, k));
    colorbar;
    caxis([0 max_density]);  
    axis equal tight;
    title(sprintf('Current Density in Conductor %d (A/m^2)', k));
    xlabel('Filament Index X');
    ylabel('Filament Index Y');
end

% Compute losses using terminal quantities
P_total_terminal = real(I' * Z_terminal * conj(I)) * 0.5;
disp('Total Losses in the Conductor (W) using terminal quantities:');
disp(P_total_terminal);

% Compute losses using filament quantities
P_total_filament = 0; 
for i = 1:length(i_lambda)
    P_total_filament = P_total_filament + Resistance(i,i) * (i_lambda(i) * conj(i_lambda(i))) * 0.5;
end
disp('Total Losses in the Conductor (W) using filament quantities:');
disp(P_total_filament);

% Error Calculation 
error_terminal = abs(P_total_terminal - reference_losses);
perc_error_terminal = (error_terminal / reference_losses) * 100;

error_filament = abs(P_total_filament - reference_losses);
perc_error_filament = (error_filament / reference_losses) * 100;

disp('Error (Terminal method):');
disp(error_terminal);
disp('Percentage Error (Terminal method) (%):');
disp(perc_error_terminal);
disp('Error (Filament method):');
disp(error_filament);
disp('Percentage Error (Filament method) (%):');
disp(perc_error_filament);

% Mesh Display 
figure;
hold on;
for k = 1:numConductors
    [X, Y] = meshgrid(0:dw:a, 0:dt:b);
    xShift = (k - 2) * (5 * a);
    X = X + xShift;
    for i = 1:size(X,1)
        plot(X(i,:), Y(i,:), 'k');
    end
    for j = 1:size(X,2)
        plot(X(:,j), Y(:,j), 'k');
    end
    CentersX = (X(1:end-1,1:end-1) + X(2:end,2:end)) / 2;
    CentersY = (Y(1:end-1,1:end-1) + Y(2:end,2:end)) / 2;
    scatter(CentersX(:), CentersY(:), 10, 'r', 'filled');
end
xlabel('x (m)');
ylabel('y (m)');
title('Mesh and Conductor Positions');
axis equal;
hold off;

% 3D Surface Plot of Current Distribution 
figure;
hold on;
for k = 1:numConductors
    % Generate grid for the filament centers for plotting
    [Xc, Yc] = meshgrid(linspace(0, a, nw-1), linspace(0, b, nt-1));
    xShift = (k - 2) * (5 * a);
    Xc = Xc + xShift;
    
    % Plot the 3D surface for the current distribution of conductor k
    surf(Xc, Yc, abs(i_lambda_matrix(:,:,k)), 'EdgeColor', 'none');
end
colorbar;
title('3D Current Distribution on All Conductors');
xlabel('x (m)');
ylabel('y (m)');
zlabel('Current (A)');
view(3);  % 3D view
hold off;

%% Convergence Test: Varying the Number of Elements

nElementsVec = 10:10:40;
numTests = length(nElementsVec);
err_terminal = zeros(numTests, 1);
err_filament = zeros(numTests, 1);

% See the converge current step
h_outer = waitbar(0, 'Running convergence test...');

for idx = 1:numTests

    waitbar(idx/numTests, h_outer, sprintf('Test %d of %d', idx, numTests));
    
    % Current discretization
    nw = nElementsVec(idx);
    nt = nElementsVec(idx);
    N = (nw - 1) * (nt - 1);  % Filaments per conductor
    
    % Aspects ratio
    dw = a / (nw - 1);
    dt = b / (nt - 1);
    
    totalPositions = [];
    totalCenters = [];
    totalAreas = [];
    
    for k = 1:numConductors
        [X, Y] = meshgrid(0:dw:a, 0:dt:b);
        xShift = (k - 2) * (5 * a);
        X = X + xShift;
        
        CentersX = (X(1:end-1, 1:end-1) + X(2:end, 2:end)) / 2;
        CentersY = (Y(1:end-1, 1:end-1) + Y(2:end, 2:end)) / 2;
        Areas = dw * dt * ones(size(CentersX));
        
        totalPositions = [totalPositions; [X(:), Y(:)]];
        totalCenters = [totalCenters; [CentersX(:), CentersY(:)]];
        totalAreas = [totalAreas; Areas(:)];
    end
    
    % Compute resistance matrix
    Resistance = diag(rho * l ./ totalAreas);
    
    % Total number of filaments
    totalFilaments = N * numConductors;
    
    D = pdist2(totalCenters, totalCenters); % See documentation.
    
    % For the off-diagonal elements, apply the mutual inductance formula.
    D_no_diag = D + eye(totalFilaments);  
    Inductance = (mu0 * l / (2 * pi)) * ( log(l ./ D_no_diag + sqrt(1 + (l ./ D_no_diag).^2)) ...
                  - sqrt(1 + (D_no_diag / l).^2) + D_no_diag / l );
    
    u_val = l / dw;
    omega_val = dt / dw;
    Lpi_val = L_pii(mu0, omega_val, u_val, l) / l; % Self-inductances are the same: uniform mesh
    
    Inductance(eye(totalFilaments)==1) = Lpi_val;
    
    % Construct the total impedance matrix
    Z_Lambda = Resistance + 1i * omega_f * Inductance;
    
    % Terminal calculations
    I = [1; 0.8; 1.2];
    C = kron(eye(numConductors), ones(N, 1));
    Z_terminal = inv(C' * inv(Z_Lambda) * C);
    
    % Compute current in each filament
    i_lambda = Z_Lambda^(-1) * C * Z_terminal * I;
    
    % Compute losses using terminal quantities
    P_total_terminal = real(I' * Z_terminal * conj(I)) * 0.5;
    
    % Compute losses using filament quantities (vectorized summation)
    P_total_filament = sum(diag(Resistance) .* (abs(i_lambda).^2)) * 0.5;
    
    % Compute the absolute error relative to the reference losses
    err_terminal(idx) = abs(P_total_terminal - reference_losses);
    err_filament(idx) = abs(P_total_filament - reference_losses);
end

close(h_outer);

% Plot the convergence results
figure;
plot(nElementsVec, err_terminal, '-o', 'LineWidth', 2);
hold on;
plot(nElementsVec, err_filament, '-s', 'LineWidth', 2);
xlabel('Number of Elements (nw = nt)');
ylabel('Absolute Error in Loss Calculation (W)');
title('Convergence Test: Error vs. Number of Elements');
legend('Terminal Method', 'Filament Method');
grid on;

