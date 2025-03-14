clear all, close all, clc

%% Constants and Parameters
mu0 = 4 * pi * 1e-7;             % Permeability of free space (H/m)
R = 0.1;                         % Radius of the cylinder (m) -> match COMSOL 
rho = 1.68e-8;                   % Resistivity of copper (Ohm*m)
f = 50e2;                          % Frequency (Hz)
omega = 2 * pi * f;              % Pulsation (rad/s)
l = 1;                           % Length of the conductor (m)

skin_depth = 1 / sqrt(mu0 * (1 / rho) * omega); % Skin depth (m)

a = 0.003;                       % Width of the conductor (m)
b = 0.002;                       % Height of the conductor (m)
nx = 50;                         % Number of filaments along x
ny = 50;                         % Number of filaments along y

reference_losses = 1.547812645914713e-3; % COMSOL see: SingleConductor.mph

% Filament dimensions
dx = a / nx;
dy = b / ny;

r_hat = sqrt(a * b / pi);        % Equivalent radius

N = nx * ny;                     % Total number of filaments
Resistance = rho * l / (dx * dy) * eye(N); % Per unit of length 
Inductance = zeros(N, N);

[X, Y] = meshgrid(linspace(-a/2 + dx/2, a/2 - dx/2, nx), linspace(-b/2 + dy/2, b/2 - dy/2, ny));
positions = [X(:), Y(:)];

% Calculate mutual and self inductances using normalized coordinates
for i = 1:N
    for j = 1:N
        xi = [positions(i,1) / R, positions(i,2) / R];
        xj = [positions(j,1) / R, positions(j,2) / R];

        distance = norm(xi - xj);  % Euclidean distance between normalized positions
        xi_dot_xj = dot(xi, xj);   % Dot product of normalized positions

        if i ~= j
            G_mutual = log(distance) - 0.5 * log(norm(xi)^2 * norm(xj)^2 - 2 * xi_dot_xj + 1);
            Inductance(i,j) = -mu0 * l / (2 * pi) * G_mutual;  % Mutual inductance
        end
    end

    % Self-inductance for filament i
    rho_i = r_hat / R; % Normalized cross section radius
    G_self = log(rho_i) - 0.5 * log((1 - norm(xi)^2))^2 + (norm(xi)^2 * rho_i)^2;
    Inductance(i, i) = -mu0 * l / (2 * pi) * G_self;  % Self-inductance
end

% Combined impedance matrix
Z_Lambda = Resistance + 1i * omega * Inductance;

% Connectivity matrix 
C = ones(N,1); 

% Terminal current and voltage
I = 1; % Terminal current (Ampere)
U = Z_Lambda * C * I; % Terminal voltage using Ohm's law

% Terminal Impedance
Z_terminal = inv(C' * inv(Z_Lambda) * C);

% Filament voltages
u_lambda = C .* U; 

% Filament currents
i_lambda = Z_Lambda^(-1) * C * Z_terminal * I;

%% Display results: resistance, inductance, Z_Lambda and Z_terminal
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
title('Impedance matrix');
xlabel('Filament Index');
ylabel('Filament Index');

disp('Terminal Impedance:');
disp(Z_terminal);

%% Display current and current density in each filament

i_lambda_matrix = reshape(i_lambda, [ny, nx]);

current_density = abs(i_lambda_matrix) / (dx * dy);

figure;
subplot(1, 2, 1); 
imagesc(abs(i_lambda_matrix));
colorbar;
title('Current Distribution (A)');
xlabel('Filament Index X');
ylabel('Filament Index Y');
axis equal tight;

subplot(1, 2, 2); 
imagesc(current_density);
colorbar;
title('Current Density (A/m^2)');
xlabel('Filament Index X');
ylabel('Filament Index Y');
axis equal tight; 

disp('Total current summing all the current contributions from the filaments:')
disp(sum(i_lambda_matrix,'all'))

%% Calculate losses in the conductor
P_total = 0.5*abs(I)^2*real(Z_terminal);

disp('Total Losses in the Conductor (W):');
disp(P_total);

%% Calculate error between estimated and reference losses
error_abs = abs(P_total - reference_losses);
percentage_error = (error_abs / reference_losses) * 100;

disp('Absolute error between estimated and reference losses (W):');
disp(error_abs);

disp('Percentage error (%):');
disp(percentage_error);