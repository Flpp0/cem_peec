clc, clear all, close all

format short e 

% Constants
mu0 = 4 * pi * 1e-7;
l = 1;  % Consistent with the function: must be in meter (m)!

% Ranges for l/W (u) and T/W (omega)
u_values = logspace(-1, 3, 100);
omega_values = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0];

% Here with the figure we want to test whether the function fo calculate
% the Lpii is correct by reproducing the:
% Fig 5. Inductance Calculations in a Complex Integrated Circuit Environment. A. E. Ruehli
figure;
hold on;

for omega = omega_values
    correction_factor = 1e-6; % nH/mm : 1e-9/1e-3 -> 1e-6
    L_pii_values = arrayfun(@(u) L_pii(mu0, omega, u, l), u_values)/correction_factor;
    plot(u_values, L_pii_values, 'DisplayName', sprintf('T/W = %.2f', omega));
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('l/W');
ylabel('L_{p_{ii}}/l (nH/mm)');
title('Partial self-inductance L_{p_{ii}}/l (nH/mm) for rectangular conductors');
legend('Location','southeast');
grid on;
hold off;

%% Apply the function to example case
t = 3.556e-5;       % Thickness (m)
w = 3.81e-4;        % Width (m)
l = 0.254;          % Length (m) 

NT = 16; 
NW = 172;

dt = t / NT; % (m)
dw = w / NW; % (m)

u = l/dw; 
omega = dt/dw;

L_pii_bar = L_pii(mu0,omega,u,l);

disp(['u = ', num2str(u)])
disp(['omega = ', num2str(omega)])
disp(['L_pii_bar = ', num2str(L_pii_bar)])

% Check approximation of round wire
L_pii_bar_round_wire = L_pii_round_wire(mu0, omega, u, l);
disp(['L_pii_bar_round_wire = ', num2str(L_pii_bar_round_wire)])


%% Single Conductor
% References
% [1]: Paul C.R.-Inductance_Loop and partial-Willey (2009)
clc; clear all; close all;

% Constants
mu0 = 4 * pi * 1e-7; % Permeability of free space (H/m)

sigma_cu = 5.8e7;  % Conductivity of copper (S/m)

inch_m = 0.0254;
t_mils = 1.4; 
w_mils = 15; 
l_inch = 10;

t = 1e-3*inch_m*t_mils; 
w = 1e-3*inch_m*w_mils; 
l = inch_m*l_inch;
l = 1; % [m] to scale the values 

% Given parameters
%t = 3.81e-5;        % Thickness (m)
%w = 3.556e-4;         % Width (m)

data = readmatrix('SingleConductor_ReferenceLosses.xlsx');  % Losses calculated with COMSOL
frequencies = data(:,1);
target_losses = data(:,2);
%frequencies = [10e1, 10e2, 10e3, 10e4, 100e3, 10e6, 100e6, 1e9];

% Skin depth at frequencies
skin_depth = sqrt(1./(sigma_cu*pi*frequencies*mu0));

NT = 16; % Number of elements in thickness
NW = 172; % Number of elements in width

dt = t / NT; % Thickness step (m)
dw = w / NW; % Width step (m)

subbar_area = dt*dw;

% Print discretization steps
disp(['dw = ', num2str(dw)]);
disp(['dt = ', num2str(dt)]);

% Mesh with conductor centered at (0, 0)
[x, y] = meshgrid(linspace(-w/2, w/2, NW+1), linspace(-t/2, t/2, NT+1));
x_center = x(1:end-1, 1:end-1) + dw / 2;
y_center = y(1:end-1, 1:end-1) + dt / 2;

% Plot the mesh and the centers
figure;
hold on;
plot(x, y, 'k'); 
plot(x', y', 'k'); 
plot(x_center(:), y_center(:), 'r.');
xlabel('Width (w)');
ylabel('Thickness (t)');
title('Uniform Mesh and Element Centers for Rectangular Conductor Centered at (0,0)');
axis equal;
hold off;

% Resistance of each element
r_i = l / (sigma_cu * dw * dt);
disp(['Resistance of each subelement r_i = ', num2str(r_i), ' ohms']);

u = l / dw; 
omega = dt / dw;

% Self-partial inductances
L_pii_temp = L_pii(mu0, omega, u, l); 
L_pii_bar = L_pii_temp*l; % Multiply by l for correct scaling to total inductance

% High-frequency resistance asymptote
r_hf = l ./ (2 * sigma_cu * skin_depth .* (dw + dt)); % (6.71a [1]). Not using the per-unit-length

% Precompute distances
[X1, X2] = meshgrid(x_center(:), x_center(:));
[Y1, Y2] = meshgrid(y_center(:), y_center(:));
distances = sqrt((X1 - X2).^2 + (Y1 - Y2).^2);

% Calculate mutual partial inductance M_pij
Mpij_values = (mu0 * l) / (2 * pi) * (log((l + sqrt(l^2 + distances.^2)) ./ distances) - sqrt(1 + (distances / l).^2) + (distances / l));

% Set the diagonal elements to zero (self-inductance is handled separately)
Mpij_values(1:size(Mpij_values, 1) + 1:end) = 0;

% Distances between mesh centers 
figure;
imagesc(distances); 
colorbar;
xlabel('Sub-element Index');
ylabel('Sub-element Index');
title('Heatmap of Distances Between Mesh Centers (dij)');
axis equal;

% Mutual inductance values 
figure;
imagesc(Mpij_values); 
colorbar;
xlabel('Sub-element Index');
ylabel('Sub-element Index');
title('Heatmap of Mutual Partial Inductances (Mpij)');
axis equal;

% Scatter plot of Mpij values versus distances 
figure;
scatter(distances(:), Mpij_values(:), 'filled');
xlabel('Distance between sub-element centers (m)');
ylabel('Mutual Partial Inductance (H)');
title('Mutual Partial Inductance vs. Distance');
grid on;

% Initialize results
net_resistance_per_unit_length = zeros(size(frequencies));
total_inductance_per_unit_length = zeros(size(frequencies));
total_currents = zeros(size(frequencies)); 
power_losses = zeros(size(frequencies)); % Store power losses

% Define the voltage 
V0 = 1; 

% Different frequencies
for idx = 1:length(frequencies)
    f = frequencies(idx);
    omega_f = 2 * pi * f;

    N = NT * NW;
    Z = zeros(N, N);

    % Self-impedances (per unit length)
    for i = 1:N
        Z(i, i) = r_i + 1i * omega_f * L_pii_bar; 
    end

    % Mutual impedances (per unit length)
    Z = Z + 1i * omega_f * Mpij_values;

    % Visualize the Impedance Matrix Z
    % figure;
    % imagesc(abs(Z)); 
    % colorbar;
    % xlabel('Sub-element Index');
    % ylabel('Sub-element Index');
    % title(['Magnitude of Impedance Matrix |Z| at ', num2str(f/1e6), ' MHz']);
    % axis equal;
    % 
    % % Plot the Resistance Matrix Re(Z)
    % figure;
    % imagesc(real(Z)); 
    % colorbar;
    % xlabel('Sub-element Index');
    % ylabel('Sub-element Index');
    % title(['Resistance Matrix Re(Z) at ', num2str(f/1e6), ' MHz']);
    % axis equal;
    % 
    % % Plot the Inductive Reactance Matrix Im(Z)
    % figure;
    % imagesc(imag(Z)); 
    % colorbar;
    % xlabel('Sub-element Index');
    % ylabel('Sub-element Index');
    % title(['Inductive Reactance Matrix Im(Z) at ', num2str(f/1e6), ' MHz']);
    % axis equal;

    % Voltage-driven approach for current distribution and power loss
    V = V0*ones(N, 1); % Voltage vector (1 V per subbar)
    I_subbars = Z \ V; % Solve for currents using matrix division: more stable

    % Calculate the total current
    total_current = sum(I_subbars);
    total_currents(idx) = total_current;

    % Calculate Z_total
    Z_total = V0 / total_current; 

    % Extract the real and imaginary parts for the net resistance and inductance
    R_total = real(Z_total);
    L_p_total = imag(Z_total) / omega_f;

    % Store net resistance per unit length
    net_resistance_per_unit_length(idx) = R_total;

    % Store the total inductance per unit length
    total_inductance_per_unit_length(idx) = L_p_total;

    % Calculate power losses using sum of real parts of impedances times current squared
    P_total = r_i * (I_subbars' * I_subbars)/2; 

    power_losses(idx) = P_total;
    
    plotCurrentDistribution = false;

    % Plot current distribution at specified frequencies
    if ismember(f, frequencies) && plotCurrentDistribution
        I_matrix = reshape(I_subbars, [NT, NW]);
        
        % Plot the current distribution
        figure;
        surf(x_center * 1e3, y_center * 1e3, abs(I_matrix) * 1e6); % Convert to mm and mA
        xlabel('Width (mm)');
        ylabel('Thickness (mm)');
        zlabel('Current Magnitude (mA)');
        title(['Current Distribution at ', num2str(f), ' Hz']);
        colorbar;

        % % Plot the absolute value of current in each subbar
        % figure;
        % imagesc(abs(I_matrix./subbar_area) * 1e6); % Convert to mA
        % colorbar;
        % xlabel('Width (subbars)');
        % ylabel('Thickness (subbars)');
        % title(['Abs(Current) in Each Subbar at ', num2str(f), ' Hz']);
    end
end

% Analyze dij and Mpij values
min_dij = min(distances(:));
max_dij = max(distances(:));
disp(['Minimum dij: ', num2str(min_dij)]);
disp(['Maximum dij: ', num2str(max_dij)]);
disp(['Minimum Mpij: ', num2str(min(Mpij_values(:)))]);
disp(['Maximum Mpij: ', num2str(max(Mpij_values(:)))]);

% Plot net resistance vs frequency (on a logarithmic scale)
figure;
loglog(frequencies / 1e6, net_resistance_per_unit_length, '-o', 'DisplayName', 'Calculated Net Resistance');
hold on;
%loglog(frequencies / 1e6, sqrt(frequencies), '--', 'DisplayName', 'High-Frequency Asymptote r_{hf}');
xlabel('Frequency (MHz)');
ylabel('Net Resistance (Ohm/m)');
title('Net Resistance vs Frequency');
legend('Location', 'best');
grid on;
% Additional settings to make the plot look like Figure 6.16(a)
set(gca, 'XScale', 'log', 'YScale', 'log'); 

% Plot total current vs frequency
figure;
loglog(frequencies / 1e6, abs(total_currents), '-o', 'DisplayName', 'Total Current');
xlabel('Frequency (MHz)');
ylabel('Total Current (A)');
title('Total Current vs Frequency');
legend('Location', 'best'); 
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log'); 

% Plot power losses vs frequency
figure;
plot(frequencies / 1e6, power_losses, '-o', 'DisplayName', 'Power Loss');
xlabel('Frequency (MHz)');
ylabel('Power Loss (W)');
title('Power Loss vs Frequency');
legend('Location', 'best'); 
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log'); 


% Plot power losses vs Target power losses LOG LOG SCALE
figure; 
subplot(1,2,1)
hold on
plot(frequencies, power_losses, '-o', 'DisplayName', 'Predicted Power Loss');
plot(frequencies, target_losses, '-o', 'DisplayName', 'Target Power Loss');
xlabel('Frequency (MHz)');
ylabel('Power Loss (W)');
title('Power Loss vs Frequency');
legend('Location', 'best'); 
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log'); 

subplot(1,2,2)
plot(frequencies, (power_losses - target_losses)*100, '-o', 'DisplayName', 'Percentage of error');
xlabel('Frequency (MHz)');
ylabel('Power Loss (W)');
title('Percentage of error');
legend('Location', 'best'); 
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log'); 

% Plot power losses vs Target power losses 
figure; 
subplot(1,2,1)
hold on
plot(frequencies, power_losses, '-o', 'DisplayName', 'Predicted Power Loss');
plot(frequencies, target_losses, '-o', 'DisplayName', 'Target Power Loss');
xlabel('Frequency (MHz)');
ylabel('Power Loss (W)');
title('Power Loss vs Frequency');
legend('Location', 'best');
grid on;
set(gca, 'XScale', 'log'); 

subplot(1,2,2)
plot(frequencies, (power_losses - target_losses)*100, '-o', 'DisplayName', 'Percentage of error');
xlabel('Frequency (MHz)');
ylabel('Power Loss (W)');
title('Percentage of error');
legend('Location', 'best'); 
grid on;
set(gca, 'XScale', 'log');

% Plot the total current in the conductor
figure
plot(frequencies, total_currents, '-o', 'DisplayName', 'Total Current')
xlabel('Frequency (Hz)');
ylabel('Total Current (A)'); 
title('Total Current vs Frequency')
legend('Location', 'Best');
grid on; 
set(gca, 'XScale', 'log')