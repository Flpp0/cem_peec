%% References for the methods:
% [1] - D. P. Morisco, S. Kurz, H. Rapp, and A. M¨ockel, “A hybrid modeling
%       approach for current diffusion in rectangular conductors,” IEEE
%       Transactions on Magnetics, vol. 55, no. 9, pp. 1–9, 2019.
% [2] - C. R. Paul, Inductance: Loop and Partial.
%       Hoboken, NJ: John Wiley & Sons, 2010.

% Before saving the results, verify the output figures; then re-run the full
% code to save the final data and images.
%
% -------------------------------------------------------------------------

clear all; close all; clc;

format short e

%% Load frequency data
data = readmatrix('ExportsCOMSOL/C_1C_Original_Dim_Loss.xlsx'); 
end_freq_idx = length(data);
frequencies = data(1:end_freq_idx, 1);
reference_losses_list = data(1:end_freq_idx, 2); 

% Constants and Parameters
mu0 = 4 * pi * 1e-7;            
R = 0.1;                        
rho = 1/5.998e7;                
l = 1;                          

a = 3.81e-4;                    
b = 3.556e-5;                   
numConductors = 1;              

%% Mesh Configurations for Method 1 and Method 2
meshConfigurationsMethod1 = [...
    40, 40;
    50, 50;
    60, 60
];  % List of (nx, ny) pairs for Method 1

meshConfigurationsMethod2 = [...
    172*0.25, 16*0.25;
    172*0.5, 16*0.5;
    172, 16; 
];  % List of (nx, ny) pairs for Method 2

% % For testing, you can use a single configuration:
% meshConfigurationsMethod1 = [...
%     40, 40;
%     50, 50
% ];  % List of (nx, ny) pairs for Method 1
% 
% meshConfigurationsMethod2 = [...
%     172*0.5, 16*0.5;
%     172, 16;
% ];  % List of (nx, ny) pairs for Method 2

numConfigurationsMethod1 = size(meshConfigurationsMethod1, 1);
numConfigurationsMethod2 = size(meshConfigurationsMethod2, 1);

% Storage
P_total_method1_direct_all = cell(length(frequencies), 1);
P_total_method2_direct_all = cell(length(frequencies), 1);
errors_vec_method1_direct_all = cell(length(frequencies), 1);
errors_vec_method2_direct_all = cell(length(frequencies), 1);
elements_vec_method1_all = cell(length(frequencies), 1);
elements_vec_method2_all = cell(length(frequencies), 1);

final_method1_direct = zeros(length(frequencies), 1);
final_method2_direct = zeros(length(frequencies), 1);

% Loop over frequencies
for freq_idx = 1:length(frequencies)
    f = frequencies(freq_idx);  
    omega = 2 * pi * f;         
    reference_losses = reference_losses_list(freq_idx);  
    
    errors_vec_method1_direct = zeros(numConfigurationsMethod1, 1);
    errors_vec_method2_direct = zeros(numConfigurationsMethod2, 1);
    
    P_total_method1_direct = zeros(numConfigurationsMethod1, 1);  
    P_total_method2_direct = zeros(numConfigurationsMethod2, 1);  
    elements_vec_method1 = zeros(numConfigurationsMethod1, 1);
    elements_vec_method2 = zeros(numConfigurationsMethod2, 1);
    
    %% Method 1 - [1]
    for config_idx = 1:numConfigurationsMethod1
        nx_method1 = meshConfigurationsMethod1(config_idx, 1); 
        ny_method1 = meshConfigurationsMethod1(config_idx, 2);  
        N_method1 = (nx_method1 - 1) * (ny_method1 - 1);             
        
        [X_method1, Y_method1, Areas_method1] = generateNonUniformMesh(a, b, nx_method1, ny_method1);
        CentersX_method1 = (X_method1(1:end-1, 1:end-1) + X_method1(2:end, 2:end)) / 2;
        CentersY_method1 = (Y_method1(1:end-1, 1:end-1) + Y_method1(2:end, 2:end)) / 2;
        totalCenters_method1 = [CentersX_method1(:), CentersY_method1(:)];
        totalAreas_method1 = Areas_method1(:);
        
        % Resistance matrix
        Resistance_method1 = diag(rho * l ./ totalAreas_method1);
        
        r_hat = sqrt(a * b / pi); % Characteristic lenght

        Inductance_method1 = zeros(N_method1, N_method1);
        for i = 1:N_method1
            xi = totalCenters_method1(i,:) / R;
            for j = 1:N_method1
                if i ~= j
                    xj = totalCenters_method1(j,:) / R;
                    distance = norm(xi - xj);
                    xi_dot_xj = dot(xi, xj);
                    G_mutual = log(distance) - 0.5 * log(norm(xi)^2 * norm(xj)^2 - 2*xi_dot_xj + 1);
                    Inductance_method1(i,j) = -mu0 * l / (2 * pi) * G_mutual;
                end
            end
            rho_i = r_hat / R;
            G_self = log(rho_i) - 0.5 * (log(1 - norm(xi)^2))^2 + (norm(xi)^2 * rho_i)^2;
            Inductance_method1(i,i) = -mu0 * l / (2 * pi) * G_self;
        end
        
        Z_Lambda_method1 = Resistance_method1 + 1i * omega * Inductance_method1;
        
        % Current stimulation: set terminal current I=1
        I_method1 = 1;  
        C_method1 = ones(N_method1,1);
        Z_terminal_method1 = inv(C_method1' * inv(Z_Lambda_method1) * C_method1);
        
        % Compute the current in each filament using direct method
        i_lambda_method1_direct = Z_Lambda_method1 \ (C_method1 * Z_terminal_method1 * I_method1);
        P_total_method1_direct(config_idx) = sum(diag(Resistance_method1) .* abs(i_lambda_method1_direct).^2 * 0.5); 
        errors_vec_method1_direct(config_idx) = abs(P_total_method1_direct(config_idx) - reference_losses);
        
        elements_vec_method1(config_idx) = N_method1 * numConductors;
        
        disp(['Frequency: ', num2str(f), ' Hz, Configuration ', num2str(config_idx), ' (Method 1): Direct Error = ', num2str(errors_vec_method1_direct(config_idx))]);
    end
    
    %% Method 2 - [2]
    for config_idx = 1:numConfigurationsMethod2
        nx_method2 = meshConfigurationsMethod2(config_idx, 1); 
        ny_method2 = meshConfigurationsMethod2(config_idx, 2);  
        N_method2 = (nx_method2 - 1) * (ny_method2 - 1);             
    
        % Uniform mesh
        [X_method2, Y_method2, Areas_method2] = generateUniformMesh(a, b, nx_method2, ny_method2);
        CentersX_method2 = (X_method2(1:end-1, 1:end-1) + X_method2(2:end, 2:end)) / 2;
        CentersY_method2 = (Y_method2(1:end-1, 1:end-1) + Y_method2(2:end, 2:end)) / 2;
        totalCenters_method2 = [CentersX_method2(:), CentersY_method2(:)];
        totalAreas_method2 = Areas_method2(:);
    
        Resistance_method2 = diag(rho * l ./ totalAreas_method2);
    
        xi_method2 = totalCenters_method2;
        Inductance_method2 = zeros(N_method2, N_method2);
    
        distances_method2 = sqrt((xi_method2(:,1) - xi_method2(:,1)').^2 + (xi_method2(:,2) - xi_method2(:,2)').^2);
        
        % Aspect ratios
        dw_method2 = a / nx_method2;
        dt_method2 = b / ny_method2; 
    
        u_method2 = l / dw_method2; 
        omega_method2 = dt_method2 / dw_method2;
    
        Mpij_values_method2 = (mu0 * l) / (2 * pi) * (log((l ./ distances_method2) + sqrt((l ./ distances_method2).^2 + 1)) ...
                       - sqrt(1 + (distances_method2 ./ l).^2) + distances_method2 ./ l);
        Mpij_values_method2(distances_method2 == 0) = 0; 
    
        Inductance_method2 = Mpij_values_method2; 
    
        L_pii_temp_method2 = L_pii(mu0, omega_method2, u_method2, l); 
        L_pii_bar_method2 = L_pii_temp_method2 * l; 
        Inductance_method2(1:N_method2 + 1:end) = L_pii_bar_method2; 
    
        % Combined impedance matrix
        Z_Lambda_method2 = Resistance_method2 + 1i * omega * Inductance_method2;
    
        % Current stimulation: set terminal current I=1
        I_method2 = 1;  
        C_method2 = ones(N_method2,1);  % Connectivity matrix
    
        % Compute terminal impedance
        Z_terminal_method2 = inv(C_method2' * inv(Z_Lambda_method2) * C_method2);
    
        % Compute the current in each filament using direct method
        i_lambda_method2_direct = Z_Lambda_method2 \ (C_method2 * Z_terminal_method2 * I_method2);
        P_total_method2_direct(config_idx) = sum(diag(Resistance_method2) .* abs(i_lambda_method2_direct).^2 * 0.5); 
        errors_vec_method2_direct(config_idx) = abs(P_total_method2_direct(config_idx) - reference_losses);
    
        elements_vec_method2(config_idx) = N_method2 * numConductors;
    
        disp(['Frequency: ', num2str(f), ' Hz, Configuration ', num2str(config_idx), ' (Method 2): Direct Error = ', num2str(errors_vec_method2_direct(config_idx))]);
    end
    
    % Results for Methods 1 and 2
    P_total_method1_direct_all{freq_idx} = P_total_method1_direct;
    P_total_method2_direct_all{freq_idx} = P_total_method2_direct;
    errors_vec_method1_direct_all{freq_idx} = errors_vec_method1_direct;
    errors_vec_method2_direct_all{freq_idx} = errors_vec_method2_direct;
    elements_vec_method1_all{freq_idx} = elements_vec_method1;
    elements_vec_method2_all{freq_idx} = elements_vec_method2;
    
    % Store final results for Methods 1 and 2 (last configuration)
    final_method1_direct(freq_idx) = P_total_method1_direct(end);
    final_method2_direct(freq_idx) = P_total_method2_direct(end);
end

%% Plotting the Convergence and Time Series
% Flag to control saving and closing of figures
saveAndClose = true; % Set to true to save and close figures, false to only display them

saveDirectory = 'C_1C_Original_Dim_Current_Distribution_Fixed_Mesh/Images_Method1_and_Method2/';

if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);
end

numFrequencies = length(frequencies);

for freq_idx = 1:numFrequencies
    f = frequencies(freq_idx);  
    
    fig = figure;
    
    subplot(1, 2, 1);
    hold on;
    plot(elements_vec_method1_all{freq_idx}, P_total_method1_direct_all{freq_idx}, '-^', 'LineWidth', 2, 'DisplayName', 'Method 1 Direct');
    plot(elements_vec_method2_all{freq_idx}, P_total_method2_direct_all{freq_idx}, '-s', 'LineWidth', 2, 'DisplayName', 'Method 2 Direct');
    yline(reference_losses_list(freq_idx), '--r', 'LineWidth', 2, 'DisplayName', 'Reference Losses');
    xlabel('Number of Elements');
    ylabel('Losses');
    title(['Losses vs. Number of Elements at ', num2str(f/1e6), ' MHz']);
    legend('Location', 'best');
    grid on;
    hold off;
    
    subplot(1, 2, 2);
    hold on;
    plot(elements_vec_method1_all{freq_idx}, errors_vec_method1_direct_all{freq_idx}, '-^', 'LineWidth', 2, 'DisplayName', 'Method 1 Direct');
    plot(elements_vec_method2_all{freq_idx}, errors_vec_method2_direct_all{freq_idx}, '-s', 'LineWidth', 2, 'DisplayName', 'Method 2 Direct');
    xlabel('Number of Elements');
    ylabel('Error in Losses');
    title(['Error in Losses vs. Number of Elements at ', num2str(f/1e6), ' MHz']);
    legend('Location', 'best');
    grid on;
    hold off;
    
    if saveAndClose
        savefig(fig, fullfile(saveDirectory, ['Losses_Convergence_', num2str(f/1e6), 'MHz.fig']));
        saveas(fig, fullfile(saveDirectory, ['Losses_Convergence_', num2str(f/1e6), 'MHz.pdf']));
        close(fig);
    end
end

%% Summary Plot for All Frequencies (All Configurations)

% Updated colors for better distinction
colors = [
    0.0000, 0.4470, 0.7410;  % Method 1 Config 1 - Blue
    0.8500, 0.3250, 0.0980;  % Method 1 Config 2 - Red
    0.9290, 0.6940, 0.1250;  % Method 1 Config 3 - Yellow
    0.4940, 0.1840, 0.5560;  % Method 2 Config 1 - Purple
    0.4660, 0.6740, 0.1880;  % Method 2 Config 2 - Green
    0.8500, 0.5250, 0.0980;  % Method 2 Config 3 - Dark Orange
];

markerSize = 7;    % Marker size
lineWidth = 2;     % Line width
transparency = 0.8; % Transparency for lines

% Adjust legend properties
legendFontSize = 3;  % Smaller font size for legends

%% Figure for Method 1 (X-axis Log Scale, Y-axis Linear Scale)
fig1 = figure;

% Plot Losses for Method 1
subplot(3, 1, 1);  
hold on;
plot(frequencies, reference_losses_list, '--r', 'LineWidth', 1.5, 'DisplayName', 'Reference Losses', 'Color', 'k'); % Reference line

for config_idx = 1:numConfigurationsMethod1
    % Plot lines with transparency
    plot(frequencies, cellfun(@(x) x(config_idx), P_total_method1_direct_all), '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
         
    % Overlay scatter plot for markers with the same color, but hide them from the legend
    scatter(frequencies, cellfun(@(x) x(config_idx), P_total_method1_direct_all), markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');  % Exclude from legend
end

xlabel('Frequency (Hz)'); 
ylabel('Losses');
title('Method 1: Losses Across All Frequencies (Y-axis Linear, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);  % Adjust legend font size
grid on;
set(gca, 'XScale', 'log');
hold off;

% Plot Errors for Method 1
subplot(3, 1, 2);  
hold on;
for config_idx = 1:numConfigurationsMethod1
    errors_final_method1_direct = abs(cellfun(@(x) x(config_idx), P_total_method1_direct_all) - reference_losses_list);
    
    % Plot lines with transparency
    plot(frequencies, errors_final_method1_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
    
    % Overlay scatter plot for markers with the same color, but hide them from the legend
    scatter(frequencies, errors_final_method1_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');  % Exclude from legend
end
xlabel('Frequency (Hz)'); 
ylabel('Error in Losses');
title('Method 1: Error in Losses Across All Frequencies (Y-axis Linear, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);  % Adjust legend font size
grid on;
set(gca, 'XScale', 'log');
hold off;

% Plot Percentage Errors for Method 1
subplot(3, 1, 3);
hold on;
for config_idx = 1:numConfigurationsMethod1
    errors_final_method1_direct = abs(cellfun(@(x) x(config_idx), P_total_method1_direct_all) - reference_losses_list);
    percentage_error_method1_direct = (errors_final_method1_direct ./ reference_losses_list) * 100;
    
    % Plot lines with transparency
    plot(frequencies, percentage_error_method1_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct % Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
    
    % Overlay scatter plot for markers with the same color, but hide them from the legend
    scatter(frequencies, percentage_error_method1_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');  % Exclude from legend
end
xlabel('Frequency (Hz)'); 
ylabel('Percentage Error (%)');
title('Method 1: Percentage Error Across All Frequencies (Y-axis Linear, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);  % Adjust legend font size
grid on;
set(gca, 'XScale', 'log');
hold off;

if saveAndClose
    savefig(fig1, fullfile(saveDirectory, 'Method1_Summary_Losses_Errors_Percentage_Linear.fig'));
    exportgraphics(fig1, fullfile(saveDirectory, 'Method1_Summary_Losses_Errors_Percentage_Linear.pdf'));
    close(fig1);
end

%% Figure for Method 1 (X-axis Log Scale, Y-axis Log Scale)
fig2 = figure;

% Plot Losses for Method 1
subplot(3, 1, 1);  
hold on;
plot(frequencies, reference_losses_list, '--r', 'LineWidth', 1.5, 'DisplayName', 'Reference Losses', 'Color', 'k');
for config_idx = 1:numConfigurationsMethod1
    plot(frequencies, cellfun(@(x) x(config_idx), P_total_method1_direct_all), '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
    scatter(frequencies, cellfun(@(x) x(config_idx), P_total_method1_direct_all), markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');  % Exclude from legend
end
xlabel('Frequency (Hz)'); 
ylabel('Losses');
title('Method 1: Losses Across All Frequencies (Y-axis Log, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);  % Adjust legend font size
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

% Plot Errors for Method 1
subplot(3, 1, 2);  
hold on;
for config_idx = 1:numConfigurationsMethod1
    errors_final_method1_direct = abs(cellfun(@(x) x(config_idx), P_total_method1_direct_all) - reference_losses_list);
    plot(frequencies, errors_final_method1_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
    scatter(frequencies, errors_final_method1_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');  % Exclude from legend
end
xlabel('Frequency (Hz)'); 
ylabel('Error in Losses');
title('Method 1: Error in Losses Across All Frequencies (Y-axis Log, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

% Plot Percentage Errors for Method 1
subplot(3, 1, 3);
hold on;
for config_idx = 1:numConfigurationsMethod1
    errors_final_method1_direct = abs(cellfun(@(x) x(config_idx), P_total_method1_direct_all) - reference_losses_list);
    percentage_error_method1_direct = (errors_final_method1_direct ./ reference_losses_list) * 100;
    plot(frequencies, percentage_error_method1_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct % Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
    scatter(frequencies, percentage_error_method1_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Percentage Error (%)');
title('Method 1: Percentage Error Across All Frequencies (Y-axis Log, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

if saveAndClose
    savefig(fig2, fullfile(saveDirectory, 'Method1_Summary_Losses_Errors_Percentage_LogLog.fig'));
    exportgraphics(fig2, fullfile(saveDirectory, 'Method1_Summary_Losses_Errors_Percentage_LogLog.pdf'));
    close(fig2);
end

%% Figure for Method 2 (X-axis Log Scale, Y-axis Linear Scale)
fig3 = figure;

% Plot Losses for Method 2
subplot(3, 1, 1);  
hold on;
plot(frequencies, reference_losses_list, '--r', 'LineWidth', 1.5, 'DisplayName', 'Reference Losses', 'Color', 'k');
for config_idx = 1:numConfigurationsMethod2
    plot(frequencies, cellfun(@(x) x(config_idx), P_total_method2_direct_all), '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]); % Offset for Method 2 colors
    scatter(frequencies, cellfun(@(x) x(config_idx), P_total_method2_direct_all), markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Losses');
title('Method 2: Losses Across All Frequencies (Y-axis Linear, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log');
hold off;

% Plot Errors for Method 2
subplot(3, 1, 2);  
hold on;
for config_idx = 1:numConfigurationsMethod2
    errors_final_method2_direct = abs(cellfun(@(x) x(config_idx), P_total_method2_direct_all) - reference_losses_list);
    plot(frequencies, errors_final_method2_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]);
    scatter(frequencies, errors_final_method2_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Error in Losses');
title('Method 2: Error in Losses Across All Frequencies (Y-axis Linear, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log');
hold off;

% Plot Percentage Errors for Method 2
subplot(3, 1, 3);
hold on;
for config_idx = 1:numConfigurationsMethod2
    errors_final_method2_direct = abs(cellfun(@(x) x(config_idx), P_total_method2_direct_all) - reference_losses_list);
    percentage_error_method2_direct = (errors_final_method2_direct ./ reference_losses_list) * 100;
    plot(frequencies, percentage_error_method2_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct % Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]);
    scatter(frequencies, percentage_error_method2_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Percentage Error (%)');
title('Method 2: Percentage Error Across All Frequencies (Y-axis Linear, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log');
hold off;

if saveAndClose
    savefig(fig3, fullfile(saveDirectory, 'Method2_Summary_Losses_Errors_Percentage_Linear.fig'));
    exportgraphics(fig3, fullfile(saveDirectory, 'Method2_Summary_Losses_Errors_Percentage_Linear.pdf'));
    close(fig3);
end

%% Figure for Method 2 (X-axis Log Scale, Y-axis Log Scale)
fig4 = figure;

% Plot Losses for Method 2
subplot(3, 1, 1);  
hold on;
plot(frequencies, reference_losses_list, '--r', 'LineWidth', 1.5, 'DisplayName', 'Reference Losses', 'Color', 'k');
for config_idx = 1:numConfigurationsMethod2
    plot(frequencies, cellfun(@(x) x(config_idx), P_total_method2_direct_all), '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]);
    scatter(frequencies, cellfun(@(x) x(config_idx), P_total_method2_direct_all), markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Losses');
title('Method 2: Losses Across All Frequencies (Y-axis Log, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

% Plot Errors for Method 2
subplot(3, 1, 2);  
hold on;
for config_idx = 1:numConfigurationsMethod2
    errors_final_method2_direct = abs(cellfun(@(x) x(config_idx), P_total_method2_direct_all) - reference_losses_list);
    plot(frequencies, errors_final_method2_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]);
    scatter(frequencies, errors_final_method2_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Error in Losses');
title('Method 2: Error in Losses Across All Frequencies (Y-axis Log, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

% Plot Percentage Errors for Method 2
subplot(3, 1, 3);
hold on;
for config_idx = 1:numConfigurationsMethod2
    errors_final_method2_direct = abs(cellfun(@(x) x(config_idx), P_total_method2_direct_all) - reference_losses_list);
    percentage_error_method2_direct = (errors_final_method2_direct ./ reference_losses_list) * 100;
    plot(frequencies, percentage_error_method2_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct % Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]);
    scatter(frequencies, percentage_error_method2_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Percentage Error (%)');
title('Method 2: Percentage Error Across All Frequencies (Y-axis Log, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

if saveAndClose
    savefig(fig4, fullfile(saveDirectory, 'Method2_Summary_Losses_Errors_Percentage_LogLog.fig'));
    exportgraphics(fig4, fullfile(saveDirectory, 'Method2_Summary_Losses_Errors_Percentage_LogLog.pdf'));
    close(fig4);
end

%% Figure for Method 1 and Method 2 Combined (X-axis Log Scale, Y-axis Linear Scale)
fig5 = figure;

% Plot Losses for Method 1 and 2
subplot(3, 1, 1);  
hold on;
plot(frequencies, reference_losses_list, '--r', 'LineWidth', 1.5, 'DisplayName', 'Reference Losses', 'Color', 'k');
for config_idx = 1:numConfigurationsMethod1
    plot(frequencies, cellfun(@(x) x(config_idx), P_total_method1_direct_all), '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
    scatter(frequencies, cellfun(@(x) x(config_idx), P_total_method1_direct_all), markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');
end
for config_idx = 1:numConfigurationsMethod2
    plot(frequencies, cellfun(@(x) x(config_idx), P_total_method2_direct_all), '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]);
    scatter(frequencies, cellfun(@(x) x(config_idx), P_total_method2_direct_all), markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Losses');
title('Method 1 and 2: Losses Across All Frequencies (Y-axis Linear, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log');
hold off;

% Plot Errors for Method 1 and 2
subplot(3, 1, 2);  
hold on;
for config_idx = 1:numConfigurationsMethod1
    errors_final_method1_direct = abs(cellfun(@(x) x(config_idx), P_total_method1_direct_all) - reference_losses_list);
    plot(frequencies, errors_final_method1_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
    scatter(frequencies, errors_final_method1_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');
end
for config_idx = 1:numConfigurationsMethod2
    errors_final_method2_direct = abs(cellfun(@(x) x(config_idx), P_total_method2_direct_all) - reference_losses_list);
    plot(frequencies, errors_final_method2_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]);
    scatter(frequencies, errors_final_method2_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Error in Losses');
title('Method 1 and 2: Error in Losses Across All Frequencies (Y-axis Linear, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log');
hold off;

% Plot Percentage Errors for Method 1 and 2
subplot(3, 1, 3);
hold on;
for config_idx = 1:numConfigurationsMethod1
    errors_final_method1_direct = abs(cellfun(@(x) x(config_idx), P_total_method1_direct_all) - reference_losses_list);
    percentage_error_method1_direct = (errors_final_method1_direct ./ reference_losses_list) * 100;
    plot(frequencies, percentage_error_method1_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct % Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
    scatter(frequencies, percentage_error_method1_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');
end
for config_idx = 1:numConfigurationsMethod2
    errors_final_method2_direct = abs(cellfun(@(x) x(config_idx), P_total_method2_direct_all) - reference_losses_list);
    percentage_error_method2_direct = (errors_final_method2_direct ./ reference_losses_list) * 100;
    plot(frequencies, percentage_error_method2_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct % Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]);
    scatter(frequencies, percentage_error_method2_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Percentage Error (%)');
title('Method 1 and 2: Percentage Error Across All Frequencies (Y-axis Linear, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log');
hold off;

if saveAndClose
    savefig(fig5, fullfile(saveDirectory, 'Method1_and_2_Summary_Losses_Errors_Percentage_Linear.fig'));
    exportgraphics(fig5, fullfile(saveDirectory, 'Method1_and_2_Summary_Losses_Errors_Percentage_Linear.pdf'));
    close(fig5);
end

%% Figure for Method 1 and Method 2 Combined (X-axis Log Scale, Y-axis Log Scale)
fig6 = figure;

% Plot Losses for Method 1 and 2
subplot(3, 1, 1);  
hold on;
plot(frequencies, reference_losses_list, '--r', 'LineWidth', 1.5, 'DisplayName', 'Reference Losses', 'Color', 'k');
for config_idx = 1:numConfigurationsMethod1
    plot(frequencies, cellfun(@(x) x(config_idx), P_total_method1_direct_all), '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
    scatter(frequencies, cellfun(@(x) x(config_idx), P_total_method1_direct_all), markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');
end
for config_idx = 1:numConfigurationsMethod2
    plot(frequencies, cellfun(@(x) x(config_idx), P_total_method2_direct_all), '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]);
    scatter(frequencies, cellfun(@(x) x(config_idx), P_total_method2_direct_all), markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Losses');
title('Method 1 and 2: Losses Across All Frequencies (Y-axis Log, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

% Plot Errors for Method 1 and 2
subplot(3, 1, 2);  
hold on;
for config_idx = 1:numConfigurationsMethod1
    errors_final_method1_direct = abs(cellfun(@(x) x(config_idx), P_total_method1_direct_all) - reference_losses_list);
    plot(frequencies, errors_final_method1_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
    scatter(frequencies, errors_final_method1_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');
end
for config_idx = 1:numConfigurationsMethod2
    errors_final_method2_direct = abs(cellfun(@(x) x(config_idx), P_total_method2_direct_all) - reference_losses_list);
    plot(frequencies, errors_final_method2_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]);
    scatter(frequencies, errors_final_method2_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Error in Losses');
title('Method 1 and 2: Error in Losses Across All Frequencies (Y-axis Log, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

% Plot Percentage Errors for Method 1 and 2
subplot(3, 1, 3);
hold on;
for config_idx = 1:numConfigurationsMethod1
    errors_final_method1_direct = abs(cellfun(@(x) x(config_idx), P_total_method1_direct_all) - reference_losses_list);
    percentage_error_method1_direct = (errors_final_method1_direct ./ reference_losses_list) * 100;
    plot(frequencies, percentage_error_method1_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 1 Direct % Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx, :), transparency]);
    scatter(frequencies, percentage_error_method1_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx, :), ...
            'HandleVisibility', 'off');
end
for config_idx = 1:numConfigurationsMethod2
    errors_final_method2_direct = abs(cellfun(@(x) x(config_idx), P_total_method2_direct_all) - reference_losses_list);
    percentage_error_method2_direct = (errors_final_method2_direct ./ reference_losses_list) * 100;
    plot(frequencies, percentage_error_method2_direct, '-', ...
         'LineWidth', lineWidth, 'DisplayName', ['Method 2 Direct % Error (Config ' num2str(config_idx) ')'], ...
         'Color', [colors(config_idx + numConfigurationsMethod1, :), transparency]);
    scatter(frequencies, percentage_error_method2_direct, markerSize * 3, ...
            'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(config_idx + numConfigurationsMethod1, :), ...
            'HandleVisibility', 'off');
end
xlabel('Frequency (Hz)'); 
ylabel('Percentage Error (%)');
title('Method 1 and 2: Percentage Error Across All Frequencies (Y-axis Log, X-axis Log Scale)');
legend('Location', 'best', 'FontSize', legendFontSize);
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');
hold off;

if saveAndClose
    savefig(fig6, fullfile(saveDirectory, 'Method1_and_2_Summary_Losses_Errors_Percentage_LogLog.fig'));
    exportgraphics(fig6, fullfile(saveDirectory, 'Method1_and_2_Summary_Losses_Errors_Percentage_LogLog.pdf'));
    close(fig6);
end
