clear all, close all, clc
% Code to generate image of source points for the presentation

%% Only one source point
R = 1;  
x_prime = [0.5, 0.5, 0.5];  % Source point inside the sphere
lambda = R / norm(x_prime);  
image_charge_point = lambda^2 * x_prime;  % Image charge point

% Sphere
[phi, theta] = meshgrid(linspace(0, 2*pi, 50), linspace(0, pi, 50));
x_sphere = R * sin(theta) .* cos(phi);
y_sphere = R * sin(theta) .* sin(phi);
z_sphere = R * cos(theta);


figure;
hold on;
surf(x_sphere, y_sphere, z_sphere, 'FaceAlpha', 0.3, 'EdgeColor', 'none');  
colormap winter;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Visualization of Sphere, Source Point, and Image Charge Point');
axis equal;
grid on;

scatter3(x_prime(1), x_prime(2), x_prime(3), 100, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'r', 'DisplayName', 'Source Point');
scatter3(image_charge_point(1), image_charge_point(2), image_charge_point(3), 100, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'b', 'DisplayName', 'Image Charge Point');
line([x_prime(1), image_charge_point(1)], [x_prime(2), image_charge_point(2)], [x_prime(3), image_charge_point(3)], 'Color', 'k', 'LineStyle', '--');


legend('Sphere', 'Source Point', 'Image Charge Point', 'Location', 'best');
view(3);  
hold off;


%% More source points
R = 1;  
source_points = [0.5, 0.5, 0.5; -0.3, -0.3, 0.2; 0.1, -0.4, 0.6];  
num_sources = size(source_points, 1);
lambda_values = R ./ vecnorm(source_points, 2, 2);  
image_charge_points = (lambda_values.^2) .* source_points;  

[phi, theta] = meshgrid(linspace(0, 2*pi, 50), linspace(0, pi, 50));
x_sphere = R * sin(theta) .* cos(phi);
y_sphere = R * sin(theta) .* sin(phi);
z_sphere = R * cos(theta);

figure;
hold on;
surf(x_sphere, y_sphere, z_sphere, 'FaceAlpha', 0.3, 'EdgeColor', 'none');  
colormap winter;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Visualization of Sphere, Multiple Source Points, and Image Charge Points');
axis equal;
grid on;

for i = 1:num_sources
    scatter3(source_points(i, 1), source_points(i, 2), source_points(i, 3), 100, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'r');
    scatter3(image_charge_points(i, 1), image_charge_points(i, 2), image_charge_points(i, 3), 100, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor', 'b');
    line([source_points(i, 1), image_charge_points(i, 1)], [source_points(i, 2), image_charge_points(i, 2)], [source_points(i, 3), image_charge_points(i, 3)], 'Color', 'k', 'LineStyle', '--');
end

legend('Sphere', 'Source Points', 'Image Charge Points', 'Location', 'best');
view(3);  
hold off;
