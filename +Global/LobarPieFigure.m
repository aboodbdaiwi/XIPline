% MATLAB code to generate multi-layered lobar pie chart with proportional defects

% Lobe volume percentages (outer pie)
lobe_data = [40, 20, 45, 45, 50]; % RUL, RML, RLL, LUL, LLL
lobe_labels = {'RUL', 'RML', 'RLL', 'LLL','LUL' };

% Defect percentages (inner pie)
scalefactor = 5;
% defect_data = [4.2, 5.0, 2.7, 3.1, 2.3].*scalefactor;
defect_data = [6, 7, 2, 5, 3].*scalefactor;

% Define outer colors
outer_colors = [200, 255, 200;   % RUL - light green
                255, 255, 150;   % RML - light yellow
                200, 230, 255;   % RLL - light blue
                255, 200, 200;   % LUL - pink
                255, 200, 255] / 255; % LLL - salmon

% Darker shades for defects
inner_colors = [0, 200, 0;   % RUL -  green
                200, 200, 0;   % RML -  yellow
                0, 200, 200;   % RLL -  blue
                255, 0, 0;   % LUL - pink
                255, 0, 255] / 255; % LLL - salmon


% Normalize total for full circle
lobe_data_norm = lobe_data / sum(lobe_data) * 100;

theta0 = pi/2; % starting angle rotated by 90 degrees
r_outer = 1.2; % radius for outer lobes
r_scale = 0.012; % defect pie scale factor (defect % * scale = radius)

figure;
hold on;
axis equal;
axis equal off;
set(gcf,'Color','w');

% Draw outer lobes
theta_accum = 0;
for i = 1:length(lobe_data)
    theta = lobe_data_norm(i) / 100 * 2 * pi;
    % Outer lobe wedge
    t = linspace(theta0 + theta_accum, theta0 + theta_accum + theta, 100);
    x = [0, r_outer * cos(t), 0];
    y = [0, r_outer * sin(t), 0];
    fill(x, y, outer_colors(i,:), 'EdgeColor', 'none');

    % Lobe label placement (midpoint angle)
    angle_mid = theta0 + theta_accum + theta / 2;
    r_text = r_outer * 1.1;
    text(r_text * cos(angle_mid), r_text * sin(angle_mid), lobe_labels{i}, ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 15);

    theta_accum = theta_accum + theta;
end

% Draw defect pies proportionally scaled
theta_accum = 0;
for i = 1:length(defect_data)
    theta = lobe_data_norm(i) / 100 * 2 * pi;
    t = linspace(theta0 + theta_accum, theta0 + theta_accum + theta, 100);
    r_defect = defect_data(i) * r_scale;
    x = [0, r_defect * cos(t), 0];
    y = [0, r_defect * sin(t), 0];
    fill(x, y, inner_colors(i,:), 'EdgeColor', 'none');

    % Place defect percentage text outside defect pie
    angle_mid = theta0 + theta_accum + theta / 2;
    r_text = r_defect * 1.4;
    text(r_text * cos(angle_mid), r_text * sin(angle_mid), sprintf('%.1f%%', defect_data(i)/scalefactor), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 10);

    theta_accum = theta_accum + theta;
end

% title('Lobar Distribution with Ventilation Defects');


