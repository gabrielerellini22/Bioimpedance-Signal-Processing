clear; clc; close all;

%% 1. ACTUAL GEOMETRIC PARAMETERS (mm)
% IWAKI 24-Well Plate Specifications
D_Well_Outer = 15.6;     
R_Well = D_Well_Outer / 2;
H_Well = 10.0;           

% Greiner ThinCert Insert Specifications
Membrane_Area = 0.336;   
D_Insert = 2 * sqrt((Membrane_Area * 100) / pi); 
R_Insert = D_Insert / 2;
Y_Membrane = 3.0;        
H_Insert = 9.5;          
Wall_Thick = 0.4;        

% STX01 Electrode Specifications
Thickness = 1.0;         
Gap_Legs = 5.0;          
X_Wall_Ref = R_Insert;   
X_Leg_In = X_Wall_Ref - 2.0;  
X_Leg_Out = X_Wall_Ref + 3.0; 

% Electrode Tips and Active Centers
Y_Tip_In = Y_Membrane + 2.0; 
Y_Tip_Out = 1.0;
Y_Active_In = Y_Tip_In + 0.5;
Y_Active_Out = Y_Tip_Out + 0.5;

% Charges
Q_Source = 1; Q_Sink = -1;

%% 2. FIELD CALCULATION (Extended to cover full background)
% NOTE: Grid spans from -15 to 15 (wider than view) to avoid white borders
[x, y] = meshgrid(linspace(-15, 15, 800), linspace(-5, 20, 800)); 
Xc_In = X_Leg_In - 0.5; Yc_In = Y_Active_In;
Xc_Out = X_Leg_Out + 0.5; Yc_Out = Y_Active_Out;

r_src = sqrt((x - Xc_In).^2 + (y - Yc_In).^2 + 0.1);
r_sink = sqrt((x - Xc_Out).^2 + (y - Yc_Out).^2 + 0.1);

Ex = Q_Source .* (x - Xc_In)./r_src.^3 + Q_Sink .* (x - Xc_Out)./r_sink.^3;
Ey = Q_Source .* (y - Yc_In)./r_src.^3 + Q_Sink .* (y - Yc_Out)./r_sink.^3;
J_mag = sqrt(Ex.^2 + Ey.^2);

%% 3. GRAPHICAL PLOT
% High-resolution figure
fig = figure('Color', 'w', 'Name', 'STX01 Actual Configuration', ...
    'Position', [50 50 1200 900], 'Visible', 'on'); 
hold on;

% A. HEATMAP (Full Background)
J_log = log10(J_mag);
contourf(x, y, J_log, 120, 'LineColor', 'none'); 
colormap(jet);
caxis([-1.5 1.5]);
cb = colorbar;
cb.Label.String = 'Norm. Current Density (Log Scale)';
cb.Label.FontSize = 12;
cb.Label.FontWeight = 'bold';

% B. STREAMLINES
hS = streamslice(x, y, Ex, Ey, 3); 
set(hS, 'Color', [1 1 1], 'LineWidth', 0.6);

% C. PHYSICAL STRUCTURES
Color_Plastic = [0.2 0.2 0.2]; 
Color_Insert = [0.6 0.6 0.6]; 
Color_Probe = [0 0 0]; 

% IWAKI Well
fill([-R_Well-0.5 -R_Well -R_Well -R_Well-0.5], [0 0 H_Well H_Well], Color_Plastic, 'EdgeColor', 'none');
fill([R_Well R_Well+0.5 R_Well+0.5 R_Well], [0 0 H_Well H_Well], Color_Plastic, 'EdgeColor', 'none');
fill([-R_Well-0.5 R_Well+0.5 R_Well+0.5 -R_Well-0.5], [-0.5 -0.5 0 0], Color_Plastic, 'EdgeColor', 'none');

% Greiner Insert
fill([-R_Insert -R_Insert+Wall_Thick -R_Insert+Wall_Thick -R_Insert], [Y_Membrane Y_Membrane H_Insert H_Insert], Color_Insert, 'EdgeColor', 'none');
fill([R_Insert R_Insert-Wall_Thick R_Insert-Wall_Thick R_Insert], [Y_Membrane Y_Membrane H_Insert H_Insert], Color_Insert, 'EdgeColor', 'none');
fill([-R_Insert-1 -R_Insert+Wall_Thick -R_Insert+Wall_Thick -R_Insert-1], [H_Insert H_Insert H_Insert+0.5 H_Insert+0.5], Color_Insert, 'EdgeColor', 'none');
fill([R_Insert-Wall_Thick R_Insert+1 R_Insert+1 R_Insert-Wall_Thick], [H_Insert H_Insert H_Insert+0.5 H_Insert+0.5], Color_Insert, 'EdgeColor', 'none');

% Membrane
line([-R_Insert+Wall_Thick R_Insert-Wall_Thick], [Y_Membrane Y_Membrane], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 3);

% STX01 Electrode
rectangle('Position', [X_Leg_In-Thickness/2, Y_Tip_In, Thickness, 14-Y_Tip_In], 'FaceColor', Color_Probe, 'EdgeColor', 'none');
rectangle('Position', [X_Leg_Out-Thickness/2, Y_Tip_Out, Thickness, 14-Y_Tip_Out], 'FaceColor', Color_Probe, 'EdgeColor', 'none');
Block_Y = 12.5; 
rectangle('Position', [X_Leg_In-Thickness, Block_Y, (X_Leg_Out-X_Leg_In)+2*Thickness, 2.5], 'FaceColor', [0.1 0.1 0.3], 'EdgeColor', 'none');

% D. MARKERS AND LABELS
plot(Xc_In, Yc_In, 'wo', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 1.2);
plot(Xc_Out, Yc_Out, 'wo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'LineWidth', 1.2);
plot(X_Leg_In+0.5, Yc_In, 'ys', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');
plot(X_Leg_Out-0.5, Yc_Out, 'ys', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');

% Text labels
text(0, Y_Membrane + 0.8, '\bf MEMBRANE', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', 10);
text(0, Y_Membrane - 1.2, {'FLOW BENEATH', 'INSERT WALL'}, 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
text(X_Leg_In, Y_Tip_In-0.6, '\bf Apical', 'Color', 'k', 'HorizontalAlignment', 'center', 'FontSize', 10);
text(X_Leg_Out, Y_Tip_Out-0.6, '\bf Basal', 'Color', 'k', 'HorizontalAlignment', 'center', 'FontSize', 10);
text(-5.5, H_Well/2, 'IWAKI Well', 'Color', 'w', 'FontSize', 11, 'Rotation', 90, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(-R_Insert+1.5, H_Insert/2+1, 'Greiner Insert', 'Color', 'k', 'FontSize', 11, 'Rotation', 90, 'HorizontalAlignment', 'center');

% E. LEGEND (Top Left - Optimized)
Leg_X = -11.0; Leg_Y = 13.0; 
rectangle('Position', [Leg_X, Leg_Y, 4.8, 2.5], 'FaceColor', [1 1 1 0.9], 'EdgeColor', 'k', 'LineWidth', 0.8);
plot(Leg_X + 0.6, Leg_Y + 1.9, 'wo', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
text(Leg_X + 1.4, Leg_Y + 1.9, 'I+ (Source)', 'FontSize', 11, 'Color', 'k', 'FontName', 'Arial');
plot(Leg_X + 0.6, Leg_Y + 1.2, 'wo', 'MarkerFaceColor', 'b', 'MarkerSize', 9);
text(Leg_X + 1.4, Leg_Y + 1.2, 'I- (Sink)', 'FontSize', 11, 'Color', 'k', 'FontName', 'Arial');
plot(Leg_X + 0.6, Leg_Y + 0.5, 'ys', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
text(Leg_X + 1.4, Leg_Y + 0.5, 'V (Voltage)', 'FontSize', 11, 'Color', 'k', 'FontName', 'Arial');

% Titles and Axes
title({'Current Density Distribution (J)', '24-Well & STX01 Actual Geometry'}, 'FontSize', 16, 'FontName', 'Arial');
xlabel('Distance (mm)', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Height (mm)', 'FontSize', 12, 'FontName', 'Arial');
axis equal;

% LOCKED VIEW (To crop grid edges and maintain full background)
xlim([-12 12]); 
ylim([-2 16]); 
set(gca, 'Layer', 'top', 'Box', 'on', 'FontSize', 11, 'LineWidth', 1.0);
hold off;

%% HIGH-RESOLUTION SAVE
try
    desktopPath = fullfile(getenv('USERPROFILE'), 'Desktop');
    filename = fullfile(desktopPath, 'STX01_Simulation_Thesis_FullBlue.jpg');
    if ~exist(desktopPath, 'dir')
        filename = 'STX01_Simulation_Thesis_FullBlue.jpg';
    end
    print(fig, filename, '-djpeg100', '-r300');
    fprintf('Image saved: %s\n', filename);
catch
    print(fig, 'STX01_Simulation_Thesis_FullBlue.jpg', '-djpeg100', '-r300');
end