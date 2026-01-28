clearvars; clc; close all;
% Time range, constant values, and initial conditions
tstart = 0;
tstop = 10;
tspan = [tstart, tstop];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s

% Defining 3 microcarriers
% r_i = [90e-6; 90e-6; 90e-6; 90e-6; 90e-6]; % Microcarrier radii, in m
% m = [2.79e-9; 2.79e-9; 2.79e-9; 2.79e-9; 2.79e-9]; % Microcarrier mass, in kg
r_i = [90e-6; 90e-6];
m = [2.79e-9; 2.79e-9];

g = [0, -9.81]; % m/s^2
gx = g(1);
gy = g(2);
rho_f = 1000; % kg/m^3
n_microcarriers = length(m);
rho_i = m ./ (4 / 3 * pi * r_i.^3); % Density in g/cm^3

% x0 = [0.01; 0.02; 0.1; 0.08; 0.08]; % Initial x position of each particle, m
% y0 = [0.01; 0.02; 0.02; 0.04; 0.025]; % Initial y positions, m
x0 = [0.09; 0.0901]; % Initial x position, m
y0 = [-0.01; -0.01]; % Initial y position, m
s0 = [x0; y0]; % Initial positions, m

% vx0 = [0; 0; 0; 0; 0]; % Initial x velocities of each particle, m/s
% vy0 = [0; 0; 0; 0; 0]; % Initial y velocities, m/s
vx0 = [0; 0]; % Initial x velocity, m/s
vy0 = [0; 0]; % Initial y velocity, m/s
v0 = [vx0; vy0]; % Initial velocity, m/s

% x0 = [0.045; 0.045]; % Initial position, m
% v0 = [0.1; -0.2]; % Initial velocity, m/s

% Initial velocity values inspired by range seen in CMMC paper, Fig. 5A

% Define the vortex, boundaries, and grid
A_TG = 0.005; % Amplitude of the TG vortex, in m^2/s
L_TG = 0.1; % Length of the TG vortex, in m

floor = -L_TG / 2; % The low y boundary
ceil = L_TG / 2; % The high y boundary (not as concerned with this in our physical setup)
l_wall = -L_TG / 2; % The left x boundary coordinate
r_wall = 3 * L_TG / 2; % The right x boundary

n_step = 199;
grid_sl = L_TG / n_step; % Side length of each square in the grid

[X,Y] = meshgrid(l_wall : grid_sl : r_wall, floor : grid_sl : ceil);

psi_TG = -A_TG*cos(pi * X / L_TG) .* cos(pi * Y / L_TG);

% Define background velocity field (x and y cpts)
vfx_disc = (A_TG * pi / L_TG) * cos(pi * X / L_TG) .* sin(pi * Y / L_TG);
vfy_disc = -(A_TG * pi / L_TG) * sin(pi * X / L_TG) .* cos(pi * Y / L_TG);

% Properly format initial conditions for ode45
z0 = [s0; v0];

% Pass the number of microcarriers into the ODE solver
odeFun = @(t, z) eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_microcarriers);
[t, z_all] = ode45(odeFun, tspan, z0);

figure
hold on
[curlz,cav] = curl(X, Y, vfx_disc, vfy_disc);
c = pcolor(X, Y, curlz);
c.FaceColor = 'interp';
c.EdgeColor = 'none';

% Define three colors for the gradient
% color_neg = [0.3010 0.7450 0.9330];   % Blue for negative vorticity
% color_zero = [0.6660 0.8740 0.3880];  % Green for zero vorticity
color_pos = [0.9290 0.6940 0.1250];   % Yellow for positive vorticity

color_pos = [0.9290 0.6940 0.1250];   % Yellow for positive vorticity
color_neg = [0.00392156862745098 0.44313725490196076 0.7372549019607844];   % Blue for negative vorticity
color_zero = [0.8480392156862745 0.8911764705882353 0.9225490196078431];  % White for zero vorticity

% Create a custom colormap with a smooth transition between these three colors
numColors = 256;  % Number of colors in the colormap (higher resolution)
cmap_neg_to_zero = [linspace(color_neg(1), color_zero(1), numColors/2)', ...
                    linspace(color_neg(2), color_zero(2), numColors/2)', ...
                    linspace(color_neg(3), color_zero(3), numColors/2)'];
                
cmap_zero_to_pos = [linspace(color_zero(1), color_pos(1), numColors/2)', ...
                    linspace(color_zero(2), color_pos(2), numColors/2)', ...
                    linspace(color_zero(3), color_pos(3), numColors/2)'];

% Combine the two halves to form the full colormap
cmap = [cmap_neg_to_zero; cmap_zero_to_pos];

colormap(cmap);
cb = colorbar;
% fontsize(18,"points")
yl = ylabel(cb,'Vorticity','Rotation',270);

n_step_q = 20;
[X_q,Y_q] = meshgrid(l_wall : L_TG / n_step_q : r_wall, floor : L_TG / n_step_q : ceil);
vfx_disc_q = (A_TG * pi / L_TG) * cos(pi * X_q / L_TG) .* sin(pi * Y_q / L_TG);
vfy_disc_q = -(A_TG * pi / L_TG) * sin(pi * X_q / L_TG) .* cos(pi * Y_q / L_TG);
q = quiver(X_q, Y_q, vfx_disc_q, vfy_disc_q, 'k');
col_list = ["black", "red", "blue", "magenta"];

x_plot = z_all(:,1:n_microcarriers);
y_plot = z_all(:,n_microcarriers+1:2*n_microcarriers);
vx_plot = z_all(:,2*n_microcarriers+1:3*n_microcarriers);
vy_plot = z_all(:,3*n_microcarriers+1:4*n_microcarriers);
v_plot = sqrt(vx_plot.^2+vy_plot.^2);
dist_plot = sqrt((x_plot(:, 2) - x_plot(:, 1)).^2 + (y_plot(:, 2) - y_plot(:, 1)).^2);

for i=1:n_microcarriers
    plot(x_plot(:,i), y_plot(:,i), 'color', col_list(i));
    plot(s0(i), s0(i+n_microcarriers), "or", "LineWidth", 1.25)
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
title('No cell bonding implemented, used as a comparison')
%title('Microcarrier trajectories with boundary reflections, no background flow or forces')

% Fully vectorized numerical solution
function dzdt = eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_microcarriers)
    dzdt = zeros(4 * n_microcarriers, 1);  % Initialize the output vector for all microcarriers
                % 4*n_microcarriers rows, 1 column
    % z vector takes the form [x1, x2, x3, y1, y2, y3, vx1, vx2, vx3, vy1, vy2, vy3]
    % Extract positions, velocities, and other parameters for each particle
    x = z(1:n_microcarriers);
    y = z(n_microcarriers+1:2*n_microcarriers);
    vx = z(2*n_microcarriers+1:3*n_microcarriers);  % x-velocity of particle i
    vy = z(3*n_microcarriers+1:end);      % y-velocity of particle i

    % Interpolating to find background flow at the particle's position
    vfx = interp2(X, Y, vfx_disc, x, y);
    vfy = interp2(X, Y, vfy_disc, x, y);

    % Position derivatives
    dxdt = vx;
    dydt = vy;

    % Drag and gravity force equations
    dvxdt = -(6 * pi * mu * r_i .* (vx - vfx)) ./ m + gx * (1 - rho_f ./ rho_i);
    dvydt = -(6 * pi * mu * r_i .* (vy - vfy)) ./ m + gy * (1 - rho_f ./ rho_i);

    % Store the position and velocity derivatives
    dzdt(1:n_microcarriers) = dxdt;
    dzdt(n_microcarriers+1:2*n_microcarriers) = dydt;
    dzdt(2*n_microcarriers+1:3*n_microcarriers) = dvxdt;
    dzdt(3*n_microcarriers+1:end) = dvydt;
end

% Plotting time step size over time of simulation
deltat = [0; diff(t)];
figure
plot(t, deltat)
grid on
xlabel('time (s)')
ylabel('time step size (s)')
title('Time step size throughout the simulation — no bonding')

% Plotting velocity over time of simulation
% figure
% for i=1:n_microcarriers
%     hold on
%     plot(t, v_plot(:, i))
%     grid on
%     xlabel('time (s)')
%     ylabel('velocity (m/s)')
% end

% Plotting distance between particles over time
figure
hold on
plot(t, dist_plot)
yline(2*r_i(1))
grid on
xlabel('time (s)')
ylabel('particle distance (m)')
title('distance between particles 1 and 2 – no bonding')