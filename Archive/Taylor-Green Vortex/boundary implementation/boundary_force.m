clearvars; close all; clc;
% Time range, constant values, and initial conditions
tspan = [0, 10];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s

% Defining 3 microcarriers
r_i = [90e-6]; % Microcarrier radii, in m
m = [2.79e-9]; % Microcarrier mass, in kg

g = [0, 0]; % m/s^2
gx = g(1);
gy = g(2);
rho_f = 1000; % kg/m^3
n_microcarriers = length(m);
rho_i = m ./ (4 / 3 * pi * r_i.^3); % Density in g/cm^3

x0 = [0.00, -0.04]; % Initial position, m
v0 = [1, -2]; % Initial velocity, m/s
% Initial velocity values inspired by range seen in CMMC paper, Fig. 5A

% Define the vortex, boundaries, and grid
A_TG = 0; % Amplitude of the TG vortex, in m^2/s
L_TG = 0.1; % Length of the TG vortex, in m

floor = -L_TG / 2; % The low y boundary
ceil = L_TG / 2; % The high y boundary (not as concerned with this in our physical setup)
l_wall = -L_TG / 2; % The left x boundary coordinate
r_wall = L_TG / 2; % The right x boundary

n_step = 199;
grid_sl = L_TG / n_step; % Side length of each square in the grid

[X,Y] = meshgrid(l_wall : grid_sl : r_wall, floor : grid_sl : ceil);

psi_TG = -A_TG*cos(pi * X / L_TG) .* cos(pi * Y / L_TG);

% Define background velocity field (x and y cpts)
vfx_disc = (A_TG * pi / L_TG) * cos(pi * X / L_TG) .* sin(pi * Y / L_TG);
vfy_disc = -(A_TG * pi / L_TG) * sin(pi * X / L_TG) .* cos(pi * Y / L_TG);

% Properly format initial conditions for ode45
x0_all = reshape(x0', [], 1);
v0_all = reshape(v0', [], 1);  % This converts the v0 matrix to a column vector
z0_all = [x0_all; v0_all];

% Pass the number of microcarriers into the ODE solver
odeFun = @(t, z) eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_microcarriers, l_wall, r_wall, floor, ceil);
[t, z_all] = ode45(odeFun, tspan, z0_all);

% Try implementing RK4 method in MATLAB – this goes sequentially in
% timesteps

figure
hold on
[curlz,cav] = curl(X, Y, vfx_disc, vfy_disc);
c = pcolor(X, Y, curlz);
c.FaceColor = 'interp';
c.EdgeColor = 'none';

% Define three colors for the gradient
color_neg = [0.3010 0.7450 0.9330];   % Blue for negative vorticity
color_zero = [0.6660 0.8740 0.3880];  % Green for zero vorticity
color_pos = [0.9290 0.6940 0.1250];   % Yellow for positive vorticity

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
yl = ylabel(cb,'Vorticity','FontSize',10,'Rotation',270);

n_step_q = 20;
[X_q,Y_q] = meshgrid(l_wall : L_TG / n_step_q : r_wall, floor : L_TG / n_step_q : ceil);
vfx_disc_q = (A_TG * pi / L_TG) * cos(pi * X_q / L_TG) .* sin(pi * Y_q / L_TG);
vfy_disc_q = -(A_TG * pi / L_TG) * sin(pi * X_q / L_TG) .* cos(pi * Y_q / L_TG);
q = quiver(X_q, Y_q, vfx_disc_q, vfy_disc_q, 'k');
%col_list = ["red", "blue", "magenta"];
for i=1:n_microcarriers
    plot(z_all(:, 2*i - 1), z_all(:, 2*i), 'color', "black");
    plot(x0_all(2*i - 1), x0_all(2*i), "or", "LineWidth", 1.25)
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
title('Microcarrier trajectories in counterflow vortex with buoyancy/gravity and drag')

% Fully vectorized numerical solution
function dzdt = eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_microcarriers, l_wall, r_wall, floor, ceil)
    dzdt = zeros(4 * n_microcarriers, 1);  % Initialize the output vector for all microcarriers
    % z vector takes the form [x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3]
    % Extract positions, velocities, and other parameters for each particle
    x = z(1:2:2*n_microcarriers);
    y = z(2:2:2*n_microcarriers);
    vx = z(2*n_microcarriers+1:2:end);
    vy = z(2*n_microcarriers+2:2:end);

    % Flipping velocities of particles that are outside the boundary and 
    % moving away from the boundary
    bound_left = (x < l_wall) & (vx < 0);
    bound_right = (x > r_wall) & (vx > 0);
    x_flip = (bound_left | bound_right);
    vx(x_flip) = -vx(x_flip);

    bound_floor = (y < floor) & (vy < 0);
    bound_ceil = (y > ceil) & (vy > 0);
    y_flip = (bound_floor | bound_ceil);
    vy(y_flip) = -vy(y_flip);

    % Putting the out-of-bounds particles back on the boundary so the 
    % derivatives can be computed
    x(bound_left) = l_wall;
    x(bound_right) = r_wall;
    y(bound_floor) = floor;
    y(bound_ceil) = ceil;

    % Updating the z vector to reflect these changes
    z(1:2:2*n_microcarriers) = x;
    z(2:2:2*n_microcarriers) = y;
    z(2*n_microcarriers+1:2:end) = vx;
    z(2*n_microcarriers+2:2:end) = vy;

    % Interpolating to find background flow at the particle's position
    vfx = interp2(X, Y, vfx_disc, x, y);
    vfy = interp2(X, Y, vfy_disc, x, y);

    % Position derivatives
    dxdt = vx;
    dydt = vy;

    % Drag and gravity force equations
    dvxdt = -(6 * pi * mu * r_i .* (vx - vfx))*0 ./ m + gx * (1 - rho_f ./ rho_i);
    dvydt = -(6 * pi * mu * r_i .* (vy - vfy))*0 ./ m + gy * (1 - rho_f ./ rho_i);

    % Store the position and velocity derivatives
    dzdt(1:2:2*n_microcarriers) = dxdt;
    dzdt(2:2:2*n_microcarriers) = dydt;
    dzdt(2*n_microcarriers+1:2:end) = dvxdt;
    dzdt(2*n_microcarriers+2:2:end) = dvydt;
end