clearvars; clc;
% Time range, constant values, and initial conditions
tstart = 0;
tstop = 12;
tspan = [tstart, tstop];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s

% Defining 3 microcarriers
r_i = [90e-6]; % Microcarrier radii, in m
m = [3.49e-9]; % Microcarrier mass, in kg

g = [0, -9.81]; % m/s^2
gx = g(1);
gy = g(2);
rho_f = 1000; % kg/m^3
n_microcarriers = length(m);
rho_i = m ./ (4 / 3 * pi * r_i.^3); % Density in g/cm^3

x0 = [0.047]; % Initial x position, m
y0 = [0.047]; % Initial y position, m
s0 = [x0; y0]; % Initial positions, m

vx0 = [0.01]; % Initial x velocity, m/s
vy0 = [-0.02]; % Initial y velocity, m/s
v0 = [vx0; vy0]; % Initial velocity, m/s
% Properly format initial conditions for ode45
z0 = [s0; v0];

% Initial velocity values inspired by range seen in CMMC paper, Fig. 5A

% Define the vortex, boundaries, and grid
A_TG = 0.00; % Amplitude of the TG vortex, in m^2/s
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

% Pass the number of microcarriers into the ODE solver
options = odeset('Events', @(t, z) boundary(t, z, floor, ceil, l_wall, r_wall), 'Refine', 100);
odeFun = @(t, z) eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_microcarriers);
[t, z_all] = ode45(odeFun, tspan, z0, options);

tplot = t;
xplot = z_all(:, 1);
yplot = z_all(:, 2);

while t(end) < tstop
    % z vector takes the form [x1, x2, x3,..., y1, y2, y3,..., vx1, vx2, vx3,..., vy1, vy2, vy3]
    % Separate into vectors x0, y0, vx0, vy0
    % Identify indices of positions that have crossed the boundaries
    % Flip the corresponding velocities and adjust those positions to the
    % boundary
    z0 = z_all(end, :);
    % If it hits either wall, replace and reverse x velocity
    if z0(1) <= l_wall
        z0(1) = l_wall;
        z0(3) = -z0(3);
    end
    if z0(1) >= r_wall
        z0(1) = r_wall;
        z0(3) = -z0(3);
    end
    % If it hits the floor or ceiling, replace and reverse y velocity
    if z0(2) <= floor
        z0(2) = floor;
        z0(4) = -z0(4);
    end
    if z0(2) >= ceil
        z0(2) = ceil;
        z0(4) = -z0(4);
    end
    tspan = [t(end), tstop];
    [t, z_all] = ode45(odeFun, tspan, z0, options);
    tplot = [tplot; t];
    xplot = [xplot; z_all(:, 1)];
    yplot = [yplot; z_all(:, 2)];
end

zplot = [tplot, xplot, yplot];

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
    % plot(z_all(:, 2*i - 1), z_all(:, 2*i), 'color', "black");
    plot(xplot, yplot, 'color', 'black')
    plot(s0(i), s0(i+n_microcarriers), "or", "LineWidth", 1.25)
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
title('Microcarrier trajectories in counterflow vortex with buoyancy/gravity and drag')

% Boundary interaction event function
function [value, isTerminal, direction] = boundary(t, z, floor, ceil, l_wall, r_wall)
value(1) = (z(1) - l_wall);
value(2) = (z(1) - r_wall);
value(3) = (z(2) - floor);
value(4) = (z(2) - ceil);
isTerminal(value ~= 0) = 1;
direction = 0;
end

% Fully vectorized numerical solution
function dzdt = eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_microcarriers)
    dzdt = zeros(4 * n_microcarriers, 1);  % Initialize the output vector for all microcarriers
    % z vector takes the form [x1, x2, x3, y1, y2, y3, vx1, vx2, vx3, vy1, vy2, vy3]
    % Extract positions, velocities, and other parameters for each particle
    x = z(1:n_microcarriers);
    y = z(n_microcarriers+1:2*n_microcarriers);
    vx = z(2*n_microcarriers+1:3*n_microcarriers);  % x-velocity of particle i
    vy = z(3*n_microcarriers+1:end);      % y-velocity of particle i

    % Interpolating to find background flow at the particle's position
    vfx = interp2(X, Y, vfx_disc, x, y, 'linear', 0);
    vfy = interp2(X, Y, vfy_disc, x, y, 'linear', 0);

    % Position derivatives
    dxdt = vx;
    dydt = vy;

    % Drag and gravity force equations
    % dvxdt = -(6 * pi * mu * r_i .* (vx - vfx)) ./ m + gx * (1 - rho_f ./ rho_i);
    % dvydt = -(6 * pi * mu * r_i .* (vy - vfy)) ./ m + gy * (1 - rho_f ./ rho_i);
    dvxdt = 0;
    dvydt = 0;

    % Store the position and velocity derivatives
    dzdt(1:n_microcarriers) = dxdt;
    dzdt(n_microcarriers+1:2*n_microcarriers) = dydt;
    dzdt(2*n_microcarriers+1:3*n_microcarriers) = dvxdt;
    dzdt(3*n_microcarriers+1:end) = dvydt;
end
deltat = [0];
for i=2:length(tplot)
    deltat(end+1) = tplot(i) - tplot(i-1);
end
figure
plot(tplot, deltat)
grid on
xlabel('time (s)')
ylabel('time step size (s)')
title('Time step size throughout the simulation')