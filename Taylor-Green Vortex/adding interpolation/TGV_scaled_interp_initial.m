clearvars; close all; clc;
% Time range, constant values, and initial conditions
tspan = [0, 10];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s

% Defining 3 microcarriers
% Take another look at these params
r_i = [90e-6, 90e-6, 90e-6]; % Microcarrier radii, in m
m = [2.79e-9, 2.79e-9, 2.79e-9]; % Microcarrier mass, in kg

g = [0, -9.81]; % m/s^2
gx = g(1);
gy = g(2);
rho_f = 1000; % kg/m^3
n_microcarriers = length(m);
rho_i = m ./ (4 / 3 * pi * r_i.^3); % Density in g/cm^3

x0 = [0.01, 0.01; 0.02, 0.02; 0.03, 0.03]; % Initial position, m
v0 = [0, 0; 0, 0; 0, 0]; % Initial velocity, m/s
% Initial velocity values inspired by range seen in CMMC paper, Fig. 5A

% Defining background velocity, x and y cpts
A_TG = 0.005; % Amplitude of the TG vortex, in m^2/s
L_TG = 0.1; % Length of the TG vortex, in m
n_step = 10;
[X,Y] = meshgrid(-L_TG / 2 : L_TG / n_step : L_TG / 2);
psi_TG = -A_TG*cos(pi * X / L_TG) .* cos(pi * Y / L_TG);

vfx_disc = (A_TG * pi / L_TG) * cos(pi * X / L_TG) .* sin(pi * Y / L_TG);
vfy_disc = -(A_TG * pi / L_TG) * sin(pi * X / L_TG) .* cos(pi * Y / L_TG);

% Properly formatting initial conditions for ode45
x0_all = reshape(x0', [], 1);
v0_all = reshape(v0', [], 1);  % This converts the v0 matrix to a column vector
z0_all = [x0_all; v0_all];

% Calling ode45
odeFun = @(t, z) eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_microcarriers);
[t, z_all] = ode45(odeFun, tspan, z0_all);

% Plotting the results
figure
hold on
[curlz,cav] = curl(X, Y, vfx_disc, vfy_disc);
c = pcolor(X, Y, curlz);
c.FaceColor = 'interp';

% Define two colors for the gradient
color1 = [0.3010 0.7450 0.9330]; %  (RGB)
color2 = [1 1 1]; %            (RGB)

% Create a custom colormap with a smooth transition between these two colors
numColors = 128; % Number of colors in the colormap
cmap = [linspace(color1(1), color2(1), numColors)' ...
        linspace(color1(2), color2(2), numColors)' ...
        linspace(color1(3), color2(3), numColors)'];

colormap(cmap);
cb = colorbar;
yl = ylabel(cb,'Vorticity','FontSize',10,'Rotation',270);

q = quiver(X, Y, vfx_disc, vfy_disc, 'k');
col_list = ["red", "blue", "magenta"];
for i=1:n_microcarriers
    plot(z_all(:, 2*i - 1), z_all(:, 2*i), 'color', col_list(mod(i, numel(col_list)) + 1));
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
title('Microcarrier trajectories with drag and gravity; velocity and vorticity of TG vortex', 'FontSize', 14)
grid on
legend({'vorticity', 'velocity field', 'microcarrier 1','microcarrier 2', 'microcarrier 3'})

% Calculating the derivatives in the eom function
function dzdt = eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_microcarriers)
    dzdt = zeros(4 * n_microcarriers, 1);  % Initialize the output vector for all microcarriers
    for i = 1:n_microcarriers
        % z vector takes the form [x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3]
        % Extract positions, velocities, and other parameters for each particle
        x = z(2*i - 1);
        y = z(2*i);
        vx = z(2*n_microcarriers + 2*i - 1);  % x-velocity of microcarrier i
        vy = z(2*n_microcarriers + 2*i);      % y-velocity of microcarrier i
        r = r_i(i);       % radius of microcarrier i
        mi = m(i);        % mass of microcarrier i
        rhoi = rho_i(i);

        % Interpolating to find background flow at the particle's position
        vfx = interp2(X, Y, vfx_disc, x, y);
        vfy = interp2(X, Y, vfy_disc, x, y);

        % Position derivatives
        dxdt = vx;
        dydt = vy;
        
        % Drag and gravity force equations
        dvxdt = -(6 * pi * mu * r * (vx - vfx)) / mi + gx * (1 - rho_f / rhoi);
        dvydt = -(6 * pi * mu * r * (vy - vfy)) / mi + gy * (1 - rho_f / rhoi);
        
        % Store the position and velocity derivatives
        dzdt(2*i - 1) = dxdt;
        dzdt(2*i) = dydt;
        dzdt(2*n_microcarriers + 2*i - 1) = dvxdt;
        dzdt(2*n_microcarriers + 2*i) = dvydt;
    end
end

% Plot vorticity field, velocity vector field, cell trajectories all on one
% plot
% (probably much larger distance than cell size)
% Instead of this analytical form, put a function for interpolation
% Define a grid of x and y, calculate using analytical form at these
% discrete points
% Interpolate between those points to find particle velocity at any
% position
% Also, have the analytical form take more realistic parameters based on
% formula