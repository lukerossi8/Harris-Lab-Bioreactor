clearvars; close all; clc;
% Time range, constant values, and initial conditions
tspan = [0, 20];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s

% Defining 3 microcarriers
r_i = [90e-6, 80e-6, 70e-6]; % Microcarrier radii, in m
m = [2.79e-9, 1.79e-9, 0.79e-9]; % Microcarrier mass, in kg

g = [0, -9.81]; % m/s^2
gx = g(1);
gy = g(2);
rho_f = 1000; % kg/m^3
n_microcarriers = length(m);
rho_i = [];
for i=1:n_microcarriers
    rho_i = [rho_i, m(i) / (4 / 3 * pi * r_i(i)^3)]; % Density in g/cm^3
end

x0 = [0.01, 0.01; -0.02, 0.02; -0.015, -0.005]; % Initial position, m
v0 = [0.02, -0.02; 0.05, 0.03; -0.06, 0.02]; % Initial velocity, m/s
% Initial velocity values inspired by range seen in CMMC paper, Fig. 5A

% Define background velocity, x and y cpts
vfx = @(x, y) sin(x) * cos(y);
vfy = @(x, y) -cos(x) * sin(y);

% Properly format initial conditions for ode45
x0_all = reshape(x0', [], 1);
v0_all = reshape(v0', [], 1);  % This converts the v0 matrix to a column vector
z0_all = [x0_all; v0_all];

% Pass the number of microcarriers into the ODE solver
odeFun = @(t, z) eom(t, z, mu, r_i, vfx, vfy, m, gx, gy, rho_f, rho_i, n_microcarriers);
[t, z_all] = ode45(odeFun, tspan, z0_all);

% Plotting x velocity numerically over time
figure
hold on
for i=1:n_microcarriers
    plot(t, z_all(:, 2*n_microcarriers + 2*i - 1));
end
xlabel('Time (s)')
ylabel('x direction velocity (m/s)')
title('x velocity vs time with buoyancy/gravity, drag, and TG vortex')
grid on
legend({'microcarrier 1','microcarrier 2', 'microcarrier 3'})

% Plotting y velocity numerically over time
figure
hold on
for i=1:n_microcarriers
    plot(t, z_all(:, 2*n_microcarriers + 2*i));
end
xlabel('Time (s)')
ylabel('y direction velocity (m/s)')
title('y velocity vs time with buoyancy/gravity, drag, and TG vortex')
grid on
legend({'microcarrier 1','microcarrier 2', 'microcarrier 3'})

figure
hold on
for i=1:n_microcarriers
    plot(z_all(:, 2*i - 1), z_all(:, 2*i));
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
title('y vs x with buoyancy/gravity, drag, and TG vortex')
grid on
legend({'microcarrier 1','microcarrier 2', 'microcarrier 3'})


function dzdt = eom(t, z, mu, r_i, vfx, vfy, m, gx, gy, rho_f, rho_i, n_microcarriers)
    dzdt = zeros(4 * n_microcarriers, 1);  % Initialize the output vector for all microcarriers
    for i = 1:n_microcarriers
        % z vector takes the form [x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3]
        % Extract positions, velocities, and other parameters for each particle
        x = z(2*i - 1);
        y = z(2*i);
        vx = z(2*n_microcarriers + 2*i - 1);  % x-velocity of particle i
        vy = z(2*n_microcarriers + 2*i);      % y-velocity of particle i
        r = r_i(i);       % radius of particle i
        mi = m(i);        % mass of particle i
        rhoi = rho_i(i);

        % Position derivatives
        dxdt = vx;
        dydt = vy;
        
        % Drag and gravity force equations
        dvxdt = -(6 * pi * mu * r * (vx - vfx(x,y))) / mi + gx * (1 - rho_f / rhoi);
        dvydt = -(6 * pi * mu * r * (vy - vfy(x,y))) / mi + gy * (1 - rho_f / rhoi);
        
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