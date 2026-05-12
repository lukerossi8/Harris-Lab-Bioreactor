clearvars; close all; clc;
% Time range, constant values, and initial conditions
tspan = [0, 0.1];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s

% Defining 3 microcarriers
r_i = [90e-4, 80e-4, 70e-4]; % Microcarrier radii, in cm
m = [2.79e-6, 1.79e-6, 0.79e-6]; % Microcarrier mass, in g

g = [0, -9.81]; % m/s^2
gx = g(1);
gy = g(2);
p_f = 1; % g/ml
n_microcarriers = length(m);
p_i = [];
for i=1:n_microcarriers
    p_i = [p_i, m(i) / (4 / 3 * pi * r_i(i)^3)]; % Density in g/cm^3
end

% In the future, calculate the correct densities, and make this a vector

v0 = [0.02, -0.02; 0.05, 0.03; -0.06, 0.02]; % Initial velocity

vf = [0.01, 0.01]; % Fluid velocity
vfx = vf(1);
vfy = vf(2);

% Combine initial velocities of all microcarriers
v0_all = reshape(v0', [], 1);  % This converts the v0 matrix to a column vector

% Pass the number of microcarriers into the ODE solver
odeFun = @(t, v) eom(t, v, mu, r_i, vf, m, gx, gy, p_f, p_i, n_microcarriers);
[t, v_all] = ode45(odeFun, tspan, v0_all);

% Plotting the numerical solution in x direction
figure
hold on
for i=1:n_microcarriers
    plot(t, v_all(:, 2*i - 1));
end
xlabel('Time (s)')
ylabel('x direction velocity (m/s)')
title('Solution of buoyancy/gravity and drag force ODE, all microcarriers, x direction')
grid on
legend show

% Plotting the numerical solution in y direction
figure
hold on
for i=1:n_microcarriers
    plot(t, v_all(:, 2*i));
end
xlabel('Time (s)')
ylabel('y direction velocity (m/s)')
title('Solution of buoyancy/gravity and drag force ODE, all microcarriers, y direction')
grid on
legend show

function dvdt = eom(t, v, mu, r_i, vf, m, gx, gy, p_f, p_i, n_microcarriers)
    dvdt = zeros(2 * n_microcarriers, 1);  % Initialize the output vector for all microcarriers
    for i = 1:n_microcarriers
        % Extract velocities for each particle
        vx = v(2*i - 1);  % x-velocity of particle i
        vy = v(2*i);      % y-velocity of particle i
        r = r_i(i);       % radius of particle i
        mi = m(i);        % mass of particle i
        vfx = vf(1);      % Fluid velocity in x
        vfy = vf(2);      % Fluid velocity in y
        
        % Drag and gravity force equations
        dvxdt = -(6 * pi * mu * r * (vx - vfx)) / mi + gx * (1 - p_f / p_i(i));
        dvydt = -(6 * pi * mu * r * (vy - vfy)) / mi + gy * (1 - p_f / p_i(i));
        
        % Store the velocity derivatives
        dvdt(2*i - 1) = dvxdt;
        dvdt(2*i) = dvydt;
    end
end

% Clean up current code by making p_i a vector and g a vector
% Implement a Rankine vortex velocity field
% Start by figuring out how to plot this field on its own
% Then integrate it into this code
% Plot vx over time, vy over time, and microcarrier trajectory