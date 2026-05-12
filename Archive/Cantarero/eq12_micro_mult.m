function eq12_cantarero_micro_mult
clearvars; close all; clc;
% Time range, constant values, and initial conditions
tspan = [800, 1000]; % Time range starts later to show exp growth clearly
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s
r_i = 90 * 10^-6; % Microcarrier radius, in m
v0 = [0, 1000, 2000]; % Initial velocity, in m/s
m = 2.79*10^-5; % Microcarrier mass, in kg
vf = 0;

% Define the anonymous function
odeFun = @(t, V) eom(V, mu, r_i, m, vf);

[t, V] = ode45(odeFun, tspan, v0);
plot(t, V)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Solution of the drag force ODE with multiple particles, in hydrostatic environment')
grid on

end

function dVdt = eom(V, mu, r_i, m, vf)
dVdt = (6 * pi * mu * r_i * (vf - V)) / m;
% v in this eqn is relative velocity, so this assumes no background flow
% (hydrostatic environment)
end

% Next steps:
% Combine eq12 and 13 in 1d, compare with analytical solution, then do this
% for multiple particles
%     Try changing dvdt for each particle, such as by changing m for the
%     different particles
% Expand the velocities to be vectors, vx and vy - solve each separately
% Define a uniform background flow (2d) - two variables, vfx and vfy
% Decompose the force equations into their x and y components


