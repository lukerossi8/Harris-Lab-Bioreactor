function eq12_13
clearvars; close all; clc;
% Time range, constant values, and initial conditions
tspan = [0, 100];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s
r_i = 90 * 10^-6; % Microcarrier radius, in m
v0 = 100; % Initial velocity
vf = 0; % Fluid velocity
g = -9.81; % m/s^2
p_m = 1; % g/ml
p_i = 1.1; %g/ml, average cell density
m = 2.79*10^-5; % Microcarrier mass, in kg

% Define the anonymous function
odeFun = @(t, v) eom(t, v, mu, r_i, vf, m, g, p_m, p_i);

% Define the analytical solution (solved by hand)
syms t_an
a = -(6 * pi * mu * r_i) / m;
b = g * (1 - p_m / p_i);
c = -a * vf + b;
v_an = (- c / a) + (v0 + c / a) * exp(a * t_an);

% Alternate way of doing it
% a = (m * g * (1 - p_m / p_i) / (6 * pi * mu * r_i));
% b = -6 * pi * mu * r_i / m;
% C = v0 - vf - a;
% v_an = vf + a + C * exp(b * t_an);

% Plotting the numerical solution
hold on
[t, v] = ode45(odeFun, tspan, v0);
plot(t, v);
axis([0, 100, -20, 110])
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Numerical solution of buoyancy/gravity and drag force ODE for one microcarrier')
grid on

% Plotting the analytical solution
figure
fplot(v_an);
axis([0, 100, -20, 110])
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Analytical solution of buoyancy/gravity and drag force ODE for one microcarrier')
grid on

end

function dvdt = eom(t, v, mu, r_i, vf, m, g, p_m, p_i)
dvdt = -(6 * pi * mu * r_i * (v - vf)) / m + g * (1 - p_m / p_i);
end