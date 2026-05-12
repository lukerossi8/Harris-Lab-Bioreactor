clearvars; close all; clc;
% Time range, constant values, and initial conditions
tspan = [0, 100];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s
r_i = 90 * 10^-6; % Microcarrier radius, in m

v0 = [100, 100]; % Initial velocity
v0x = v0(1);
v0y = v0(2);

vf = [50, 0]; % Fluid velocity
vfx = vf(1);
vfy = vf(2);

g = [0, -9.81]; % m/s^2
gx = g(1);
gy = g(2);
p_f = 1; % g/ml
p_i = 1.1; %g/ml, average cell density
m = 2.79*10^-5; % Microcarrier mass, in kg

% Define the anonymous function
odeFun = @(t, v) eom(t, v, mu, r_i, vf, m, gx, gy, p_f, p_i);
[t, v] = ode45(odeFun, tspan, v0);

% Define the analytical solution in x direction (solved by hand)
syms t_an
a = -(6 * pi * mu * r_i) / m;
vx_an = vfx + (v0x - vfx) * exp(a * t_an);
vx_an_func = matlabFunction(vx_an, 'Vars', {t_an});

% Define analytical solution in y direction
b = gy * (1 - p_f / p_i);
c = -a * vfy + b;
vy_an = -c / a + (v0y + c / a) * exp(a * t_an);
vy_an_func = matlabFunction(vy_an, 'Vars', {t_an});

% Plotting the numerical solution in x direction
figure
hold on
plot(t, v(:, 1));
xlabel('Time (s)')
ylabel('x direction velocity (m/s)')
title('Numerical solution of buoyancy/gravity and drag force ODE, one microcarrier, x direction')
grid on

% Plotting the analytical solution in x direction
plot(t, vx_an_func(t));
xlabel('Time (s)')
ylabel('x direction velocity (m/s)')
title('Analytical and numerical solution of buoyancy/gravity and drag force ODE, one microcarrier, x direction')
grid on

% Plotting the numerical solution in y direction
figure
hold on
plot(t, v(:, 2));
xlabel('Time (s)')
ylabel('x direction velocity (m/s)')
title('Numerical solution of buoyancy/gravity and drag force ODE for one microcarrier, in y direction')
grid on

% Plotting the analytical solution in y direction
plot(t, vy_an_func(t));
xlabel('Time (s)')
ylabel('y direction velocity (m/s)')
title('Analytical and numerical solution of buoyancy/gravity and drag force ODE, one microcarrier, y direction')
grid on
legend show

function dvdt = eom(t, v, mu, r_i, vf, m, gx, gy, p_f, p_i)
vx = v(1);
vy = v(2);
vfx = vf(1);
vfy = vf(2);
dvxdt = -(6 * pi * mu * r_i * (vx - vfx)) / m + gx * (1 - p_f / p_i);
dvydt = -(6 * pi * mu * r_i * (vy - vfy)) / m + gy * (1 - p_f / p_i);
dvdt = [dvxdt; dvydt];
end