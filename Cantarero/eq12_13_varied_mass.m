%function eq12_13_varied_mass
clear all; close all; clc;
% Time range, constant values, and initial conditions
tspan = [0, 200];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s
r_i = 90 * 10^-6; % Microcarrier radius, in m
v0 = [-30, 20, 100, 150]; % Initial particle velocity
vf = 0; % Fluid velocity
g = -9.81; % m/s^2
p_m = 1; % g/ml
p_i = 1.1; %g/ml, average cell density
m = [2.79*10^-5, 3.79*10^-5, 4.79*10^-5, 5.79*10^-5]; % Microcarrier masses, in kg

% Define the anonymous function
options = odeset('RelTol',1e-10, 'AbsTol', 1e-10);
for i = 1:length(v0)
    odeFun = @(t, v) eom(t, v, mu, r_i, vf, m(i), g, p_m, p_i);
    [t{i}, v{i}] = ode45(odeFun, tspan, v0(i), options);
end

% Define the analytical solution (solved by hand)
for i = 1:length(v0)
    m_current = m(i);
    syms t_an
    a = -(6 * pi * mu * r_i) / m_current;
    b = g * (1 - p_m / p_i);
    c = -a * vf + b;
    v_an(:, i) = -c / a + (v0(i) + c / a) * exp(a * t_an);
end

% Plotting the numerical solution
figure
hold on
for i = 1:length(v0)
    plot(t{i}, v{i});
end
axis([tspan, -40, 170]);
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Numerical solution of buoyancy/gravity and drag force ODE for multiple microcarriers')
grid on

% Plotting the analytical solution
fplot(v_an);
axis([tspan, -40, 170]);
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Analytical solution of buoyancy/gravity and drag force ODE for multiple microcarriers')
grid on

%end

% Defining the ODE
function dvdt = eom(t, v, mu, r_i, vf, m_current, g, p_m, p_i)
dvdt = -(6 * pi * mu * r_i * (v - vf)) / m_current + g * (1 - p_m / p_i);
end

% function eq12_13_varied_mass
% clearvars; close all; clc;
% % Time range, constant values, and initial conditions
% tspan = [0, 80];
% mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s
% r_i = 90 * 10^-6; % Microcarrier radius, in m
% v0 = [-30, 20, 100, 130]; % Initial particle velocity
% vf = 0; % Fluid velocity
% g = -9.81; % m/s^2
% p_m = 1; % g/ml
% p_i = 1.1; %g/ml, average cell density
% m = [2.79*10^-5, 3.79*10^-5, 4.79*10^-5, 5.79*10^-5]; % Microcarrier masses, in kg
% 
% % Define the anonymous function for the numerical solution
% for i = 1:length(v0)
%     m_current = m(i);
%     odeFun = @(t, v) eom(t, v, mu, r_i, vf, m_current, g, p_m, p_i);
%     [t, v(:, i)] = ode45(odeFun, tspan, v0(i));
% end
% 
% % Define the analytical solution (solved by hand)
% for i = 1:length(v0)
%     m_current = m(i);
%     syms t_an
%     a = (m_current * g * (1 - p_m / p_i) / (6 * pi * mu * r_i));
%     b = -6 * pi * mu * r_i / m_current;
%     C = v0 - vf - a;
%     v_an(:, 1) = vf + a + C * exp(b * t_an);
%     %v_an(:, 1) = -c / a + (v0 + c / a) * exp(a * t_an);
% end
% 
% % Plotting the numerical solution
% figure
% hold on
% plot(t, v);
% axis([tspan, -40, 170]);
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% title('Numerical solution of buoyancy/gravity and drag force ODE for multiple microcarriers')
% grid on
% 
% % Plotting the analytical solution
% fplot(v_an);
% axis([tspan, -40, 170]);
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% title('Analytical solution of buoyancy/gravity and drag force ODE for multiple microcarriers')
% grid on
% 
% end
% 
% % Defining the ODE
% function dvdt = eom(t, v, mu, r_i, vf, m_current, g, p_m, p_i)
% dvdt = (6 * pi * mu * r_i * (vf - v)) / m_current + g * (1 - p_m / p_i);
% end