function eq12_cantarero_micro
clearvars; close all; clc;
% Time range, constant values, and initial conditions
tspan = [0, 1000000];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s
r_i = 90 * 10^-6; % Microcarrier radius, in m
v0 = 1000; % Initial velocity
vf = 0;

% Define the anonymous function
odeFun = @(t, v) eom(v, mu, r_i, vf);

[t, v] = ode45(odeFun, tspan, v0);
plot(t, v)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Solution of the drag force ODE, in hydrostatic environment')
grid on

end

function dvdt = eom(v, mu, r_i, vf)
dvdt = 6 * pi * mu * r_i * (vf - v); 
% v in this eqn is relative velocity, so this assumes no background flow
% (hydrostatic environment)
end