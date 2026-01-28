clear all; close all; clc; format compact;

g = -9.81; % m/s^2
pm = 1; % g/ml
pi = 1.1; %g/ml, average cell density

tspan = [0,100];
v0 = 0;
eom_with_params = @(t, v) eom(t, v, g, pm, pi);
[t, v] = ode45(eom_with_params, tspan, v0);


plot(t, v)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Solution of the buoyancy/gravity force ODE')
grid on

function dvdt = eom(t, v, g, pm, pi)
dvdt = g*(1 - pm/pi);
end