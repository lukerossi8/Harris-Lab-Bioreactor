function simpleODE
clearvars; clc; close all;
tspan = [0,5]; % Time range
y0 = 0; % Initial condition
[t, y] = ode45(@eom, tspan, y0);
plot(t, y);
xlabel('Time (s)')
ylabel('Solution value')
title('Solution of the ODE dydt = -2*y + 1')
grid on
end

function dydt = eom(t, y)
dydt = -2*y +1;
end