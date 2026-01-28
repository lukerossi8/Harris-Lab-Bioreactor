function simpleODE2
clearvars; close all; clc;
tspan = [0,10];
y0 = 0;
[t, y] = ode45(@eom, tspan, y0);
plot(t, y);
xlabel('Time (s)')
ylabel('Solution value')
title('Solution of the ODE dydt = -2')
grid on
end

function dydt = eom(~, ~)
dydt = -2;
end