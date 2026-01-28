clearvars; close all; clc;
% Time range, constant values, and initial conditions
tstart = 0;
tstop = 20;
tspan = [tstart; tstop];
x0 = 0;
y0 = 1;
vx0 = 0.1;
vy0 = -0.2;
z0 = [x0; y0; vx0; vy0];

options = odeset('Events', @boundary, 'Refine', 100);
odeFun = @(t, z) eom(t, z);
[t, z_all] = ode45(odeFun, tspan, z0, options);
xplot = z_all(:, 1);
yplot = z_all(:, 2);

for i=1:9
    z0 = z_all(end, :);
    z0(2) = 0;
    z0(4) = -z0(4) * 0.9;
    [t, z_all] = ode45(odeFun, tspan, z0, options);
    xplot = [xplot; z_all(:, 1)];
    yplot = [yplot; z_all(:, 2)];
end

plot(xplot, yplot)

function [value, isTerminal, direction] = boundary(t, z)
value = z(2);
isTerminal = 1;
direction = 0;
end

function dzdt = eom(t, z)
dzdt = [z(3); z(4); 0; -9.8];
end