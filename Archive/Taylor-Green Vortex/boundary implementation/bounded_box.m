clearvars; close all; clc;
% Time range, constant values, and initial conditions
tstart = 0;
tstop = 400;
tspan = [tstart; tstop];
x0 = 0;
y0 = 1;
vx0 = 0.1;
vy0 = -1;
z0 = [x0; y0; vx0; vy0];

% Boundary parameters
floor = -5;
ceil = 5;
l_wall = -15;
r_wall = 5;

options = odeset('Events', @(t, z) boundary(t, z, floor, ceil, l_wall, r_wall), 'Refine', 100);
odeFun = @(t, z) eom(t, z);
[t, z_all] = ode45(odeFun, tspan, z0, options);
tplot = t;
xplot = z_all(:, 1);
yplot = z_all(:, 2);

while t(end) < tstop
    z0 = z_all(end, :);
    % If it hits either wall, replace and reverse y velocity
    if z0(1) <= l_wall
        z0(1) = l_wall;
        z0(3) = -z0(3);
    end
    if z0(1) >= r_wall
        z0(1) = r_wall;
        z0(3) = -z0(3);
    end
    % If it hits the floor or ceiling, replace and reverse y velocity
    if z0(2) <= floor
        z0(2) = floor;
        z0(4) = -z0(4);
    end
    if z0(2) >= ceil
        z0(2) = ceil;
        z0(4) = -z0(4);
    end
    tspan = [t(end), tstop];
    [t, z_all] = ode45(odeFun, tspan, z0, options);
    tplot = [tplot; t];
    xplot = [xplot; z_all(:, 1)];
    yplot = [yplot; z_all(:, 2)];
end


plot(xplot, yplot)

function [value, isTerminal, direction] = boundary(t, z, floor, ceil, l_wall, r_wall)
value(1) = (z(1) - l_wall);
value(2) = (z(1) - r_wall);
value(3) = (z(2) - floor);
value(4) = (z(2) - ceil);
isTerminal(value ~= 0) = 1;
direction = 0;
end

function dzdt = eom(t, z)
dzdt = [z(3); z(4); 0; 0];
end