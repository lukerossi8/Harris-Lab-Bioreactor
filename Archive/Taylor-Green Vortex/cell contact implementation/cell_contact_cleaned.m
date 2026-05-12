% Cell contact assumptions: 2 agents only, both w/same (constant) radius
clearvars; clc;
%% Time range, constant values, and initial conditions
tstart = 0;
tstop = 1;
tspan = [tstart, tstop];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s

%% Defining 2 microcarriers
r_i = [90e-6; 90e-6]; % Microcarrier radii, in m
m = [2.79e-9; 2.79e-9]; % Microcarrier mass, in kg

g = [0, -9.81]; % gravitational acceleration, m/s^2
gx = g(1);
gy = g(2);
rho_f = 1000; % kg/m^3
n_microcarriers = length(m); % storing the # of microcarriers modeled
rho_i = m ./ (4 / 3 * pi * r_i.^3); % Density in g/cm^3

x0 = [0.09; 0.02]; % Initial x position, m
y0 = [-0.01; -0.02]; % Initial y position, m
s0 = [x0; y0]; % Initial positions, m

vx0 = [0; 0]; % Initial x velocity, m/s
vy0 = [0; 0]; % Initial y velocity, m/s
v0 = [vx0; vy0]; % Initial velocity, m/s

% Properly format initial conditions for ode45
z0 = [s0; v0];

% Initial velocity values inspired by range seen in CMMC paper, Fig. 5A

%% Define the vortex, bounds, bonded pair list, & background velocity field
A_TG = 0.005; % Amplitude of the TG vortex, in m^2/s
L_TG = 0.1; % Length of the TG vortex, in m

floor = -L_TG / 2; % low y boundary
ceil = L_TG / 2; % high y boundary
l_wall = -L_TG / 2; % left x boundary coordinate
r_wall = 3 * L_TG / 2; % right x boundary

n_step = 200;
grid_sl = L_TG / n_step; % Side length of each square in the grid
[X,Y] = meshgrid(l_wall : grid_sl : r_wall, floor : grid_sl : ceil);
psi_TG = -A_TG*cos(pi * X / L_TG) .* cos(pi * Y / L_TG);

bonded_pairs = [0,0]; % List of indices representing pairs of bonded agents

% Define a discrete background velocity field (x and y cpts)
vfx_disc = (A_TG * pi / L_TG) * cos(pi * X / L_TG) .* sin(pi * Y / L_TG);
vfy_disc = -(A_TG * pi / L_TG) * sin(pi * X / L_TG) .* cos(pi * Y / L_TG);

%% ode45 set up and call
RelTol = 1e-6;
AbsTol = 1e-9;
options = odeset('Events', @(t, z) boundary(t, z, floor, ceil, l_wall, r_wall, n_microcarriers), ...
    'Refine', 100, ...
    'RelTol', RelTol, ...
    'AbsTol', AbsTol);

%options = odeset('Events', @(t, z) boundary(t, z, floor, ceil, l_wall, r_wall, n_microcarriers), 'Refine', 100);

% Would like to better understand the difference between these two options

odeFun = @(t, z) eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_microcarriers, bonded_pairs);
[t, z_all] = ode45(odeFun, tspan, z0, options);

t_plot = t; % storing vector with all values of t throughout the ode45 call
z_all_plot = z_all; % storing vector with all values of z

%% Model boundary interactions
% Define boundary logic in the while loop
while t(end) < tstop
    % z vector takes the form [x1, x2, x3,..., y1, y2, y3,..., vx1, vx2, vx3,..., vy1, vy2, vy3,...]
    % Separate into vectors x0, y0, vx0, vy0
    % Identify indices of positions that have crossed the boundaries
    % Flip the corresponding velocities and adjust those positions to the
    % boundary
    z1 = z_all(end, :);
    x1 = z1(1:n_microcarriers)';
    y1 = z1(n_microcarriers+1:2*n_microcarriers)';
    vx1 = z1(2*n_microcarriers+1:3*n_microcarriers)';
    vy1 = z1(3*n_microcarriers+1:end)';

    left = find(x1 <= l_wall);
    right = find(x1 >= r_wall);
    low = find(y1 <= floor);
    high = find(y1 >= ceil);

    % If any agent hits either wall, replace & reverse x velocity
    if ~isempty(left)
        x1(left) = l_wall;
        vx1(left) = -vx1(left);
    end

    if ~isempty(right)
        x1(right) = r_wall;
        vx1(right) = -vx1(right);
    end

    % If any agent hits the floor or ceiling, replace & reverse y velocity
    if ~isempty(low)
        y1(low) = floor;
        vy1(low) = -vy1(low);
    end

    if ~isempty(high)
        y1(high) = ceil;
        vy1(high) = -vy1(high);
    end

    z1 = [x1; y1; vx1; vy1]; % New initial state

    % Reset time range to what's remaining
    tspan = [t(end), tstop];

    % Call ode45 with new adjusted conditions
    [t, z_all] = ode45(odeFun, tspan, z1, options);

    % Accumulate values of t, x, and y across all the calls
    t_plot = [t_plot; t];
    z_all_plot = [z_all_plot; z_all];
end

% Boundary interaction event function
function [value, isTerminal, direction] = boundary(t, z, floor, ceil, l_wall, r_wall, n_microcarriers)
% z vector takes the form [x1, x2, x3, y1, y2, y3, vx1, vx2, vx3, vy1, vy2, vy3]

x = z(1:n_microcarriers);
y = z(n_microcarriers+1:2*n_microcarriers);

value = [(x - l_wall); (x - r_wall); (y - floor); (y - ceil)];
isTerminal(value ~= 0) = 1;
direction = 0;
end

%% Numerical solution
function dzdt = eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_microcarriers, bonded_pairs)
dzdt = zeros(4 * n_microcarriers, 1);  % Initialize the output vector for all microcarriers
% z vector takes the form [x1, x2, x3, y1, y2, y3, vx1, vx2, vx3, vy1, vy2, vy3]

% Extracting positions, velocities, and other parameters for each agent
x = z(1:n_microcarriers);  % x positions
y = z(n_microcarriers+1:2*n_microcarriers); % y positions
vx = z(2*n_microcarriers+1:3*n_microcarriers);  % x velocities
vy = z(3*n_microcarriers+1:end);      % y velocities

% Interpolating to find background flow at the agent's position
vfx = interp2(X, Y, vfx_disc, x, y, 'linear', 0);
vfy = interp2(X, Y, vfy_disc, x, y, 'linear', 0);

% Position derivatives
dxdt = vx;
dydt = vy;

% Gravity/buoyancy and drag forces in x and y dir
Fgx = m .* gx .* (1 - rho_f ./ rho_i); % N
Fdx = -(6 * pi * mu * r_i .* (vx - vfx)); % N
% Fgx = 0;
% Fdx = 0;

Fgy = m .* gy .* (1 - rho_f ./ rho_i); % N
Fdy = -(6 * pi * mu * r_i .* (vy - vfy)); % N
% Fgy = 0;
% Fdy = 0;

% Calculating agent distances
delta_c = 2*r_i(1); % Bond formation threshold, m
% ^Highly simplified, principle is threshold = sum of radii of the two agents
delta_d = 1.4*delta_c; % Bond breaking threshold, m
dist_x = z(2)-z(1);
dist_y = z(4)-z(3);
dist = sqrt((dist_x)^2 + (dist_y)^2); % Distance between the two agents, m
delta_ij = delta_c - dist; % Degree of separation or overlap, m
theta_1 = atan2(-dist_y, -dist_x); % Angle of agent 1 relative to agent 2, radians
theta_2 = atan2(dist_y, dist_x); % Angle of agent 2 relative to agent 1, radians
theta = [theta_1; theta_2];

bond = ismember(bonded_pairs, [1, 2], "rows");
if any(bond) % If the agents are bonded
    if dist > delta_d % If the bond is broken
        bonded_pairs(bond, :) = []; % Remove the pair (row) from bonded_pairs
    end
else % the agents are not already bonded
    if dist <= delta_c % the agents have bonded in the previous timestep
        % Add the pair to bonded_pairs
        bonded_pairs(end+1, :) = [min(1,2), max(1,2)];
    end
end

% Solve with the bond force (plus drag and buoyancy)
delta_ij = delta_c - dist; % Degree of separation or overlap, m
K_ij = 1e-3; % Spring constant, N/m (SA)
s_ij = 0.2; % Bond sensitivity (SA)
F_ij = K_ij * delta_ij * tanh(s_ij/abs(delta_ij));
F_ijx = F_ij .* cos(theta);
F_ijy = F_ij .* sin(theta);

dvxdt = (Fgx + Fdx + F_ijx) ./ m; % Acceleration in x dir, m/s^2
dvydt = (Fgy + Fdy + F_ijy) ./ m; % Acceleration in y dir, m/s^2

% Store the position and velocity derivatives
dzdt(1:n_microcarriers) = dxdt;
dzdt(n_microcarriers+1:2*n_microcarriers) = dydt;
dzdt(2*n_microcarriers+1:3*n_microcarriers) = dvxdt;
dzdt(3*n_microcarriers+1:end) = dvydt;
end

%% Plot agent trajectories w/ background velocity & vorticity in TG vortex
figure
hold on
[curlz,cav] = curl(X, Y, vfx_disc, vfy_disc);
c = pcolor(X, Y, curlz);
c.FaceColor = 'interp';
c.EdgeColor = 'none';

% Define three colors for the gradient
color_pos = [0.9290 0.6940 0.1250];   % Yellow for positive vorticity
color_neg = [0.00392156862745098 0.44313725490196076 0.7372549019607844];   % Blue for negative vorticity
color_zero = [0.8480392156862745 0.8911764705882353 0.9225490196078431];  % White for zero vorticity

% Create a custom colormap with a smooth transition between these three colors
numColors = 256;  % Number of colors in the colormap (higher resolution)
cmap_neg_to_zero = [linspace(color_neg(1), color_zero(1), numColors/2)', ...
    linspace(color_neg(2), color_zero(2), numColors/2)', ...
    linspace(color_neg(3), color_zero(3), numColors/2)'];

cmap_zero_to_pos = [linspace(color_zero(1), color_pos(1), numColors/2)', ...
    linspace(color_zero(2), color_pos(2), numColors/2)', ...
    linspace(color_zero(3), color_pos(3), numColors/2)'];

% Combine the two halves to form the full colormap
cmap = [cmap_neg_to_zero; cmap_zero_to_pos];

colormap(cmap);
cb = colorbar;
yl = ylabel(cb,'Vorticity','FontSize',10,'Rotation',270);

n_step_q = 20;
[X_q,Y_q] = meshgrid(l_wall : L_TG / n_step_q : r_wall, floor : L_TG / n_step_q : ceil);
vfx_disc_q = (A_TG * pi / L_TG) * cos(pi * X_q / L_TG) .* sin(pi * Y_q / L_TG);
vfy_disc_q = -(A_TG * pi / L_TG) * sin(pi * X_q / L_TG) .* cos(pi * Y_q / L_TG);
q = quiver(X_q, Y_q, vfx_disc_q, vfy_disc_q, 'k');
col_list = ["black", "red", "blue", "magenta"];

x_plot = z_all_plot(:,1:n_microcarriers);
y_plot = z_all_plot(:,n_microcarriers+1:2*n_microcarriers);
vx_plot = z_all_plot(:,2*n_microcarriers+1:3*n_microcarriers);
vy_plot = z_all_plot(:,3*n_microcarriers+1:4*n_microcarriers);
v_plot = sqrt(vx_plot.^2+vy_plot.^2);
dist_plot = sqrt((x_plot(:, 2) - x_plot(:, 1)).^2 + (y_plot(:, 2) - y_plot(:, 1)).^2);

for i=1:n_microcarriers
    plot(x_plot(:,i), y_plot(:,i), 'color', col_list(i));
    plot(s0(i), s0(i+n_microcarriers), "or", "LineWidth", 1.25)
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
title('Cell-cell bonding proof-of-concept in TG Vortex')

%% Plotting time step size over time of simulation

deltat = [0];
for i=2:length(t_plot)
    deltat(end+1) = t_plot(i) - t_plot(i-1);
end

figure
plot(t_plot, deltat)
grid on
xlabel('time (s)')
ylabel('time step size (s)')
title('Time step size throughout the simulation — with bonding')

%% Plotting velocity over time of simulation
% figure
% for i=1:n_microcarriers
%     hold on
%     plot(t_plot, v_plot(:, i))
%     grid on
%     xlabel('time (s)')
%     ylabel('velocity (m/s)')
% end

%% Plotting distance between agents over time
figure
hold on
plot(t_plot, dist_plot)
yline(2*r_i(1))
grid on
xlabel('time (s)')
ylabel('particle distance (m)')
title('distance between agents 1 and 2 — with bonding')