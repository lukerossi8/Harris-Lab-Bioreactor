%% Progress 9/24
% Implemented way to initialize arbitrary number of paired (bonded) agents,
% tested varying number of bonds while keeping number of agents constant

% Updated grid square tracking system, was using an empty cell array, but
% now using a dictionary (more streamlined, space efficient)

tic
% Cell contact assumptions: 2 agents only, both w/same (constant) radius
clearvars; clc;
%% Time range, constant values, and initial conditions
tstart = 0;
tstop = 1;
tspan = [tstart, tstop];
mu = 10^-3; % Water dynamic viscosity at 20°C, in Pa*s

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

global bonded_pairs;
bonded_pairs = [0, 0]; % List of indices representing pairs of bonded agents (initialize as [0,0])

% Define a discrete background velocity field (x and y cpts)
vfx_disc = (A_TG * pi / L_TG) * cos(pi * X / L_TG) .* sin(pi * Y / L_TG);
vfy_disc = -(A_TG * pi / L_TG) * sin(pi * X / L_TG) .* cos(pi * Y / L_TG);

%% Defining agents
n_agents = 5; % number of agents modeled
r_i = 90e-6 + zeros(n_agents, 1); % Agent radii, in m
m = 2.79e-9 + zeros(n_agents, 1); % Agent mass, in kg

g = [0, -9.81]; % gravitational acceleration, m/s^2
gx = g(1);
gy = g(2);
rho_f = 1000; % kg/m^3
rho_i = m ./ (4 / 3 * pi * r_i.^3); % Density in g/cm^3

% Random position setting
x0 = l_wall + 0.05*(r_wall - l_wall) + 0.9*(r_wall - l_wall).*rand(n_agents, 1);
y0 = floor + 0.05*(ceil - floor) + 0.9*(ceil - floor).*rand(n_agents, 1);

% Random bonded position setting
% n_bonds = 100;
% n_unbonded = n_agents - 2*n_bonds;
% 
% x01 = l_wall + 0.05*(r_wall - l_wall) + 0.9*(r_wall - l_wall).*rand(n_bonds, 1);
% x02 = x01 + 14e-5;
% x0_unbonded = l_wall + 0.05*(r_wall - l_wall) + 0.9*(r_wall - l_wall).*rand(n_unbonded, 1);
% x0 = [x01;x02;x0_unbonded];
% 
% y01 = floor + 0.05*(ceil - floor) + 0.9*(ceil - floor).*rand(n_bonds, 1);
% y0_unbonded = floor + 0.05*(ceil - floor) + 0.9*(ceil - floor).*rand(n_unbonded, 1);
% y0 = [y01;y01;y0_unbonded];

% Manual position setting
%x0 = [0.09; 0.0901; 0.09015; 0.088; 0.0881; 0.09]; % Initial x position, m
%y0 = [0; 0; -0.0001; -0.001; -0.001; -0.002]; % Initial y position, m

s0 = [x0; y0]; % Initial positions, m

vx0 = zeros(n_agents, 1); % Initial x velocity, m/s
vy0 = zeros(n_agents, 1); % Initial y velocity, m/s
v0 = [vx0; vy0]; % Initial velocity, m/s

% Properly format initial conditions for ode45
z0 = [s0; v0];

%% ode45 set up and call
RelTol = 1e-6;
AbsTol = 1e-9;
options = odeset('Events', @(t, z) boundary(t, z, floor, ceil, l_wall, r_wall, n_agents), ...
    'Refine', 100, ...
    'RelTol', RelTol, ...
    'AbsTol', AbsTol);

%options = odeset('Events', @(t, z) boundary(t, z, floor, ceil, l_wall, r_wall, n_agents), 'Refine', 100);

% Would like to better understand the difference between these two options

odeFun = @(t, z) eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_agents);
[t, z_all] = ode45(odeFun, tspan, z0, options);

t_plot = t; % storing vector with all values of t throughout the ode45 call
z_all_plot = z_all; % storing vector with all values of z

%% Model boundary interactions
% Define boundary logic in the while loop
while t(end) < tstop
    % z vector takes the form [x1; x2; x3;...; y1; y2; y3;...; vx1; vx2; vx3;...; vy1; vy2; vy3;...]
    % Separate into vectors x0, y0, vx0, vy0
    % Identify indices of positions that have crossed the boundaries
    % Flip the corresponding velocities and adjust those positions to the
    % boundary
    z1 = z_all(end, :);
    x1 = z1(1:n_agents)';
    y1 = z1(n_agents+1:2*n_agents)';
    vx1 = z1(2*n_agents+1:3*n_agents)';
    vy1 = z1(3*n_agents+1:end)';

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
function [value, isTerminal, direction] = boundary(t, z, floor, ceil, l_wall, r_wall, n_agents)
% z vector takes the form [x1; x2; x3;...; y1; y2; y3;...; vx1; vx2; vx3;...; vy1; vy2; vy3;...]

x = z(1:n_agents);
y = z(n_agents+1:2*n_agents);

value = [(x - l_wall); (x - r_wall); (y - floor); (y - ceil)];
isTerminal(value ~= 0) = 1;
direction = 0;
end

%% Numerical solution
function dzdt = eom(t, z, mu, r_i, vfx_disc, vfy_disc, X, Y, m, gx, gy, rho_f, rho_i, n_agents)
dzdt = zeros(4 * n_agents, 1);  % Initialize the output vector for all agents
% z vector takes the form [x1; x2; x3;...; y1; y2; y3;...; vx1; vx2; vx3;...; vy1; vy2; vy3;...]

% Extracting positions, velocities, and other parameters for each agent
x = z(1:n_agents);  % x positions
y = z(n_agents+1:2*n_agents); % y positions
vx = z(2*n_agents+1:3*n_agents);  % x velocities
vy = z(3*n_agents+1:end);      % y velocities

% Interpolating to find background flow at the agent's position
vfx = interp2(X, Y, vfx_disc, x, y, 'linear', 0);
vfy = interp2(X, Y, vfy_disc, x, y, 'linear', 0);

% Position derivatives
dxdt = vx;
dydt = vy;

% Gravity/buoyancy and drag forces in x and y dir
Fgx = m .* gx .* (1 - rho_f ./ rho_i); % N
Fdx = -(6 * pi * mu * r_i .* (vx - vfx)); % N

Fgy = m .* gy .* (1 - rho_f ./ rho_i); % N
Fdy = -(6 * pi * mu * r_i .* (vy - vfy)); % N

% Bonding thresholds
delta_c = 2*r_i(1); % Bond formation threshold, m
% ^Threshold = sum of radii of the two agents. This needs to become more
% complex once we have agents with changing or different radii. Really, it
% can just be brought inside the for loop and consider the specific radii
% of the two agents at that moment in time
delta_d = 1.4*delta_c; % Bond breaking threshold, m

global bonded_pairs;

% Defining a new grid for bond calculations

% Calculating which grid square each agent is in, storing in dictionary
grid_pts_x = round(x, 3); % To the nearest 0.001
grid_pts_y = round(y, 3);
grid_pts = [grid_pts_x, grid_pts_y];

agent_grid_pts = dictionary();

for i=1:size(grid_pts, 1)
    current_pt = grid_pts(i, :);
    if not(numEntries(agent_grid_pts) == 0)
        if not(isKey(agent_grid_pts, current_pt))
            agent_grid_pts = insert(agent_grid_pts, current_pt, i);
        else
            point = agent_grid_pts(current_pt);
            point(end+1, :) = i;
            agent_grid_pts(current_pt) = point;
        end
    else
        agent_grid_pts = insert(agent_grid_pts, current_pt, i);
    end
end

% Iterating through each occupied grid square to find bonds
occupied_pts = keys(agent_grid_pts);
for i=1:length(occupied_pts) % looping through the occupied grid points
    current_pt = occupied_pts(i);
    if length(agent_grid_pts(current_pt)) > 1
        pairs = nchoosek(agent_grid_pts(current_pt), 2);
        for k = 1:size(pairs, 1) % for each pair of agents identified
                current_pair = pairs(k, :);
                agent_i = current_pair(:, 1); % Lower indexed agent in the pair
                agent_j = current_pair(:, 2); % Higher indexed agent in the pair

                % Measure distances and angle between the two agents
                dist_x = z(agent_j)-z(agent_i);
                dist_y = z(n_agents+agent_j)-z(n_agents+agent_i);
                dist = sqrt((dist_x)^2 + (dist_y)^2); % Distance between the two agents, m

                % Check if the agents are already bonded
                is_bonded = ismember(current_pair, bonded_pairs, "rows");
                % If the agents are already bonded, check if the bond has broken
                if is_bonded
                    if dist > delta_d % If the bond has broken in the previous timestep
                        bonded_pairs(bond, :) = []; % Remove the pair (row) from bonded_pairs
                    end
                else % the agents are not already bonded
                    if dist <= delta_c % the agents have bonded in the previous timestep
                        if bonded_pairs == [0,0]
                            bonded_pairs = current_pair;
                        else
                            % Add the pair to bonded_pairs
                            bonded_pairs(end+1, :) = current_pair;
                        end
                    end
                end
            end
    end
end

% Loop through each pair in bonded_pairs, calculating force b/w them
bond_forces_x = zeros(n_agents, 1);
bond_forces_y = zeros(n_agents, 1);

if bonded_pairs == [0,0] % If there are no bonds
    % Then no updates are needed to the bond forces vectors, they can stay
    % as 0 vectors of size n_agents
else % There are bonds, so we need to solve for the bond forces
    for i=1:size(bonded_pairs, 1)
        bond = bonded_pairs(i, :);
        agent_i = bond(:, 1); % Lower indexed agent in the bond
        agent_j = bond(:, 2); % Higher indexed agent in the bond

        dist_x = z(agent_j)-z(agent_i);
        dist_y = z(n_agents+agent_j)-z(n_agents+agent_i);
        dist = sqrt((dist_x)^2 + (dist_y)^2); % Distance between the two agents, m
        theta_1 = atan2(-dist_y, -dist_x); % Angle of agent 1 relative to agent 2, radians
        theta_2 = atan2(dist_y, dist_x); % Angle of agent 2 relative to agent 1, radians

        % Calculate the bond force on each agent
        K_ij = 1e-3; % Spring constant, N/m (SA)
        s_ij = 0.2; % Bond sensitivity (SA)
        delta_ij = delta_c - dist; % Degree of separation or overlap, m
        F_ij = K_ij * delta_ij * tanh(s_ij/abs(delta_ij));
        F_ix = F_ij .* cos(theta_1); % Bond force on i in x dir
        F_iy = F_ij .* sin(theta_1); % Bond force on i in y dir
        F_jx = F_ij .* cos(theta_2); % Bond force on j in x dir
        F_jy = F_ij .* sin(theta_2); % Bond force on j in y dir

        % Updating bond_forces vector
        bond_forces_x(agent_i, 1) = bond_forces_x(agent_i, 1) + F_ix;
        bond_forces_y(agent_i, 1) = bond_forces_y(agent_i, 1) + F_iy;
        bond_forces_x(agent_j, 1) = bond_forces_x(agent_j, 1) + F_jx;
        bond_forces_y(agent_j, 1) = bond_forces_y(agent_j, 1) + F_jy;
    end
end

% Velocity derivatives represented by sum of forces div by masses
dvxdt = (Fgx + Fdx + bond_forces_x) ./ m;
dvydt = (Fgy + Fdy + bond_forces_y) ./ m;
% A way to turn forces off for testing purposes
% dvxdt = 0;
% dvydt = 0;

% Store the position and velocity derivatives
dzdt(1:n_agents) = dxdt;
dzdt(n_agents+1:2*n_agents) = dydt;
dzdt(2*n_agents+1:3*n_agents) = dvxdt;
dzdt(3*n_agents+1:end) = dvydt;
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

x_plot = z_all_plot(:,1:n_agents);
y_plot = z_all_plot(:,n_agents+1:2*n_agents);
vx_plot = z_all_plot(:,2*n_agents+1:3*n_agents);
vy_plot = z_all_plot(:,3*n_agents+1:4*n_agents);
v_plot = sqrt(vx_plot.^2+vy_plot.^2);

for i=1:n_agents
    plot(x_plot(:,i), y_plot(:,i), 'color', col_list(mod(i, 4)+1));
    plot(s0(i), s0(i+n_agents), "or", "LineWidth", 1.25)
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
title('Cell-cell bonding in TG Vortex')

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
% for i=1:n_agents
%     hold on
%     plot(t_plot, v_plot(:, i))
%     grid on
%     xlabel('time (s)')
%     ylabel('velocity (m/s)')
% end

%% Plotting distance between agents over time
% dist_plot_4_5 = sqrt((x_plot(:, 5) - x_plot(:, 4)).^2 + (y_plot(:, 5) - y_plot(:, 4)).^2);
% dist_plot_2_3 = sqrt((x_plot(:, 3) - x_plot(:, 2)).^2 + (y_plot(:, 3) - y_plot(:, 2)).^2);

% figure
% hold on
% plot(t_plot, dist_plot_2_3, 'LineWidth',2)
% yline(2*r_i(1))
% grid on
% xlabel('time (s)')
% ylabel('particle distance (m)')
% title('distance between agents 2 and 3 — with bonding')
% 
% figure
% hold on
% plot(t_plot, dist_plot_4_5, 'LineWidth',2)
% yline(2*r_i(1))
% grid on
% xlabel('time (s)')
% ylabel('particle distance (m)')
% title('distance between agents 4 and 5 — with bonding')

toc