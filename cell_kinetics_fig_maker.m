function cell_kinetics_fig_maker(tstop, n_agents, position)
% cell_kinetics_fig_maker creates visualizations of agent behavior in a
% bioreactor environment, by calling the cell_kinetics function
% Syntax:
%   cell_kinetics_fig_maker(tstop, n_agents, position)
% Inputs:
%   tstop       Double representing the length of the simulation (seconds)
%   n_agents    Double representing the number of agents in the simulation
%   position    Either a string in ["random", "random bonded"], or a matrix 
%               of size [n_agents, 2] of doubles (column 1 = x vals, 
%               column 2 = y vals)
% Outputs
%   Figure plotting the trajectories of the cells on a plain background
%   Figure plotting the trajectories of the cells with the background 
%   vorticity
%   Figure plotting the trajectories of the cells with the background 
%   volume fraction 
%   Histogram of the distribution of cell average velocities over the
%   simulation

arguments
    tstop (1,1) double = 3;
    n_agents (1,1) double = 6;
    position {mustBeRandomRandomBondedOrPosMatrix(position, n_agents)} = ...
        [0, -0.01;
        -0.002, -0.012;
        -0.004, -0.014;
        -0.006, -0.016;
        -0.008, -0.018;
        -0.01, -0.02];
end

% Calling the simulator function

my_output = cell_kinetics(tstop, n_agents, position);

% Defining the required variables as the relevant fields from the output
% struct
X_grid = my_output.X_grid;
Y_grid = my_output.Y_grid;
vfx_disc = my_output.vfx_disc;
vfy_disc = my_output.vfy_disc;
z_all_plot = my_output.z_all_plot;
n_agents = my_output.n_agents;
s0 = my_output.s0;
flow_data = my_output.flow_data;
t_plot = my_output.t_plot;
period = my_output.period;

grid_length = size(X_grid, 2);
grid_height = size(X_grid, 1);
floor = min(Y_grid(:,1)); % low y boundary, m
ceil = max(Y_grid(:,1)); % high y boundary, m
l_wall = min(X_grid(1,:)); % left x boundary, m
r_wall = max(X_grid(1,:)); % right x boundary, m

% Defining position and velocity data for all agents
x_plot = z_all_plot(:,1:n_agents);
y_plot = z_all_plot(:,n_agents+1:2*n_agents);
vx_plot = z_all_plot(:,2*n_agents+1:3*n_agents);
vy_plot = z_all_plot(:,3*n_agents+1:4*n_agents);
v_plot = sqrt(vx_plot.^2+vy_plot.^2);

% A list of colors to be used for plotting trajectories
col_list = ["black", "red", "blue", "magenta"];

%% Plot agent trajectories
figure
hold on
rectangle('Position', [l_wall, floor, r_wall-l_wall, ceil-floor],...
    'FaceColor',[0.55, 0.86, 0.92], 'LineWidth',2)
for i=1:n_agents
    plot(x_plot(:,i), y_plot(:,i), 'color', col_list(mod(i, 4)+1),...
        "LineWidth", 1.25);
    plot(s0(i), s0(i+n_agents), "or", "LineWidth", 1.25)
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
if r_wall-l_wall > ceil-floor % If the geometry is a wide rectangle
    xlim([l_wall*1.1 r_wall*1.1])
    ylim([l_wall*1.1 r_wall*1.1])
else % If the geometry is a square or a narrow rectangle
    xlim([floor*1.1 ceil*1.1])
    ylim([floor*1.1 ceil*1.1])
end
title(['Cell trajectories in time-dependent bioreactor simulation ' ...
    'flow field - plain background'])
subtitle('simulation time: ' + string(tstop) + [' seconds — rocking ' ...
    'period: '] + period + ' seconds')

%% Plot agent trajectories w/ initial background velocity & vorticity
figure
hold on
[curlz,cav] = curl(X_grid, Y_grid, vfx_disc, vfy_disc);
c = pcolor(X_grid, Y_grid, curlz);
c.FaceColor = 'interp';
c.EdgeColor = 'none';

% Define three colors for the gradient
color_pos = [0.9590 0.7240 0.1550];   % Yellow for positive vorticity
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
yl = ylabel(cb,'Initial vorticity','FontSize',10,'Rotation',270);

q = quiver(X_grid, Y_grid, vfx_disc, vfy_disc, 'k');

for i=1:n_agents
    plot(x_plot(:,i), y_plot(:,i), 'color', col_list(mod(i, 4)+1), ...
        "LineWidth", 1.25);
    plot(s0(i), s0(i+n_agents), "or", "LineWidth", 1.25)
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
if r_wall-l_wall > ceil-floor % If the geometry is a wide rectangle
    xlim([l_wall*1.1 r_wall*1.1])
    ylim([l_wall*1.1 r_wall*1.1])
else % If the geometry is a square or a narrow rectangle
    xlim([floor*1.1 ceil*1.1])
    ylim([floor*1.1 ceil*1.1])
end
title('Cell trajectories in time-dependent bioreactor simulation flow field')
subtitle('simulation time: ' + string(tstop) + [' seconds — rocking ' ...
    'period: '] + period + ' seconds')

%% Plot agent trajectories w/ background velocity & volume fraction in TG vortex
figure
hold on
vol_frac = flow_data(:,5,1); % Initial volume fraction data
vol_frac_grid = reshape(vol_frac, grid_length, grid_height)';
c = pcolor(X_grid, Y_grid, vol_frac_grid);
c.FaceColor = 'interp';
c.EdgeColor = 'none';

% Define three colors for the gradient
[0.9290 0.6940 0.1250];   % Yellow for 1 vf
color_pos = [0.9590 0.7240 0.1550];   % Yellow for positive vf
color_neg = [0.894, 0.914, 0.984];   % Light blue for 0 vf
% Create a custom colormap with a smooth transition between these three colors
numColors = 256;  % Number of colors in the colormap (higher resolution)
cmap = [linspace(color_neg(1), color_pos(1), numColors)', ...
    linspace(color_neg(2), color_pos(2), numColors)', ...
    linspace(color_neg(3), color_pos(3), numColors)'];

colormap(cmap);
cb = colorbar;
yl = ylabel(cb,'Initial volume fraction','FontSize',10,'Rotation',270);

q = quiver(X_grid, Y_grid, vfx_disc, vfy_disc, 'k');

for i=1:n_agents
    plot(x_plot(:,i), y_plot(:,i), 'color', col_list(mod(i, 4)+1), ...
        "LineWidth", 1.25);
    plot(s0(i), s0(i+n_agents), "or", "LineWidth", 1.25)
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
if r_wall-l_wall > ceil-floor % If the geometry is a wide rectangle
    xlim([l_wall*1.1 r_wall*1.1])
    ylim([l_wall*1.1 r_wall*1.1])
else % If the geometry is a square or a narrow rectangle
    xlim([floor*1.1 ceil*1.1])
    ylim([floor*1.1 ceil*1.1])
end
title('Cell trajectories in time-dependent bioreactor simulation flow field')
subtitle('background: initial volume fraction')
subtitle('simulation time: ' + string(tstop) + [' seconds — rocking ' ...
    'period: '] + period + ' seconds')

%% Plotting histogram of average velocities
v_avg_plot = mean(v_plot); % Average speed of each particle, m/s
figure
histogram(v_avg_plot, 20)
xlabel('velocity (m/s)')
ylabel('number of particles')
title('Average velocity distribution')
subtitle(append('simulation time: ', num2str(tstop), [' s - total ' ...
    'agents: '], num2str(n_agents)))
end
function mustBeRandomRandomBondedOrPosMatrix(position, n_agents)
eidType = 'mustBeRandomRandomBondedOrPosMatrix:notAllowedValue';
msgType = ['Input must be a string in ["random", "random bonded"], ' ...
    'or a matrix of size [n_agents, 2] column.'];
if isa(position, "string")
    if ~ismember(position, ["random", "random bonded"])
        error(eidType, msgType)
    end
elseif ~isequal(size(position), [n_agents, 2])
    error(eidType, msgType)
end
end