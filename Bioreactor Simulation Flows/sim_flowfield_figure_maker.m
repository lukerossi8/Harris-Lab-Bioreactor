% Probably should put a function wrapper here too so that users can call
% this from the command window with their own inputs, rather than having to
% enter their inputs into this script (thus still needing to edit a script)

function sim_flowfield_figure_maker(tstop, n_agents, position)
% sim_flowfield_figure_maker creates visualizations of agent behavior in a
% bioreactor environment, by calling the cell_contact_sim_flowfield
% function
% Syntax:
%   sim_flowfield_figure_maker(tstop, n_agents, position)
% Inputs:
%   tstop       Double representing the length of the simulation (seconds)
%   n_agents    Double representing the number of agents in the simulation
%   position    Either a string in ["random", "random bonded"], or a matrix 
%               of size [n_agents, 2] of doubles (column 1 = x vals, 
%               column 2 = y vals)
% Outputs
%   Figures representing the trajectories of the agents with the background 
%   vorticity, the trajectories with the background volume fraction, and 
%   the timestep size over the course of the simulation.

arguments
    tstop (1,1) double = 0.5;
    n_agents (1,1) double = 6;
    position {mustBeRandomRandomBondedOrPosMatrix(position, n_agents)} = ...
        [0, -0.05;
        -0.01, -0.06;
        -0.02, -0.07;
        -0.03, -0.08;
        -0.04, -0.09;
        -0.05, -0.10];
end

% Calling the simulator function

my_output = cell_contact_sim_flowfield(tstop, n_agents, position);

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
rectangle('Position', [-0.5, -0.15, 1, 0.3], 'FaceColor',[0.55, 0.86, 0.92],...
    'LineWidth',2)
for i=1:n_agents
    plot(x_plot(:,i), y_plot(:,i), 'color', col_list(mod(i, 4)+1), "LineWidth", 1.25);
    plot(s0(i), s0(i+n_agents), "or", "LineWidth", 1.25)
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
xlim([-0.6 0.6])
ylim([-0.2 0.2])
title('Cell trajectories in time-dependent bioreactor simulation flow field - plain background')
subtitle('simulation time: ' + string(tstop) + ' seconds — period time: 0.595 seconds')

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
yl = ylabel(cb,'Vorticity','FontSize',10,'Rotation',270);

q = quiver(X_grid, Y_grid, vfx_disc, vfy_disc, 'k');

for i=1:n_agents
    plot(x_plot(:,i), y_plot(:,i), 'color', col_list(mod(i, 4)+1), "LineWidth", 1.25);
    plot(s0(i), s0(i+n_agents), "or", "LineWidth", 1.25)
end
grid on
xlabel('x position (m)')
ylabel('y position (m)')
ylim([-0.15 0.15])
title('Cell trajectories in time-dependent bioreactor simulation flow field')
subtitle(append('simulation time: ', num2str(tstop), ' seconds — period time: 0.595 seconds'))

%% Plot agent trajectories w/ background velocity & volume fraction in TG vortex
figure
hold on
vol_frac = flow_data(:,5,1); % INITIAL volume fraction data
vol_frac_grid = reshape(vol_frac, 64, 64)';
c = pcolor(X_grid, Y_grid, vol_frac_grid);
c.FaceColor = 'interp';
c.EdgeColor = 'none';

% Define three colors for the gradient
[0.9290 0.6940 0.1250];   % Yellow for 1 vf
color_pos = [0.9590 0.7240 0.1550];   % Yellow for positive vorticity
color_neg = [1 1 1];   % White for 0 vf
% Create a custom colormap with a smooth transition between these three colors
numColors = 256;  % Number of colors in the colormap (higher resolution)
cmap = [linspace(color_neg(1), color_pos(1), numColors)', ...
    linspace(color_neg(2), color_pos(2), numColors)', ...
    linspace(color_neg(3), color_pos(3), numColors)'];

colormap(cmap);
cb = colorbar;
yl = ylabel(cb,'Volume fraction','FontSize',10,'Rotation',270);

q = quiver(X_grid, Y_grid, vfx_disc, vfy_disc, 'k');

for i=1:n_agents
    plot(x_plot(:,i), y_plot(:,i), 'color', col_list(mod(i, 4)+1), "LineWidth", 1.25);
    plot(s0(i), s0(i+n_agents), "or", "LineWidth", 1.25)
end
grid on
ylim([-0.15 0.15])
xlabel('x position (m)')
ylabel('y position (m)')
title('Cell kinetics model in bioreactor simulation flow field')


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
title('Time step size throughout the simulation')

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
% dist_plot_1_3 = sqrt((x_plot(:, 3) - x_plot(:, 1)).^2 + (y_plot(:, 3) - y_plot(:, 1)).^2);
% dist_plot_2_4 = sqrt((x_plot(:, 4) - x_plot(:, 2)).^2 + (y_plot(:, 4) - y_plot(:, 2)).^2);

% figure
% hold on
% plot(t_plot, dist_plot_2_4, 'LineWidth',2)
% yline(2*r_i(1))
% grid on
% xlabel('time (s)')
% ylabel('particle distance (m)')
% title('distance between agents 2 and 4 — with bonding')
% 
% figure
% hold on
% plot(t_plot, dist_plot_1_3, 'LineWidth',2)
% yline(2*r_i(1))
% grid on
% xlabel('time (s)')
% ylabel('particle distance (m)')
% title('distance between agents 1 and 3 — with bonding')

function mustBeRandomRandomBondedOrPosMatrix(position, n_agents)
eidType = 'mustBeRandomRandomBondedOrPosMatrix:notAllowedValue';
msgType = 'Input must be a string in ["random", "random bonded"], or a matrix of size [n_agents, 2] column.';
if isa(position, "string")
    if ~ismember(position, ["random", "random bonded"])
        error(eidType, msgType)
    end
elseif ~isequal(size(position), [n_agents, 2])
    error(eidType, msgType)
end