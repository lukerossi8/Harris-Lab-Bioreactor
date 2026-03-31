function video_maker_vol_frac(tstop, n_agents, position)
% figure_maker_vol_frac creates a video of agent trajectories in a
% bioreactor environment with the volume fraction in the background, by
% calling the cell_contact_sim_flowfield function
% Syntax:
%   sim_flowfield_figure_maker(tstop, n_agents, position)
% Inputs:
%   tstop       Double representing the length of the simulation (seconds)
%   n_agents    Double representing the number of agents in the simulation
%   position    Either a string in ["random", "random bonded"], or a matrix
%               of size [n_agents, 2] of doubles (column 1 = x vals,
%               column 2 = y vals)
% Outputs
%   A video modeling the trajectories of the agents throughout the
%   simulation.

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
times = my_output.times;

% Defining position and velocity data for all agents
x_plot = z_all_plot(:,1:n_agents);
y_plot = z_all_plot(:,n_agents+1:2*n_agents);
vx_plot = z_all_plot(:,2*n_agents+1:3*n_agents);
vy_plot = z_all_plot(:,3*n_agents+1:4*n_agents);
v_plot = sqrt(vx_plot.^2+vy_plot.^2);

% A list of colors to be used for plotting trajectories
col_list = ["black", "red", "blue", "magenta"];

% Plot agent trajectories over time and create the video

% Video writer setup
vid = VideoWriter(append('cell_traj_vf_', num2str(tstop), 's.mp4'), 'MPEG-4');
vid.Quality = 100;
vid.FrameRate = 30;
open(vid);

% Creating the video loop
fig = figure;
step = 2000;
for i=1:step:length(t_plot)
    cla
    hold on

    t_eff = mod(t_plot(i), times(end)); % effective time inside a single period

    % Find i such that times(i) <= t_eff < times(i+1)
    layer = find(times <= t_eff, 1, 'last');

    vol_frac = flow_data(:,5, layer); % volume fraction data at current time
    vol_frac_grid = reshape(vol_frac, 64, 64)';
    c = pcolor(X_grid, Y_grid, vol_frac_grid);
    c.FaceColor = 'interp';
    c.EdgeColor = 'none';

    % Define three colors for the gradient
    color_pos = [0.9290 0.6940 0.1250];   % Yellow for 1 vf
    color_neg = [1 1 1];   % White for 0 vf
    % Create a custom colormap with a smooth transition between these three colors
    numColors = 256;  % Number of colors in the colormap (higher resolution)
    cmap = [linspace(color_neg(1), color_pos(1), numColors)', ...
        linspace(color_neg(2), color_pos(2), numColors)', ...
        linspace(color_neg(3), color_pos(3), numColors)'];

    colormap(cmap);
    cb = colorbar;
    yl = ylabel(cb,'Volume fraction','FontSize',10,'Rotation',270);
    
    for j=1:n_agents
        plot(x_plot(1:i,j), y_plot(1:i,j), 'color', col_list(mod(j, 4)+1), "LineWidth", 1.25);
        plot(x_plot(i,j), y_plot(i,j), "or", "LineWidth", 1.25)
    end
    grid on
    xlabel('x position (m)')
    ylabel('y position (m)')
    xlim([-0.6 0.6])
    ylim([-0.2 0.2])
    title('Cell trajectories in time-dependent bioreactor simulation flow field - plain background')
    subtitle('simulation time: ' + string(tstop) + ' seconds — period time: 0.595 seconds')

    % Capture the frame
    frame = getframe(gcf);
    writeVideo(vid, frame);
end

% Close video writer
close(vid);
close(fig);

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