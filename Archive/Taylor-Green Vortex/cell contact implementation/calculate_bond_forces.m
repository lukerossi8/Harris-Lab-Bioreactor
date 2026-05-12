% For bond in bonded pairs,
% calculate the distance and angle between the agents
% calculate the forces on each agent
% add to a matrix containing sum of bond forces on each agent 
clearvars; clc;
r_i = [90e-6; 90e-6; 90e-6; 90e-6; 90e-6]; % Agent radii, in m
m = [2.79e-9; 2.79e-9; 2.79e-9; 2.79e-9; 2.79e-9]; % Agent mass, in kg
n_agents = length(m);
delta_c = 2*r_i(1); % Bond formation threshold, m
% ^Highly simplified, principle is threshold = sum of radii of the two agents
delta_d = 1.4*delta_c; % Bond breaking threshold, m
bonded_pairs = [1, 2; 2, 3; 3, 5];
z= [0; sqrt(2)*r_i(1); sqrt(2)*r_i(1)+2*r_i(1); 0.1; sqrt(2)*r_i(1)+2*r_i(1);
    0; sqrt(2)*r_i(1); sqrt(2)*r_i(1); sqrt(2)*r_i(1); sqrt(2)*r_i(1)+2*r_i(1); 
    0; 0; 0; 0; 0;
    0; 0; 0; 0; 0];

bond_forces_x = zeros(n_agents, 1);
bond_forces_y = zeros(n_agents, 1);

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


