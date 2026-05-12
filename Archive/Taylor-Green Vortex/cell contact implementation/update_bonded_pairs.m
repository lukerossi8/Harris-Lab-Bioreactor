% clearvars; clc;
% r_i = [90e-6; 90e-6; 90e-6]; % Agent radii, in m
% m = [2.79e-9; 2.79e-9; 2.79e-9]; % Agent mass, in kg
% n_agents = length(m);
% bonded_pairs = [0, 0];
% delta_c = 2*r_i(1); % Bond formation threshold, m
% % ^Highly simplified, principle is threshold = sum of radii of the two agents
% delta_d = 1.4*delta_c; % Bond breaking threshold, m
% 
% agent_indices = (1:size(m, 1))'; % Column vector with elts 1 to n_agents
% pairs = nchoosek(agent_indices, 2); % All possible pairs of agents
% for i = 1:size(pairs, 1)
%     current_pair = pairs(i, :);
%     agent_i = current_pair(:, 1); % Lower indexed agent in the pair
%     agent_j = current_pair(:, 2); % Higher indexed agent in the pair
% 
%     % Measure distances and angle between the two agents
%     % z vector takes the form [x1; x2; x3;...; y1; y2; y3;...; vx1; vx2; vx3;...; vy1; vy2; vy3;...]
%     z = [0.1; 0.1; 0.1; 0.2; 0.2+91e-6; 0.2-91e-6; 0; 0; 0; 0; 0; 0];
%     dist_x = z(agent_j)-z(agent_i);
%     dist_y = z(n_agents+agent_j)-z(n_agents+agent_i);
%     dist = sqrt((dist_x)^2 + (dist_y)^2); % Distance between the two agents, m
% 
%     % Check if the agents are already bonded
%     is_bonded = ismember(current_pair, bonded_pairs, "rows");
%     % If the agents are already bonded, check if the bond has broken
%     if is_bonded
%         if dist > delta_d % If the bond has broken in the previous timestep
%             bonded_pairs(bond, :) = []; % Remove the pair (row) from bonded_pairs
%         end
%     else % the agents are not already bonded
%         if dist <= delta_c % the agents have bonded in the previous timestep
%             if bonded_pairs == [0,0]
%                 bonded_pairs = current_pair;
%             else
%                 % Add the pair to bonded_pairs
%                 bonded_pairs(end+1, :) = current_pair;
%             end
%         end
%     end
% end

a = randi([1, 201], 1, 100);
a(a > 100) = a(a > 100) + 99;
a