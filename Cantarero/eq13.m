function eq13_cantarero
clearvars; close all; clc;
% Time range, constant values, and initial conditions
tspan = [0,100];
g = -9.81; % m/s^2
p_m = 1; % g/ml
p_i = 1.1; %g/ml, average cell density
v0 = 0;

% Define the anonymous function
odeFun = @(t, v) eom(t, v, g, p_m, p_i);

[t, v] = ode45(odeFun, tspan, v0);
plot(t, v)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Solution of the buoyancy/gravity force ODE')
grid on

end

function dvdt = eom(t, v, g, p_m, p_i)
dvdt = g*(1 - p_m/p_i);
end


% function eq13_cantarero
% clearvars; close all; clc;
% % Time range, constant values, and initial conditions
% tspan = [0,100];
% v0 = 0;
% [t, v] = ode45(@eom, tspan, v0);
% plot(t, v)
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% title('Solution of the buoyancy/gravity force ODE')
% grid on
% end
% 
% function dvdt = eom(t, v)
% g = -9.81; % m/s^2
% pm = 1; % g/ml
% pi = 1.1; %g/ml, average cell density
% dvdt = g*(1 - pm/pi);
% end
