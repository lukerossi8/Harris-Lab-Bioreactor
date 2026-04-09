%% Plotting runtime vs. number of agents
% Holding simulated time constant at 5 seconds
agent_counts = [1; 50; 100; 150; 200; 250; 300];
runtimes = zeros(1, size(agent_counts, 1));

for i=1:length(agent_counts)    
    my_output = cell_contact_sim_flowfield(5, agent_counts(i), "random");
    runtimes(i) = my_output.runtime/60;
end

figure
hold on
plot(agent_counts, runtimes, '-o');
xlabel('number of agents')
ylabel('elapsed time (min)')
title('Runtime vs. number of agents')
subtitle('random position setting, simulated time = 5 sec')
