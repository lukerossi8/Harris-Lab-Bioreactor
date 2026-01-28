clearvars; clc; close all;
t_sim = 1; % 1-second simulation

%runtimes for different numbers of agents
trials10 = [4.366009; 4.632792; 4.209334; 4.093835; 4.433517];
t10 = mean(trials10, "all");
trials20 = [10.323133; 9.487368; 9.571022; 9.615910; 9.425754];
t20 = mean(trials20, "all");
trials30 = [18.541262; 17.930895; 18.358983; 18.172831; 17.823959];
t30 = mean(trials30, "all");
t50 = 43.194785;
t70 = 84.825573;
t100 = 169.708278;
t150 = 386.332793;
t200 = 693.065624;

agent_counts = [10; 20; 30; 50; 70; 100; 150; 200];
times = [t10; t20; t30; t50; t70; t100; t150; t200];

% loglog(agent_counts, times, '-o', 'Color',[0.5 0 0.8], 'LineWidth',1.5)
% grid on
% xlabel('number of agents')
% ylabel('run time (s)')
% title('run times vs. number of agents')
% subtitle('simulated time: 1 second')
% xlim([8 300])
% ylim([1 10^3])

t10grid = 7.564287;
t20grid = 8.054310;
t30grid = 8.172512;
t50grid = 8.291364;
t70grid = 10.725395;
t100grid = 9.426697;
t150grid = 13.683561;
t200grid = 13.896612;

times_grid = [t10grid, t20grid, t30grid, t50grid, t70grid, t100grid, t150grid, t200grid];

% figure
% loglog(agent_counts, times_grid, '-o', 'Color',[0.5 0 0.8], 'LineWidth',1.5)
% grid on
% xlabel('number of agents')
% ylabel('run time (s)')
% title('run times vs. number of agents with grid optimization')
% subtitle('simulated time: 1 second')
% xlim([8 300])
% ylim([1 10^3])


% varying number of bonds total agent count = 20
t0bonds = 7.095221;
t2bonds = 9.129771;
t4bonds = 9.962552;
t6bonds = 11.708535;
t8bonds = 11.933175;
t10bonds = 12.914513;

bond_counts_20 = [0,2,4,6,8,10];
times_bond_20 = [t0bonds, t2bonds, t4bonds, t6bonds, t8bonds, t10bonds];

figure
plot(bond_counts_20, times_bond_20, '-o', 'Color',[0.5 0 0.8], 'LineWidth',1.5)
grid on
xlabel('number of bonds')
ylabel('run time (s)')
title('run times vs. number of bonds in 20-agent simulation with grid optimization')
subtitle('simulated time: 1 second')
xlim([0,10])
ylim([0,100])


% varying number of bonds, total agent count = 200
t0_200 = 11.615815;
t10_200 = 16.191390;
t20_200 = 17.992354;
t30_200 = 19.512927;
t40_200 = 22.213033;
t50_200 = 25.517081;
t60_200 = 25.529975;
t70_200 = 26.653953;
t80_200 = 37.581334;
t90_200 = 38.490351;
t100_200 = 37.160540;

bond_counts_200 = [0,10,20,30,40,50,60,70,80,90,100];
times_bond_200 = [t0_200,t10_200,t20_200,t30_200,t40_200,t50_200,t60_200,t70_200,t80_200,t90_200,t100_200];

figure
plot(bond_counts_200, times_bond_200, '-o', 'Color',[0.5 0 0.8], 'LineWidth',1.5)
grid on
xlabel('number of bonds')
ylabel('run time (s)')
title('run times vs. number of bonds in 200-agent simulation with grid optimization')
subtitle('simulated time: 1 second')
xlim([0,100])
ylim([0,100])

