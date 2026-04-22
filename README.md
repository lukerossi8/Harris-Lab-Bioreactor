# Harris-Lab-Bioreactor
This repository contains a collection of MATLAB functions that model cell dynamics in a rocking bioreactor environment. It models cell motion in a periodic background flow field based on prior CFD simulation data. It also models bonding interactions between cells. Users can customize the length of the simulation, the number of cells, and the cells' initial positions.

## Overview
Cell motion is modeled by implementing four forces on each cell: Gravity/buoyancy, drag, a boundary interactive (repulsion) force, and a cell bonding force. Using these forces and a periodic fluid flow field imported from previous simulation data, the repository is able to compute each cell's position throughout the simulation. This data can be used to plot graphs and videos of the cells' trajectories. It can also be used to evaluate the shear stress experienced by the cells throughout the life of the simulation.

## Background
 Rocking bioreactors are a promising candidate in  cultivated meat bioprocess optimization. They impart relatively low shear stresses on cells, along with being disposable and cost-effective. 

## Quickstart guide
Follow the steps outlined in this section to gain familiarity with the core functionality of the repository, and to prepare to run simulations of your own design.

### Cloning the repository
### Preparing to run a simulation
### Creating a basic simulation
Let's start by modeling a simple case — the repository's default simulation — to ensure everything is working smoothly. This simulation will model the motion of 6 cells in the rocking bioreactor for 0.5 seconds. To run this simulation, simply call `sim_flowfield_figure_maker()` in your command window or terminal. This will produce a figure of the particle trajectories, as well as plots of the cell velocities over time and the shear stress experienced by each cell over time.

To generate a video of the same default simulation, call `vol_frac_video_maker()` in your command window or terminal. This will produce a video of the cells moving through the bioreactor environment.

### A more complex simulation
Now, let's model a more complex scenario. This time, we'll model 100 cells, randomly placed, for 10 seconds. To run this simulation, call `sim_flowfield_figure_maker(10, 100, 'random')` in your command window or terminal. The same types of figures will be produced — the cell trajectories, velocities over time, and shear stress over time. To generate the video of this simulation, call `vol_frac_video_maker()` in your command window or terminal. 

### Creating your own simulations
Now that you've run a few simulations, you can create new ones by varying the parameters in the figure maker and/or the video maker, calling either one or both as best suits your needs. For more information on the parameters that can be changed to create unique simulations, see the "Parameters" sub-section of the "Using the repository" section.

## Using the repository
### Files
### Parameters 
You can modify three parameters in these simulations: the length of the simulation (`tstop`), the number of cells in the simulation (`n_agents`), and the initial positions of the cells (`position`). So, the syntax for calling either the figure maker is `sim_flowfield_figure_maker(tstop, n_agents, position)`, and for the video maker, it's the same: `vol_frac_video_maker(tstop, n_agents, position)`.

Let's talk about each parameter.

`tstop` is the amount of time being simulated in the simulation, in seconds. The default value is 0.5 seconds.

`n_agents` is the number of cells being modeled in the simulation. The default value is 6 cells.

`position` is where you define the initial positions of the cells. There are three types of initial positions: random, random bonded, or manual. 
- Random positions mean the cells are allocated randomly throughout the domain, with a buffer zone to prevent them from being initialized right next to the walls. Set `position = "random"`
- Random bonded positions mean that the cells are randomly positioned as well, but under the condition that each cell will begin right next to another cell, such that it will bond with that cell immediately. If there are an odd number of cells in the simulation, the last one will not be bonded to any other cell. Set `position = "random bonded"`
- Manual positions allow you to define exactly where each cell will be initialized. To do this, set `position` equal to a matrix of size [n_agents, 2]. Each row represents the position of each agent; the first column is the x position (range = [-0.5, 0.5]), and the second is the y position (range = [-0.15, 0.15]). Make sure to create a matrix that is the correct size - that is, that you define a position for all the cells, and no more. For example, here is the position matrix for the default case, where `n_agents = 6`:
> `position = [0, -0.05;
        -0.01, -0.06;
        -0.02, -0.07;
        -0.03, -0.08;
        -0.04, -0.09;
        -0.05, -0.10];`

### Using your own flow data

## Interpreting the results

## Equations

## Contributing
