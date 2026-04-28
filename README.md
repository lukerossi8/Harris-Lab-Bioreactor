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
The simulation solves the equation of motion for each cell as it moves through the simulated bioreactor environment. 

### Overall equation of motion
Equation (1) defines the forces experienced by each cell and modeled in this simulation.

1. $dvdt = \vec{F_g} + \vec{F_d} + ∑\vec{F_{bond}} + \vec{F_{boundary}}$

Note that in the code base, the above equation, as well as the following equations for each force, are broken into their $x$ and $y$ components. For simplicity, here, they are presented in vector form. In the following sub-sections, each term in the above equation is explained.

### Gravity and buoyancy equations
$\vec{F_g}$ is the force due to gravity and buoyancy, modeled by equation (2):

2. $\vec{F_g} = m\vec{g}(1-\frac{\rho_{local}}{\rho_c})$

Where:
- $m$ is the mass of the cell;
- $g$ is the gravity vector — see equations (3) and (4);
- $\rho_{local}$ is the local density of the fluid; and
- $\rho_c$ is the density of the cell.

Because the bioreactor environment is rotating, the gravity vector must evolve in time to reflect this. The gravity vector's components, $g_x$ and $g_y$, are modeled by equations (3) and (4), respectively:

3. $g_x = 9.81\sin(\theta)$
4. $g_y = -9.81\cos(\theta)$

Where $\theta$ is given by equation (5):

5. $\theta = A\sin(\frac{2\pi t_{eff}}{T})$;

Where:
- A is the rocking amplitude of $7^{\circ}$;
- $t_{eff}$ is the effective time, meaning the simulated time modded by the period; and
- $T$ is the period of ~ $0.5953$ seconds.

### Drag equation
The drag force is based upon Stokes' flow past a sphere. It is implemented using equation (6):

6. $\vec{F_d} = -6\pi \mu_{local} r_c \vec{v_{rel}}$

Where:
- $\mu_{local}$ is the local density of the background fluid;
- $r_c$ is the radius of the cell; and
- $\vec{v_{rel}}$ is the relative velocity of the cell to the background fluid, $\vec{v_{rel}} = \vec{v_{cell}} - \vec{v_{fluid}}$.

### Cell bonding equations
When cells come into contact with one another, they bond together. For bonded pairs of cells, a bonding force is applied between them according to equation (7):

7. $|\vec{F_{bond}}| = K_{ij} \delta_{ij} \tanh(s_{ij}|\delta_{ij}|)$

Where:
- $K_{ij}$ is the bond spring constant of $10^{-3}$ Nm;
- $\delta_{ij}$ is the overlap of the two cells — see equation (8); and
- $s_{ij}$ is the bond stiffness parameter set at a value of $0.2$.
- To find the x and y components, $|\vec{F_{bond}}|$ is multiplied by $\cos{\theta_b}$ and $\sin{\theta_b}$, respectively, where $\theta_b$ is the angle between the two bonded cells.

The overlap of the two bonded cells, $\delta_{ij}$, is given by equation (8):

8. $\delta_{ij} = D_{ij} - (R_i + R_j)$

Where $D_{ij}$ is the distance between the centers of the two cells, and $(R_i + R_j)$ is the sum of the radii of the cells.

The cell bonding equation is selectively applied to pairs of cells depending on whether an active bond exists between them. A bond is formed between cells when their overlap $D_{ij}$ falls below 0, meaning they have come into contact. This is modulated by the parameter $\delta_c$, set equal to $(R_i + R_j)$. Once a bond forms, it remains active until the cells move far enough apart to break the bond. This is modulated by the parameter $\delta_d$, set equal to $1.4\delta_c$. A matrix tracks the bonded pairs of cells throughout the simulation.

### Boundary repulsion equations
Cells are kept within the bounds of the simulation, i.e. the walls of the bioreactor, by a repulsive force active at each wall of the bioreactor. The force is formulated slightly differently for each wall, as shown in equations (9)–(12):

9. $F_{wall-left} = A_{wall}e^{-B_{wall}(x-l_{wall})}+K_{wall}(l_{wall}-x)$
10. $F_{wall-right} = -A_{wall}e^{-B_{wall}(r_{wall}-x)}+K_{wall}(x-r_{wall})$
11. $F_{wall-low} = A_{wall}e^{-B_{wall}(y-f)}+K_{wall}(f-y)$
12. $F_{wall-high} = -A_{wall}e^{-B_{wall}(c-y)}+K_{wall}(y-c)$

Where:
- $A_{wall}$ is the base soft wall force at the boundary, with a value of $10^{-4}$ N;
- $B_{wall}$ is the decay constant of $10^5$ 1/m;
- $K_{wall}$ is the spring constant of $1.0$ N/m; and
- $l_{wall}$, $r_{wall}$, $f$, and $c$ are the one-dimensional coordinates of the left wall ($-0.15$ m), right wall ($0.15$ m), floor ($-0.05$ m), and ceiling ($0.05$ m), respectively, of the bioreactor.

The forces are modulated so as to only take effect when an agent breaches one of the boundaries.

## Contributing
