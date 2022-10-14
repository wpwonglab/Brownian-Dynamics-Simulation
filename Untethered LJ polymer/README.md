# Brownian-Dynamics-Simulation

Brownian dynamic polymer simulation under flow for a tethered polymer 
by Hans Bergal, Darren Yang and Wesley Wong

Written in Matlab and Fortran,
Run in Matlab 


Main function is "bd_main_RPY_free_ID_v11(name,np,u1,flow,rep)"

%Filename, number of particles, LJ parameter, flow rate, trial number 
% name  = root name of save files
% parameter_out = parameters for the simulation.
% np = number of beads 
% u1 = lennard jones constant in kb/t. Values <0.314 uncollapsed polymer. 
% Values = 0.314 are ideal, Values >0.314 are collapsed polymers  
% flow = integer value which controls the flow rate based on exponential base 2
% 1 corresponds to 2000 s^-1, 8 corresponds to 128,000 s^-1. Dependent on
% real radius of bead and should be modfiied as needed below

% rep = used for seeding rng in cluster

% Flow veloctiy can be modified to look at flow extension
% or flow relaxation. 
% Equilibration and production times should be modififed as needed 
% Detailed methods can be found in supplemental of 
% "https://ashpublications.org/blood/article/doi/10.1182/blood.2022016969/486450/Conformation-of-von-Willebrand-factor-in-shear"


Flow is applied from negative to postive x direction and velocity is dependent on height from the surface in Z

Requires compiling of two fortran mex files ("bd_RPY_hyd_v1.F" and "bd_force_free_ID.F")

Example compiling in matlab with proper fortran compiler setup

"mex bd_RPY_hyd_v1.F"

Mex fortran compiler set up help:
https://www.mathworks.com/matlabcentral/answers/472860-mex-setup-fortran


Trajectories can be visualized in VMD file type XYZ


