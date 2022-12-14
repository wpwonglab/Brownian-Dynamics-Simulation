function bd_main_blake_equi_ID_v11(name,np,u1,flow,rep)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Filename, number of particles, LJ parameter, flow rate, trial number 
% name  = root name of save files
% parameter_out = parameters for the simulation.
% np = number of beads 
% u1 = lennard jones constant in kb/t. Values <0.314 -> uncollapsed polymer. 
% Values = 0.314 are ideal, Values >0.314 are collapsed polymers  
% flow = integer value which controls the flow rate based on exponential base 2
% 1 corresponds to 2000 s^-1, 8 corresponds to 128,000 s^-1. Dependent on
% real radius of bead and should be modfiied as needed below

% rep = used for seeding rng in cluster

% Flow veloctiy can be modified to look at flow extension
% or flow relaxation. 
% Equilibration and production times should be modififed as needed 

% by Hans Bergal, Darren Yang and Wesley Wong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up
k  = 1;             % Boltzmann Constant (nm2*kg/s2/K)
T  = 1;             % Temperature (K) 
kT = k*T;           % Thermal Energy (nm2*kg/s2)
mu = 1;             % Viscosity (kg/(s*nm), Pa*s)
u  = u1*kT;         % Lennard-Jones potential Strength (kT)
a  = 1;             % Sphere Radius (nm)
L  = 2*a;           % Bond Length (nm)
n3 = 3*np;          % 3 x Number of Particle
h  = L;  % Distance from surface (a)  
K  = 800.*kT/a^2; % Bond Spring Constant (J/nm^2)
tK = K;    % Tethered Spring Constant (J/nm^2)

real_a = 3.7; %real radius in nm%
real_tau = 4.0522*10^-9 .*real_a.^3; % real tau in seconds
flow_const = 128000.*real_tau./2^8; % Converts 128000s^-1 to correspond to input flow 8
Wi = (flow_const.*2.^(flow)); % Define Wiesenberg number for shear rate 



%make file names
init=strcat('init_',name);
para=strcat('para_',name);
traj=strcat('traj_',name);
siza=strcat('siza_',name);
endt=strcat('endt_',name);
gyra=strcat('gyra_',name);

% Simulation Parameters
n_tau  = 100*np;           % Number of tau in Production
n_step = 5000;             % Number of step per tau
n_equi = 4*np*n_step;     % Number of Step in Equilibration
n_prod = n_tau*n_step/8;     % Number of Step in production
w_step = 500;       % Number of Step between each output
tau    = 6*pi*mu*a^3/(kT); % Relaxation Time of a single Sphere
g      = Wi ./tau;       %Convert unitless shear parameter to shear rate 
dt     = tau/n_step;       % delta-Time per step
period = 100;              % steps between HD update 

%% Set Up random number generator for cluster 
rng('shuffle')
for i = 1:rep
    rand(100000,1);
end
rng(randi(2^32-1));


%% load init_file if exists 
rep_num = num2str(rep);
intial_conf = strcat('init_',num2str(np),'_',rep_num);

if isfile(intial_conf )
[~, x, y, z] = bd_read(intial_conf);
x = x(1:np);
y = y(1:np);
z = z(1:np);

else
% or  Randomly Place Chain
 [x, y, z] = bd_conf_tether(np, 2*a, init, 0, 'linear');
 z = z + h;	
end 

%% Write initial configuration to a Initial.xyz
conf_out = fopen(init, 'w');
fprintf(conf_out,'%d\n', np);
fprintf(conf_out,'Initial Position generated by bd_conf.m \n');
for I = 1:np
    fprintf(conf_out,'D %+6.4f %+6.4f %+6.4f \n', x(I), y(I), z(I));
end

%% Create output file
p_out = fopen(para, 'a');
c_out = fopen(traj, 'a');
b_out = fopen(siza, 'a');
e_out = fopen(endt, 'a');
g_out = fopen(gyra, 'a');
% Write out the parameter to parameter_out for analysis
fprintf(p_out, 'Boltzmann Constant : %6.4e (m2*kg/s2/K) \n', k);
fprintf(p_out, 'Temperature        : %6.4e (K) \n', T);
fprintf(p_out, 'Solvent Viscosity  : %6.4e (N*s/m2) \n', mu);
fprintf(p_out, 'LJ Strength        : %6.4e (J) \n', u);
fprintf(p_out, 'Sphere Radius      : %6.4e (m) \n', a);
fprintf(p_out, 'Bond Strength      : %6.4e (m) \n', K);
fprintf(p_out, 'Bond Length        : %6.4e (m) \n', L);
fprintf(p_out, 'Relaxation Time    : %6.4e (s) \n', tau);
fprintf(p_out, 'Steps per tau      : %6.4e (m) \n', n_step);
fprintf(p_out, 'Equilibration Steps: %6.4e     \n', n_equi);
fprintf(p_out, 'Production Steps   : %6.4e     \n', n_prod);
fprintf(p_out, 'Save Step          : %d        \n', w_step);
fprintf(p_out, 'Equilibration (tau): %6.4e     \n', n_equi/n_step);
fprintf(p_out, 'Production (tau)   : %6.4e     \n', n_prod/n_step);
fprintf(p_out, 'Save Step (tau)    : %d        \n', w_step/n_step);
fprintf(p_out, 'delta Time Step    : %6.4e (tau) \n', dt);
fprintf(p_out, 'Number of Particle : %d        \n', np);
fprintf(p_out, 'Ideal Gyration     : %6.4e (m) \n', np*L^2/6);
fprintf(p_out, 'Wi                 : %6.4e     \n', Wi);
fprintf(p_out, 'shear rate: %6.4e \n', g);
fprintf(p_out, 'HD tensor update step    : %6.4e \n', period);
%% Start!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equilibrate with no Flow
fprintf('\n Equilibration Step 1 \n');
Vx=0; Vy=0; Vz = 0; % set flow to zero
D = bd_hyd_Blake_v1(np,n3,a,L,mu,x,y,z); % calculate inital hydrodynamic tensor
D = D/(6*pi*mu*a);
[B,~] = chol(kT*D,'lower');
  for i = 1:n_equi
     [~, TFX, TFY, TFZ, ~] = bd_force_surf_ID(np,n3,u,a,K,L,x,y,z,mu);
          if mod(i,period) == 0
            D = bd_hyd_Blake_v1(np,n3,a,L,mu,x,y,z);
            D = D/(6*pi*mu*a);
            [B,~] = chol(kT*D,'lower');
          end
      [x, y, z] = bd_move_blake_v2(n3, x, y, z, Vx, Vy, Vz, TFX, TFY, TFZ, D,B, dt, kT,mu,a);
  end

 
%% Extension
fprintf('\n Production Step \n');

D = bd_hyd_Blake_v1(np,n3,a,L,mu,x,y,z);
D = D/(6*pi*mu*a);
[B,~] = chol(kT*D,'lower');
bd_write(c_out,b_out,e_out,g_out,np,x,y,z,kT);

    for i = 1:(n_prod)
        [Vx, Vy, Vz] = bd_fluid_v2(x, y, z, g);
        [~, TFX, TFY, TFZ, ~] = bd_force_surf_ID(np,n3,u,a,K,L,x,y,z,mu);

       if mod(i,period) == 0
            D = bd_hyd_Blake_v1(np,n3,a,L,mu,x,y,z);
            D = D/(6*pi*mu*a);
           [B,~] = chol(kT*D,'lower');
       end
       
        [x, y, z] = bd_move_blake_v2(n3, x, y, z, Vx, Vy, Vz, TFX, TFY, TFZ, D,B, dt, kT,mu,a);
 
        if mod(i,w_step) == 0
            bd_write(c_out,b_out,e_out,g_out,np,x,y,z,kT);
        end   
    end
	fprintf('\n Done! \n');
end