function [x, y, z] = bd_move_blake_v2(n3, x, y, z, Vx, Vy, Vz, TFX, TFY, TFZ,D,B, dt, kT,mu,a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x, y, z] = bd_move(n3, x, y, z, Vx, Vy, Vz, TFX, TFY, TFZ, D, dt, kT)
%  
% Subroutine to move particles once the following variables are calculated
% np = number of particle x 3
% x, y, and z = position 
% Vx, Vy, and Vz = unperturbed velocity 
% TFX, TFY, and TFZ = total force 
% D = mobility matrix (not multiply by kT)
% dt = time step
% kT = Thermal Energy
%
% by Darren Yang (ydarren@gmail.com) and Wesley Wong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = zeros(n3, 1);
R(1:3:end) = x;
R(2:3:end) = y;
R(3:3:end) = z;


V = zeros(n3, 1);
V(1:3:end) = Vx;
V(2:3:end) = Vy;
V(3:3:end) = Vz;


F = zeros(n3, 1);
F(1:3:end) = TFX;
F(2:3:end) = TFY;
F(3:3:end) = TFZ;


h = (2^0.5) * dt^0.5;
flow = (diag(D)).*V.*(dt.*6*pi*mu*a);
force = D*(F*dt); 
therm = B*randn(n3, 1);
therm = h * therm;
R = R + flow+force+therm;

%  R = R +(diag(D)).*V.*(dt.*6*pi*mu*a) + D*F*dt + 2^0.5*B*randn(n3, 1)*dt^0.5;
% Langevin update step, breaking apart was found to be computationally
% faster 

x = R(1:3:end);
y = R(2:3:end);
z = R(3:3:end)+kT*dt/(6*pi*mu*a)*(9/8*a./(z.^2)-3/2*a^3./(z.^4));

end
