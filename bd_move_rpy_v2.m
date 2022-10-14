function [x, y, z] = bd_move_rpy_v2(n3, x, y, z, Vx, Vy, Vz, TFX, TFY, TFZ, D,B,dt, kT)
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

% R = R +(diag(D)).*V.*(dt.*6*pi*mu*a) + D*F*dt + 2^0.5*B*randn(n3, 1)*dt^0.5;
% update equation 

% Calculating in parts if faster than a single line 
vel = (diag(D)).*V.*(dt.*6*pi); 
force = D*F*dt; 
brown = B*randn(n3, 1);
bmod = h * brown;


R = R + vel+force+bmod;

x = R(1:3:end);
y = R(2:3:end);
z = R(3:3:end);

end
