function [VX, VY, VZ] = bd_fluid_v2(~, ~, z, g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [VX, VY, VZ] = bd_fluid_v2(x, y, z, g)
% Hans Bergal, Darren Yang, and Wesley Wong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Shear flow
VX = g .*z; 
VY = zeros(size(z));
VZ = zeros(size(z));

end