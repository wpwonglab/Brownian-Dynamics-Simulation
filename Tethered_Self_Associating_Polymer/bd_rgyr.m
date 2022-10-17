function [Rgsq] = bd_rgyr(np, x, y, z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rgsq] = bd_rgyr(np, x, y, z)
% Determine radius of gyration squared for the polymer
% Hans Bergal,Darren Yang, and Wesley Wong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
avg_x = mean(x);
avg_y = mean(y);
avg_z = mean(z);

mysum = sum((x - avg_x).^2 + (y - avg_y).^2 + (z - avg_z).^2);

Rgsq = mysum/np;

end



