function [ete, etev] = bd_endt_x(~, x,y,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ete] = bd_endt_x(np, x, y, z)
% Calculate End-to-End Distant squared (ete)
% by Hans Bergal, Darren Yang, and Wesley Wong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a=(x(end)-x(1)).^2;
b=(y(end)-y(1)).^2;
c=(z(end)-z(1)).^2;
d=a+b+c;

ete = d;
etev = [a, b, c];

end