function [mx, my, mz] = bd_size(x, y, z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [mx, my, mz] = bd_maxi(x, y, z)
% Determine the box size enclose the polymer
% Darren Yang (ydarren@gmail.com) with Wesley Wong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mx = max(x) - min(x);
my = max(y) - min(y);
mz = max(z) - min(z);
end