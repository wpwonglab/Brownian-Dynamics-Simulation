function [np, x, y, z]= bd_read(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[np, x, y, z]= bd_conf_read(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in all the position into a matrix
position_matrix = dlmread(filename);

% Measure the number of particle (np), and save postion (xi, yi, zi)
[np, ~] = size(position_matrix);
x = position_matrix(:,1);
y = position_matrix(:,2);
z = position_matrix(:,3);

end