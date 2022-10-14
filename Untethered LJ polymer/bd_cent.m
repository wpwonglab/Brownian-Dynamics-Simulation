function [x, y ,z] = bd_cent(x, y, z,h)

x = x - mean(x);
y = y - mean(y);
z = z - mean(z)+h;

end