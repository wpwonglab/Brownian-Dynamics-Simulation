function bd_write(c, b, e, g,np, x, y, z,kT)

% Write to Trajectory
fprintf(c,'%d \n', np);
fprintf(c,'Unit: (nm) \n');
for j = 1:np
    fprintf(c, 'D %+6.4f %+6.4f %+6.4f \n', x(j), y(j), z(j));
end

% Calculate

[ ete, etev] = bd_endt_x(np, x, y, z);
[   Rgsq   ] = bd_rgyr(np, x, y, z);
[mx, my, mz] = bd_size(x, y, z);

% Write outputs
fprintf(b,'%+6.4f %+6.4f %+6.4f \n', mx, my, mz);
fprintf(e,'%+6.4f %+6.4f %+6.4f %+6.4f \n',etev(1),etev(2),etev(3),ete);
fprintf(g,'%+6.4f \n',Rgsq);


end