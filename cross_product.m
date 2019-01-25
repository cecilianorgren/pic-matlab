function out = cross_product(ax,ay,az,bx,by,bz);
% Calculates cross product
% out = cross_product(ax,ay,az,bx,by,bz);

out_x = ay.*bz - az.*by;
out_y = az.*bz - ax.*bz;
out_z = ax.*by - ay.*bx;

out.x = out_x;
out.y = out_y;
out.z = out_z;