function out = cross_product(ax,ay,az,bx,by,bz)
% Calculates cross product
% out = cross_product(ax,ay,az,bx,by,bz);

out_x_yz = ay.*bz;
out_x_zy = - az.*by;
out_y_xz = az.*bx;
out_y_zx = - ax.*bz;
out_z_xy = ax.*by;
out_z_yx = - ay.*bx;

out_x = ay.*bz - az.*by;
out_y = az.*bx - ax.*bz;
out_z = ax.*by - ay.*bx;

out.x = out_x;
out.y = out_y;
out.z = out_z;

out.x_yz = out_x_yz;
out.x_zy = out_x_zy;
out.y_xz = out_y_xz;
out.y_zx = out_y_zx;
out.z_xy = out_z_xy;
out.z_yx = out_z_yx;