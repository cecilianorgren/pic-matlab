function [new_x,new_y,new_z] = rotate_xyz(old_x,old_y,old_z,rx,ry,rz)

size_data = size(old_x);

new_x = old_x*rx(1) + old_y*rx(2) + old_z*rx(3);
new_y = old_x*ry(1) + old_y*ry(2) + old_z*ry(3);
new_z = old_x*rz(1) + old_y*rz(2) + old_z*rz(3);