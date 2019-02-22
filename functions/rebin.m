function fe = rebin(f,old_grid,new_grid)

method = 'cart_3D_to_radial_1D';
% cart_3D_to_radial_1D
% bins can 

switch method
  case 'cart_3D_to_radial_1D'
    old_grid_x = old_grid;
    old_grid_y = old_grid;
    old_grid_z = old_grid;
    old_grid_abs = sqrt(old_grid_x.^2 + old_grid_y.^2 + old_grid_z.^2);
    [N,X] = hist(old_grid_abs,new_grid);
end
