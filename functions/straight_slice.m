function out = straight_slice(x_grid,z_grid,A,x_line,z_line)

method = 'linear';
[nx,nz] = size(A);


if all(size(x_grid) == [nx,nz]) && all(size(z_grid) == [nx,nz])
  Ai = interp2(x_grid,z_grid,A,tocolumn(z_line),torow(x_line),method);
else
  %[X_grid,Z_grid] = meshgrid(x_grid,z_grid);
  [X_grid,Z_grid] = ndgrid(x_grid,z_grid);  
  [X_inter,Z_inter] = ndgrid(x_grid(1):0.5:x_grid(end),z_grid(1):0.5:z_grid(end));  
  Ai = interp2(X_grid,X_grid,A',tocolumn(z_line),tocolumn(x_line),method);  
  Ai = interp2(X_grid,X_grid,A,X_inter,Z_inter,method);  
end

%[x_line,z_line] = meshgrid(x_line,z_line);


out = Ai;
