function out = div_scalar(x,z,T)

nx = numel(x);
nz = numel(z);

dx = x(2)-x(1);
dy = Inf;
dz = z(2)-z(1);


% Txx = reshape(T.xx,[nx,ny,nz]);
% Txy = reshape(T.xy,[nx,ny,nz]);
% Txz = reshape(T.xz,[nx,ny,nz]);
% Tyx = reshape(T.xy,[nx,ny,nz]);
% Tyy = reshape(T.yy,[nx,ny,nz]);
% Tyz = reshape(T.yz,[nx,ny,nz]);
% Tzx = reshape(T.xz,[nx,ny,nz]);
% Tzy = reshape(T.yz,[nx,ny,nz]);
% Tzz = reshape(T.zz,[nx,ny,nz]);


diff_order = 1;
if diff_order == 1
  dxT = [diff(T(1:2,:),1,1); T(2:end,:)-T(1:(end-1),:)];
  dyT = T*0;
  %dzT = [diff(T(:,1:2),1,2), diff(T,diff_order,2)];
  dzT = [diff(T(:,1:2),1,2), T(:,2:end)-T(:,1:(end-1))];
else diff_order == 2
  dxT = [2*diff(T(1:2,:),1,1); diff(T,diff_order,1); 2*diff(T((end-1):end,:),1,1)];
  dyT = T*0;
  dzT = [2*diff(T(:,1:2),1,2), diff(T,diff_order,2), 2*diff(T(:,(end-1):end),1,2)];
end

div_x = dxT/dx/diff_order;
div_y = dyT/dy/diff_order;
div_z = dzT/dz/diff_order;

out.x = div_x;
out.y = div_y;
out.z = div_z;