function out = div_tensor(x,z,T)

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

Txx = T.xx;
Txy = T.xy;
Txz = T.xz;
Tyx = T.xy;
Tyy = T.yy;
Tyz = T.yz;
Tzx = T.xz;
Tzy = T.yz;
Tzz = T.zz;

diff_order = 2;
dxTxx = [2*diff(Txx(1:2,:),1,1); diff(Txx,diff_order,1); 2*diff(Txx((end-1):end,:),1,1)];
dxTxy = [2*diff(Txy(1:2,:),1,1); diff(Txy,diff_order,1); 2*diff(Txy((end-1):end,:),1,1)];
dxTxz = [2*diff(Txz(1:2,:),1,1); diff(Txz,diff_order,1); 2*diff(Txz((end-1):end,:),1,1)];
dyTyx = T.xx*0;
dyTyy = T.xx*0;
dyTyz = T.xx*0;
dzTzx = [2*diff(Tzx(:,1:2),1,2), diff(Tzx,diff_order,2), 2*diff(Tzx(:,(end-1):end),1,2)];
dzTzy = [2*diff(Tzy(:,1:2),1,2), diff(Tzy,diff_order,2), 2*diff(Tzy(:,(end-1):end),1,2)];
dzTzz = [2*diff(Tzz(:,1:2),1,2), diff(Tzz,diff_order,2), 2*diff(Tzz(:,(end-1):end),1,2)];

div_x = (dxTxx/dx + dyTyx/dy + dzTzx/dz)/diff_order;
div_y = (dxTxy/dx + dyTyy/dy + dzTzy/dz)/diff_order;
div_z = (dxTxz/dx + dyTyz/dy + dzTzz/dz)/diff_order;

out.x = div_x;
out.y = div_y;
out.z = div_z;