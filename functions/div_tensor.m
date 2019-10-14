function out = div_tensor(x,z,T,varargin)

returnComponents = 0;

nargs = numel(varargin);
if nargs == 2
  if strcmp(varargin{1},'comp')
    returnComponents = varargin{2};
  end
end

nx = numel(x);
nz = numel(z);

dx = x(2)-x(1);
dy = Inf;
dz = z(2)-z(1);

Txx = T.xx;
Txy = T.xy;
Txz = T.xz;
Tyx = T.xy;
Tyy = T.yy;
Tyz = T.yz;
Tzx = T.xz;
Tzy = T.yz;
Tzz = T.zz;

% dxTxx = zeros(nx,nz);
% dxTxy = zeros(nx,nz);
% dxTxz = zeros(nx,nz);
% dyTyx = zeros(nx,nz);
% dyTyy = zeros(nx,nz);
% dyTyz = zeros(nx,nz);
% dzTzx = zeros(nx,nz);
% dzTzy = zeros(nx,nz);
% dzTzz = zeros(nx,nz);

diff_order = 1;
dxTxx = [1*diff(Txx(1:2,:),1,1); diff(Txx,diff_order,1)]/dx/diff_order;
dxTxy = [1*diff(Txy(1:2,:),1,1); diff(Txy,diff_order,1)]/dx/diff_order;
dxTxz = [1*diff(Txz(1:2,:),1,1); diff(Txz,diff_order,1)]/dx/diff_order;
dyTyx = Tyx*0/dy/diff_order;
dyTyy = Tyy*0/dy/diff_order;
dyTyz = Tyz*0/dy/diff_order;
dzTzx = [1*diff(Tzx(:,1:2),1,2), diff(Tzx,diff_order,2)]/dz/diff_order;
dzTzy = [1*diff(Tzy(:,1:2),1,2), diff(Tzy,diff_order,2)]/dz/diff_order;
dzTzz = [1*diff(Tzz(:,1:2),1,2), diff(Tzz,diff_order,2)]/dz/diff_order;

div_x = (dxTxx + dyTyx + dzTzx);
div_y = (dxTxy + dyTyy + dzTzy);
div_z = (dxTxz + dyTyz + dzTzz);

out.x = div_x;
out.y = div_y;
out.z = div_z;
  
if returnComponents % is this the correct placements ???
  out.x_xx = dxTxx;
  out.x_yy = dyTyx;
  out.x_zz = dzTzx;
  out.y_xx = dxTxy;
  out.y_yy = dyTyy;
  out.y_zz = dzTzy;
  out.z_xx = dxTxz;
  out.z_yy = dyTyz;
  out.z_zz = dzTzz;
end