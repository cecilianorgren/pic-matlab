function out = div_vector(x,z,T)

nx = numel(x);
nz = numel(z);

dx = x(2)-x(1);
dy = Inf;
dz = z(2)-z(1);

Tx = T.x;
Ty = T.y;
Tz = T.z;

if 1
  %%
  dxTx = zeros(nx,nz);
  dxTy = zeros(nx,nz);
  dxTz = zeros(nx,nz);
  
  diff_order = 1;
  dxTx = [1*diff(Tx(1:2,:),1,1); diff(Tx,diff_order,1)]/dx;  
  dyTy = Ty*0/dy;
  dzTz = [1*diff(Tz(:,1:2),1,2), diff(Tz,diff_order,2)]/dz;
end

out = dxTx + dyTy + dzTz;