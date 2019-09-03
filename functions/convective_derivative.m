function out = convective_derivative(x,z,T)

nx = numel(x);
nz = numel(z);

dx = x(2)-x(1);
dy = Inf;
dz = z(2)-z(1);

Tx = T.x;
Ty = T.y;
Tz = T.z;

% dxTx = zeros(nx,nz);
% dxTy = zeros(nx,nz);
% dxTz = zeros(nx,nz);
% dyTx = zeros(nx,nz);
% dyTy = zeros(nx,nz);
% dyTz = zeros(nx,nz);
% dzTx = zeros(nx,nz);
% dzTy = zeros(nx,nz);
% dzTz = zeros(nx,nz);
  
diff_order = 1;
dxTx = [1*diff(Tx(1:2,:),1,1); diff(Tx,diff_order,1)]/dx/diff_order;
dxTy = [1*diff(Ty(1:2,:),1,1); diff(Ty,diff_order,1)]/dx/diff_order;
dxTz = [1*diff(Tz(1:2,:),1,1); diff(Tz,diff_order,1)]/dx/diff_order;
dyTx = Tx*0;
dyTy = Tx*0;
dyTz = Tx*0;
dzTx = [1*diff(Tx(:,1:2),1,2), diff(Tx,diff_order,2)]/dz/diff_order;
dzTy = [1*diff(Ty(:,1:2),1,2), diff(Ty,diff_order,2)]/dz/diff_order;
dzTz = [1*diff(Tz(:,1:2),1,2), diff(Tz,diff_order,2)]/dz/diff_order;

% Wo ing
der_x = (Tx.*dxTx + Ty.*dyTx + Tz.*dzTx);
der_y = (Tx.*dxTy + Ty.*dyTy + Tz.*dzTy);
der_z = (Tx.*dxTz + Ty.*dyTz + Tz.*dzTz);
  
out.x = der_x;
out.y = der_y;
out.z = der_z;