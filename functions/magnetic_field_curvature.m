function bcurv = magnetic_field_curvature(x,z,bx,by,bz)
% bcurv = magnetic_field_curvature(x,z,bx,by,bz);
% bcurv = dot(B,nabla) B
% bcurv_x = (Bxdx + Bydy + Bzdz) Bx
% bcurv_y = (Bxdx + Bydy + Bzdz) By
% bcurv_z = (Bxdx + Bydy + Bzdz) Bz

% normalize b
babs = sqrt(bx.^2 + by.^2 + bz.^2);
bx = bx./babs;
by = by./babs;
bz = bz./babs;

dx = x(2)-x(1);
dz = z(2)-z(1);
dy = Inf;

dxbx = diff(bx,1,1); dxbx(end+1,:) = dxbx(end,:);
dybx = 0;
dzbx = diff(bx,1,2); dzbx(:,end+1) = dzbx(:,end);
dxby = diff(by,1,1); dxby(end+1,:) = dxby(end,:);
dyby = 0;
dzby = diff(by,1,2); dzby(:,end+1) = dzby(:,end);
dxbz = diff(bz,1,1); dxbz(end+1,:) = dxbz(end,:);
dybz = 0;
dzbz = diff(bz,1,2); dzbz(:,end+1) = dzbz(:,end);

bcurv_x = bx.*dxbx/dx + by.*dybx/dy + bz.*dzbx/dz;
bcurv_y = bx.*dxby/dx + by.*dyby/dy + bz.*dzby/dz;
bcurv_z = bx.*dxbz/dx + by.*dybz/dy + bz.*dzbz/dz;

bcurv.units = 'wpi/c';
bcurv.x = bcurv_x;
bcurv.y = bcurv_y;
bcurv.z = bcurv_z;
bcurv.abs = sqrt(bcurv_x.^2 + bcurv_y.^2 + bcurv_z.^2);