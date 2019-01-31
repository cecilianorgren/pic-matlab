function out = grad_scalar(x,z,T)
% GRAD_SCALAR calculates the gradient of a scalar
%   grad_s = GRAD_SCALAR(x,z,s);
%
%   If data seems very noisy, consider smoothening it first, using smooth2.
nx = numel(x);
nz = numel(z);

dx = x(2)-x(1);
dy = Inf;
dz = z(2)-z(1);

diff_order = 1; % using diff_order = 2 seems to not work as well
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
out.abs = sqrt(div_x.^2 + div_y.^2 + div_z.^2);