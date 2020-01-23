function out = curl_vector(x,z,T)

nx = numel(x);
nz = numel(z);

dx = x(2)-x(1);
dy = Inf;
dz = z(2)-z(1);

%Tx = T.x;
Ty = T;
%Tz = T.z;

if 1
  %%
  dxTy = zeros(nx,nz);
  dzTy = zeros(nx,nz);    
  diff_order = 1;
  
  dzTy = [1*diff(Ty(1:2,:),1,1); diff(Ty,diff_order,1)]/dz;
  dxTy = [1*diff(Ty(:,1:2),1,2), diff(Ty,diff_order,2)]/dz;  
end

out = -dzTy;