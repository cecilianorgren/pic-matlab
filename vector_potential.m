function A = vector_potential(x,z,bx,bz)
% A Calculates the xz-plane magnetic vector potential.
%   A = A(x,z,bx,bz)

% Grid
dx = x(2)-x(1);
dz = z(2)-z(1);
[X,Z] = meshgrid(x,z);

nnz = numel(z);
nnx = numel(x);
% Vector potential, for plotting xz-plane magnetic field lines
%Ay1 = cumsum(bz,1)*dx; - cumsum(bx,2)*dz;

% Pauls, looks better
ixm = 10;
Ay2(ixm,1) = 0;
for j=2:nnz 
  Ay2(ixm,j) = Ay2(ixm,j-1) + dz*bx(ixm,j-1);
end
for ix=ixm+1:nnx
  Ay2(ix,:) = Ay2(ix-1,:) - bz(ix-1,:)*dx; % ;     advance to the right
end
for ix=ixm-1:-1:1
  Ay2(ix,:) = Ay2(ix+1,:)+ bz(ix,:)*dx; % ;     advance to the left
end

A = Ay2;


