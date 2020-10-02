function A = vector_potential(x,z,bx,bz)
% Calculates the xz-plane magnetic vector potential.
%   A = vector_potential(x,z,bx,bz)
%   B = rot(A) = xhat(dA/dz) - zhat(dA/dx)
%   Bx = dA/dz
%   Bz = -dA/dx

% Grid
dx = x(2)-x(1);
dz = z(2)-z(1);
nz = numel(z);
nx = numel(x);

% Check input, can be timeseries..
if ndims(bx) == 3 % time series, with [time,x,z]
  ntimes = size(bx,1);
  A = bx*0;
  for itime = 1:ntimes
    bx_tmp = squeeze(bx(itime,:,:));
    bz_tmp = squeeze(bz(itime,:,:));
    A_tmp = zeros(nx,nz);
    % Dont put zero right at the edge 
    ixm = 10;    
    % Advance up
    A_tmp(ixm,:) = cumsum(-bx_tmp(ixm,:),2)*dz;
    % Advance to the right
    A_tmp((ixm+1):end,:) = repmat(A_tmp(ixm,:),nx-ixm,1) + cumsum(bz_tmp((ixm+1):end,:),1)*dx;
    % Advance to the left
    A_tmp((ixm-1):-1:1,:) = repmat(A_tmp(ixm,:),ixm-1,1) + cumsum(bz_tmp((ixm-1):-1:1,:),1)*dx;
    A(itime,:,:) = A_tmp;
    1;
  end
else
  % Dont put zero right at the edge 
  ixm = 10;
  A = zeros(nx,nz);
  % Advance up
  A(ixm,:) = cumsum(-bx(ixm,:),2)*dz;
  % Advance to the right
  A((ixm+1):end,:) = repmat(A(ixm,:),nx-ixm,1) + cumsum(bz((ixm+1):end,:),1)*dx;
  % Advance to the left
  A((ixm-1):-1:1,:) = repmat(A(ixm,:),ixm-1,1) + cumsum(bz((ixm-1):-1:1,:),1)*dx;
end