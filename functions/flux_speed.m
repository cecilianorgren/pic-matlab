function varargout = flux_speed(t,x,z,A,varargin)
% FLUX_SPEED Calculates (perpendicular) speed of flux.
%   [vAx,vAz] = FLUX_SPEED(t,x,z,A,options);
%   A(t,x,z) - (nt x nx x nz) matrix
%   options - 
if not(isempty(varargin)) && strcmp(varargin{1},'plot'); doPlot = 1; else doPlot = 0; end

%B = varargin{2};

% Initialize variables
matsize = size(A);
dAdt = zeros(matsize); % partial time derivative
Bx = zeros(matsize);
Bz = zeros(matsize);
vx = zeros(matsize);
vz = zeros(matsize);

% Grid size
dt = t(2)-t(1); % need to implement uneven timestep
dx = x(2)-x(1);
dz = z(2)-z(1);

% Calculate Bz
% Bx = -dAdz;
% Bz = dAdx;

Bx(:,:,1) = -(A(:,:,2)-A(:,:,1))/dz;
Bx(:,:,2:end-1) = -(A(:,:,3:end)-A(:,:,1:end-2))/(2*dz);
Bx(:,:,end) = -(A(:,:,end)-A(:,:,end-1))/dz;

Bz(:,1,:) = (A(:,2,:)-A(:,1,:))/dx;
Bz(:,2:end-1,:) =(A(:,3:end,:)-A(:,1:end-2,:))/(2*dx);
Bz(:,end,:) = (A(:,end,:)-A(:,end-1,:))/dx;

% partial time derivative
dAdt(1,:,:) = (A(2,:,:)-A(1,:,:))/dt;
dAdt(2:end-1,:,:) = (A(3:end,:,:)-A(1:end-2,:,:))/(2*dt);
dAdt(end,:,:) = (A(end,:,:)-A(end-1,:,:))/dt;

% Bx = interp(z(2:end)-0.5*dz,diff(A,3),z);
% Bz = -dAdz;

% Using total dAdt (laplacian) = 0, and vA.dot(B) = 0;
vx = -Bz.*dAdt./(Bx.^2+Bz.^2);
vz = Bx.*dAdt./(Bx.^2+Bz.^2);

if nargout == 1
  varargout = cat(vx,vz);  
elseif nargout == 2
  varargout{1} = vx;
  varargout{2} = vz;
end