function [ax,az] = flux_speed(t,x,z,vA,varargin)
% vA(t,x,z)

%if not(isempty(varargin)) && strcmp(varargin,'plot'); doPlot = 1; else doPlot = 0; end

vx = vA.x;
vz = vA.z;

% Grid size
dt = t(2)-t(1);
dx = x(2)-x(1);
dz = z(2)-z(1);

% Initialize variables
matsize = size(vz);
dVxdt = zeros(matsize); % total time derivative
dVzdt = zeros(matsize); % total time derivative

pdVxdt = zeros(matsize); % partial time derivative
pdVzdt = zeros(matsize); % partial time derivative
dVxdx = zeros(matsize); % spatial derivative
dVxdz = zeros(matsize); % spatial derivative
dVzdx = zeros(matsize); % spatial derivative
dVzdz = zeros(matsize); % spatial derivative

% spatial derivative
dVxdz(:,:,1) = -(vx(:,:,2)-vx(:,:,1))/dz;
dVxdz(:,:,2:end-1) = -(vx(:,:,3:end)-vx(:,:,1:end-2))/(2*dz);
dVxdz(:,:,end) = -(vx(:,:,end)-vx(:,:,end-1))/dz;
dVzdz(:,:,1) = -(vz(:,:,2)-vz(:,:,1))/dz;
dVzdz(:,:,2:end-1) = -(vz(:,:,3:end)-vz(:,:,1:end-2))/(2*dz);
dVzdz(:,:,end) = -(vz(:,:,end)-vz(:,:,end-1))/dz;

dVxdx(:,1,:) = (vx(:,2,:)-vx(:,1,:))/dx;
dVxdx(:,2:end-1,:) =(vx(:,3:end,:)-vx(:,1:end-2,:))/(2*dx);
dVxdx(:,end,:) = (vx(:,end,:)-vx(:,end-1,:))/dx;
dVzdx(:,1,:) = (vz(:,2,:)-vz(:,1,:))/dx;
dVzdx(:,2:end-1,:) =(vz(:,3:end,:)-vz(:,1:end-2,:))/(2*dx);
dVzdx(:,end,:) = (vz(:,end,:)-vz(:,end-1,:))/dx;

% partial time derivative
pdVxdt(1,:,:) = (vx(2,:,:)-vx(1,:,:))/dt;
pdVxdt(2:end-1,:,:) = (vx(3:end,:,:)-vx(1:end-2,:,:))/(2*dt);
pdVxdt(end,:,:) = (vx(end,:,:)-vx(end-1,:,:))/dt;
pdVzdt(1,:,:) = (vz(2,:,:)-vz(1,:,:))/dt;
pdVzdt(2:end-1,:,:) = (vz(3:end,:,:)-vz(1:end-2,:,:))/(2*dt);
pdVzdt(end,:,:) = (vz(end,:,:)-vz(end-1,:,:))/dt;

% Collect terms
ax = pdVxdt + vx.*dVxdx + vz.*dVxdz;
az = pdVzdt + vx.*dVzdx + vz.*dVzdz;

 