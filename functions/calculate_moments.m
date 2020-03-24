function varargout = calculate_moments(fin)
% CALCULATE_MOMENTS Calculate moments of distribution f.

v = fin.v;
f = fin.f;
dv = v(2)-v(1);
%dx = fin.x(2)-fin.x(1);
%dz = fin.z(2)-fin.z(1);
[VX,VY,VZ] = ndgrid(v,v,v);

% density: n = int(f)dv3
n = sum(f(:))*dv^3;

% flux/velocity: int(fv)dv3
FVX = f.*VX;
FVY = f.*VY;
FVZ = f.*VZ;

jx = sum(FVX(:))*dv^3;
jy = sum(FVY(:))*dv^3;
jz = sum(FVZ(:))*dv^3;

vx = jx/n;
vy = jy/n;
vz = jz/n;

vout.x = vx;
vout.y = vy;
vout.z = vz;

% pressure/temperature, only do scalar pressure
VX_cent = VX - vx;
VY_cent = VY - vy;
VZ_cent = VZ - vz; 
VABS_cent = sqrt(VX_cent.^2+VY_cent.^2+VZ_cent.^2);
%FVX_cent2 = f.*VX_cent2;
%FVY_cent2 = f.*VY_cent2;
%FVZ_cent2 = f.*VZ_cent2;

p = sum(sum(sum(f.*VABS_cent.^2)))*dv^3;
t = p/n;



varargout{1} = n;
varargout{2} = vout;
varargout{3} = p;
1;
