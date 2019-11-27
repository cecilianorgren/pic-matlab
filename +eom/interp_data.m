function  x_res = general(t,x_vect,x_,y_,z_,Bx_,By_,Bz_,Ex_,Ey_,Ez_)
% x_res = eom_general(t,x_vect,Bx_,By_,Bz_,Ex_,Ey_,Ez_)

% physical constants
e = 1.6022e-19;
me = 9.10939999e-31;

x = x_vect(1);
y = x_vect(2);
z = x_vect(3);
vx = x_vect(4);
vy = x_vect(5);
vz = x_vect(6);

%Vq = interp1(X,V,Xq)
% Bx = Bx_(x,y,z);
% By = By_(x,y,z);
% Bz = Bz_(x,y,z);
% Ex = Ex_(x,y,z);
% Ey = Ey_(x,y,z);
% Ez = Ez_(x,y,z);
ind = find(abs(z_*1e3-z) == min(abs(z_*1e3-z)));
BE = [Bx_,By_,Bz_,Ex_,Ey_,Ez_];
BE = BE(ind,:);
%BE = interp1(z_*1e3,[Bx_,By_,Bz_,Ex_,Ey_,Ez_],z,'nearest')*1e-9;
% Bx = interp1(z_*1e3,Bx_,z)*1e-9;
% By = interp1(z_*1e3,By_,z)*1e-9;
% Bz = interp1(z_*1e3,Bz_,z)*1e-9;
% Ex = interp1(z_*1e3,Ex_,z)*1e-3;
% Ey = interp1(z_*1e3,Ey_,z)*1e-3;
% Ez = interp1(z_*1e3,Ez_,z)*1e-3;
Bx = BE(:,1)*1e-9;
By = BE(:,2)*1e-9;
Bz = BE(:,3)*1e-9;
Ex = BE(:,4)*1e-3;
Ey = BE(:,5)*1e-3;
Ez = BE(:,6)*1e-3;

x_res = zeros(6,1);
x_res(1) = vx; % dx/dt = vx;
x_res(2) = vy; % dy/dt = vy;
x_res(3) = vz; % dz/dt = vz;
x_res(4) = (-e/me)*(Ex + vy*Bz - vz*By);
x_res(5) = (-e/me)*(Ey + vz*Bx - vx*Bz);
x_res(6) = (-e/me)*(Ez + vx*By - vy*Bx);                                              
                