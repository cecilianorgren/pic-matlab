function  x_res = eom_pic_back(t,x_vect,pic,m,q)

x = x_vect(1);
y = x_vect(2);
z = x_vect(3);
vx = x_vect(4);
vy = x_vect(5);
vz = x_vect(6);

%disp(sprintf('t = %g, x = %g, y = %g, z = %g',t,x,y,z))

[Ex,Ey,Ez,Bx,By,Bz] = pic.interp_EB(x,z,t);

% Equations to be solved
x_res = zeros(6,1);
x_res(1) = -vx; % dx/dt = vx;
x_res(2) = -vy; % dy/dt = vy;
x_res(3) = -vz; % dz/dt = vz;
x_res(4) = -(q/m)*(Ex + vy*Bz - vz*By);
x_res(5) = -(q/m)*(Ey + vz*Bx - vx*Bz);
x_res(6) = -(q/m)*(Ez + vx*By - vy*Bx);                                              
                