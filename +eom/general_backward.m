function  x_res = general_backward(t,x_vect,Bx_,By_,Bz_,Ex_,Ey_,Ez_)
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

Bx = Bx_(x,y,z);
By = By_(x,y,z);
Bz = Bz_(x,y,z);
Ex = Ex_(x,y,z);
Ey = Ey_(x,y,z);
Ez = Ez_(x,y,z);

x_res = zeros(6,1);
x_res(1) = -vx; % dx/dt = vx;
x_res(2) = -vy; % dy/dt = vy;
x_res(3) = -vz; % dz/dt = vz;
x_res(4) = -(-e/me)*(Ex + vy*Bz - vz*By);
x_res(5) = -(-e/me)*(Ey + vz*Bx - vx*Bz);
x_res(6) = -(-e/me)*(Ez + vx*By - vy*Bx);                                              
                