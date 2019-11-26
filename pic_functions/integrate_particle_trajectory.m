% integrate_particle_trajectory
pic = df04;

% Particle position
x = 9.1;
y = 0;
z = 1.03;
r = [x,y,z];
t = 10.234; % wci-1

% Get closest field points
tic
nClosest = 2;
tmppic = pic.xlim(x,'closest',nClosest).zlim(z,'closest',nClosest).twcilim(t,'closest',nClosest);
tmpt =  tmppic.twci;
tmpx =  tmppic.xi;
tmpz =  tmppic.zi;
tmpEx = tmppic.Ex;
tmpEy = tmppic.Ey;
tmpEz = tmppic.Ez;
tmpBx = tmppic.Bx;
tmpBy = tmppic.By;
tmpBz = tmppic.Bz;
toc
tic
% Interpolate to particle position
[X,Z,T] = meshgrid(tmpx,tmpz,tmpt);
intEx = interp3(tmpx,tmpz,tmpt,tmpEx,x,z,t);
%[x,y,z] = integrate_particle_trajectory(pic,r0,v0)
toc
