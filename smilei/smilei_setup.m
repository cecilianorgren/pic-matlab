n0 = 1;
nc = 0.2;
L = 1;
z0 = 5;
n1 = @(z,L) n0*cosh((z-z0)/L).^(-2);
n2 = @(z,L) nc*0.5*(1+tanh((abs(z-z0)-2*L)/(0.5*L)));
n3 = @(z,L) nc*0.5*(1+tanh(((z-z0)-2*L)/(0.5*L)));
n4 = @(z,L) nc*0.5*(1+tanh((-(z-z0)-2*L)/(0.5*L)));
%n3 = @(z,L) nc*tanh(z/(3*L)).^2;

z = linspace(0,10,100);

plot(z,n1(z,L),z,n3(z,L),z,n4(z,L))

%%
n0 = 1;
nc = 0.1;
L = 1;
Lf = 1;
z0 = 40;
n1_ = @(z,L) n0*cosh((z-z0)/L).^(-2);
n2_ = @(z,L) nc*0.5*(1+tanh((abs(z-z0)-2*L)/(0.5*L)));
n3_ = @(z,L) nc*0.5*(1+tanh(((z-z0)-2*L)/(0.5*L))).*(3+cos((2*pi/(Lf))*(z-z0)))/4;
%            nb*0.5*(1+tanh(((y-y0)-2*L)/(0.5*L)))*(3+cos((2*pi/(Lf))*(y-y0)))/4
n4_ = @(z,L) 2*nc*0.5*(1+tanh((-(z-z0)-2*L)/(0.5*L))).*(3+cos((2*pi/(Lf))*(z-z0)))/4;

n5_ = @(z,L) cos((2*pi/2/L)*(z-z0));
%n3 = @(z,L) nc*tanh(z/(3*L)).^2;

z = linspace(0,80,2000);

plot(z,n1_(z,L),z,n3_(z,L),z,n4_(z,L))

%% 3D sinusoidal grid
mime = 100;
resx = 4; dx = 1/resx;
resy = 4; dy = 1/resy;
resz = 4; dz = 1/resz;
box_size = [4*sqrt(mime) 5*sqrt(mime) 6*sqrt(mime)];
nx = box_size(1)*resx + 1;
ny = box_size(2)*resy + 1;
nz = box_size(3)*resz + 1;
ncell = nx*ny*nz;
x = 0:dx:box_size(1);
y = 0:dy:box_size(2);
z = 0:dz:box_size(3);
[X,Y,Z] = meshgrid(x,y,z);


nwx = 2; % number of wavelengths
nwy = 2; % number of wavelengths
nwz = 2; % number of wavelengths
lx = box_size(1)/nwx;
ly = box_size(2)/nwy;
lz = box_size(3)/nwz;

fBx = @(x,y,z) sin(2*pi*x/lx);
fBy = @(x,y,z) sin(2*pi*y/ly);
fBz = @(x,y,z) sin(2*pi*z/lz);

Bmag = sqrt(fBx(X,Y,Z).^2 + fBy(X,Y,Z).^2 + fBz(X,Y,Z).^2);
%%
h = setup_subplots(1,1);
isub = 1;
% 3D density plot
hca = h(isub); isub = isub + 1;

isovalue = 0.6;
surf1 = isosurface(X,Y,Z,Bmag,isovalue);
p1 = patch(hca,surf1);
%isonormals(x,y,z,normalized_Free_Energy_map,p1);
set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
daspect([1,1,1])
view(hca,3); axis tight
camlight; lighting gouraud
isonormals(X,Y,Z,Bmag,p1);
set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
daspect([1,1,1])
view(hca,3); axis tight
camlight; lighting gouraud

hold(hca,'off')
%%
h = setup_subplots(2,2);
isub = 1;
hca = h(isub); isub = isub + 1;
plot(hca,x,fBx(x,0,0))
hca = h(isub); isub = isub + 1;
plot(hca,y,fBy(0,y,0))
hca = h(isub); isub = isub + 1;
plot(hca,z,fBz(0,0,z))

%%
% 3D density plot
hca = h(isub); isub = isub + 1;

isovalue = 1;
surf1 = isosurface(x,y,z,Bmag,isovalue);
p1 = patch(surf1);
%isonormals(x,y,z,normalized_Free_Energy_map,p1);
set(p1,'FaceColor','red','EdgeColor','none','FaceAlpha',0.1); % set the color, mesh and transparency level of the surface
daspect([1,1,1])
view(3); axis tight
camlight; lighting gouraud





