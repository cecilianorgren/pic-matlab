mass = 25;
wpewce = 2;

xmax = 200;
zmax = 100;
lz = 2*zmax;
lx = xmax;
a0 = 1.5*sqrt(mass);
a0 = 0.2;
xm = 12.8;
zm = 12.8;

x = linspace(-xmax,xmax,1000);
z = linspace(-zmax,zmax,1000);
[X,Z] = meshgrid(x,z);
%xx = x - 0.5*xmax;

jy = @(x,z) a0/xm*(pi*xm^2/lz + 1 - ((x - 0.5*xmax)/xm).^2).*exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5).*cos(pi*z/lz);
jy = @(x,z) 0 + exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5*0 - 0.5*(z).^2/xm^2);

%pcolor(X,Z,jy(X,Z)); shading flat; colorbar;

xm = 10;
zm = 10;
syms x y z
J = [0, -exp(-0.5*x^2/xm^2-0.5*z^2/zm^2), 0];
X = [x y z];
B = curl(J,X);

Bx = B(1);
By = B(2);
Bz = B(3);
Jx = J(1);
Jy = J(2);
Jz = J(3);

fB = symfun(B,[x y z]);
fBx = symfun(Bx,[x y z]);
fBy = symfun(By,[x y z]);
fBz = symfun(Bz,[x y z]);
mfBx = matlabFunction(fBx);
mfBy = matlabFunction(fBy);
mfBz = matlabFunction(fBz);
fJ = symfun(J,[x y z]);
fJx = symfun(Jx,[x y z]);
fJy = symfun(Jy,[x y z]);
fJz = symfun(Jz,[x y z]);
mfJx = matlabFunction(fJx);
mfJy = matlabFunction(fJy);
mfJz = matlabFunction(fJz);

xvec = linspace(-xmax,xmax,1000);
zvec = linspace(-zmax,zmax,1000);
[Xvec,Zvec] = meshgrid(xvec,zvec);


nrows = 2;
ncols = 2;
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);  
end

isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfJy(Xvec,0,Zvec));
shading(hca,'flat')

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfBz(Xvec,0,Zvec));
shading(hca,'flat')

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfBx(Xvec,0,Zvec));
shading(hca,'flat')

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfJy(Xvec,0,Zvec));
shading(hca,'flat')


%%
mass = 25;
wpewce = 2;

xmax = 200;
zmax = 100;
lz = 2*zmax;
lx = xmax;
a0 = 1.5*sqrt(mass);
a0 = 0.2;
xm = 5;
zm = 5;

x = linspace(-xmax,xmax,1000);
z = linspace(-zmax,zmax,1000);
[X,Z] = meshgrid(x,z);
%xx = x - 0.5*xmax;

jy = @(x,z) a0/xm*(pi*xm^2/lz + 1 - ((x - 0.5*xmax)/xm).^2).*exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5).*cos(pi*z/lz);
jy = @(x,z) 0 + exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5*0 - 0.5*(z).^2/xm^2);

%pcolor(X,Z,jy(X,Z)); shading flat; colorbar;

xm = 10;
zm = 10;
syms x y z
B = [-z/xm*exp(-0.5*x^2/xm^2-0.5*z^2/zm^2), 0, x/zm*exp(-0.5*x^2/xm^2-0.5*z^2/zm^2)];
X = [x y z];
J = curl(B,X);

Bx = B(1);
By = B(2);
Bz = B(3);
Jx = J(1);
Jy = J(2);
Jz = J(3);

fB = symfun(B,[x y z]);
fBx = symfun(Bx,[x y z]);
fBy = symfun(By,[x y z]);
fBz = symfun(Bz,[x y z]);
mfBx = matlabFunction(fBx);
mfBy = matlabFunction(fBy);
mfBz = matlabFunction(fBz);
fJ = symfun(J,[x y z]);
fJx = symfun(Jx,[x y z]);
fJy = symfun(Jy,[x y z]);
fJz = symfun(Jz,[x y z]);
mfJx = matlabFunction(fJx);
mfJy = matlabFunction(fJy);
mfJz = matlabFunction(fJz);

xvec = linspace(-xmax,xmax,1000);
zvec = linspace(-zmax,zmax,1000);
[Xvec,Zvec] = meshgrid(xvec,zvec);


nrows = 2;
ncols = 2;
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);  
end

isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfJy(Xvec,0,Zvec));
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfBz(Xvec,0,Zvec));
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfBx(Xvec,0,Zvec));
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfJy(Xvec,0,Zvec));
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)

colormap(cn.cmap('blue_red'))

%%
mass = 25;
wpewce = 2;

xmax = 200;
zmax = 100;
lz = 2*zmax;
lx = xmax;
a0 = 1.5*sqrt(mass);
a0 = 0.2;
xm = 12.8;
zm = 12.8;

x = linspace(-xmax,xmax,1000);
z = linspace(-zmax,zmax,1000);
[X,Z] = meshgrid(x,z);
%xx = x - 0.5*xmax;

jy = @(x,z) a0/xm*(pi*xm^2/lz + 1 - ((x - 0.5*xmax)/xm).^2).*exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5).*cos(pi*z/lz);
jy = @(x,z) 0 + exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5*0 - 0.5*(z).^2/xm^2);

%pcolor(X,Z,jy(X,Z)); shading flat; colorbar;

xm = 10;
zm = 10;
syms x y z
J = [0, -exp(-0.5*x^2/xm^2-0.5*z^2/zm^2), 0];
X = [x y z];
B = curl(J,X);

Bx = B(1);
By = B(2);
Bz = B(3);
Jx = J(1);
Jy = J(2);
Jz = J(3);

fB = symfun(B,[x y z]);
fBx = symfun(Bx,[x y z]);
fBy = symfun(By,[x y z]);
fBz = symfun(Bz,[x y z]);
mfBx = matlabFunction(fBx);
mfBy = matlabFunction(fBy);
mfBz = matlabFunction(fBz);
fJ = symfun(J,[x y z]);
fJx = symfun(Jx,[x y z]);
fJy = symfun(Jy,[x y z]);
fJz = symfun(Jz,[x y z]);
mfJx = matlabFunction(fJx);
mfJy = matlabFunction(fJy);
mfJz = matlabFunction(fJz);

xvec = linspace(-xmax,xmax,1000);
zvec = linspace(-zmax,zmax,1000);
[Xvec,Zvec] = meshgrid(xvec,zvec);


nrows = 2;
ncols = 2;
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);  
end

isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfJy(Xvec,0,Zvec));
shading(hca,'flat')

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfBz(Xvec,0,Zvec));
shading(hca,'flat')

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfBx(Xvec,0,Zvec));
shading(hca,'flat')

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfJy(Xvec,0,Zvec));
shading(hca,'flat')

%%
mass = 25;
wpewce = 2;

xmax = 200;
zmax = 100;
lz = 2*zmax;
lx = xmax;
a0 = 1.5*sqrt(mass);
a0 = 0.2;
xm = 5;
zm = 5;

x = linspace(-xmax,xmax,1000);
z = linspace(-zmax,zmax,1000);
[X,Z] = meshgrid(x,z);
%xx = x - 0.5*xmax;

jy = @(x,z) a0/xm*(pi*xm^2/lz + 1 - ((x - 0.5*xmax)/xm).^2).*exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5).*cos(pi*z/lz);
jy = @(x,z) 0 + exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5*0 - 0.5*(z).^2/xm^2);

%pcolor(X,Z,jy(X,Z)); shading flat; colorbar;

xm = 10;
zm = 10;
syms x y z
B = [xm*exp(-0.5*x^2/xm^2-0.5*z^2/zm^2), 0, x/zm*exp(-0.5*x^2/xm^2-0.5*z^2/zm^2)];
X = [x y z];
J = curl(B,X);

Bx = B(1);
By = B(2);
Bz = B(3);
Jx = J(1);
Jy = J(2);
Jz = J(3);

fB = symfun(B,[x y z]);
fBx = symfun(Bx,[x y z]);
fBy = symfun(By,[x y z]);
fBz = symfun(Bz,[x y z]);
mfBx = matlabFunction(fBx);
mfBy = matlabFunction(fBy);
mfBz = matlabFunction(fBz);
fJ = symfun(J,[x y z]);
fJx = symfun(Jx,[x y z]);
fJy = symfun(Jy,[x y z]);
fJz = symfun(Jz,[x y z]);
mfJx = matlabFunction(fJx);
mfJy = matlabFunction(fJy);
mfJz = matlabFunction(fJz);

xvec = linspace(-xmax,xmax,1000);
zvec = linspace(-zmax,zmax,1000);
[Xvec,Zvec] = meshgrid(xvec,zvec);


nrows = 2;
ncols = 2;
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);  
end

isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfJy(Xvec,0,Zvec));
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)
hca.Title.String = 'Jy';

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfBz(Xvec,0,Zvec));
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)
hca.Title.String = 'Bz';

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfBx(Xvec,0,Zvec));
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)
hca.Title.String = 'Jy';

hca = h(isub); isub = isub + 1;
pcolor(hca,Xvec,Zvec,mfJy(Xvec,0,Zvec));
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)
hca.Title.String = 'Jy';

colormap(cn.cmap('blue_red'))