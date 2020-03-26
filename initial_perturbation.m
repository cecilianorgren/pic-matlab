%
%
%

mass = 25; % not used
wpewce = 2; % not used
xmax = 200; % sixe of box
zmax = 25; % size of box
B0 = 1; % harris sheet asymptotic amplitude
BG = 0.5;
zh = 1; % harris sheet width
xp = 1*zh; % scale length of perturbation
zp = 1*zh; % scale length of perturbation
ah = B0*zh; % amplitude of harris sheet
ap = ah*0.25; % amplitude of perturbation
TeTi = 5; % not used
Ttot = 0.5; % not used

% If you want to get the expressions in terms of all parameters (zh,zp,xp,
% ah,ap, etc), you need to add these as syms here.
% syms x y z zh xp zp ah ap
% and
% add them as variables below, ex
% fA = symfun(A,[x y z]); 
% -->
% fA = symfun(A,[x y z zh xp zp ah ap]);
% Otherwise they are mulitplied together as numbers
syms x y z %xp zp ah ap

R = [x y z];

doLocal = 0;
% Harris current sheet
AH = [0; -ah*log(cosh(z/zh)); 0]; 
% Perturbation
AP = [0; -ap*exp(-x^2/2/xp^2 -z^2/2/zp^2 + 0.0); 0]; % circular perturbation
AG = [BG*z; 0; 0]; % guide field
%AP = [0; ah*log(cosh(z/zh))*exp(-x^2/2/xp^2 -z^2/2/zp^2); 0]; 
%AP = [0 a0*sin(0.5*pi*(1+x/xm))*sin(0.5*pi*(1+z/zm)) 0];
%AP = [0 ap*exp(-x^2/2/xp^2 + 0.5)*cos(pi*z/lz) 0];
%AP = [0 a0*sech(pi*x/xm)*sech(pi*z/zm) 0];
%AP = [0 0+a0*cos(0.5*pi*(x/xm))*cos(0.5*pi*(z/zm)) 0];
  
A = AH + AP + AG;
%A = AP;
B = curl(A,R);
J = curl(B,R);

Ax = A(1);
Ay = A(2);
Az = A(3);
Bx = B(1);
By = B(2);
Bz = B(3);
Jx = J(1);
Jy = J(2);
Jz = J(3);

nx = 1000;
nz = 500;
xvec = linspace(-xmax,xmax,nx);
zvec = linspace(-zmax,zmax,nz);
[X,Z] = meshgrid(xvec,zvec);

fA = symfun(A,[x y z]);
fAx = symfun(Ax,[x y z]);
fAy = symfun(Ay,[x y z]);
fAz = symfun(Az,[x y z]);
mfAx = matlabFunction(fAx);
mfAy = matlabFunction(fAy);
mfAz = matlabFunction(fAz);
matAx = mfAx(X,0,Z); if isscalar(matAx); matAx = repmat(matAx,nz,nx); end
matAy = mfAy(X,0,Z); if isscalar(matAy); matAy = repmat(matAy,nz,nx); end
matAz = mfAz(X,0,Z); if isscalar(matAz); matAz = repmat(matAz,nz,nx); end

fB = symfun(B,[x y z]);
fBx = symfun(Bx,[x y z]);
fBy = symfun(By,[x y z]);
fBz = symfun(Bz,[x y z]);
mfBx = matlabFunction(fBx);
mfBy = matlabFunction(fBy);
mfBz = matlabFunction(fBz);
matBx = mfBx(X,0,Z); if isscalar(matBx); matBx = repmat(matBx,nz,nx); end
matBy = mfBy(X,0,Z); if isscalar(matBy); matBy = repmat(matBy,nz,nx); end
matBz = mfBz(X,0,Z); if isscalar(matBz); matBz = repmat(matBz,nz,nx); end

fJ = symfun(J,[x y z]);
fJx = symfun(Jx,[x y z]);
fJy = symfun(Jy,[x y z]);
fJz = symfun(Jz,[x y z]);
mfJx = matlabFunction(fJx);
mfJy = matlabFunction(fJy);
mfJz = matlabFunction(fJz);
matJx = mfJx(X,0,Z); if isscalar(matJx); matJx = repmat(matJx,nz,nx); end
matJy = mfJy(X,0,Z); if isscalar(matJy); matJy = repmat(matJy,nz,nx); end
matJz = mfJz(X,0,Z); if isscalar(matJz); matJz = repmat(matJz,nz,nx); end

if doLocal
  matAx(any([abs(X(:))>xp abs(Z(:))>xp],2)) = 0;
  matAy(any([abs(X(:))>xp abs(Z(:))>xp],2)) = 0;
  matAz(any([abs(X(:))>xp abs(Z(:))>xp],2)) = 0;
  matBx(any([abs(X(:))>xp abs(Z(:))>xp],2)) = 0;
  matBy(any([abs(X(:))>xp abs(Z(:))>xp],2)) = 0;
  matBz(any([abs(X(:))>xp abs(Z(:))>xp],2)) = 0;
  matJx(any([abs(X(:))>xp abs(Z(:))>xp],2)) = 0;
  matJy(any([abs(X(:))>xp abs(Z(:))>xp],2)) = 0;
  matJz(any([abs(X(:))>xp abs(Z(:))>xp],2)) = 0;
end

% Particle velocities, NOT IMPLEMENTED
velocities = 'diamagnetic only';
switch velocities
  case 'diamagnetic only' % no electric field
    
  case 'stationary ions' % finite electric field
    
end

% Plot
nrows = 3;
ncols = 2;
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);  
end
isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,matAy);
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)
hca.Title.String = 'Ay';

hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,matBx);
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)
hca.Title.String = 'Bx';

hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,matBy);
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)
hca.Title.String = 'By';

hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,matBz);
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)
hca.Title.String = 'Bz';

hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,matJx);
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)
hca.Title.String = 'Jx';

hca = h(isub); isub = isub + 1;
pcolor(hca,X,Z,matJy);
shading(hca,'flat')
hca.CLim = max(abs(hca.CLim))*[-1 1];
colorbar('peer',hca)
hca.Title.String = 'Jy';

if 0
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matJz);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Jz';
end
colormap(cn.cmap('blue_red'))

hlinks = linkprop(h,{'XLim','YLim'});
hlinks.Targets(1).XLim = [-10 10];
hlinks.Targets(1).YLim = [-10 10];

%% Old
if 0
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
[R,Z] = meshgrid(x,z);
%xx = x - 0.5*xmax;

jy = @(x,z) a0/xm*(pi*xm^2/lz + 1 - ((x - 0.5*xmax)/xm).^2).*exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5).*cos(pi*z/lz);
jy = @(x,z) 0 + exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5*0 - 0.5*(z).^2/xm^2);

%pcolor(X,Z,jy(X,Z)); shading flat; colorbar;

xm = 10;
zm = 10;
syms x y z
J = [0, -exp(-0.5*x^2/xm^2-0.5*z^2/zm^2), 0];
R = [x y z];
B = curl(J,R);

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
[R,Z] = meshgrid(x,z);
%xx = x - 0.5*xmax;

jy = @(x,z) a0/xm*(pi*xm^2/lz + 1 - ((x - 0.5*xmax)/xm).^2).*exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5).*cos(pi*z/lz);
jy = @(x,z) 0 + exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5*0 - 0.5*(z).^2/xm^2);

%pcolor(X,Z,jy(X,Z)); shading flat; colorbar;

xm = 10;
zm = 10;
syms x y z
B = [-z/xm*exp(-0.5*x^2/xm^2-0.5*z^2/zm^2), 0, x/zm*exp(-0.5*x^2/xm^2-0.5*z^2/zm^2)];
R = [x y z];
J = curl(B,R);

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
[R,Z] = meshgrid(x,z);
%xx = x - 0.5*xmax;

jy = @(x,z) a0/xm*(pi*xm^2/lz + 1 - ((x - 0.5*xmax)/xm).^2).*exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5).*cos(pi*z/lz);
jy = @(x,z) 0 + exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5*0 - 0.5*(z).^2/xm^2);

%pcolor(X,Z,jy(X,Z)); shading flat; colorbar;

xm = 10;
zm = 10;
syms x y z
J = [0, -exp(-0.5*x^2/xm^2-0.5*z^2/zm^2), 0];
R = [x y z];
B = curl(J,R);

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
[R,Z] = meshgrid(x,z);
%xx = x - 0.5*xmax;

jy = @(x,z) a0/xm*(pi*xm^2/lz + 1 - ((x - 0.5*xmax)/xm).^2).*exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5).*cos(pi*z/lz);
jy = @(x,z) 0 + exp(-0.5*(x - 0.5*xmax).^2/xm^2 + 0.5*0 - 0.5*(z).^2/xm^2);

%pcolor(X,Z,jy(X,Z)); shading flat; colorbar;

xm = 10;
zm = 10;
syms x y z
B = [xm*exp(-0.5*x^2/xm^2-0.5*z^2/zm^2), 0, x/zm*exp(-0.5*x^2/xm^2-0.5*z^2/zm^2)];
R = [x y z];
J = curl(B,R);

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

%%
end