%
nSimp = 50;
wpewce = 2;
mime = 100;
% de, sixe of box
xemin = 0; 
xemax = 0.25*2048; 
zemin = -128;
zemax = 128;
%zemax = 56;
%zemin = 0;
% di/de = (c/wpi)/(c/wpe) = wpe/wpi = sqrt(mi/me) 
ximin = xemin/sqrt(mime);
ximax = xemax/sqrt(mime);
zimin = zemin/sqrt(mime);
zimax = zemax/sqrt(mime);

xmin = ximin;
xmax = ximax;
zmin = zimin;
zmax = zimax;

Lhe = 20; % Harris sheet width in de
Lhi = Lhe/sqrt(mime); % Harris sheet width in di
l = Lhi; % harris sheet width
zh = l;
B0 = 1; % harris sheet asymptotic amplitude
BG = 0*0.25;
ah = B0*l; % rot(A) = B -> A0/l = B0 ....... rot(B) = J -> B0/l = Jy0

% no_hot_ng_n02, very violent, but quick
xp = 3*l; % scale length of perturbation
zp = 1.5*l; % scale length of perturbation
ap = ah*1; % amplitude of perturbation


% no_hot_bg_n02_m100, calm start, also slow
xp = 2*l; % scale length of perturbation
zp = 1*l; % scale length of perturbation
ap = 0.5*ah; % amplitude of perturbation

% no_hot_bg_test, calm start, but kind of slow
xp = 1*l; % scale length of perturbation
zp = 1*l; % scale length of perturbation
ap = 0.25*ah; % amplitude of perturbation

%
xp = 1.5*l; % scale length of perturbation
zp = 1.5*l; % scale length of perturbation
ap = 0.6*ah; % amplitude of perturbation

TeTi = 1/5; % not used, implement to get densities
Ttot = 0.5; % not used, implement to get densities

nx = 640;
nz = 160;
x0 = 0.5*(xmax-xmin);
xvec = linspace(xmin,xmax,nx)-x0;
zvec = linspace(zmin,zmax,nz);
[X,Z] = meshgrid(xvec,zvec);

% If you want to get the expressions in terms of all parameters (zh,zp,xp,
% ah,ap, etc), you need to add these as syms here.
% syms x y z zh xp zp ah ap
% and
% add them as variables below, ex
% fA = symfun(A,[x y z]); 
% -->
% fA = symfun(A,[x y z zh xp zp ah ap]);
% Otherwise they are mulitplied together as numbers
syms x y z% xp zp ah ap

R = [x y z];

z0_BG = 5;
zL_BG = 1;

doLocal = 0;
% Harris current sheet
AH = [0; -ah*log(cosh(z/zh)); 0]; 
% Perturbation
AP = [0; -ap*exp(-x^2/2/xp^2 -z^2/2/zp^2 + 0.0); 0]; % circular perturbation
AP = [0; -ap*exp(-x^2/2/xp^2 -z^2/2/zp^2 + 0.0); 0]; % circular perturbation
%AP = [0; -ap*cos(2*pi*(x)/xmax)*cos(pi*(z)/zmax)*exp(-x^2/2/xp^2 -z^2/2/zp^2 + 0.0); 0]; % circular perturbation
AG = [BG*z + 1*BG*log(cosh((z-z0_BG)/zL_BG)); 0; 0]; % guide field

%AP = [0; ah*log(cosh(z/zh))*exp(-x^2/2/xp^2 -z^2/2/zp^2); 0]; 
APsin = [0; -0.2*ah*cos(2*pi*x/(xmax-xmin)/1)*cos(1*pi*abs(z)/(zmax-zmin)); 0];
%AP = [0 ap*exp(-x^2/2/xp^2 + 0.5)*cos(pi*z/lz) 0];
%AP = [0 a0*sech(pi*x/xm)*sech(pi*z/zm) 0];
%AP = [0 0+a0*cos(0.5*pi*(x/xm))*cos(0.5*pi*(z/zm)) 0];
  
A = AH + 1*AP + 0*APsin;% + AP + AG + 1*APsin;
%A = AP;
B = curl(A,R);
J = curl(B,R);
divB = divergence(B,R);
PB = 0.5*(B(1).^2+B(2).^2+B(3).^2);
% Given T, n should be calculated such that PB+PT=const
% 
PT = 0.5*B0^2-PB;
% PT = PTe + PTi = PTi/(Ti/Te) + PTi = PTi*(Te/Ti + 1)
%                = PTe + PTe/(Te/Ti) = PTe*(1 + Ti/Te)
PTi = PT/(TeTi+1);
PTe = PT/(1+1/TeTi);
n = PT/Ttot;


% How to add a polarized inflow current sheet, such that the electrons are
% drifting, and not the ions?
% Particle velocities, NOT IMPLEMENTED
velocities = 'momentum conservation';
switch velocities
  case 'diamagnetic only' % no electric field
    
  case 'stationary ions' % finite electric field
    
  case 'momentum conservation'
    % mi*ni*vi + me*ne*ve = 0
    % ni*vi = - (me/mi)*ne*ve
    % J = ni*vi-ne*ve 
    %   = ni*vi + (mi/me)*ni*vi = (1 + mi/me)*ni*vi
    %   = - (me/mi)*ne*ve - ne*ve = - (me/mi + 1)*ne*ve
    % vi = J/[ni*(1 + mi/me)]
    % ve = -J/[ne*(me/mi + 1)]
    Vi = J./(n*(mime+1));
    Ve = -J./(n*(1/mime+1));
    VixB = simplify(cross(Vi,B),nSimp);
    VexB = cross(Ve,B);
    % E+vixB - div(Pti)/ne = 0
    GradPi = gradient(PTi,R);
    GradPe = gradient(PTe,R);    
    GradPin = GradPi/n;
    GradPen = GradPe/n;
    E = -VixB + GradPi./n;
    % E becomes a tanh function, which is not what I want
end

Ax = A(1);
Ay = A(2);
Az = A(3);
Bx = B(1);
By = B(2);
Bz = B(3);
Jx = J(1);
Jy = J(2);
Jz = J(3);
Viy = Vi(2);
Vey = Ve(2);
VixBx = VixB(1);
VixBy = VixB(2);
VixBz = VixB(3);
VexBx = VixB(1);
VexBy = VixB(2);
VexBz = VixB(3);
Ex = E(1);
Ey = E(2);
Ez = E(3);
GradPix = GradPi(1);
GradPiy = GradPi(2);
GradPiz = GradPi(3);

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

fVi = symfun(Vi,[x y z]);
fViy = symfun(Viy,[x y z]);
mfViy = matlabFunction(fViy);
matViy = mfViy(X,0,Z); if isscalar(matViy); matViy = repmat(matViy,nz,nx); end
fVe = symfun(Ve,[x y z]);
fVey = symfun(Vey,[x y z]);
mfVey = matlabFunction(fVey);
matVey = mfVey(X,0,Z); if isscalar(matVey); matVey = repmat(matVey,nz,nx); end

fVixB = symfun(VixB,[x y z]);
fVixBx = symfun(VixBx,[x y z]);
fVixBy = symfun(VixBy,[x y z]);
fVixBz = symfun(VixBz,[x y z]);
mfVixBx = matlabFunction(fVixBx);
mfVixBy = matlabFunction(fVixBy);
mfVixBz = matlabFunction(fVixBz);
matVixBx = mfVixBx(X,0,Z); if isscalar(matVixBx); matVixBx = repmat(matVixBx,nz,nx); end
matVixBy = mfVixBy(X,0,Z); if isscalar(matVixBy); matVixBy = repmat(matVixBy,nz,nx); end
matVixBz = mfVixBz(X,0,Z); if isscalar(matVixBz); matVixBz = repmat(matVixBz,nz,nx); end

fVexB = symfun(VexB,[x y z]);
fVexBx = symfun(VexBx,[x y z]);
fVexBy = symfun(VexBy,[x y z]);
fVexBz = symfun(VexBz,[x y z]);
mfVexBx = matlabFunction(fVexBx);
mfVexBy = matlabFunction(fVexBy);
mfVexBz = matlabFunction(fVexBz);
matVexBx = mfVexBx(X,0,Z); if isscalar(matVexBx); matVexBx = repmat(matVexBx,nz,nx); end
matVexBy = mfVexBy(X,0,Z); if isscalar(matVexBy); matVexBy = repmat(matVexBy,nz,nx); end
matVexBz = mfVexBz(X,0,Z); if isscalar(matVexBz); matVexBz = repmat(matVexBz,nz,nx); end

fN = symfun(n,[x y z]);
mfN = matlabFunction(fN);
matN = mfN(X,0,Z); if isscalar(matN); matN = repmat(matN,nz,nx); end

fPi = symfun(PTi,[x y z]);
mfPi = matlabFunction(fPi);
matPi = mfPi(X,0,Z); if isscalar(matPi); matPi = repmat(matPi,nz,nx); end

fGradPi = symfun(GradPi,[x y z]);
fGradPix = symfun(GradPix,[x y z]);
fGradPiy = symfun(GradPiy,[x y z]);
fGradPiz = symfun(GradPiz,[x y z]);
mfGradPix = matlabFunction(fGradPix);
mfGradPiy = matlabFunction(fGradPiy);
mfGradPiz = matlabFunction(fGradPiz);
matGradPix = mfGradPix(X,0,Z); if isscalar(matGradPix); matGradPix = repmat(matGradPix,nz,nx); end
matGradPiy = mfGradPiy(X,0,Z); if isscalar(matGradPiy); matGradPiy = repmat(matGradPiy,nz,nx); end
matGradPiz = mfGradPiz(X,0,Z); if isscalar(matGradPiz); matGradPiz = repmat(matGradPiz,nz,nx); end

fE = symfun(E,[x y z]);
fEx = symfun(Ex,[x y z]);
fEy = symfun(Ey,[x y z]);
fEz = symfun(Ez,[x y z]);
mfEx = matlabFunction(fEx);
mfEy = matlabFunction(fEy);
mfEz = matlabFunction(fEz);
matEx = mfEx(X,0,Z); if isscalar(matEx); matEx = repmat(matEx,nz,nx); end
matEy = mfEy(X,0,Z); if isscalar(matEy); matEy = repmat(matEy,nz,nx); end
matEz = mfEz(X,0,Z); if isscalar(matEz); matEz = repmat(matEz,nz,nx); end

fE = symfun(E,[x y z]);
fEx = symfun(Ex,[x y z]);
fEy = symfun(Ey,[x y z]);
fEz = symfun(Ez,[x y z]);
mfEx = matlabFunction(fEx);
mfEy = matlabFunction(fEy);
mfEz = matlabFunction(fEz);
matEx = mfEx(X,0,Z); if isscalar(matEx); matEx = repmat(matEx,nz,nx); end
matEy = mfEy(X,0,Z); if isscalar(matEy); matEy = repmat(matEy,nz,nx); end
matEz = mfEz(X,0,Z); if isscalar(matEz); matEz = repmat(matEz,nz,nx); end



% Densitites and fluxes for ions and electrons respectively.
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

% Plot
nrows = 5;
ncols = 1;
npanels = nrows*ncols;
%for ipanel = 1:npanels
%  h(ipanel) = subplot(nrows,ncols,ipanel);  
%end
h = setup_subplots(nrows,ncols);
isub = 1;

if 0 % Ax
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matAx);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Ax';
end
if 1 % Ay
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matAy);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Ay';
  hold(hca,'on')
  Alev = linspace(floor(min(matAy(:))),ceil(max(matAy(:))),15);
  Alev = floor(min(matAy(:))):0.5:ceil(max(matAy(:)));
  contour(hca,X,Z,matAy,Alev,'k')
  hold(hca,'off')
end
if 0 % Az
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matAz);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Az';
end
if 1 % Bx
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matBx);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Bx';
end
if 0 % By
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matBy);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'By';
end
if 1 % Bz
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matBz);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Bz';
end
if 0 % Jx
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matJx);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Jx';
end
if 1 % Jy
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matJy);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Jy';
end
if 1 % Jy
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matJy);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Jy';
  hca.CLim = 0.01*[-1 1];
end
if 0 % log10(abs(Jy))
  hca = h(isub); isub = isub + 1;
  plmatJy = matJy; plmatJy(plmatJy<0) = NaN;
  pcolor(hca,X,Z,log10((plmatJy)));
  shading(hca,'flat')
  %hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'log10(Jy)';
end
if 0 % Jz
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matJz);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Jz';
end

if 0
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matJz);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Jz';
end
if 0 % Ez
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matEz);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Ez';
end
colormap(cn.cmap('blue_red'))
compact_panels(0.03,0.08)

hlinks = linkprop(h,{'XLim','YLim'});
%h(1).YLim = [0 6];
%hlinks.Targets(1).XLim = [-10 10];
%hlinks.Targets(1).YLim = [-10 10];

%% Keeping constants
nSimp = 50;
wpewce = 2;
mime = 100;
% de, sixe of box
xemin = 0; 
xemax = 2048; 
zemin = -128;
zemax = 128;
%zemax = 56;
%zemin = 0;
% di/de = (c/wpi)/(c/wpe) = wpe/wpi = sqrt(mi/me) 
ximin = xemin/sqrt(mime);
ximax = xemax/sqrt(mime);
zimin = zemin/sqrt(mime);
zimax = zemax/sqrt(mime);

xmin = ximin;
xmax = ximax;
zmin = zimin;
zmax = zimax;

Lhe = 20; % Harris sheet width in de
Lhi = Lhe/sqrt(mime); % Harris sheet width in di
l = Lhi; % harris sheet width
zh = l;
B0 = 1; % harris sheet asymptotic amplitude
BG = 0*0.25;
ah = B0*l; % rot(A) = B -> A0/l = B0 ....... rot(B) = J -> B0/l = Jy0
xp = 2.0*l; % scale length of perturbation
zp = 1.5*l; % scale length of perturbation
ap = ah*0.7; % amplitude of perturbation

TeTi = 1/5; % not used, implement to get densities
Ttot = 0.5; % not used, implement to get densities

nx = 640;
nz = 160;
x0 = 0.5*(xmax-xmin);
xvec = linspace(xmin,xmax,nx)-x0;
zvec = linspace(zmin,zmax,nz);
[X,Z] = meshgrid(xvec,zvec);

% If you want to get the expressions in terms of all parameters (zh,zp,xp,
% ah,ap, etc), you need to add these as syms here.
% syms x y z zh xp zp ah ap
% and
% add them as variables below, ex
% fA = symfun(A,[x y z]); 
% -->
% fA = symfun(A,[x y z zh xp zp ah ap]);
% Otherwise they are mulitplied together as numbers
syms x y z %ah zh ap xp zp

R = [x y z];

z0_BG = 5;
zL_BG = 1;

doLocal = 0;
% Harris current sheet
AH = [0; -ah*log(cosh(z/zh)); 0]; 
% Perturbation
AP = [0; -ap*exp(-x^2/2/xp^2 -z^2/2/zp^2 + 0.0); 0]; % circular perturbation
AP = [0; -ap*exp(-x^2/2/xp^2 -z^2/2/zp^2 + 0.0); 0]; % circular perturbation
%AP = [0; -ap*cos(2*pi*(x)/xmax)*cos(pi*(z)/zmax)*exp(-x^2/2/xp^2 -z^2/2/zp^2 + 0.0); 0]; % circular perturbation
AG = [BG*z + 1*BG*log(cosh((z-z0_BG)/zL_BG)); 0; 0]; % guide field

%AP = [0; ah*log(cosh(z/zh))*exp(-x^2/2/xp^2 -z^2/2/zp^2); 0]; 
APsin = [0; -0.2*ap*cos(2*pi*x/(xmax-xmin))*cos(1*pi*abs(z)/(zmax-zmin)); 0];
%AP = [0 ap*exp(-x^2/2/xp^2 + 0.5)*cos(pi*z/lz) 0];
%AP = [0 a0*sech(pi*x/xm)*sech(pi*z/zm) 0];
%AP = [0 0+a0*cos(0.5*pi*(x/xm))*cos(0.5*pi*(z/zm)) 0];
  
A = AH;% + AP + AG + 1*APsin;
%A = AP;
B = curl(A,R);
J = curl(B,R);
divB = divergence(B,R);
PB = 0.5*(B(1).^2+B(2).^2+B(3).^2);
% Given T, n should be calculated such that PB+PT=const
% 
PT = 0.5*B0^2-PB;
% PT = PTe + PTi = PTi/(Ti/Te) + PTi = PTi*(Te/Ti + 1)
%                = PTe + PTe/(Te/Ti) = PTe*(1 + Ti/Te)
PTi = PT/(TeTi+1);
PTe = PT/(1+1/TeTi);
n = PT/Ttot;
nH = cosh(z/zh)^-2;

% How to add a polarized inflow current sheet, such that the electrons are
% drifting, and not the ions?
% Particle velocities, NOT IMPLEMENTED
velocities = 'momentum conservation';
switch velocities
  case 'diamagnetic only' % no electric field
    
  case 'stationary ions' % finite electric field
    
  case 'momentum conservation'
    % mi*ni*vi + me*ne*ve = 0
    % ni*vi = - (me/mi)*ne*ve
    % J = ni*vi-ne*ve 
    %   = ni*vi + (mi/me)*ni*vi = (1 + mi/me)*ni*vi
    %   = - (me/mi)*ne*ve - ne*ve = - (me/mi + 1)*ne*ve
    % vi = J/[ni*(1 + mi/me)]
    % ve = -J/[ne*(me/mi + 1)]
    Vi = J./(n*(mime+1));
    Ve = -J./(n*(1/mime+1));
    VixB = simplify(cross(Vi,B),nSimp);
    VexB = cross(Ve,B);
    % E+vixB - div(Pti)/ne = 0
    GradPi = gradient(PTi,R);
    GradPe = gradient(PTe,R);    
    GradPin = GradPi/n;
    GradPen = GradPe/n;
    E = -VixB + GradPi./n;
    % E becomes a tanh function, which is not what I want
end

Ax = A(1);
Ay = A(2);
Az = A(3);
Bx = B(1);
By = B(2);
Bz = B(3);
Jx = J(1);
Jy = J(2);
Jz = J(3);
Viy = Vi(2);
Vey = Ve(2);
VixBx = VixB(1);
VixBy = VixB(2);
VixBz = VixB(3);
VexBx = VixB(1);
VexBy = VixB(2);
VexBz = VixB(3);
Ex = E(1);
Ey = E(2);
Ez = E(3);
GradPix = GradPi(1);
GradPiy = GradPi(2);
GradPiz = GradPi(3);

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

fVi = symfun(Vi,[x y z]);
fViy = symfun(Viy,[x y z]);
mfViy = matlabFunction(fViy);
matViy = mfViy(X,0,Z); if isscalar(matViy); matViy = repmat(matViy,nz,nx); end
fVe = symfun(Ve,[x y z]);
fVey = symfun(Vey,[x y z]);
mfVey = matlabFunction(fVey);
matVey = mfVey(X,0,Z); if isscalar(matVey); matVey = repmat(matVey,nz,nx); end

fVixB = symfun(VixB,[x y z]);
fVixBx = symfun(VixBx,[x y z]);
fVixBy = symfun(VixBy,[x y z]);
fVixBz = symfun(VixBz,[x y z]);
mfVixBx = matlabFunction(fVixBx);
mfVixBy = matlabFunction(fVixBy);
mfVixBz = matlabFunction(fVixBz);
matVixBx = mfVixBx(X,0,Z); if isscalar(matVixBx); matVixBx = repmat(matVixBx,nz,nx); end
matVixBy = mfVixBy(X,0,Z); if isscalar(matVixBy); matVixBy = repmat(matVixBy,nz,nx); end
matVixBz = mfVixBz(X,0,Z); if isscalar(matVixBz); matVixBz = repmat(matVixBz,nz,nx); end

fVexB = symfun(VexB,[x y z]);
fVexBx = symfun(VexBx,[x y z]);
fVexBy = symfun(VexBy,[x y z]);
fVexBz = symfun(VexBz,[x y z]);
mfVexBx = matlabFunction(fVexBx);
mfVexBy = matlabFunction(fVexBy);
mfVexBz = matlabFunction(fVexBz);
matVexBx = mfVexBx(X,0,Z); if isscalar(matVexBx); matVexBx = repmat(matVexBx,nz,nx); end
matVexBy = mfVexBy(X,0,Z); if isscalar(matVexBy); matVexBy = repmat(matVexBy,nz,nx); end
matVexBz = mfVexBz(X,0,Z); if isscalar(matVexBz); matVexBz = repmat(matVexBz,nz,nx); end

fNH = symfun(nH,[x y z]);
mfNH = matlabFunction(fNH);
matNH = mfNH(X,0,Z); if isscalar(matNH); matNH = repmat(matNH,nz,nx); end

fN = symfun(n,[x y z]);
mfN = matlabFunction(fN);
matN = mfN(X,0,Z); if isscalar(matN); matN = repmat(matN,nz,nx); end

fPB = symfun(PB,[x y z]);
mfPB = matlabFunction(fPB);
matPB = mfPB(X,0,Z); if isscalar(matPB); matPB = repmat(matPB,nz,nx); end
fPi = symfun(PTi,[x y z]);
mfPi = matlabFunction(fPi);
matPi = mfPi(X,0,Z); if isscalar(matPi); matPi = repmat(matPi,nz,nx); end
fPe = symfun(PTe,[x y z]);
mfPe = matlabFunction(fPe);
matPe = mfPe(X,0,Z); if isscalar(matPe); matPi = repmat(matPe,nz,nx); end

fGradPi = symfun(GradPi,[x y z]);
fGradPix = symfun(GradPix,[x y z]);
fGradPiy = symfun(GradPiy,[x y z]);
fGradPiz = symfun(GradPiz,[x y z]);
mfGradPix = matlabFunction(fGradPix);
mfGradPiy = matlabFunction(fGradPiy);
mfGradPiz = matlabFunction(fGradPiz);
matGradPix = mfGradPix(X,0,Z); if isscalar(matGradPix); matGradPix = repmat(matGradPix,nz,nx); end
matGradPiy = mfGradPiy(X,0,Z); if isscalar(matGradPiy); matGradPiy = repmat(matGradPiy,nz,nx); end
matGradPiz = mfGradPiz(X,0,Z); if isscalar(matGradPiz); matGradPiz = repmat(matGradPiz,nz,nx); end

fE = symfun(E,[x y z]);
fEx = symfun(Ex,[x y z]);
fEy = symfun(Ey,[x y z]);
fEz = symfun(Ez,[x y z]);
mfEx = matlabFunction(fEx);
mfEy = matlabFunction(fEy);
mfEz = matlabFunction(fEz);
matEx = mfEx(X,0,Z); if isscalar(matEx); matEx = repmat(matEx,nz,nx); end
matEy = mfEy(X,0,Z); if isscalar(matEy); matEy = repmat(matEy,nz,nx); end
matEz = mfEz(X,0,Z); if isscalar(matEz); matEz = repmat(matEz,nz,nx); end

fE = symfun(E,[x y z]);
fEx = symfun(Ex,[x y z]);
fEy = symfun(Ey,[x y z]);
fEz = symfun(Ez,[x y z]);
mfEx = matlabFunction(fEx);
mfEy = matlabFunction(fEy);
mfEz = matlabFunction(fEz);
matEx = mfEx(X,0,Z); if isscalar(matEx); matEx = repmat(matEx,nz,nx); end
matEy = mfEy(X,0,Z); if isscalar(matEy); matEy = repmat(matEy,nz,nx); end
matEz = mfEz(X,0,Z); if isscalar(matEz); matEz = repmat(matEz,nz,nx); end



% Densitites and fluxes for ions and electrons respectively.
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