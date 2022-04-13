%
mtot = 1.4e10/40;
M_target = [1e9 3e9 3e9]/4;
M_target = mtot*[0.05 0.23 0.21];

nSimp = 50;
wpewce = 4;
mime = 100;
resx = 2;
resz = 3;
% de, sixe of box
xemin = 0; 
xemax = 3000;256*sqrt(mime); 
zemin = 0;
zemax = 2*256;48*sqrt(mime);
zemean = (zemin+zemax)/2;
zemin = zemin-zemean;
zemax = zemax-zemean;
zemean = (zemin+zemax)/2;

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

nb = 0.1;
Lhe = 10*2; % Harris sheet width in de, 5de = 1di
Lhi = Lhe/sqrt(mime); % Harris sheet width in di
l = Lhi; % harris sheet width
zh = l;
B0 = 1; % harris sheet asymptotic amplitude
BG = 0*0.25;
ah = B0*l; % rot(A) = B -> A0/l = B0 ....... rot(B) = J -> B0/l = Jy0

% no_hot_ng_n02, very violent, but quick
xp = 1*3*l; % scale length of perturbation
zp = 1*1.5*l; % scale length of perturbation
ap = ah*1; % amplitude of perturbation

% no_hot_bg_n02_m100, calm start, also slow
xp = 2*l; % scale length of perturbation
zp = 1*l; % scale length of perturbation
ap = 1.0*ah; % amplitude of perturbation

% adapted to LH=1di
xp = 2*l; % scale length of perturbation
zp = 1*l; % scale length of perturbation
ap = 0.5*ah; % amplitude of perturbation

% no_hot_bg_test, calm start, but kind of slow
% xp = 1*l; % scale length of perturbation
% zp = 1*l; % scale length of perturbation
% ap = 0.25*ah; % amplitude of perturbation
% 
% %
% xp = 1.5*l; % scale length of perturbation
% zp = 1.5*l; % scale length of perturbation
% ap = 0.6*ah; % amplitude of perturbation

TeTi = 1/5; % not used, implement to get densities
Ttot = 0.5; % not used, implement to get densities

nx = (xmax-xmin)*sqrt(mime)*resx;
nz = (zmax-zmin)*sqrt(mime)*resz;
x0 = 0.5*(xmax-xmin);
xvec = linspace(xmin,xmax,nx)-x0;
zvec = linspace(zmin,zmax,nz);
[X,Z] = meshgrid(xvec,zvec);
%
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

z0_BG = 5;
zL_BG = 2;

doLocal = 0;
% Harris current sheet
AH = [0; -ah*log(cosh(z/zh)); 0]; 
% Guide field
%AG = [0; 0; ah];
AG = [BG*z + 1*BG*(cosh((z-z0_BG)/zL_BG)); 0; 0]; % guide field
AG = [ah*z; 0; 0]; % guide field
AG = [ah*log(cosh((z-0)/zL_BG)); 0; 0]; % guide field
% Perturbation
AP = [0; -ap*exp(-x^2/2/xp^2 -z^2/2/zp^2 + 0.0); 0]; % circular perturbation
AP = [0; -ap*exp(-x^2/2/xp^2 -z^2/2/zp^2 + 0.0); 0]; % circular perturbation
%AP = [0; -ap*cos(2*pi*(x)/xmax)*cos(pi*(z)/zmax)*exp(-x^2/2/xp^2 -z^2/2/zp^2 + 0.0); 0]; % circular perturbation


%AP = [0; ah*log(cosh(z/zh))*exp(-x^2/2/xp^2 -z^2/2/zp^2); 0]; 
APsin = [0; -0.2*ah*cos(2*pi*x/(xmax-xmin)/1)*cos(1*pi*abs(z)/(zmax-zmin)); 0];
%AP = [0 ap*exp(-x^2/2/xp^2 + 0.5)*cos(pi*z/lz) 0];
%AP = [0 a0*sech(pi*x/xm)*sech(pi*z/zm) 0];
%AP = [0 0+a0*cos(0.5*pi*(x/xm))*cos(0.5*pi*(z/zm)) 0];
  
A = AH +0*AG + 1*AP + 0*APsin;% + AP + AG + 1*APsin;
%A = AP;
B = curl(A,R);
BH = curl(AH,R);
J = curl(B,R);
divB = divergence(B,R);
PB = 0.5*(B(1).^2+B(2).^2+B(3).^2);
PBH = 0.5*(BH(1).^2+BH(2).^2+BH(3).^2);
% Given T, n should be calculated such that PB+PT=const
% 
PT = 0.5*B0^2-PB; % Plasma pressure to balance Harris sheet pressure
% PT = PTe + PTi = PTi/(Ti/Te) + PTi = PTi*(Te/Ti + 1)
%                = PTe + PTe/(Te/Ti) = PTe*(1 + Ti/Te)
PTi = PT/(TeTi+1);
PTe = PT/(1+1/TeTi);
nh = PT/Ttot; % Hot plasma density (mainly Harris, but also some density due to B perturbation)


% Specify some variables
dn = 0.05; % amplitude of transition
% locations of transitions
z1 = 5;
z2 = 8;
z3 = 12;
z4 = 15;

%z1 = 1;
z1 = 4;
z2 = 7;
z3 = 10;
z4 = 13;
z5 = 16;
z6 = 19;
lin = 0.25;
l = 0.5;

B0 = 1;
R = [x; y; z];
% You can also try it with a constant B0, and later add the tangent field.
%Btot = B0*tanh(z/l);
BH = B0*tanh(z/l);
Btot = B0;
iPert = 6;
switch iPert
  case 1
    f1 = +0.5*dn*(1 + tanh((abs(z)-z1)/lin));
    f2 = -0.5*dn*(1 + tanh((abs(z)-z2)/lin));
    f3 = +0.5*dn*(1 + tanh((abs(z)-z3)/lin));
    f4 = -0.5*dn*(1 + tanh((abs(z)-z4)/lin));
    f5 = +0.5*dn*(1 + tanh((abs(z)-z5)/lin));
    f6 = -0.5*dn*(1 + tanh((abs(z)-z6)/lin));
    f = f1 + f2 + f3 + f4 + f5 + f6;
  case 2 % ok for L = 2di;
    f1a = +0.5*dn*(1 - tanh((abs(z)-6)/lin).^1);
    f1b = -0.5*dn*(1 - tanh((abs(z)-4)/(lin*4)).^1);
    f2a = +0.5*dn*(1 - tanh((abs(z)-9)/lin).^1);
    f2b = -0.5*dn*(1 - tanh((abs(z)-4)/(lin*4)).^1);
    f = f1a + f1b + f2a + f2b;
  case 3 % ok for L = 1di;
    f1a = +0.5*dn*(1 - tanh((abs(z)-5)/lin).^1);
    f1b = -0.5*dn*(1 - tanh((abs(z)-2)/(lin*4)).^1);
    f2a = +0.5*dn*(1 - tanh((abs(z)-8)/lin).^1);
    f2b = -0.5*dn*(1 - tanh((abs(z)-2)/(lin*4)).^1);
    f = f1a + f1b + f2a + f2b;
  case 4 % 
    % From the top (z>0)
    f1a = +4*0.5*dn*(1 + tanh(((z)-3.5)/1.4).^1);
    f1b = -1.0*0.5*dn*(1 + tanh(((z)-6)/0.5).^1);
    f1c = -1.0*0.5*dn*(1 + tanh(((z)-9)/0.5).^1);
    
    
    f2a = +2*0.5*dn*(1 + tanh(((-z)-3.8)/1.4).^1);
    f2b = -1*0.5*dn*(1 + tanh(((-z)-6)/0.5).^1);
    f2c = -1*0.5*dn*(1 + tanh(((-z)-9)/0.5).^1);
    %f1a = +0.5*dn*(1 - tanh((abs(z)-5)/lin).^1);
    %f1b = -0.5*dn*(1 - tanh((abs(z)-2)/(lin*4)).^1);
    %f2a = +0.5*dn*(1 - tanh((abs(z)-8)/lin).^1);
    %f2b = -0.5*dn*(1 - tanh((abs(z)-2)/(lin*4)).^1);
    %f = f1a + f1b + f2a + f2b;
    %f = f1a + f1b + f2a + f2b;
    ftop = f1a + f1b + f1c + 0.001;
    fbot = f2a + 0.001; % + f2b + f2c;
    f = ftop + fbot;
  case 5 
    % From t top (z>0)
    f1a = +3*0.5*dn*(1 + tanh(((z)-4)/1.4).^1);
    %f1b = -1*0.5*dn*(1 + tanh(((z)-6)/0.5).^1);
    %f1c = -1*0.5*dn*(1 + tanh(((z)-9)/0.5).^1);
    
    
    f2a = +1*0.5*dn*(1 + tanh(((-z)-6)/1.4).^1);
    %f2b = -1*0.5*dn*(1 + tanh(((-z)-6)/0.5).^1);
    %f2c = -1*0.5*dn*(1 + tanh(((-z)-9)/0.5).^1);
    %f1a = +0.5*dn*(1 - tanh((abs(z)-5)/lin).^1);
    %f1b = -0.5*dn*(1 - tanh((abs(z)-2)/(lin*4)).^1);
    %f2a = +0.5*dn*(1 - tanh((abs(z)-8)/lin).^1);
    %f2b = -0.5*dn*(1 - tanh((abs(z)-2)/(lin*4)).^1);
    %f = f1a + f1b + f2a + f2b;
    %f = f1a + f1b + f2a + f2b;
    ftop = f1a + 0.001;
    fbot = f2a + 0.001; % + f2b + f2c;
    f = ftop + fbot;
  case 6 % 
    % From the top (z>0)
    f1a = +6*0.5*dn*(1 + tanh(((z)-3)/1.0).^1);
    f1b = -3*0.5*dn*(1 + tanh(((z)-6)/0.5).^1);
    f1c = -3*0.5*dn*(1 + tanh(((z)-9)/0.5).^1);
    
    % From the bottom (z>0)
    f2a = +2*0.5*dn*(1 + tanh(((-z)-3.8)/1.4).^1);
    
    ftop = f1a + f1b + f1c + 0.001;
    fbot = f2a + 0.001; % + f2b + f2c;
    f = ftop + fbot;
end
n_pert = symfun(f,[x z]);
n_top = symfun(ftop,[x z]);
n_bot = symfun(fbot,[x z]);

nb = 0.05*tanh(abs(z)/4)^2;

n = nh + 0*nb + 1*n_pert;
n1 = nh;
%n2 = nh;
n2 = n_top;
n3 = n_bot;
%n5 =
%n6 =
% Construct function for the number of particles given density profile and
% target number of particles for each species.
% for i in np.arange(0,int(nx)-1):
%         for j in np.arange(0,int(ny)-1):
% #            n_center[i,j] = nfun(x_center[i],y_center[j]) # make this a sum directly to save memory
%             Ntot = Ntot + nfun(x_center[i],y_center[j]) # make sum directly to save memory    
% #    Ntot = sum(n_center.flatten())*dx*dy # approximate total number of micro particles
%     def mfun(x,y):
%         return m.ceil(nfun(x,y))*(Mtot_target/Ntot) # at least one per cell, unless it's absolutely zero
%     return mfun

%mfun = @(Mtot_target,nfun) nfun(X,Z)*(Mtot_target/);
c_eval('mf_n? = matlabFunction(n?);',1:numel(M_target))
c_eval('sum_n? = sum(sum(mf_n?(X,Z)));',1:numel(M_target))
c_eval('mf_mfun? = matlabFunction(n?*M_target(?)/sum_n?);',1:numel(M_target))

% How to add a polarized inflow current sheet, such that the electrons are
% drifting, and not the ions?

Ax = A(1);
Ay = A(2);
Az = A(3);
Bx = B(1);
By = B(2);
Bz = B(3);
Jx = J(1);
Jy = J(2);
Jz = J(3);
% Viy = Vi(2);
% Vey = Ve(2);
% VixBx = VixB(1);
% VixBy = VixB(2);
% VixBz = VixB(3);
% VexBx = VixB(1);
% VexBy = VixB(2);
% VexBz = VixB(3);
% Ex = E(1);
% Ey = E(2);
% Ez = E(3);
% GradPix = GradPi(1);
% GradPiy = GradPi(2);
% GradPiz = GradPi(3);

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

% fVi = symfun(Vi,[x y z]);
% fViy = symfun(Viy,[x y z]);
% mfViy = matlabFunction(fViy);
% matViy = mfViy(X,0,Z); if isscalar(matViy); matViy = repmat(matViy,nz,nx); end
% fVe = symfun(Ve,[x y z]);
% fVey = symfun(Vey,[x y z]);
% mfVey = matlabFunction(fVey);
% matVey = mfVey(X,0,Z); if isscalar(matVey); matVey = repmat(matVey,nz,nx); end

% fVixB = symfun(VixB,[x y z]);
% fVixBx = symfun(VixBx,[x y z]);
% fVixBy = symfun(VixBy,[x y z]);
% fVixBz = symfun(VixBz,[x y z]);
% mfVixBx = matlabFunction(fVixBx);
% mfVixBy = matlabFunction(fVixBy);
% mfVixBz = matlabFunction(fVixBz);
% matVixBx = mfVixBx(X,0,Z); if isscalar(matVixBx); matVixBx = repmat(matVixBx,nz,nx); end
% matVixBy = mfVixBy(X,0,Z); if isscalar(matVixBy); matVixBy = repmat(matVixBy,nz,nx); end
% matVixBz = mfVixBz(X,0,Z); if isscalar(matVixBz); matVixBz = repmat(matVixBz,nz,nx); end

% fVexB = symfun(VexB,[x y z]);
% fVexBx = symfun(VexBx,[x y z]);
% fVexBy = symfun(VexBy,[x y z]);
% fVexBz = symfun(VexBz,[x y z]);
% mfVexBx = matlabFunction(fVexBx);
% mfVexBy = matlabFunction(fVexBy);
% mfVexBz = matlabFunction(fVexBz);
% matVexBx = mfVexBx(X,0,Z); if isscalar(matVexBx); matVexBx = repmat(matVexBx,nz,nx); end
% matVexBy = mfVexBy(X,0,Z); if isscalar(matVexBy); matVexBy = repmat(matVexBy,nz,nx); end
% matVexBz = mfVexBz(X,0,Z); if isscalar(matVexBz); matVexBz = repmat(matVexBz,nz,nx); end

fN = symfun(n,[x y z]);
mfN = matlabFunction(fN);
matN = mfN(X,0,Z); if isscalar(matN); matN = repmat(matN,nz,nx); end

fNh = symfun(nh,[x y z]);
mfNh = matlabFunction(fNh);
matNh = mfNh(X,0,Z); if isscalar(matNh); matNh = repmat(matNh,nz,nx); end

fNtop = symfun(n_top,[x y z]);
mfNtop = matlabFunction(fNtop);
matNtop = mfNtop(X,0,Z); if isscalar(matNtop); matNtop = repmat(matNtop,nz,nx); end

fNbot = symfun(n_bot,[x y z]);
mfNbot = matlabFunction(fNbot);
matNbot = mfNbot(X,0,Z); if isscalar(matNbot); matNbot = repmat(matNbot,nz,nx); end

% fNp = symfun(n_pert,[x y z]);
% mfNp = matlabFunction(fNp);
% matNp = mfNp(X,0,Z); if isscalar(matNp); matNp = repmat(matNp,nz,nx); end
% 
% fNb = symfun(nb,[x y z]);
% mfNb = matlabFunction(fNb);
% matNb = mfNb(X,0,Z); if isscalar(matNb); matNb = repmat(matNb,nz,nx); end

fPi = symfun(PTi,[x y z]);
mfPi = matlabFunction(fPi);
matPi = mfPi(X,0,Z); if isscalar(matPi); matPi = repmat(matPi,nz,nx); end

% fGradPi = symfun(GradPi,[x y z]);
% fGradPix = symfun(GradPix,[x y z]);
% fGradPiy = symfun(GradPiy,[x y z]);
% fGradPiz = symfun(GradPiz,[x y z]);
% mfGradPix = matlabFunction(fGradPix);
% mfGradPiy = matlabFunction(fGradPiy);
% mfGradPiz = matlabFunction(fGradPiz);
% matGradPix = mfGradPix(X,0,Z); if isscalar(matGradPix); matGradPix = repmat(matGradPix,nz,nx); end
% matGradPiy = mfGradPiy(X,0,Z); if isscalar(matGradPiy); matGradPiy = repmat(matGradPiy,nz,nx); end
% matGradPiz = mfGradPiz(X,0,Z); if isscalar(matGradPiz); matGradPiz = repmat(matGradPiz,nz,nx); end

% fE = symfun(E,[x y z]);
% fEx = symfun(Ex,[x y z]);
% fEy = symfun(Ey,[x y z]);
% fEz = symfun(Ez,[x y z]);
% mfEx = matlabFunction(fEx);
% mfEy = matlabFunction(fEy);
% mfEz = matlabFunction(fEz);
% matEx = mfEx(X,0,Z); if isscalar(matEx); matEx = repmat(matEx,nz,nx); end
% matEy = mfEy(X,0,Z); if isscalar(matEy); matEy = repmat(matEy,nz,nx); end
% matEz = mfEz(X,0,Z); if isscalar(matEz); matEz = repmat(matEz,nz,nx); end

% fE = symfun(E,[x y z]);
% fEx = symfun(Ex,[x y z]);
% fEy = symfun(Ey,[x y z]);
% fEz = symfun(Ez,[x y z]);
% mfEx = matlabFunction(fEx);
% mfEy = matlabFunction(fEy);
% mfEz = matlabFunction(fEz);
% matEx = mfEx(X,0,Z); if isscalar(matEx); matEx = repmat(matEx,nz,nx); end
% matEy = mfEy(X,0,Z); if isscalar(matEy); matEy = repmat(matEy,nz,nx); end
% matEz = mfEz(X,0,Z); if isscalar(matEz); matEz = repmat(matEz,nz,nx); end



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
figure(101)
nrows = 1;
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
if 1*0 % Ay
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
if 0 % Bx
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matBx);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Bx';
end
if 1*0 % By
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matBy);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'By';
end
if 0 % Bz
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
if 1*0 % Jy
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matJy);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'Jy';
end
if 0 % Jy
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
if 1*0 % n
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,matN);
  shading(hca,'flat')
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colorbar('peer',hca)
  hca.Title.String = 'n';
  hca.CLim = [0 1.2];
end


if 1 % n_cut
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,fN(20,0,zvec),zvec,fNh(20,0,zvec),zvec,fNtop(20,0,zvec),zvec,fNbot(20,0,zvec))%,zvec,fNb(20,0,zvec))  
  legend(hca,{'n_{tot}','n_{hot}','n_{cold,top}','n_{cold,bot}'})
  %hca.YLim = hca.YLim+diff(hca.YLim)*0.05*[-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XTick = -30:2:30;
  hca.XLabel.String = 'z/d_{i0}';
  hca.YLabel.String = 'n/n_{0}';
  hca.XLim = zvec([1 end]);
end
if 1*0 % n_cut (abs(zvec))
  hca = h(isub); isub = isub + 1;  
  plot(hca,abs(zvec),fN(20,0,zvec),abs(zvec),fNh(20,0,zvec),abs(zvec),fNtop(20,0,zvec),abs(zvec),fNbot(20,0,zvec))%,zvec,fNb(20,0,zvec))  
  legend(hca,{'n_{tot}','n_{hot}','n_{cold,top}','n_{cold,bot}'})
  %hca.YLim = hca.YLim+diff(hca.YLim)*0.05*[-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XTick = -30:2:30;
  hca.XLabel.String = '|z|/d_{i0}';
  hca.YLabel.String = 'n/n_{0}';
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

%hlinks = linkprop(h(1:5),{'XLim','YLim'});
%h(1).YLim = [0 6];
%hlinks.Targets(1).XLim = [-10 10];
%hlinks.Targets(1).YLim = [-10 10];

%% Plot macroparticles
figure(102)
nrows = 5;
ncols = 1;
npanels = nrows*ncols;
%for ipanel = 1:npanels
%  h(ipanel) = subplot(nrows,ncols,ipanel);  
%end
nPop = 3;
h = setup_subplots(nPop,1);
isub = 1;

for iPop = 1:nPop  
    hca = h(isub); isub = isub + 1;
    c_eval('mfun = mf_mfun?;',iPop);
    pcolor(hca,X,Z,mfun(X,Z));
    shading(hca,'flat')
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = sprintf('Species %g',iPop);    
    hca.Title.String = sprintf('# Macro particles');
end
colormap(pic_colors('candy4'))
compact_panels(0.03,0.08)

hlinks = linkprop(h,{'XLim','YLim','CLim'});
%h(1).CLim = [0 500];
%h(1).YLim = [0 6];
%hlinks.Targets(1).XLim = [-10 10];
%hlinks.Targets(1).YLim = [-10 10];


disp('Done.')