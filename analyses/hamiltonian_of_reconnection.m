% Hamiltonian vector potential of magnetic reconenction
%df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
twci = 100;
xlim = [185 225];
zlim = [-3 3];
pic = df04.twcilim(twci).xlim(xlim).zlim(zlim);

m = 1;
q = 1;
A = pic.A;
B0 = 1;
BH = -0.5;
EH = 1;
L = .5;
LE = 0.2*L;
l = L*30;
x0 = 205;
f_Ay = @(x,z) B0*z.^2/2/L + -B0*(x-x0).^2/2/l;
%f_Ay = @(x,z) B0*z.^2/2/L + exp(-(x-x0).^2/2/l);
%f_Ay = @(x,z) B0*z.^2/2/L.*(x-x0).^2/2;
%f_Ay = @(x,z) B0*z.^2/2/L.*cos(-(x-x0).^2/(40*2*pi));
f_Ax = @(x,z) BH*abs(z).*(x-x0).^2; 
f_Ax = @(x,z) BH*cos(z).*(x-x0).^1; 
f_Ax = @(x,z) BH*exp(-z.^2/2/L).*(x-x0).^1; 
f_Phi = @(x,z) EH*exp(-z.^2/2/(LE)); 
f_Ez = @(x,z) EH*(z).^2/2/L; 

f_px = @(x,z,vx) vx*m + q*f_Ax(x,z);
f_py = @(x,z,vy) vy*m + q*f_Ay(x,z);
f_pz = @(x,z,vz) vz*m;

f_vx = @(x,z,px) px/m - q*f_Ax(x,z)/m;
f_vy = @(x,z,py) py/m - q*f_Ay(x,z)/m;
f_vz = @(x,z,pz) pz/m;

f_Ek = @(x,z,px,py,pz) 0.5*m*(f_vx(x,z,px).^2 + f_vy(x,z,py).^2 + f_vz(x,z,pz).^2);

px0 = f_px(x0,0,0.2);
py0 = f_py(x0,0,-0.2);
pz0 = f_pz(x0,0,0.5);

%f_Ax = @(x,z) BH*exp(-z.^2).*exp(-(x-x0).^2); 

%f_Ay = @(x,z) +B0*z.^2/2/L - B0*(x-x0).^2/2/l;
%f_Ay = @(x,z) B0*z.^2/2/L + -abs((x-x0)/100);
[XI,ZI] = ndgrid(pic.xi,pic.zi);
dx = pic.xi(2)-pic.xi(1);
dz = pic.zi(2)-pic.zi(1);

EK = f_Ek(XI,ZI,px0,py0,pz0);

PHI = f_Phi(XI,ZI);
dorder = 1;
ddim = 2;
EZ = [1*diff(PHI(1:2,:),dorder,ddim); diff(PHI,dorder,ddim)]/dz;


AY = f_Ay(XI,ZI);
AX = f_Ax(XI,ZI);
dorder = 1;
ddim = 2;
BX = [1*diff(AY(1:2,:),dorder,ddim); diff(AY,dorder,ddim)]/dz;
ddim = 1;
BZ = -[1*diff(AY(1:2,:),dorder,ddim); diff(AY,dorder,ddim)]/dx;
ddim = 2;
BY = [1*diff(AX(1:2,:),dorder,ddim); diff(AX,dorder,ddim)]/dz;




%B=Bzxˆ and B=Byˆ, 0Lgg

%Ax =Bgz, Ay =−B0z2, Az =0.


nrows = 7;
ncols = 2;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % Ay data
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.A')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'A_y^{data}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % Ay model
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,AY'-20)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'A_y^{mod}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % Bx from data
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.Bx')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'B_x^{data}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % Bx from Ay model
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,BX')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'B_x^{mod}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % Bz from Ay data
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.Bz')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'B_z^{data}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % Bz from Ay model
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,BZ')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'B_z^{mod}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % empty
  hca = h(isub); isub = isub + 1;  
end
if 1 % Ax model
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,AX')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'A_x^{mod}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % By from Ax data
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.By')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'B_y^{data}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % By from Ax model
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,BY')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'B_y^{mod}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % Ez from model
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,EZ')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_z^{mod}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % Ez from data
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.Ez')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_z^{data}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % EK from data
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,EK')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_k^{mod}';
  hca.YDir = 'normal';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
colormap(pic_colors('blue_red'))
