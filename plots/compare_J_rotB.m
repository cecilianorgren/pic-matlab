% turb = PIC('/Volumes/Fountain/Data/PIC/turbulencerun/data_h5/fields.h5');
% gf05 = PIC('/Volumes/Fountain/Data/PIC/guide_field_05_cold_ions/data_h5/fields.h5');
twpe = 5987;
twpe_minus = twpe - 1;
twpe_plus  = twpe + 1;
xlim = [00 50];
zlim = [0 10];
pic = turb.twpelim(twpe).xlim(xlim).zlim(zlim);
pic_minus = turb.twpelim(twpe_minus).xlim(xlim).zlim(zlim);
pic_plus  = turb.twpelim(twpe_plus).xlim(xlim).zlim(zlim);

% dEdt
tic
dtwpe = (twpe_plus-twpe_minus);
dtwci = (twpe_plus-twpe_minus)/(pic.mime*pic.wpewce);
c2 = pic.mime*pic.wpewce^2; % c^2-vA^2
Ex_minus = pic_minus.Ex;
Ey_minus = pic_minus.Ey;
Ez_minus = pic_minus.Ez;
Ex_plus = pic_plus.Ex;
Ey_plus = pic_plus.Ey;
Ez_plus = pic_plus.Ez;
dEdt.x = (1/c2)*(Ex_plus - Ex_minus)/dtwci;
dEdt.y = (1/c2)*(Ey_plus - Ey_minus)/dtwci;
dEdt.z = (1/c2)*(Ez_plus - Ez_minus)/dtwci;
toc
% rot(B)
tic
rotB = pic.curlb;
toc
% J
tic
J.x = pic.Jx;
J.y = pic.Jy;
J.z = pic.Jz;
toc
clim = 0.2*[-1 1];

% Plot
nrows = 2;
ncols = 3;
npanels = nrows*ncols;
% for ipanel = 1:npanels
%   h(ipanel) = subplot(nrows,ncols,ipanel);  
% end
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;

if 0 % rot(B)_x
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,rotB.x');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'rot(B)_x';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 0 % rot(B)_y
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,rotB.y');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'rot(B)_y';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 0 % rot(B)_z
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,rotB.z');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'rot(B)_z';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 0 % J_x
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,J.x');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'J_x';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 0 % J_y
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,J.y');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'J_y';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 0 % J_z
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,J.z');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'J_z';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % rot(B)_x-Jx
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,rotB.x'-J.x');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'rot(B)_x-J_x';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % rot(B)_y -Jy
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,rotB.y'-J.y');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'rot(B)_y-J_y';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % rot(B)_z - Jz
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,rotB.z'-J.z');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'rot(B)_z-J_z';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % dEdt_x
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,dEdt.x');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'dE/dt_x';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % dEdt_y
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,dEdt.y');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'dE/dt_y';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
if 1 % dEdt_z
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,dEdt.z');
  shading(hca,'flat')  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'dE/dt_z';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end
  
colormap(pic_colors('blue_red'))
hlinks = linkprop(h,{'XLim','YLim','CLim'});
hlinks.Targets(1).CLim = clim;
compact_panels(0.01,0.08)
%compact_panels(0.01,0.07)

for ipanel = 1:npanels
  irf_legend(h(ipanel),sprintf('twpe = %g',twpe),[0.98 0.98])
end
