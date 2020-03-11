%% Plot timeline of E + vxB for all species
% turb = PIC('/Volumes/Fountain/Data/PIC/turbulencerun/data_h5/fields.h5');
% df04n = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/data_h5/fields.h5');
% df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
x0 = mean(pic.xi);
pic = df04n;
dx = 0.2; % width of band to average over
dz = 0.2; % width of band to average over
xlim = x0 + [-50 50]; zval = 0;
zlim = pic.zi([1 end]); xval = x0;

x_pic = pic.xlim(xlim).zlim(zval+0.5*dz*[-1 1]);
z_pic = pic.zlim(zlim).xlim(xval+0.5*dx*[-1 1]);

%meandim = 2; % z
iSpeciesAll = [1 3 5];
eSpeciesAll = [2 4 6];
iSpeciesCold = [3 5];
eSpeciesCold = [4 6];
iSpeciesHot = [1];
eSpeciesHot = [2];
x_vex_all = x_pic.vx(eSpeciesAll);
x_vey_all = x_pic.vy(eSpeciesAll);
x_vez_all = x_pic.vz(eSpeciesAll);
x_vix_all = x_pic.vx(iSpeciesAll);
x_viy_all = x_pic.vy(iSpeciesAll);
x_viz_all = x_pic.vz(iSpeciesAll);
x_vex_cold = x_pic.vx(eSpeciesCold);
x_vey_cold = x_pic.vy(eSpeciesCold);
x_vez_cold = x_pic.vz(eSpeciesCold);
x_vix_cold = x_pic.vx(iSpeciesCold);
x_viy_cold = x_pic.vy(iSpeciesCold);
x_viz_cold = x_pic.vz(iSpeciesCold);
x_vex_hot = x_pic.vx(eSpeciesHot);
x_vey_hot = x_pic.vy(eSpeciesHot);
x_vez_hot = x_pic.vz(eSpeciesHot);
x_vix_hot = x_pic.vx(iSpeciesHot);
x_viy_hot = x_pic.vy(iSpeciesHot);
x_viz_hot = x_pic.vz(iSpeciesHot);
x_Ex = x_pic.Ex;
x_Ey = x_pic.Ey;
x_Ez = x_pic.Ez;
x_Bx = x_pic.Bx;
x_By = x_pic.By;
x_Bz = x_pic.Bz;

x_vexB_hot = cross_product(x_vex_hot,x_vey_hot,x_vez_hot,x_Bx,x_By,x_Bz);
x_vexB_cold = cross_product(x_vex_cold,x_vey_cold,x_vez_cold,x_Bx,x_By,x_Bz);
x_vexB_all = cross_product(x_vex_all,x_vey_all,x_vez_all,x_Bx,x_By,x_Bz);
x_vixB_hot = cross_product(x_vix_hot,x_viy_hot,x_viz_hot,x_Bx,x_By,x_Bz);
x_vixB_cold = cross_product(x_vix_cold,x_viy_cold,x_viz_cold,x_Bx,x_By,x_Bz);
x_vixB_all = cross_product(x_vix_all,x_viy_all,x_viz_all,x_Bx,x_By,x_Bz);

x_E_vexB_hot.x = squeeze(mean(x_vexB_hot.x + x_Ex,2));
x_E_vexB_cold.x = squeeze(mean(x_vexB_cold.x + x_Ex,2));
x_E_vexB_all.x = squeeze(mean(x_vexB_all.x + x_Ex,2));
x_E_vixB_hot.x = squeeze(mean(x_vixB_hot.x + x_Ex,2));
x_E_vixB_cold.x = squeeze(mean(x_vixB_cold.x + x_Ex,2));
x_E_vixB_all.x = squeeze(mean(x_vixB_all.x + x_Ex,2));
x_E_vexB_hot.y = squeeze(mean(x_vexB_hot.y + x_Ey,2));
x_E_vexB_cold.y = squeeze(mean(x_vexB_cold.y + x_Ey,2));
x_E_vexB_all.y = squeeze(mean(x_vexB_all.y + x_Ey,2));
x_E_vixB_hot.y = squeeze(mean(x_vixB_hot.y + x_Ey,2));
x_E_vixB_cold.y = squeeze(mean(x_vixB_cold.y + x_Ey,2));
x_E_vixB_all.y = squeeze(mean(x_vixB_all.y + x_Ey,2));
x_B.x = squeeze(mean(x_Bx,2));
x_B.y = squeeze(mean(x_By,2));
x_B.z = squeeze(mean(x_Bz,2));

%%
clim = 0.5*[-1 1];
nrows = 7;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % Bz
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_B.z');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'B_x (B_0)'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end

% x
if 0 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vexB_hot.x');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_x (v_{A0}B_0)','hot electrons'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
if 0 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vexB_cold.x');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_x (v_{A0}B_0)','cold electrons'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
if 0 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vexB_all.x');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_x (v_{A0}B_0)','all electrons'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
if 0 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vixB_hot.x');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_x (v_{A0}B_0)','hot ions'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
if 0 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vixB_cold.x');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_x (v_{A0}B_0)','cold ions'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
if 0 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vixB_all.x');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_x (v_{A0}B_0)','all ions'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
% y
if 1 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vexB_hot.y');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_y (v_{A0}B_0)','hot electrons'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
if 1 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vexB_cold.y');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_y (v_{A0}B_0)','cold electrons'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
if 1 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vexB_all.y');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_y (v_{A0}B_0)','all electrons'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
if 1 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vixB_hot.y');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_y (v_{A0}B_0)','hot ions'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
if 1 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vixB_cold.y');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_y (v_{A0}B_0)','cold ions'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
if 1 % x_E_vexB_hot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,x_pic.xi,x_pic.twci,x_E_vixB_all.y');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = {'(E+vxB)_y (v_{A0}B_0)','all ions'};
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end

colormap(pic_colors('blue_red'))
h(1).Title.String = sprintf('z = [%g, %g]',zval+0.5*dz,zval-0.5*dz);
hlinks = linkprop(h,{'XLim','YLim'});
compact_panels(0.01)

for ip = 1:npanels
  h(ip).XTick = 0:10:max(pic.xi);
  h(ip).YTick = 0:20:500;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).Layer = 'top';
  h(ip).CLim = clim;
  h(ip).FontSize = 14;
end
