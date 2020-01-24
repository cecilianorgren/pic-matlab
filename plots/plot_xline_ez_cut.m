df08 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data_h5/fields.h5');
df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
df04n = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/data_h5/fields.h5');

%% Plot
xlim = 204 + [-0.5 0.5];
zlim = 1.5*[-1 1];

nrows = 2;
ncols = 3;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

pic = df08.xlim(xlim).zlim(zlim);
if 1 % Ez(z,t)
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.twci,pic.zi,squeeze(mean(pic.Ez,1)))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_z';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'z/d_i';
  hca.Title.String = sprintf('n_c = 0.8 n_0');
end
if 1 % n3(z,t)
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.twci,pic.zi,squeeze(mean(pic.n(3),1)))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_c';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'z/d_i';
end

pic = df04.xlim(xlim).zlim(zlim);
if 1 % Ez(z,t)
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.twci,pic.zi,squeeze(mean(pic.Ez,1)))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_z';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'z/d_i';
  hca.Title.String = sprintf('n_c = 0.4 n_0');
end
if 1 % n35(z,t)
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.twci,pic.zi,squeeze(mean(pic.n([3 5]),1)))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_c';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'z/d_i';
end

pic = df04n.xlim(xlim).zlim(zlim);
if 1 % Ez(z,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.twci,pic.zi,squeeze(mean(pic.Ez,1)))
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_z';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'z/d_i';
  hca.Title.String = sprintf('n_c = 0.4 n_0');
end
if 1 % n3(z,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.twci,pic.zi,squeeze(mean(pic.n([3 5]),1)))
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_c';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'z/d_i';
end

colormap(h(1),pic_colors('blue_red'))
colormap(h(3),pic_colors('blue_red'))
colormap(h(5),pic_colors('blue_red'))
hlinks_all = linkprop(h,{'XLim','YLim'});
hlinks_Ez = linkprop(h([1 3 5]),{'CLim'});
hlinks_Ez.Targets(1).CLim = max(abs(hlinks_Ez.Targets(1).CLim))*[-1 1];
hlinks_nc = linkprop(h([2 4 6]),{'CLim'});

