df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
df08 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data_h5/fields.h5');

%% Compare density at equivalent stage

zlim = [-0.1 0.1];
nh04 = squeeze(mean(df04.zlim(zlim).n(1),2));
nc04 = squeeze(mean(df04.zlim(zlim).n(3),2));
nh08 = squeeze(mean(df08.zlim(zlim).n(1),2));
nc08 = squeeze(mean(df08.zlim(zlim).n(3),2));

%% Plot
nrows = 6;
ncold = 1;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % nh04 + nc04
  hca = h(isub); isub = isub + 1;
  pcolor(hca,df04.twci,df04.xi,nh04+nc04)
  shading(hca,'flat')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'n_{i} (df04)';
end
if 1 % nh04
  hca = h(isub); isub = isub + 1;
  pcolor(hca,df04.twci,df04.xi,nh04)
  shading(hca,'flat')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'n_{ih} (df04)';
end
if 1 % nc04
  hca = h(isub); isub = isub + 1;
  pcolor(hca,df04.twci,df04.xi,nc04)
  shading(hca,'flat')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'n_{ic} (df04)';
end
if 1 % nh08 + nc08
  hca = h(isub); isub = isub + 1;
  pcolor(hca,df08.twci,df08.xi,nh08+nc08)
  shading(hca,'flat')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'n_{i} (df08)';
end
if 1 % nh08
  hca = h(isub); isub = isub + 1;
  pcolor(hca,df08.twci,df08.xi,nh08)
  shading(hca,'flat')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'n_{ih} (df08)';
end
if 1 % nc08
  hca = h(isub); isub = isub + 1;
  pcolor(hca,df08.twci,df08.xi,nc08)
  shading(hca,'flat')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'n_{ic} (df08)';
end

hlinks = linkprop(h,{'CLim','XLim','YLim'});
colormap(pic_colors('thermal'))
h(1).XLim = [0 max([df04.twci(end) df08.twci(end)])]