%no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

doLoad = 0;
if doLoad
twpe_map = 24000;
xlim = mean(no02m.xi) + [-50 50];
zlim = [-10 10];
pic = no02m.twpelim(twpe_map).xlim(xlim).zlim(zlim);
P = pic.p([3 5]); % pressure of cold ions
A = pic.A; % magnetic vector potential
end

h(1) = subplot(1,3,[1 2]); % color map
h(2) = subplot(1,3,3); % energy spectra
isub = 1;

if 1 % color map of temperature
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.zi,squeeze(P)')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  colormap(hca,pic_colors('thermal'))
  
  hold(hca,'on')
  clim = hca.CLim;
  contour(hca,pic.xi,pic.zi,squeeze(A)',[-30:30],'color',0.8*[1 1 1])
  hca.CLim = clim;
  hold(hca,'off')
  
  hca.CLim = [0 0.2];
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
end