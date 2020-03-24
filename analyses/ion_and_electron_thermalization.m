
twpe = 7000;
xlim = [100 200];
zlim = [0 15];
pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);
pic_inflow = df04.twpelim(twpe).xlim([10 20 ]).zlim([20 24]);

iIon = [1 3 5];
iEle = [2 4 6];
iIonCold = [3 5];
iEleCold = [4 6];
iIonHot = [1];
iEleHot = [2];


nicold = pic.n(iIonCold);
necold = pic.n(iEleCold);
picold = pic.p(iIonCold);
pecold = pic.p(iEleCold);

ticold = picold./nicold;
tecold = pecold./necold;


nicold_inflow = pic_inflow.n(iIonCold);
necold_inflow = pic_inflow.n(iEleCold);
picold_inflow = pic_inflow.p(iIonCold);
pecold_inflow = pic_inflow.p(iEleCold);

ticold_inflow = mean(mean(picold_inflow./nicold_inflow));
tecold_inflow = mean(mean(pecold_inflow./necold_inflow));

disp('Done.')

%% Figure
% Figure
nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
isMap = [];
isDist = [];
isDistIon = [];
isDistEle = [];


if 1 % Ti cold
  hca = h(isub); isub = isub + 1;
  isMap = [isMap isub-1];
  pcolor(hca,pic.xi,pic.zi,ticold')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 'z (d_{i0})';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'T_{i,cold} (...)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = 0.5*[-1 1];
  ds.plot_boxes(hca);
end
if 1 % Te cold
  hca = h(isub); isub = isub + 1;
  isMap = [isMap isub-1];
  pcolor(hca,pic.xi,pic.zi,tecold'-tecold_inflow)
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 'z (d_{i0})';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'T_{e,cold}-T_{e,cold}^{inflow} (...)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = 0.5*[-1 1];
  ds.plot_boxes(hca);
end
if 0 % ni cold
  hca = h(isub); isub = isub + 1;
  isMap = [isMap isub-1];
  pcolor(hca,pic.xi,pic.zi,pic.n(iIonCold)')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 'z (d_{i0})';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_{i,cold} (v_{A0}B_0)';
  colormap(hca,pic_colors('waterfall'))
  ds.plot_boxes(hca);
  clim = hca.CLim;
  if 1 % A    
    hold(hca,'on')
    A = pic.twpelim(twpe).A;
    iplx = 1:5:pic.nx;
    iplz = 1:5:pic.nz;
    contour(hca,pic.xi(iplx),pic.zi(iplz),A(iplx,iplz)',-25:1:0,'k')
    hold(hca,'off')
  end
  hca.CLim = clim;
end