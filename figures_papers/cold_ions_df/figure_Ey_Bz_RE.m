% figures
df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
df04n = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/data_h5/fields.h5');
bs = PIC('/Volumes/Fountain/Data/PIC/baselinev4/data_h5/fields.h5');

df04n_nxline_z1 = interp(df04n,df04n.x_xline,df04n.z_xline+1,df04n.twci,'ne');
df04n_nxline_z0 = interp(df04n,df04n.x_xline,df04n.z_xline+0,df04n.twci,'ne');

df04_nxline_z1 = interp(df04,df04.x_xline,df04.z_xline+1,df04.twci,'ne');
df04_nxline_z0 = interp(df04,df04.x_xline,df04.z_xline+0,df04.twci,'ne');

bs_nxline_z1 = interp(bs,bs.x_xline,bs.z_xline+1,bs.twci,'ne');
bs_nxline_z0 = interp(bs,bs.x_xline,bs.z_xline+0,bs.twci,'ne');

%% Figure
zlim = [-0.5 0.5];
xlim = [100 300];
pic = bs.zlim(zlim).xlim(xlim);
pic_nxline_z1 = bs_nxline_z1;
pic_nxline_z0 = bs_nxline_z0;

%pic1 = df04.zlim(zlim).xlim(xlim);
%pic2 = df04n.zlim(zlim).xlim(xlim);
%A1 = squeeze(mean(pic1.A,2));
%A2 = squeeze(mean(pic2.A,2));
%R1 = reconnection_rate(timesteps/200,'A',A_ts,'E',E_ts(:,:,:,2));
Alev = -25:1:25;
doA = 1;
istepA = 3;

nrows = 1;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
hb = gobjects(0);


if 1 % RE2
  hca = h(isub); isub = isub + 1;  
  hRE = plot(hca,pic.RE,pic.twci);
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'E_R (B_0v_{A0})';
  hca.XLim = [0 0.15];
  
  ax1_pos = hca.Position; % position of first axes
  ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
  %ax2.YTick = [];
  ax2.XLim = [0 0.75];
  ax2.XTick = [0 0.25 0.5 0.75];
  ax2.YTick = [];
  ax2.YLim = hca.YLim;
  ax2.XLabel.String = 'n (n_0)';
  
  hl = line(pic_nxline_z1,pic.twci,'Parent',ax2,'Color','k');
  legend([hRE hl],{'E_R','n(x_x,z_x+1)'},'location','bestoutside')
end
if 0 % Ey1
  hca = h(isub); isub = isub + 1;
  varplot = squeeze(mean(pic.Ey,2));
  pcolor(hca,pic.xi,pic.twci,varplot')
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'E_y';
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    plA = pic.A;
    %contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',Alev,'k');        
    contour(hca,pic.xi(1:istepA:end),pic.twci,plA(1:istepA:end,:)',Alev,'color',[0 0 0])
    hold(hca,'off')
  end
end
if 1 % Bz1
  hca = h(isub); isub = isub + 1;
  varplot = squeeze(mean(pic.Bz,2));
  pcolor(hca,pic.xi,pic.twci,varplot')
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'B_z';
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    plA = squeeze(mean(pic.A,2));
    %contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',Alev,'k');        
    contour(hca,pic.xi(1:istepA:end),pic.twci,plA(1:istepA:end,:)',Alev,'color',[0 0 0])
    hold(hca,'off')
  end
  hca.CLim = [-0.7 0.7];
  colormap(hca,pic_colors('blue_red'))
end
if 1 % nall
  hca = h(isub); isub = isub + 1;
  varplot = squeeze(mean(pic.n([1 3]),2));
  pcolor(hca,pic.xi,pic.twci,varplot')
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'n_{ic}';
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    plA = squeeze(mean(pic.A,2));
    %contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',Alev,'k');        
    contour(hca,pic.xi(1:istepA:end),pic.twci,plA(1:istepA:end,:)',Alev,'color',[0 0 0])
    hold(hca,'off')
  end
  hca.CLim = [0 0.6];
  colormap(hca,pic_colors('candy'))
end
if 0 % nhot
  hca = h(isub); isub = isub + 1;
  varplot = squeeze(mean(pic.n(1),2));
  pcolor(hca,pic.xi,pic.twci,varplot')
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'n_{ih}';
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    plA = squeeze(mean(pic.A,2));
    %contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',Alev,'k');        
    contour(hca,pic.xi(1:istepA:end),pic.twci,plA(1:istepA:end,:)',Alev,'color',[0 0 0])
    hold(hca,'off')
  end
  hca.CLim = [0 1.6];
  colormap(hca,pic_colors('candy'))
end
if 0 % ncold
  hca = h(isub); isub = isub + 1;
  varplot = squeeze(mean(pic.n(3),2));
  pcolor(hca,pic.xi,pic.twci,varplot')
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'n_{ic}';
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    plA = squeeze(mean(pic.A,2));
    %contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',Alev,'k');        
    contour(hca,pic.xi(1:istepA:end),pic.twci,plA(1:istepA:end,:)',Alev,'color',[0 0 0])
    hold(hca,'off')
  end
  hca.CLim = [0 0.6];
  colormap(hca,pic_colors('candy'))
end
if 0 % Ey2
  hca = h(isub); isub = isub + 1;
  varplot = squeeze(mean(pic.Ey,2));
  pcolor(hca,pic.xi,pic.twci,varplot')
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'E_y';
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';
end
if 0 % Bz2
  hca = h(isub); isub = isub + 1;
  varplot = squeeze(mean(pic.Bz,2));
  pcolor(hca,pic.xi,pic.twci,varplot')
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'B_z';
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';
end



for ip = 2:npanels
  h(ip).YTickLabel = [];
  h(ip).YLabel = [];
end

fig = gcf;
allh = findobj(fig.Children,'type','axes');
for ip = 1:numel(allh)
  allh(ip).Position(2) = 0.21;
  allh(ip).Position(4) = 0.5;
  allh(ip).Position(3) = 0.20;
  %h(ip).XTick = 0:50:500;
  allh(ip).YLim = [0 240];
  allh(ip).XGrid = 'on';
  allh(ip).YGrid = 'on';
  allh(ip).FontSize = 14;
  %hb(ip).FontSize = 14;
  
end