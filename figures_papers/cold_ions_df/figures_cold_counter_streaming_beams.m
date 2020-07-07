%% Figure 1, reconnection rate vs density and temperature
pic = nobg;

zlim = [-0.1 0.1];
xlim = diff(pic.xi([1 end]))/2 + [-130 130];
twcilim = [100 300];
pic_Bz = pic.zlim(zlim).xlim(xlim).twcilim(twcilim); 
Bz_tx = pic_Bz.Bz; 
Bz_tx = squeeze(mean(Bz_tx,2));

n_tx = pic_Bz.ni; 
n_tx = squeeze(mean(n_tx,2));

if 0% Calculate density at and above xline
  % reconnection_rate_vs_density
  pic_Bxline_z1 = interp(pic,pic.x_xline,pic.z_xline+1,pic.twci,'Bx');
  pic_nxline_z1 = interp(pic,pic.x_xline,pic.z_xline+1,pic.twci,'ne');
  pic_nxline_z0 = interp(pic,pic.x_xline,pic.z_xline+0,pic.twci,'ne');
end

nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % bz(x,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic_Bz.twci,pic_Bz.xi,Bz_tx)
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('B_z(x,t,%.1f<z<%.1f)',pic_Bz.zi(1),pic_Bz.zi(end));
  hca.XLabel.String = 't (\omega_{ci}^{-1})';
  hca.YLabel.String = 'x (d_i)';
end
if 1 % n(x,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic_Bz.twci,pic_Bz.xi,n_tx)
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('n(x,t,%.1f<z<%.1f)',pic_Bz.zi(1),pic_Bz.zi(end));
  hca.XLabel.String = 't (\omega_{ci}^{-1})';
  hca.YLabel.String = 'x (d_i)';
end
%% Figure 1
pic = nobg;
% Calculate density at and above xline
pic_Bxline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'Bx');
pic_nxline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'ne');
pic_nxline_z0 = pic.interp(pic.x_xline,pic.z_xline+0,pic.twci,'ne');
pic_tixline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'ti');
% Figure
zlim = [-0.5 0.5];
xlim = [100 300];
pic = pic.zlim(zlim).xlim(xlim);

%pic1 = df04.zlim(zlim).xlim(xlim);
%pic2 = df04n.zlim(zlim).xlim(xlim);
%A1 = squeeze(mean(pic1.A,2));
%A2 = squeeze(mean(pic2.A,2));
%R1 = reconnection_rate(timesteps/200,'A',A_ts,'E',E_ts(:,:,:,2));
Alev = -25:1:25;
doA = 0;
istepA = 3;

nrows = 1;
ncols = 4;
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
  
  %ax2.XTick = [0:0.2:2];
  
  ax2.YTick = [];
  ax2.YLim = hca.YLim;
  ax2.XLabel.String = 'n (n_0)';
  
  hl = line(pic_nxline_z1,pic.twci,'Parent',ax2,'Color','k');
  legend([hRE hl],{'E_R','n(x_x,z_x+1)'},'location','best','box','off')
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
  varplot = squeeze(mean(pic.ni,2));
  pcolor(hca,pic.xi,pic.twci,varplot')
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'n_{i}';
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
if 1 % tall
  hca = h(isub); isub = isub + 1;
  varplot = squeeze(mean(pic.ti,2));
  pcolor(hca,pic.xi,pic.twci,varplot')
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'T_{i}';
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

fig = gcf; h = findobj(fig.Children,'type','axes');
hlinks = linkprop(h,{'YLim'});
h(1).YLim = pic.twci([1 end]);

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

%% Figure 2. Reduced distributions f(x,vx), f(x,vy), f(x,vz) for a few different times