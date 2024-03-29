%% Figure 1, change of inflow parameters, compression of current sheet
%% Figure 1, prepare data
pic = no02m;

x0 = no02m.xi(end)/2;
twpe1 = 1000; twpe1 = no02m.twpelim(twpe1).twpe; twci1 = no02m.twpelim(twpe1).twci;
twpe2 = 24000; twpe2 = no02m.twpelim(twpe2).twpe; twci2 = no02m.twpelim(twpe2).twci;

%pic_Bxline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'Babs');
pic_Bxline_z1_ = pic.get_points(pic.x_xline,pic.z_xline+1,pic.twci,[-0.1 0.1],'Babs');
%pic_nxline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'ni');
pic_nxline_z1_ = pic.get_points(pic.x_xline,pic.z_xline+1,pic.twci,[-0.1 0.1],'ni');
%pic_nhot_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'n(1)');
pic_nhot_z1_ = pic.get_points(pic.x_xline,pic.z_xline+1,pic.twci,[-0.1 0.1],'n(1)');
%pic_ncold_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'n(3)');
pic_ncold_z1_ = pic.get_points(pic.x_xline,pic.z_xline+1,pic.twci,[-0.1 0.1],'n(3)');
%pic_tixline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'ti');
%pic_vA_z1 = squeeze(pic_Bxline_z1./sqrt(pic_nxline_z1));
pic_vA_z1_ = squeeze(pic_Bxline_z1_./sqrt(pic_nxline_z1_));



zlim = [-0.5 0.5];
xlim = x0 + 40*[-1 1];
pic = no02m.zlim(zlim).xlim(xlim);
pic_Bz_tx = squeeze(mean(pic.Bz,2));
pic_A_tx = squeeze(mean(pic.A,2));
%% Figure 1, plot
colors = pic_colors('matlab');
Alev = -25:1:25;
doA = 0;
istepA = 3;

nrows = 5;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
hb = gobjects(0);


if 1 % zcut at early times
  hca = h(isub); isub = isub + 1;
  pic = no02m.twpelim(twpe1).xlim(x0+[-1 1]);
  hn1 = plot(hca,pic.zi,mean(pic.n(1),1),'color',colors(2,:));
  hold(hca,'on')
  hn2 = plot(hca,pic.zi,smooth(mean(pic.n([3 5]),1),5),'color',colors(1,:));
  hB = plot(hca,pic.zi,mean(pic.Babs,1),'k');
  hold(hca,'off')
  hca.XLim = [-10 10];
  if 0 % Add axes for A
    %hold(hca,'on')
    ax1 = hca;
    ax1_pos = ax1.Position; % position of first axes
    ax2 = axes('Position',ax1_pos,...
      'XAxisLocation','top',...
      'YAxisLocation','right',...
      'Color','none');
    ax2all(isub-1) = ax2;
    hA = line(ax2,pic.zi,mean(pic.A,1),'color',[0 0 0],'linestyle','--');
    ax2.XLim = hca.XLim;
    ax2.XTick = [];
    ax2.YLabel.String = 'A (B_0d_i)';
  end
  hca.XLabel.String = 'z (d_i)';
  hca.YLabel.String = 'n, B';  
  %legend([hn1,hn2,hB,hA],{'n_{hot}','n_{cold}','|B|','A_y'},'location','northwest','Box','off')
  %pos = hca.Position;
  legend([hn1,hn2,hB],{'n_{hot}','n_{cold}','|B|'},'location','east','Box','off')
  %legend([hn1,hn2,hB],{'n_{hot}','n_{cold}','|B|'},'location','northoutside','Box','off','Orientation','horizontal')
  %hca.Position = pos;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  irf_legend(hca,{['inflow at t\omega_{ci}= ' num2str(twci1)]},[0.02 0.9],'fontsize',14,'color',[0 0 0])
end

if 0 % zcut at, vs A 
  hca = h(isub); isub = isub + 1;
  pic = no02m.twpelim(twpe1).xlim(x0+[-1 1]);
  A = mean(pic.A,1);
  plot(hca,A,mean(pic.n(1),1),'color',colors(2,:))
  hold(hca,'on')
  plot(hca,A,smooth(mean(pic.n([3 5]),1),5),'color',colors(1,:))
  plot(hca,A,mean(pic.Babs,1),'k')
  hold(hca,'off')
  %hca.XLim = [-10 10];  
  hca.XLim = [min(A) max(A)];
  hca.XLabel.String = 'A (B_0d_i)';
  hca.YLabel.String = 'n, B';
  %legend(hca,{'n_{hot}','n_{cold}','|B|'},'location','northwest','box','off')
end
if 0 % xcut vs A
  hca = h(isub); isub = isub + 1;
  pic = no02m.twpelim(twpe2).xgrid(2:no02m.nx-1).zlim(0+[-1 1]);
  A = mean(pic.A,2);
  nh = mean(pic.n(1),2);
  plot(hca,A,nh,'color',colors(2,:))
  hold(hca,'on')
  plot(hca,A,mean(pic.n([3 5]),2),'color',colors(1,:))
  plot(hca,A,mean(pic.Babs,2),'k')
  hold(hca,'off')  
  %hca.XLim = x0 + [-50 50];
  hca.XLim = [min(A) max(A)];
  hca.XLabel.String = 'A (B_0d_i)';
  hca.YLabel.String = 'n, B';
  %legend(hca,{'n_{hot}','n_{cold}','|B|'},'location','northwest','box','off')
end
if 0 % ni, Bx, vA, RE
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic_nxline_z1,pic.twci,pic_Bxline_z1,pic.twci,pic_vA_z1,pic.twci,pic.RE);
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'n, B_x';
  %hca.XLim = [0 0.4];
  legend(hca,{'n(x_x,z_x+1)','B_x(x_x,z_x+1)'},'location','southwest','box','off')
end
if 0 % ni, Bx, vA
  hca = h(isub); isub = isub + 1;  
  [AX,H1,H2] = plotyy(hca,pic.twci,[pic_nxline_z1,pic_Bxline_z1]',pic.twci,smooth(pic_vA_z1,1));  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'n, B_x';
  AX(2).YLabel.String = 'v_A';
  AX(2).YColor = [0 0 0];
  H2.LineStyle = '-.';
  H2.Color = [0 0 0];
  %hca.XLim = [0 0.4];
  legend(hca,{'n(x_X,z_X+1)','B_x(x_X,z_X+1)','v_A(x_X,z_X+1)'},'location','northeast','box','off')
end
if 0 % ni, tot/hot/cold, Bx
  hca = h(isub); isub = isub + 1;  
  H1 = plot(hca,pic.twci,[pic_nxline_z1,pic_nhot_z1,pic_ncold_z1,pic_Bxline_z1]');  
  H1(1).Color = colors(4,:);
  H1(2).Color = colors(2,:);
  H1(3).Color = colors(1,:);
  H1(4).Color = colors(5,:);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'n, |B|';
  AX(2).YLabel.String = 'v_A';
  AX(2).YColor = [0 0 0];
  H2.LineStyle = '-.';
  H2.Color = [0 0 0];
  %hca.XLim = [0 0.4];
  legend(hca,{'n(x_X,z_X+1)','n_{hot}(x_X,z_X+1)','n_{cold}(x_X,z_X+1)','B_x(x_X,z_X+1)'},'location','northeast','box','off')
end
pic = no02m.zlim(zlim).xlim(xlim);
if 1 % ni, tot/hot/cold, Bx, vA
  hca = h(isub); isub = isub + 1; 
  
  [AX,H1,H2] = plotyy(hca,pic.twci,[pic_nxline_z1_,pic_nhot_z1_,pic_ncold_z1_,pic_Bxline_z1_]',pic.twci,smooth(pic_vA_z1_,1));  
%  [AX,H1,H2] = plotyy(hca,pic.twci,[pic_nxline_z1,pic_nhot_z1,pic_ncold_z1,pic_Bxline_z1]',pic.twci,smooth(pic_vA_z1_,1));  
  H1(1).Color = colors(4,:);
  H1(2).Color = colors(2,:);
  H1(3).Color = colors(1,:);
  H1(4).Color = [0 0 0];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'n, |B|';
  AX(2).YLabel.String = 'v_A';
  AX(2).YColor = [0 0 0];
  H2.LineStyle = '-.';
  H2.Color = colors(5,:);
  AX(1).YLim = [0 0.75];
  AX(1).YTick = [0 0.25 0.5 0.75];
  AX(2).YLim = [0 1.5];
  AX(2).YTick = [0 0.25 0.5 0.75]*2;
  AX(1).XLim = [5 125];
  AX(2).XLim = [5 125];
  %hca.XLim = [0 0.4];
  %legend(hca,{'n(x_X,z_X+1)','B_x(x_X,z_X+1)','v_A(x_X,z_X+1)'},'location','northeast','box','off')
  %hleg = legend(hca,{'n(x_X,z_X+1)','n_{hot}(x_X,z_X+1)','n_{cold}(x_X,z_X+1)','B_x(x_X,z_X+1)','v_A(x_X,z_X+1)'},'location','northwest','box','off','orientation','horizontal');
  %hleg = legend(hca,{'n','n_{hot}','n_{cold}','|B|','v_A'},'location','northwest','box','off','orientation','horizontal');
  hleg = legend([H1(1) H2],{'n','v_A'},'location','southwest','box','off');
  irf_legend(hca,{'inflow values 1d_i above X line'},[0.02 0.98],'fontsize',14,'color',[0 0 0])
end
if 0 % ni, Bx
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic_nxline_z1,pic.twci,pic_Bxline_z1);
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'n, B_x';
  %hca.XLim = [0 0.4];
  legend(hca,{'n(x_x,z_x+1)','B_x(x_x,z_x+1)'},'location','southwest','box','off')
end
if 0 % vA
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic_vA_z1);
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'v_x';
  %hca.XLim = [0 0.4];
  legend(hca,{'v_A(x_x,z_x+1)'},'location','southwest','box','off')
end
if 1 % RE
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic.RE);
  h_.Color = colors(3,:);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'E_R';
  %hca.XLim = [0 0.4];
  %legend(hca,{'v_A(x_x,z_x+1)'},'location','southwest','box','off')
  hca.XLim = [5 125];
  hca.YLim = [0 0.17];
end
if 1 % Bz(x,t)
  hca = h(isub); isub = isub + 1;  
  pic = no02m.zlim(zlim).xlim(xlim);
  h_ = pcolor(hca,pic.twci,pic.xi,pic_Bz_tx);
  shading(hca,'flat')
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x (d_i)';
  colormap(hca,pic_colors('blue_red'))
  clim = hca.CLim;
  pos = hca.Position;
  hcb = colorbar('peer',hca);
  hcb.Title.String = 'B_z';
  if 1 % A
    hold(hca,'on')
    contour(hca,pic.twci,pic.xi,pic_A_tx,0:0.5:25,'color',[0 0 0]);
    plot(hca,pic.twci,pic.x_xline,'linewidth',2,'color',[0 0 0])
    hold(hca,'off')
  end
  colormap(hca,pic_colors('blue_red'))    
  hca.CLim = clim;
  hca.Position = pos;
  %legend(hca,{'v_A(x_x,z_x+1)'},'location','southwest','box','off')
  hca.XLim = [5 125];
  %hca.YLim = [0 0.17];
end
if 1 % xcut at late times
  hca = h(isub); isub = isub + 1;
  pic = no02m.twpelim(twpe2).xgrid(2:no02m.nx-1).zlim(0+[-1 1]);
  hn1 = plot(hca,pic.xi,mean(pic.n(1),2),'color',colors(2,:));
  hold(hca,'on')
  hn2 = plot(hca,pic.xi,mean(pic.n([3 5]),2),'color',colors(1,:));
  hB = plot(hca,pic.xi,mean(pic.Babs,2),'k');
  hold(hca,'off')  
  hca.XLim = x0 + [-50 50];
  if 0 % Add axes for A
    %hold(hca,'on')
    ax1 = hca;
    ax1_pos = ax1.Position; % position of first axes
    ax2 = axes('Position',ax1_pos,...
      'XAxisLocation','top',...
      'YAxisLocation','right',...
      'Color','none');
    ax2all(isub-1) = ax2;
    hA = line(ax2,pic.xi,mean(pic.A,2),'color',[0 0 0],'linestyle','--');
    ax2.XLim = hca.XLim;
    ax2.XTick = [];
    ax2.YLabel.String = 'A (B_0d_i)';
    legend([hn1,hn2,hB],{'n_{hot}','n_{cold}','|B|'},'location','west','Box','off')
  end
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'n, B';
  %legend([hn1,hn2,hB,hA],{'n_{hot}','n_{cold}','|B|','A_y'},'location','north','Box','off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  irf_legend(hca,{['outflow at t\omega_{ci}= ' num2str(twci2)]},[0.20 0.9],'fontsize',14,'color',[0 0 0])
end

legends = {'a)','b)','c)','d)','e)','f)'};
for ip = 1:numel(h)
  irf_legend(h(ip),legends{ip},[-0.1 1.00],'fontsize',14)  
end

hall = findobj(gcf,'type','axes'); hall = hall(end:-1:1);
for ip = 1:numel(hall)
  hall(ip).FontSize = 12;
  hall(ip).XGrid = 'on';
  hall(ip).YGrid = 'on';
end

compact_panels(h([2 3 4]))

% Set up panels
if 1
  %%
  h(1).Position([2 4]) = [0.83 0.15];
  h(2).Position([2 4]) = [0.61 0.15];
  h(3).Position([2 4]) = [0.505 0.10];
  h(4).Position([2 4]) = [0.3 0.20];
  h(5).Position([2 4]) = [0.08 0.15];
else
  h(2).Position(4) = h(1).Position(4);
  h(2).Position(2) = h(2).Position(2)+0.015;
  h(1).Position(2) = h(1).Position(2)+0.030;
end

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

%% Figure 2, overview at t=120
% what do we want to show here?
varstrs = {'n([3 5])','n(3)';'t([3 5])','t(3)'}';
%varstrs = {'n([3 5])','n(3)','t([3 5])','t(3)'}';
clims = {[0 0.5],[0 0.5],[0 0.5],[0 0.5]};

xlim = [40 165];
xlim = [61 119];
zlim = 0.99*[-10 10];
cmapth = pic_colors('thermal');
cmaps = {cmapth,cmapth,cmapth,cmapth};
cbarlabels = {'n_{i,cold}','n_{i,cold}^{top}','T_{i,cold}','T_{i,cold}^{top}'};
pic = no02m.twpelim(24000).xlim(xlim).zlim(zlim);
h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps,'cbarlabels',cbarlabels);

for ip = 1:numel(h)
  h(ip).FontSize = 12;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end
hl = findobj(gcf,'Type','Contour');
c_eval('hl(?).Color = 0.5 + [0 0 0];',1:numel(hl))

legends = {'a)','b)','c)','d)','e)','f)'};
if 1 % 2x2
  %%
  compact_panels(h,0.002,0.002)
  for ip = 1:numel(h)
    irf_legend(h(ip),{[legends{ip} ' ' cbarlabels{ip}]},[0.8 0.98],'color',[1 1 1],'fontweight','bold','fontsize',12)
  end
  hcbar = findobj(gcf,'Type','ColorBar');
  delete(hcbar(2:end))
  hcbar(1).Position = [0.85 0.11 0.02 0.815];
  hcbar(1).YLabel.String = 'n, T';
  h(3).YTickLabels = [];
  h(3).YLabel.String = [];
  h(4).YTickLabels = [];
  h(4).YLabel.String = [];  
  for ip = 1:numel(h)
    h(ip).Position(2) = h(ip).Position(2) + 0.05;
  end
  hcbar(1).Position(2) = hcbar(1).Position(2) + 0.05;
  hcbar(1).FontSize = 12;
  h(1).Title.String = '';
end
%% Figure 2, ALT2, overview at t=120, n,p
% what do we want to show here?
varstrs = {'n(1)','n([3 5])';'p(1)','p([3 5])'};
%varstrs = {'n([3 5])','n(3)','t([3 5])','t(3)'}';
clims = {[0 0.5],[0 0.5],[0 0.5],[0 0.5]};

xlim = [40 165];
xlim = [61 119];
zlim = 0.99*[-10 10];
cmapth = pic_colors('thermal');
cmaps = {cmapth,cmapth,cmapth,cmapth};
cbarlabels = {'n_{ih}','n_{ic}';'P_{ih}','P_{ic}'};
%cbarlabels = {'n_{i}','n_{i}';'P_{i}','P_{i}'};
pic = no02m.twpelim(24000).xlim(xlim).zlim(zlim);
h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps,'cbarlabels',cbarlabels);

for ip = 1:numel(h)
  h(ip).FontSize = 12;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end
hl = findobj(gcf,'Type','Contour');
c_eval('hl(?).Color = 0.5 + [0 0 0];',1:numel(hl))

legends = {'a)','b)','c)','d)','e)','f)'};
if 1 % 2x2
  %%
  compact_panels(h,0.002,0.002)
  for ip = 1:numel(h)
    irf_legend(h(ip),{[legends{ip} ' ' cbarlabels{ip}]},[0.8 0.98],'color',[1 1 1],'fontweight','bold','fontsize',12)
  end
  delete(findobj(gcf,'Type','ColorBar'))
  hcbar(1) = colorbar('peer',h(3),'location','eastoutside');
  hcbar(2) = colorbar('peer',h(4),'location','eastoutside');
  hcbar(1).YLabel.String = 'n';
  hcbar(2).YLabel.String = 'P';
  hcbar = findobj(gcf,'Type','ColorBar');
  hcbar(1).FontSize = 12;
  hcbar(2).FontSize = 12;
  
  h(3).YTickLabels = [];
  h(3).YLabel.String = [];
  h(4).YTickLabels = [];
  h(4).YLabel.String = [];    
  %c_eval('h(?).Position(2) = h(?).Position(2) + 0.05;',[2 4])
  c_eval('h(?).Position(2) = h(?).Position(2) + 0.04;',[1 2 3 4])
  c_eval('h(?).CLim = [0 0.499];',[1 3])
  c_eval('h(?).CLim = [0 0.099];',[2 4])
  %c_eval('h(?).Position(4) = h(?).Position(4) - 0.05;',[1 3])
%   for ip = 1:numel(h)
%     %h(ip).Position(2) = h(ip).Position(2) + 0.05;
%     h(ip).Position(4) = h(ip).Position(4) - 0.05;
%   end
  %hcbar(1).Position(2) = hcbar(1).Position(2) + 0.05;
  h(1).Title.String = '';
  compact_panels(h,0.002,0.002)
end
%% Figure 2, ALT3, overview at t=120, n,p,t
% what do we want to show here?
varstrs = {'n(1)','n([3 5])';'p(1)','p([3 5])';'t(1)','t([3 5])'};
%varstrs = {'n([3 5])','n(3)','t([3 5])','t(3)'}';
clims = {[0 0.5],[0 0.5],[0 0.5],[0 0.5],[0 0.5],[0 0.5]};

xlim = [40 165];
xlim = [61 119];
zlim = 0.99*[-10 10];
cmapth = pic_colors('thermal');
cmaps = {cmapth,cmapth,cmapth,cmapth,cmapth,cmapth};
cbarlabels = {'n_{ih}','n_{ic}';'P_{ih}','P_{ic}';'T_{ih}','T_{ic}'};
%cbarlabels = {'n_{i}','n_{i}';'P_{i}','P_{i}'};
pic = no02m.twpelim(24000).xlim(xlim).zlim(zlim);
h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps,'cbarlabels',cbarlabels);

for ip = 1:numel(h)
  h(ip).FontSize = 12;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end
hl = findobj(gcf,'Type','Contour');
c_eval('hl(?).Color = 0.5 + [0 0 0];',1:numel(hl))

legends = {'a)','b)','c)','d)','e)','f)'};
legends = {'a)','c)','e)','b)','d)','f)'};
if 1 % 3x2
  %%
  compact_panels(h,0.002,0.002)
  for ip = 1:numel(h)
    %irf_legend(h(ip),{[legends{ip} ' ' cbarlabels{ip}]},[0.8 0.98],'color',[1 1 1],'fontweight','bold','fontsize',12)
    irf_legend(h(ip),{[' ' cbarlabels{ip}]},[0.98 0.98],'color',[1 1 1],'fontweight','bold','fontsize',12)
  end
  delete(findobj(gcf,'Type','ColorBar'))
  hcbar(1) = colorbar('peer',h(1,2),'location','eastoutside');
  hcbar(2) = colorbar('peer',h(2,2),'location','eastoutside');
  hcbar(3) = colorbar('peer',h(3,2),'location','eastoutside');
  hcbar(1).YLabel.String = 'n';
  hcbar(2).YLabel.String = 'P';
  hcbar(3).YLabel.String = 'T';
  %hcbar = findobj(gcf,'Type','ColorBar');
  hcbar(1).FontSize = 12;
  hcbar(2).FontSize = 12;
  hcbar(3).FontSize = 12;
  
  h(1,2).YTickLabels = [];
  h(1,2).YLabel.String = [];
  h(2,2).YTickLabels = [];
  h(2,2).YLabel.String = [];  
  h(3,2).YTickLabels = [];
  h(3,2).YLabel.String = [];    
  %c_eval('h(?).Position(2) = h(?).Position(2) + 0.05;',[2 4])
  %c_eval('h(?).Position(2) = h(?).Position(2) + 0.04;',[1 2 3 4])
  c_eval('h(1,?).CLim = [0 0.499];',1:2)
  c_eval('h(2,?).CLim = [0 0.099];',1:2)
  c_eval('h(3,?).CLim = [0 0.699];',1:2)
  %c_eval('h(?).Position(4) = h(?).Position(4) - 0.05;',[1 3])
%   for ip = 1:numel(h)
%     %h(ip).Position(2) = h(ip).Position(2) + 0.05;
%     h(ip).Position(4) = h(ip).Position(4) - 0.05;
%   end
  %hcbar(1).Position(2) = hcbar(1).Position(2) + 0.05;
  h(1).Title.String = '';
  compact_panels(h,0.002,0.002)
end

%% Figure 3, reduced distributions
%ds100 = PICDist('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');
%% Figure 3, prepare data
twpe = 24000; xlim = [50 155]; zlim = [-15 15];
sep = no02m.twpelim(twpe).separatrix_location;
for zpick = [0 2 4]
  ds = ds100.twpelim(twpe).zfind(zpick).xlim(xlim).findtag({'line horizontal'});
  
  xdist = (ds.xi1{1}+ds.xi2{1})/2;
  zdist = (ds.zi1{1}+ds.zi2{1})/2;
  tdist = repmat(twpe,size(xdist));
  vExBx_tmp = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,[-0.25 0.25],'vExBx'); eval(sprintf('vExBx_z%g = vExBx_tmp;',zpick))
  vExBy_tmp = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,[-0.25 0.25],'vExBy'); eval(sprintf('vExBy_z%g = vExBy_tmp;',zpick))
  vExBz_tmp = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,[-0.25 0.25],'vExBz'); eval(sprintf('vExBz_z%g = vExBz_tmp;',zpick))
  pic_lim = no02m.xlim(xlim).zlim(zlim).twpelim(twpe);
  pic = no02m.twpelim(twpe);
  Bx_ = pic.Bx;
  By_ = pic.By;
  Bz_ = pic.Bz;
  Bx = interpfield(pic.xi,pic.zi,Bx_,xdist,zdist); 
  By = interpfield(pic.xi,pic.zi,By_,xdist,zdist); 
  Bz = interpfield(pic.xi,pic.zi,Bz_,xdist,zdist); 
  fred5_tmp = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred5_z%g = fred5_tmp;',zpick))
  fred3_tmp = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred3_z%g = fred3_tmp;',zpick))
  fred35_tmp = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred35_z%g = fred35_tmp;',zpick))
  %fred46_tmp = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred46_z%g = fred46_tmp;',zpick))    
end
%[xDF,vDF,aDF,BDF] = no02m.xva_df;
%% Figure 3, plot
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 3;
ncols = 3;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 0; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.0;
doPhi = 0; colorPhi = [0.5 0.5 0];
doSep = 1;

cmap_dist = pic_colors('waterfall');

freds = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
%freds = {fred35_z2,fred35_z2,fred35_z2};
%freds = {fred5_z4,fred5_z4,fred5_z4;fred5_z2,fred5_z2,fred5_z2;fred5_z0,fred5_z0,fred5_z0}';
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)'}; 

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.fv' labstr]);
    pcolor(hca,fred.x,fred.v,log10(fredplot)')
    shading(hca,'flat')
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy4')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %irf_legend(hca,{sprintf('%s z = %g',legends{ifred},unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('%s f(v_%s,z=%g)',legends{ifred},labstr,unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 12;
    if 0*doE
      hold(hca,'on')
      plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
      hold(hca,'off')
    end
    if 0*doV
      hold(hca,'on')
      plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
      hold(hca,'off')
    end
    if doExB
      hold(hca,'on')
      %xx = eval(['x_z' num2str(unique(fred.z))]);
      xx = fred.x;
      vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
      hExB = plot(hca,xx,vv,'color',colorExB,'linewidth',0.5);
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end    
    if doSep
      hold(hca,'on')   
      if 1           
        zz = unique(fred.z);
        [PKS,LOCS] = findpeaks(-abs(sep.z-zz),'sort','descend');
        [~,iz] = min(abs(sep.z-zz));
        if zz == 0
          xx = sep.x(LOCS(1));
          hSep = plot(hca,xx(1)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
        else
          xx = sep.x(LOCS(1:2));
          hSep = plot(hca,xx(1)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
          hSep = plot(hca,xx(2)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
        end
      elseif 0
        xx = unique(fred.x);
        [~,ix] = min(abs(sep.x-xx));
        zz = sep.z(ix);
      end      
      hold(hca,'off')      
    end
  end
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
hl = findobj(h(1),'type','line');
legend(hl([3 2]),{'v_{ExB}','separatrix'},'box','off','location','south')
%legend([hExB,hSep],{'v_{ExB}','separatrix'})
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
hl = findobj(h(8),'type','line'); delete(hl(2));
hl = findobj(h(9),'type','line'); delete(hl(2));

hcb = colorbar('peer',h(9));
hcb.Position = [0.91 0.115 0.01 0.810];
hcb.YLabel.String = 'log_{10} f_{ic}(x,v_{x,y,z})';
c_eval('h(?).YLabel.String = []; h(?).YTickLabel = [];',[2 3 5 6 8 9])
c_eval('h(?).YLabel.String = ''v'';',[1 4 7])
compact_panels(h,0.002,0.002)
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
h(1).CLim = 0.99*[-4 2];
h(1).CLim = 0.99*[-6 0];
h(1).YLim = 0.99*4*[-1 1];

if 1 % extra markings and annotations 
  %%
  % Add DF locations
  itt = no02m.twpelim(twpe).it;  
  for ip = 7:9
    hold(h(ip),'on')
    hDF = plot(h(ip),xDF(1,itt)*[1 1],h(ip).YLim,'k--');
    plot(h(ip),xDF(2,itt)*[1 1],h(ip).YLim,'k--');
    hold(h(ip),'off')
    if ip == 7
      hlegdf = legend(hDF,{'DF'},'location','southeast','box','off');
      hlegdf.Position(1) = 0.32;
    end    
  end
  % inflow/outflow
  text(h(1),102,2.5,'inflow','fontsize',13,'horizontalalignment','center')
  text(h(1),128,2.5,'outflow','fontsize',13,'horizontalalignment','center')
  text(h(1),75,2.5,'outflow','fontsize',13,'horizontalalignment','center')
  
  annotation('textarrow',[0.56 0.544],[0.14 0.14],'String','main X line','fontsize',13,'fontweight','light')
  annotation('textarrow',[0.83 0.84],[0.74 0.76],'String',{'thermalization of','cold inbound ions','at separatrices'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.81 0.825],[0.46 0.48],'String',{'inside separatrix,','inbound ions','are thermalized'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.86 0.87],[0.46 0.50],'String',{'deeper within','exhaust, inbound','ions are still cold'},'fontsize',13,'fontweight','light','horizontalalignment','left')
  annotation('textarrow',[0.86 0.87],[0.87 0.84],'String',{'heated exhaust ions'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.78 0.78],[0.82 0.79],'String',{'cold inflow'},'fontsize',13,'fontweight','light')  
  annotation('textarrow',[0.69 0.675],[0.58 0.54],'String',{'cold outbound','beam close to DF'},'fontsize',13,'fontweight','light','horizontalalignment','left')
  annotation('textarrow',[0.46 0.46],[0.20 0.23],'String',{'inbound ions','with v_y<0'},'fontsize',13,'fontweight','light')
  %annotation('textarrow',[0.43 0.43],[0.29 0.27],'String',{'inbound ions','with v_y>0'},'fontsize',13,'fontweight','light')
  %annotation('textarrow',[0.62 0.63],[0.84 0.81],'String',{'inbound ions','with v_y>0'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.61 0.63],[0.86 0.81],'String',{'inbound ions','with v_y>0'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.60 0.60],[0.73 0.77],'String',{'inbound ions','with v_y<0'},'fontsize',13,'fontweight','light')
  % islands
  text(h(7),99,2.9,'islands','fontsize',13,'horizontalalignment','center')
  annotation('arrow',[0.260 0.260],[0.33 0.28])
  annotation('arrow',[0.235 0.235],[0.33 0.29]) 
  % secondary X lines
  text(h(8),95,3.5,'secondary X lines','fontsize',13,'horizontalalignment','center')
  annotation('arrow',[0.4795 0.4795],[0.33 0.29])
  annotation('arrow',[0.5117 0.5117],[0.33 0.29])
  % counterstreaming beams
  % alt 1
%   text(h(9),119,-3.25,'counter-streaming beams','fontsize',13,'horizontalalignment','center')
%   annotation('textarrow',0.79*[1 1],[0.17 0.21],'string','cold','fontsize',13) 
%   annotation('textarrow',0.825*[1 1],[0.17 0.21],'string','warmer','fontsize',13)
%   annotation('textarrow',0.875*[1 1],[0.17 0.22],'string','cold','fontsize',13)
  % alt 2
  text(h(9),119,-3.35,'counter-streaming beams','fontsize',13,'horizontalalignment','center')
  annotation('textarrow',[0.78 0.79],[0.17 0.21],'string','cold','fontsize',13) 
  annotation('textarrow',[0.82 0.825],[0.17 0.21],'string','warmer','fontsize',13)
  annotation('textarrow',[0.86 0.870],[0.17 0.19],'string','hot','fontsize',13)
  annotation('textarrow',[0.885 0.875],[0.17 0.22],'string','cold','fontsize',13)  
end
%% Figure 3, plot, origin fraction
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 3;
ncols = 3;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 0; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.0;
doPhi = 0; colorPhi = [0.5 0.5 0];
doSep = 1;

cmap_dist = pic_colors('waterfall');

freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
freds_tot = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
%freds = {fred35_z2,fred35_z2,fred35_z2};
%freds = {fred5_z4,fred5_z4,fred5_z4;fred5_z2,fred5_z2,fred5_z2;fred5_z0,fred5_z0,fred5_z0}';
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)'}; 

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;    
    labstr = labstrs{ifred};
    fred = freds{ifred};
    fred_tot = freds_tot{ifred};    
    %['fred_one.fv' labstr]
    fred_one_plot = eval(['fred.fv' labstr]);
    fred_tot_plot = eval(['fred_tot.fv' labstr]);
    fredplot = fred_one_plot./fred_tot_plot;
    pcolor(hca,fred.x,fred.v,fredplot')
    shading(hca,'flat')
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('pasteljet')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %irf_legend(hca,{sprintf('%s z = %g',legends{ifred},unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('%s f(v_%s,z=%g)',legends{ifred},labstr,unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.CLim = [0 1];
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 12;
    if doExB
      hold(hca,'on')
      %xx = eval(['x_z' num2str(unique(fred.z))]);
      xx = fred.x;
      vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
      hExB = plot(hca,xx,vv,'color',colorExB,'linewidth',0.5);
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end    
    if doSep
      hold(hca,'on')   
      if 1           
        zz = unique(fred.z);
        [PKS,LOCS] = findpeaks(-abs(sep.z-zz),'sort','descend');
        [~,iz] = min(abs(sep.z-zz));
        if zz == 0
          xx = sep.x(LOCS(1));
          hSep = plot(hca,xx(1)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
        else
          xx = sep.x(LOCS(1:2));
          hSep = plot(hca,xx(1)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
          hSep = plot(hca,xx(2)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
        end
      elseif 0
        xx = unique(fred.x);
        [~,ix] = min(abs(sep.x-xx));
        zz = sep.z(ix);
      end      
      hold(hca,'off')      
    end
  end
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
hl = findobj(h(1),'type','line');
legend(hl([3 2]),{'v_{ExB}','separatrix'},'box','off','location','north')
%legend([hExB,hSep],{'v_{ExB}','separatrix'})
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
hl = findobj(h(8),'type','line'); delete(hl(2));
hl = findobj(h(9),'type','line'); delete(hl(2));

hcb = colorbar('peer',h(9));
hcb.Position = [0.91 0.13 0.015 0.795];
hcb.YLabel.String = 'f^{top}_{ic}/f^{tot}_{ic}';
c_eval('h(?).YLabel.String = []; h(?).YTickLabel = [];',[2 3 5 6 8 9])
c_eval('h(?).YLabel.String = ''v'';',[1 4 7])
compact_panels(h,0.002,0.002)
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
%h(1).CLim = 0.99*[-4 2];
h(1).YLim = 0.99*4*[-1 1];

if 1 % extra markings and annotations 
  %%
  % Add DF locations
  itt = no02m.twpelim(twpe).it;  
  for ip = 7:9
    hold(h(ip),'on')
    hDF = plot(h(ip),xDF(1,itt)*[1 1],h(ip).YLim,'k--');
    plot(h(ip),xDF(2,itt)*[1 1],h(ip).YLim,'k--');
    hold(h(ip),'off')
    if ip == 7
      hlegdf = legend(hDF,{'DF'},'location','southeast','box','off');
      hlegdf.Position(1) = 0.32;
    end    
  end
  annotation('textarrow',[0.56 0.544],[0.14 0.14],'String','main X line','fontsize',13,'fontweight','light')
  annotation('textarrow',[0.83 0.84],[0.74 0.76],'String',{'thermalization of','cold inbound ions','at separatrices'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.81 0.825],[0.46 0.48],'String',{'inside separatrix,','inbound ions','are thermalized'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.86 0.87],[0.46 0.50],'String',{'deeper within','exhaust, inbound','ions are still cold'},'fontsize',13,'fontweight','light','horizontalalignment','left')
  annotation('textarrow',[0.86 0.87],[0.87 0.84],'String',{'heated exhaust ions'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.78 0.78],[0.82 0.79],'String',{'cold inflow'},'fontsize',13,'fontweight','light')  
  annotation('textarrow',[0.69 0.675],[0.58 0.54],'String',{'cold outbound','beam close to DF'},'fontsize',13,'fontweight','light','horizontalalignment','left')
  annotation('textarrow',[0.46 0.46],[0.20 0.23],'String',{'inbound ions','with v_y<0'},'fontsize',13,'fontweight','light')
  % islands
  text(h(7),99,2.9,'islands','fontsize',13,'horizontalalignment','center')
  annotation('arrow',[0.260 0.260],[0.33 0.28])
  annotation('arrow',[0.235 0.235],[0.33 0.29]) 
  % secondary X lines
  text(h(8),95,3.5,'secondary X lines','fontsize',13,'horizontalalignment','center')
  annotation('arrow',[0.4795 0.4795],[0.33 0.29])
  annotation('arrow',[0.5117 0.5117],[0.33 0.29])
  % counterstreaming beams
  % alt 1
%   text(h(9),119,-3.25,'counter-streaming beams','fontsize',13,'horizontalalignment','center')
%   annotation('textarrow',0.79*[1 1],[0.17 0.21],'string','cold','fontsize',13) 
%   annotation('textarrow',0.825*[1 1],[0.17 0.21],'string','warmer','fontsize',13)
%   annotation('textarrow',0.875*[1 1],[0.17 0.22],'string','cold','fontsize',13)
  % alt 2
  text(h(9),119,-3.25,'counter-streaming beams','fontsize',13,'horizontalalignment','center')
  annotation('textarrow',[0.78 0.79],[0.17 0.21],'string','cold','fontsize',13) 
  annotation('textarrow',[0.82 0.825],[0.17 0.21],'string','warmer','fontsize',13)
  annotation('textarrow',[0.86 0.870],[0.17 0.19],'string','hot','fontsize',13)
  annotation('textarrow',[0.885 0.875],[0.17 0.22],'string','cold','fontsize',13)
  
end

%% Figure 3.5 % counter-streaming beams in f(vz,z=0)
doE = 0; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.0;
doPhi = 0; colorPhi = [0.5 0.5 0];
doSep = 1;
labstrs = {'z'};
freds = {fred3_z0};
freds_tot = {fred35_z0};
ifred = 1;
xlim = freds{1}.x([1 end]);
zlim = 0.05*[-1 1];
tlim = [75 120];
doFCont = 1; fcontlev = [-0:1:1];

nrows = 4;
ncols = 1;
npanels =nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;

if 1 % Bz(x,t)
  hca = h(isub); isub = isub + 1;   
  pic_tmp = no02m.twcilim(tlim).xlim(xlim).zlim(zlim);
  %Bz_tmp = pic_tmp.Bz;
  %A_tmp = pic_tmp.A;
  Bz_mean = squeeze(mean(Bz_tmp,2));
  A_mean = squeeze(mean(A_tmp,2));
  pcolor(hca,pic_tmp.xi,pic_tmp.twci,Bz_mean')  
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = sprintf('B_z(x,t)',labstr);    
  hca.YLabel.String = 't\omega_{ci}';
  if 1 % A
    hold(hca,'on')
    contour(hca,pic_tmp.xi,pic_tmp.twci,A_mean',[0:1:25],'k')
    hold(hca,'off')
  end
  hca.YDir = 'reverse';
  hca.CLim = [-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 0 % Ey(x,t)
  hca = h(isub); isub = isub + 1;   
  pic_tmp = no02m.twcilim(tlim).xlim(xlim).zlim(zlim);
  %Ey_tmp = pic_tmp.Ey;
  %A_tmp = pic_tmp.A;
  Ey_mean = squeeze(mean(Ey_tmp,2));
  A_mean = squeeze(mean(A_tmp,2));
  pcolor(hca,pic_tmp.xi,pic_tmp.twci,Ey_mean')  
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = sprintf('E_{y}(x,t)',labstr);    
  hca.YLabel.String = 't\omega_{ci}';
  if 1 % A
    hold(hca,'on')
    contour(hca,pic_tmp.xi,pic_tmp.twci,A_mean',[0:1:25],'k')
    hold(hca,'off')
  end
  hca.YDir = 'reverse';
  hca.CLim = [-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 1 % Tic(x,t)
  hca = h(isub); isub = isub + 1;   
  pic_tmp = no02m.twcilim(tlim).xlim(xlim).zlim(zlim);
  %Ti_tmp = pic_tmp.t([3 5]);
  %A_tmp = pic_tmp.A;
  Ti_mean = squeeze(mean(Ti_tmp,2));
  A_mean = squeeze(mean(A_tmp,2));
  pcolor(hca,pic_tmp.xi,pic_tmp.twci,Ti_mean')  
  shading(hca,'flat')
  colormap(hca,pic_colors('thermal'))
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = sprintf('T_{ic}(x,t)',labstr);    
  hca.YLabel.String = 't\omega_{ci}';
  if 1 % A
    hold(hca,'on')
    contour(hca,pic_tmp.xi,pic_tmp.twci,A_mean',[0:1:25],'k')
    hold(hca,'off')
  end
  hca.YDir = 'reverse';
  hca.CLim = [0 0.5];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 1 % fic(vz,z=0)
  hca = h(isub); isub = isub + 1;    
  labstr = labstrs{ifred};
  fred = freds_tot{ifred};  
  fredplot = eval(['fred.fv' labstr]);
  pcolor(hca,fred.x,fred.v,log10(fredplot)')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = sprintf('v_{%s}',labstr);
  colormap(hca,pic_colors('candy4')) 
  %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  %irf_legend(hca,{sprintf('%s z = %g',legends{ifred},unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  irf_legend(hca,{sprintf('%s f(v_%s,z=%g)',legends{ifred},labstr,unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = sprintf('f_{ic}(x,v_{%s})',labstr);    
  hca.CLim = 0.99*[-4 2];
  hca.YLim = 0.99*4*[-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  if doSep
    hold(hca,'on')   
    if 1           
      zz = unique(fred.z);
      [PKS,LOCS] = findpeaks(-abs(sep.z-zz),'sort','descend');
      [~,iz] = min(abs(sep.z-zz));
      if zz == 0
        xx = sep.x(LOCS(1));
        hSep = plot(hca,xx(1)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
      else
        xx = sep.x(LOCS(1:2));
        hSep = plot(hca,xx(1)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
        hSep = plot(hca,xx(2)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
      end
    elseif 0
      xx = unique(fred.x);
      [~,ix] = min(abs(sep.x-xx));
      zz = sep.z(ix);
    end      
    hold(hca,'off')      
  end
  if doFCont
    hold(hca,'on')
    contour(hca,fred.x,fred.v,log10(fredplot)',fcontlev,'k')
    hold(hca,'off')
  end
end
if 1 % fictop/fictot(vz,z=0)
  hca = h(isub); isub = isub + 1;    
  labstr = labstrs{ifred};
  fred = freds{ifred};
  fred_tot = freds_tot{ifred};    
  %['fred_one.fv' labstr]
  fred_one_plot = eval(['fred.fv' labstr]);
  fred_tot_plot = eval(['fred_tot.fv' labstr]);
  fredplot = fred_one_plot./fred_tot_plot;
  pcolor(hca,fred.x,fred.v,fredplot')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = sprintf('v_{%s}',labstr);
  colormap(hca,pic_colors('pasteljet')) 
  %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  %irf_legend(hca,{sprintf('%s z = %g',legends{ifred},unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  irf_legend(hca,{sprintf('%s f(v_%s,z=%g)',legends{ifred},labstr,unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = sprintf('f_{ic}^{top}/f_{ic}^{tot}(x,v_{%s})',labstr);  
  hca.CLim = [0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  if doSep
    hold(hca,'on')   
    if 1           
      zz = unique(fred.z);
      [PKS,LOCS] = findpeaks(-abs(sep.z-zz),'sort','descend');
      [~,iz] = min(abs(sep.z-zz));
      if zz == 0
        xx = sep.x(LOCS(1));
        hSep = plot(hca,xx(1)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
      else
        xx = sep.x(LOCS(1:2));
        hSep = plot(hca,xx(1)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
        hSep = plot(hca,xx(2)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
      end
    elseif 0
      xx = unique(fred.x);
      [~,ix] = min(abs(sep.x-xx));
      zz = sep.z(ix);
    end      
    hold(hca,'off')      
  end  
  if doFCont
    hold(hca,'on')
    contour(hca,fred_tot.x,fred_tot.v,log10(fred_tot_plot)',fcontlev,'k','linewidth',1)
    hold(hca,'off')
  end
end

compact_panels(h,0.01,0.01)
c_eval('h(?).Position(3) = [0.7];',1:npanels)
%% Figure 4, reduced vertical distribution 
%ds100 = PICDist('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');
%% Figure 4, prepare data
for xpick = 75:10:85  
  ds = ds100.twpelim(twpe).xfind(xpick).findtag({'line vertical'});
  xdist = (ds.xi1{1}+ds.xi2{1})/2;
  zdist = (ds.zi1{1}+ds.zi2{1})/2;
  
  zlim_fred = [min(zdist) max(zdist)];
  xlim_fred = xdist(1) + 0.25*[-1 1];
  eval(sprintf('z_x%g = pic.xlim(xlim_fred).zlim(zlim_fred).zi;',xpick))
  %fred3_tmp = ds.reduce_1d_new('x',[3],[]); eval(sprintf('fred3_x%g = fred3_tmp;',xpick))
  %fred5_tmp = ds.reduce_1d_new('x',[5],[]); eval(sprintf('fred5_x%g = fred5_tmp;',xpick))
  fred35_tmp = ds.reduce_1d_new('x',[3 5],[]); eval(sprintf('fred35_x%g = fred35_tmp;',xpick))

  
  eval(sprintf('vExBx_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBx,1));',xpick))
  eval(sprintf('vExBy_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBy,1));',xpick))
  eval(sprintf('vExBz_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBz,1));',xpick))
  eval(sprintf('Bx_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Bx,1));',xpick))
  eval(sprintf('By_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).By,1));',xpick))
  eval(sprintf('Bz_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Bz,1));',xpick))
%   eval(sprintf('vix_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vx([3 5]),1));',xpick))
%   eval(sprintf('viy_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vy([3 5]),1));',xpick))
%   eval(sprintf('viz_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vz([3 5]),1));',xpick))
%   eval(sprintf('vex_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vx([4 6]),1));',xpick))
%   eval(sprintf('vey_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vy([4 6]),1));',xpick))
%   eval(sprintf('vez_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vz([4 6]),1));',xpick))
%   eval(sprintf('Ex_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Ex,1));',xpick))
%   eval(sprintf('Ey_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Ey,1));',xpick))
  eval(sprintf('Ez_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Ez,1));',xpick))
  eval(sprintf('phiz_x%g = -cumtrapz(z_x%g,Ez_x%g);',xpick,xpick,xpick))
end
%sep = no02m.twpelim(twpe).separatrix_location;
disp('Done.')
%% Figure 4, plot
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 3;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1; 
doE = 0; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doVe = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];
doSep = 1;
% if doSep
%   xline_pos = no02m.twpelim(twpe).xline_position;
%   sep = no02m.twpelim(twpe).separatrix_location;
% end
hleg = gobjects(0);

cmap_dist = pic_colors('waterfall');

%freds = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
%freds = {fred35_x75,fred35_x80,fred35_x85,fred35_x90,fred35_x95};
% freds = {fred35_x75,fred35_x85,fred35_x95,fred35_x75,fred35_x85,fred35_x95,fred35_x75,fred35_x85,fred35_x95};
% freds = {fred35_x75,fred35_x80,fred35_x85,fred35_x90,fred35_x95,...
%          fred35_x75,fred35_x80,fred35_x85,fred35_x90,fred35_x95,...
%          fred35_x75,fred35_x80,fred35_x85,fred35_x90,fred35_x95};
freds = {fred35_x75,fred35_x85,...
         fred35_x75,fred35_x85,...
         fred35_x75,fred35_x85};
xpicks = [75 85];
%freds = {fred3_x75,fred3_x85,fred3_x95};
%freds = {fred3_x85};
%freds = {fred5_x75,fred5_x85,fred5_x95};
%freds = {fred46_x75,fred46_x85,fred46_x95};
%freds = {fred46_x75,fred46_x85,fred46_x95,fred46_x75,fred46_x85,fred46_x95,fred46_x75,fred46_x85,fred46_x95};
%freds = {fred46_x75,fred46_x80,fred46_x85,fred46_x90,fred46_x95};
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','x','x','x','x','x','x','x','x'};
labstrs = {'y','y','y','y','y','y','y','y','y'};
labstrs = {'z','z','z','z','z','z','z','z','z'};
labstrs = {'z','z','z','y','y','y','x','x','x'};
labstrs = {'z','z','z','z','z','y','y','y','y','y','x','x','x','x','x'};
labstrs = {'z','z','y','y','x','x'};
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)'};

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.fv' labstr]);
    pcolor(hca,fred.v,fred.z,log10(fredplot))
    shading(hca,'flat')
    hca.YLabel.String = 'z (d_i)';
    hca.XLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy4')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %irf_legend(hca,{sprintf('%s x = %g',legends{isub-1},unique(fred.x))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('%s',legends{isub-1})},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 12;
    if 1*doE
      hold(hca,'on')
      %plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['E' labstr '_x' num2str(unique(fred.x))]);
      hE = plot(hca,smooth(vv,50),zz,'color',0*colorE,'linewidth',1,'linestyle',':');      
      hold(hca,'off')
    end
    if doSep
      hold(hca,'on')      
      xx = unique(fred.x);
      [~,ix] = min(abs(sep.x-xx));
      zz = sep.z(ix);
      hSep = plot(hca,hca.XLim,zz*[1 1],'color',0*colorExB,'linewidth',1,'linestyle',':');      
      hold(hca,'off')      
    end
    if doV
      hold(hca,'on')
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['vi' labstr '_x' num2str(unique(fred.x))]);      
      hv = plot(hca,vv,zz,'color',colorV,'linewidth',1.5);
      hold(hca,'off')
      hleg()
    end
    if doVe
      hold(hca,'on')
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['ve' labstr '_x' num2str(unique(fred.x))]);      
      hv = plot(hca,vv,zz,'color',colorV,'linewidth',1.5);
      hold(hca,'off')
      hleg()
    end
    if doExB
      hold(hca,'on')
      unique(fred.x);
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['vExB' labstr '_x' num2str(unique(fred.x))]);
      hExB = plot(hca,smooth(vv,50),zz,'color',0*colorExB,'linewidth',1,'linestyle','-');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end  
    if doPhi && ifred == 2
      hold(hca,'on')
      unique(fred.x);
      %if 1
      try
      zz = eval(['z_x' num2str(unique(fred.x))]);
      pp = eval(['phi' labstr '_x' num2str(unique(fred.x))]);
      vv0 = eval(['vExB' labstr '_x' num2str(unique(fred.x))]);
      vv0 = vv0(end);
      vv_phi = sqrt(2*(abs(pp-pp(end-150)))).*sign(pp-pp(end-150));
      %(smooth(vv,10)+vv0-vv(end));
      hphi = plot(hca,vv_phi,zz,'color',0*[1 1 1],'linewidth',1,'linestyle','--');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      end
      hold(hca,'off')
    end   
  end
end
drawnow
compact_panels(h(1:end),0.002,0.002)

%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
% colorbars
if 1 % 3 colorbars
  ip = ncols;
  position_hcb_peer = h(ip).Position;
  hcb = colorbar('peer',h(ip));
  h(ip).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})','z');
  hcb.YLabel.FontSize = 14;
  hcb.Position(3)=0.02;
  
  ip = 2*ncols;
  position_hcb_peer = h(ip).Position;
  hcb = colorbar('peer',h(ip));
  h(ip).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})','y');
  hcb.YLabel.FontSize = 14;
  hcb.Position(3)=0.02;
  
  ip = 3*ncols;
  position_hcb_peer = h(ip).Position;
  hcb = colorbar('peer',h(ip));
  h(ip).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})','x');
  hcb.YLabel.FontSize = 14;
  hcb.Position(3)=0.02;
  
else % 1 colorbar
  position_hcb_peer = h(end).Position;
  hcb = colorbar('peer',h(end));
  h(end).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})',labstr);
  hcb.YLabel.FontSize = 14;
end
hl = findobj(h(ncols),'type','line');
hleg = legend(hl,{'v_{ExB}','separatrix'},'edgecolor',[1 1 1],'location','northeast');
hleg = legend(hl([2 1 3]),{'v_{ExB}','\phi_z','separatrix'},'edgecolor',[1 1 1],'location','northeast');
% if all([doV doExB])
%   legend([hv,hExB],{sprintf('v_{i,%s}',labstr),sprintf('v_{ExB,%s}',labstr)})
% elseif doV
%   legend([hv],{sprintf('v_{i,%s}',labstr)})
% elseif doExB
%   legend([hExB],{sprintf('v_{ExB,%s}',labstr)})
% end
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
%compact_panels(0.005,0.01)
h(1).CLim = 0.99*[-4 2];
h(1).CLim = [-3.5 1.8];
h(1).CLim = [-6 0];
h(1).YLim = [-1 8];
%h(1).CLim = 0.99*[-6 0]; % electrons
h(1).XLim = 0.99*3*[-1 1];
for ip = [2:ncols (ncols+2):(2*ncols) (2*ncols+2):(3*ncols)]
  h(ip).YTickLabels = '';
  h(ip).YLabel.String = '';  
end
for ip = 1:npanels
  h(ip).FontSize = 12;
  %h(ip).Position(2) = 0.17;
end
for ip = [2 4 6]
  h(ip).Position(1) = 0.5;
end
%xpicks = xpicks;
for ip = 1:numel(xpicks)
  h(ip).Title.String = sprintf('x = %g',xpicks(ip));
  h(ip).Title.FontWeight = 'normal';
end
for ip = ((nrows-1)*ncols+1):nrows*ncols
  h(ip).XLabel.String = 'v';  
end
for ip = [2 4 6]
  h(ip).Position(1) = 0.4;
end
compact_panels(h,0.002,0.002)
hleg.FontSize = 11;
hcb_ = findobj(gcf,'Type','ColorBar');
if 1 % notations
  annotation('textarrow',[0.56 0.544],[0.14 0.14],'String','main X line','fontsize',13,'fontweight','light')
  annotation('textarrow',[0.83 0.84],[0.74 0.76],'String',{'thermalization of','cold inbound ions','at separatrices'},'fontsize',13,'fontweight','light')  
end
%% Figure 4, plot, only f(vz)
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 1;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1; 
doE = 0; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doVe = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0;
doPhi = 1; colorPhi = [0.5 0.5 0];
doSep = 1;
doB = 1;
% if doSep
%   xline_pos = no02m.twpelim(twpe).xline_position;
%   sep = no02m.twpelim(twpe).separatrix_location;
% end
hleg = gobjects(0);

cmap_dist = pic_colors('waterfall');
freds = {fred35_x75,fred35_x85};
%freds = {fred3_x75,fred3_x85};
xpicks = [75 85];
labstrs = {'z','z'};
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)'};

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.fv' labstr]);
    pcolor(hca,fred.v,fred.z,log10(fredplot))
    shading(hca,'flat')
    hca.YLabel.String = 'z (d_i)';
    hca.XLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy4')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('%s x = %g',legends{isub-1},unique(fred.x))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %irf_legend(hca,{sprintf('%s',legends{isub-1})},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 12;
    if doE
      hold(hca,'on')
      %plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['E' labstr '_x' num2str(unique(fred.x))]);
      hE = plot(hca,smooth(vv,50),zz,'color',colorE,'linewidth',1,'linestyle',':');      
      hold(hca,'off')
    end
    if doB
      hold(hca,'on')
      %plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vvx = eval(['Bx_x' num2str(unique(fred.x))]);
      vvy = eval(['By_x' num2str(unique(fred.x))]);
      vvz = eval(['Bz_x' num2str(unique(fred.x))]);
      hB = plot(hca,smooth(vvx,50),zz,smooth(vvy,50),zz,smooth(vvz,50),zz,'color',colorE,'linewidth',1,'linestyle',':');      
      hold(hca,'off')
    end
    if doSep
      hold(hca,'on')      
      xx = unique(fred.x);
      [~,ix] = min(abs(sep.x-xx));
      zz = sep.z(ix);
      hSep = plot(hca,hca.XLim,zz*[1 1],'color',0*colorExB,'linewidth',1,'linestyle',':');      
      hold(hca,'off')      
    end
    if doV
      hold(hca,'on')
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['vi' labstr '_x' num2str(unique(fred.x))]);      
      hv = plot(hca,vv,zz,'color',colorV,'linewidth',1.5);
      hold(hca,'off')
      hleg()
    end
    if doVe
      hold(hca,'on')
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['ve' labstr '_x' num2str(unique(fred.x))]);      
      hv = plot(hca,vv,zz,'color',colorV,'linewidth',1.5);
      hold(hca,'off')
      hleg()
    end
    if doExB
      hold(hca,'on')
      unique(fred.x);
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['vExB' labstr '_x' num2str(unique(fred.x))]);
      hExB = plot(hca,smooth(vv,50),zz,'color',colorExB,'linewidth',1,'linestyle','-');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end  
    if doPhi && ifred == 2
      hold(hca,'on')
      unique(fred.x);
      %if 1
      try
      zz = eval(['z_x' num2str(unique(fred.x))]);
      pp = eval(['phi' labstr '_x' num2str(unique(fred.x))]);
      vv0 = eval(['vExB' labstr '_x' num2str(unique(fred.x))]);
      vv0 = vv0(end);
      refind = 150; % old
      refind = 400;
      vv_phi = sqrt(2*(abs(pp-pp(end-refind)))).*sign(pp-pp(end-refind));
      %(smooth(vv,10)+vv0-vv(end));
      hphi = plot(hca,vv_phi,zz,'color',0*[1 1 1],'linewidth',1,'linestyle','--');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      end
      hold(hca,'off')
    end   
  end
end
drawnow
compact_panels(h(1:end),0.002,0.002)

%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
% Colorbar
ip = ncols;
position_hcb_peer = h(ip).Position;
hcb = colorbar('peer',h(ip));
h(ip).Position = position_hcb_peer;
hcb.YLabel.String = sprintf('f_{i,c}(z,v_{%s})','z');
hcb.YLabel.FontSize = 14;
%hcb.Position(3)=0.02;
 
hl = findobj(h(ncols),'type','line');
hleg = legend(hl,{'v_{ExB}','separatrix'},'edgecolor',[1 1 1],'location','northeast');
hleg = legend(hl([2 1 3]),{'v_{ExB}','v_{\phi_z}','separatrix'},'edgecolor',[1 1 1],'location','northeast');

h(1).CLim = [-6 -0.2];
h(1).YLim = [-1 8];
h(1).XLim = 0.99*3*[-1 1];

for ip = [2:ncols]
  h(ip).YTickLabels = '';
  h(ip).YLabel.String = '';  
end
for ip = 1:npanels
  h(ip).FontSize = 14;
  %h(ip).Position(2) = 0.17;
end
c_eval('h(?).Position(1) = h(?).Position(1)-0.05;',1:npanels)
c_eval('h(?).Position(2) = h(?).Position(2)+0.05;',1:npanels)
% %xpicks = xpicks;
% for ip = ((nrows-1)*ncols+1):nrows*ncols
%   h(ip).XLabel.String = 'v';  
% end
% for ip = [2 4 6]
%   h(ip).Position(1) = 0.4;
% end
compact_panels(h,0.002,0.002)
hleg.FontSize = 11;
hcb_ = findobj(gcf,'Type','ColorBar');
if 1 % notations
  annotation('textarrow',[0.55 0.6]+0.025,[0.5 0.5]+0.05,'String',{'thermalization','begin at','separatrices'},'fontsize',13,'fontweight','light','horizontalalignment','center')
  annotation('textarrow',[0.55 0.6]+0.02,[0.35 0.35]+0.05,'String',{'inbound ions','follow v_{\phi_z}'},'fontsize',13,'fontweight','light','horizontalalignment','center')
  %annotation('textarrow',[0.55 0.6]+0.02,[0.35 0.35],'String',{'inbound ions','follow v_{\phi_z} and','are gradually','thermalized'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.56 0.6]-0.35,[0.4 0.4]+0.05,'String',{'inbound ions','remain','cold and','follow v_{ExB,z}','deep within','the exhaust'},'fontsize',13,'fontweight','light','horizontalalignment','center')
  %annotation('textarrow',[0.56 0.6]-0.35,[0.4 0.4]+0.05,'String',{'thermal ions','follow v_{ExB,z}','deep within','the exhaust'},'fontsize',13,'fontweight','light','horizontalalignment','center')
  %annotation('textarrow',[0.83 0.84],[0.74 0.76],'String',{'thermalization of','cold inbound ions','at separatrices'},'fontsize',13,'fontweight','light')  
end

%% Figure 4, ALT 2, prepare data
twpe = 24000;
pic = no02m.twpelim(twpe);
for xpick = [75 85]
  ds = ds100.twpelim(twpe).xfind(xpick).findtag({'line vertical'});
  fred35_tmp = ds.reduce_1d_new('x',[3 5],[]); 
  
  xdist = (ds.xi1{1}+ds.xi2{1})/2;
  zdist = (ds.zi1{1}+ds.zi2{1})/2;  
  dxdist = ds.dxi{1}; % all
  dzdist = ds.dzi{1}; % all
  zlim_fred = [min(zdist) max(zdist)];
  xlim_fred = [ds.xi1{1}(1) ds.xi2{1}(1)];
  eval(sprintf('z_x%g = pic.xlim(xlim_fred).zlim(zlim_fred).zi;',xpick))
   
  eval(sprintf('Bx_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Bx,1));',xpick))
  eval(sprintf('By_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).By,1));',xpick))
  eval(sprintf('Bz_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Bz,1));',xpick))
  eval(sprintf('vExBx_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBx,1));',xpick))
  eval(sprintf('vExBy_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBy,1));',xpick))
  eval(sprintf('vExBz_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBz,1));',xpick))
  %eval(sprintf('Ez_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Ez,1));',xpick))
  %eval(sprintf('phiz_x%g = -cumtrapz(z_x%g,Ez_x%g);',xpick,xpick,xpick))
  
  fred35_tmp.x_all = pic.zi;
  eval(sprintf('fred35_tmp.Bx_all = Bx_x%g;',xpick))
  eval(sprintf('fred35_tmp.By_all = By_x%g;',xpick))
  eval(sprintf('fred35_tmp.Bz_all = Bz_x%g;',xpick))
  eval(sprintf('fred35_tmp.vExBx_all = vExBx_x%g;',xpick))
  eval(sprintf('fred35_tmp.vExBy_all = vExBy_x%g;',xpick))
  eval(sprintf('fred35_tmp.vExBz_all = vExBz_x%g;',xpick))
  
  pic_tmp = no02m.twpelim(twpe).xlim(xlim_fred).zlim(zlim_fred);
  Bx = pic_tmp.get_points(xdist,zdist,pic_tmp.twci,[dxdist(1) dzdist(1)],'Bx');
  By = pic_tmp.get_points(xdist,zdist,pic_tmp.twci,[dxdist(1) dzdist(1)],'By');
  Bz = pic_tmp.get_points(xdist,zdist,pic_tmp.twci,[dxdist(1) dzdist(1)],'Bz');
    
  fred35_tmp.Bx = Bx;
  fred35_tmp.By = By;
  fred35_tmp.Bz = Bz;
  eval(sprintf('fred35_x%g = fred35_tmp;',xpick))
  
%   Bx_ = pic_tmp.Bx;
%   By_ = pic_tmp.By;
%   Bz_ = pic_tmp.Bz;
%   
%   Bx = interpfield(pic.xi,pic.zi,Bx_,xdist,zdist); 
%   By = interpfield(pic.xi,pic.zi,By_,xdist,zdist); 
%   Bz = interpfield(pic.xi,pic.zi,Bz_,xdist,zdist);
end
%sep = no02m.twpelim(twpe).separatrix_location;
disp('Done.')

%% Figure 4, ALT 2, plot, alternative with forces discussing the density cut-off
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 2;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1; 
doE = 0; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doVe = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];
doSep = 1;
% if doSep
%   xline_pos = no02m.twpelim(twpe).xline_position;
%   sep = no02m.twpelim(twpe).separatrix_location;
% end
hleg = gobjects(0);

cmap_dist = pic_colors('waterfall');

freds = {fred35_x75,fred35_x75,fred35_x75};
xpicks = [75];
  labstrs = {'x','y','z'};
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)'};
f_cont = [-2:1:2]; 
xlim_map = [65 100];
zlim_map = [-8 8];

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.fv' labstr]);
    pcolor(hca,fred.v,fred.z,log10(fredplot))
    shading(hca,'flat')
    hca.YLabel.String = 'z (d_i)';
    hca.XLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy4')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %irf_legend(hca,{sprintf('%s x = %g',legends{isub-1},unique(fred.x))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('%s',legends{isub-1})},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 12;
    if 1*doE
      hold(hca,'on')
      %plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['E' labstr '_x' num2str(unique(fred.x))]);
      hE = plot(hca,smooth(vv,50),zz,'color',0*colorE,'linewidth',1,'linestyle',':');      
      hold(hca,'off')
    end
    if doSep
      hold(hca,'on')      
      xx = unique(fred.x);
      [~,ix] = min(abs(sep.x-xx));
      zz = sep.z(ix);
      hSep = plot(hca,hca.XLim,zz*[1 1],'color',0*colorExB,'linewidth',1,'linestyle',':');      
      hold(hca,'off')      
    end
    if doV
      hold(hca,'on')
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['vi' labstr '_x' num2str(unique(fred.x))]);      
      hv = plot(hca,vv,zz,'color',colorV,'linewidth',1.5);
      hold(hca,'off')
      hleg()
    end
    if doVe
      hold(hca,'on')
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['ve' labstr '_x' num2str(unique(fred.x))]);      
      hv = plot(hca,vv,zz,'color',colorV,'linewidth',1.5);
      hold(hca,'off')
      hleg()
    end
    if doExB
      hold(hca,'on')
      unique(fred.x);
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['vExB' labstr '_x' num2str(unique(fred.x))]);
      hExB = plot(hca,smooth(vv,50),zz,'color',0*colorExB,'linewidth',1,'linestyle','-');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end  
    if doPhi && ifred == 2
      hold(hca,'on')
      unique(fred.x);
      %if 1
      try
      zz = eval(['z_x' num2str(unique(fred.x))]);
      pp = eval(['phi' labstr '_x' num2str(unique(fred.x))]);
      vv0 = eval(['vExB' labstr '_x' num2str(unique(fred.x))]);
      vv0 = vv0(end);
      vv_phi = sqrt(2*(abs(pp-pp(end-150)))).*sign(pp-pp(end-150));
      %(smooth(vv,10)+vv0-vv(end));
      hphi = plot(hca,vv_phi,zz,'color',0*[1 1 1],'linewidth',1,'linestyle','--');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      end
      hold(hca,'off')
    end   
  end
end

fred = fred35_x75;
% Bx = Bx_x75;
% By = By_x75;
% Bz = Bz_x75;
if 1 % Fz = vxBy
    hca = h(isub); isub = isub + 1;
    [Z,V] = meshgrid(fred.z,fred.v);
    [B,V] = meshgrid(By,fred.v);
    F = V.*B;
    F(fred.fvx'==0) = NaN;
    pcolor(hca,V,Z,F) 
    hold(hca,'on')
    contour(hca,V,Z,log10(fred.fvx'),f_cont,'k')
    hold(hca,'off')
    irf_legend(hca,'F_z = v_xB_y',[0.02 0.98],'color',[0 0 0])
end
if 1 % Fz = -vyBx
    hca = h(isub); isub = isub + 1;
    [Z,V] = meshgrid(fred.z,fred.v);
    [B,V] = meshgrid(Bx,fred.v);
    F = -V.*B;
    F(fred.fvy'==0) = NaN;
    pcolor(hca,V,Z,F) 
    hold(hca,'on')
    contour(hca,V,Z,log10(fred.fvy'),f_cont,'k')
    hold(hca,'off')
    irf_legend(hca,'F_z = -v_yB_x',[0.02 0.98],'color',[0 0 0])
end

if 1 % By
    hca = h(isub); isub = isub + 1;
    no02m.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map(hca,{'By'},'A',1)
    colormap(hca,pic_colors('blue_red'))
    irf_legend(hca,'B_y',[0.02 0.98],'color',[0 0 0])
    hca.CLim = 0.99*[-0.2 0.2];
end

for ip = 4:5 % force panels
  shading(h(ip),'flat')
  h(ip).YLabel.String = 'z';
  colormap(h(ip),pic_colors('blue_red'))
  h(ip).CLim = 0.99*[-0.4 0.4];
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).Layer = 'top';
end
%% 
drawnow
compact_panels(h(1:end),0.002,0.002)

%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
% colorbars
if 1 % 3 colorbars
  ip = ncols;
  position_hcb_peer = h(ip).Position;
  hcb = colorbar('peer',h(ip));
  h(ip).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})','z');
  hcb.YLabel.FontSize = 14;
  hcb.Position(3)=0.02;
  
  ip = 2*ncols;
  position_hcb_peer = h(ip).Position;
  hcb = colorbar('peer',h(ip));
  h(ip).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})','y');
  hcb.YLabel.FontSize = 14;
  hcb.Position(3)=0.02;
  
  ip = 3*ncols;
  position_hcb_peer = h(ip).Position;
  hcb = colorbar('peer',h(ip));
  h(ip).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})','x');
  hcb.YLabel.FontSize = 14;
  hcb.Position(3)=0.02;
  
else % 1 colorbar
  position_hcb_peer = h(end).Position;
  hcb = colorbar('peer',h(end));
  h(end).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})',labstr);
  hcb.YLabel.FontSize = 14;
end
hl = findobj(h(ncols),'type','line');
hleg = legend(hl,{'v_{ExB}','separatrix'},'edgecolor',[1 1 1],'location','northeast');
hleg = legend(hl([2 1 3]),{'v_{ExB}','\phi_z','separatrix'},'edgecolor',[1 1 1],'location','northeast');
% if all([doV doExB])
%   legend([hv,hExB],{sprintf('v_{i,%s}',labstr),sprintf('v_{ExB,%s}',labstr)})
% elseif doV
%   legend([hv],{sprintf('v_{i,%s}',labstr)})
% elseif doExB
%   legend([hExB],{sprintf('v_{ExB,%s}',labstr)})
% end
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
%compact_panels(0.005,0.01)
h(1).CLim = 0.99*[-4 2];
h(1).CLim = [-3.5 1.8];
%h(1).CLim = 0.99*[-6 0]; % electrons
h(1).XLim = 0.99*3*[-1 1];
for ip = [2:ncols (ncols+2):(2*ncols) (2*ncols+2):(3*ncols)]
  h(ip).YTickLabels = '';
  h(ip).YLabel.String = '';  
end
for ip = 1:npanels
  h(ip).FontSize = 12;
  %h(ip).Position(2) = 0.17;
end
for ip = [2 4 6]
  h(ip).Position(1) = 0.5;
end
%xpicks = xpicks;
for ip = 1:numel(xpicks)
  h(ip).Title.String = sprintf('x = %g',xpicks(ip));
  h(ip).Title.FontWeight = 'normal';
end
for ip = ((nrows-1)*ncols+1):nrows*ncols
  h(ip).XLabel.String = 'v';  
end
for ip = [2 4 6]
  h(ip).Position(1) = 0.4;
end
compact_panels(h,0.002,0.002)
hleg.FontSize = 11;
hcb_ = findobj(gcf,'Type','ColorBar');

%% Figure 5, test particle trajectories, location of integrated particles at different times
% Set up reduced distrbutions
twpe = 24000;
%pic = no02m.twpelim(twpe).xlim([60 90]).zlim([-2 8]);
ds = ds100.twpelim(twpe).dxlim([0.5 1]).findtag({'line horizontal'}).zfind(0).xlim([50 100]);
fy3 = ds.fy(1,:,3);
fz3 = ds.fz(1,:,3);
fy5 = ds.fy(1,:,5);
fz5 = ds.fz(1,:,5);
fy = fy3; fy.f = fy3.f + fy5.f;
fz = fz3; fz.f = fz3.f + fz5.f;

fy = fy3;
fz = fz3;

% Set up trajectories
times = [90 100 110 120 125];
flim = [-6 -0.5];
xlim = [65 99.99];
zlim = [-6 9];
fontsize = 14;

colors_matlab = pic_colors('matlab');
colors_tr = [colors_matlab(5,:); 1 1 1;colors_matlab(3,:); 0 0 0];
xlims = [75 85];

tr = tr100(783:917);
tr1 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]>85);
tr = tr100(783:917);
tr2 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]<0);
tr = tr100(783:917);
tr3 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]<=75,[tr.x0]>=71);
tr = tr100(783:917);
tr4 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]>75,[tr.x0]<=85);
tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0);

trajcolordot = nan(tr100.ntr,3);
trajcolordot([tr1.id],:) = repmat(colors_tr(1,:),tr1.ntr,1);
trajcolordot([tr2.id],:) = repmat(colors_tr(2,:),tr2.ntr,1);
trajcolordot([tr3.id],:) = repmat(colors_tr(3,:),tr3.ntr,1);
trajcolordot([tr4.id],:) = repmat(colors_tr(4,:),tr4.ntr,1);
trajcolordot(isnan(trajcolordot)) = [];
trajcolordot = reshape(trajcolordot,numel(trajcolordot)/3,3);


% Plot
climEz = [-1 1];

nrows = 6;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;


tr = tr100(783:917);
tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0);
doTraj = 1;

if 1 % f(x,vy,z=0,t=120), with t0 location of particles
  hca = h(isub); isub = isub + 1;
  labstr = 'y';
  fred = fy;    
  fredplot = fred.f;
  pcolor(hca,fred.x,fred.v,log10(fredplot))
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = sprintf('v_{%s}',labstr);
  colormap(hca,pic_colors('candy4')) 
  %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'location','northoutside');  
  hcb.YLabel.String = sprintf('log_{10}f_{ic}(x,v_{%s})',labstr);
  hcb.FontSize = fontsize;
  %hca.CLim(2) = prctile(fred.fvx(:),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = fontsize;
  if doTraj
    hold(hca,'on')
    for itr_ = 1:4
      xx = eval(sprintf('[tr%g.x0]',itr_));
      vv = eval(sprintf('[tr%g.v%s0]',itr_,labstr));
      hs = scatter(hca,xx,vv,20,[0 0 0],'o');
      hs.MarkerFaceColor = colors_tr(itr_,:);
      %hs. = scatter(hca,xx,vv,20,'MarkerFaceColor',trajcolordot,'Marker','o','MarkerEdgeColor',[0 0 0]);
    end
    hold(hca,'off')
  end
  hca.YLim = 0.99*[-3 3];
  hca.CLim = 0.99*flim;  
end

for itime = 1:numel(times) % Ez
  hca = h(isub); isub = isub + 1;
  time = times(itime);
  pic = no02m.twcilim(time).xlim(xlim).zlim(zlim);
  imagesc(hca,pic.xi,pic.zi,pic.Ez')
  colormap(hca,pic_colors('blue_red'));  
  hca.YDir = 'normal';
  if itime == 1
    hcb = colorbar('peer',hca,'location','northoutside');  
    hcb.YLabel.String = 'E_z';
  end
  if 1 % A
    A = pic.A;
    hold(hca,'on')
    clim = hca.CLim;
    hc = contour(hca,pic.xi,pic.zi,A',[0:1:25],'k');
    %hc.Layer = 'top';
    hca.CLim = clim;  
    hold(hca,'off')
  end
  if doTraj
    hold(hca,'on')
    for itr_ = 1:4
      tr_tmp = eval(sprintf('tr%g',itr_));
      tr_tmp = tr_tmp.resample(time);
      xx = [tr_tmp.x];
      zz = [tr_tmp.z];
      %vv = eval(sprintf('[tr_tmp.v%s0]',labstr));
      hs = scatter(hca,xx,zz,20,[0 0 0],'o');
      hs.MarkerFaceColor = colors_tr(itr_,:);
      %hs. = scatter(hca,xx,vv,20,'MarkerFaceColor',trajcolordot,'Marker','o','MarkerEdgeColor',[0 0 0]);
    end
    hold(hca,'off')
  end
  hca.CLim = climEz;[0 0.5];
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  irf_legend(hca,{sprintf('t = %g',time)},[0.02 0.98],'color','k','fontsize',fontsize)
end

if 1 % f(x,vy,z=0,t=120), with t0 location of particles
  hca = h(isub); isub = isub + 1;
  labstr = 'z';
  fred = fz;    
  fredplot = fred.f;
  pcolor(hca,fred.x,fred.v,log10(fredplot))
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = sprintf('v_{%s}',labstr);
  colormap(hca,pic_colors('candy4')) 
  %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'location','northoutside');  
  hcb.YLabel.String = sprintf('log_{10}f_{ic}(x,v_{%s})',labstr);
  hcb.FontSize = fontsize;
  %hca.CLim(2) = prctile(fred.fvx(:),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = fontsize;
  if doTraj
    hold(hca,'on')
    for itr_ = 1:4
      xx = eval(sprintf('[tr%g.x0]',itr_));
      vv = eval(sprintf('[tr%g.v%s0]',itr_,labstr));
      hs = scatter(hca,xx,vv,20,[0 0 0],'o');
      hs.MarkerFaceColor = colors_tr(itr_,:);
      %hs. = scatter(hca,xx,vv,20,'MarkerFaceColor',trajcolordot,'Marker','o','MarkerEdgeColor',[0 0 0]);
    end
    hold(hca,'off')
  end
  hca.YLim = 0.99*[-3 3];
  hca.CLim = 0.99*flim;  
end

for itime = 1:numel(times)
  hca = h(isub); isub = isub + 1;
  time = times(itime);
  pic = no02m.twcilim(time).xlim(xlim).zlim(zlim);
  imagesc(hca,pic.xi,pic.zi,pic.n(3)')  
  colormap(hca,pic_colors('thermal'));
  hca.YDir = 'normal';
  if itime == 1
    hcb = colorbar('peer',hca,'location','northoutside');  
    hcb.YLabel.String = 'n_{ic}^{top}';
  end
  if 1 % A
    A = pic.A;
    hold(hca,'on')
    clim = hca.CLim;
    hc = contour(hca,pic.xi(1:2:end),pic.zi(1:2:end),A(1:2:end,1:2:end)',[0:1:25],'k');
    %hc.Layer = 'top';
    hca.CLim = clim;  
    hold(hca,'off')
  end
  if doTraj
    hold(hca,'on')
    for itr_ = 1:4
      tr_tmp = eval(sprintf('tr%g',itr_));
      tr_tmp = tr_tmp.resample(time);
      xx = [tr_tmp.x];
      zz = [tr_tmp.z];
      %vv = eval(sprintf('[tr_tmp.v%s0]',labstr));
      hs = scatter(hca,xx,zz,20,[0 0 0],'o');
      hs.MarkerFaceColor = colors_tr(itr_,:);
      %hs. = scatter(hca,xx,vv,20,'MarkerFaceColor',trajcolordot,'Marker','o','MarkerEdgeColor',[0 0 0]);
    end
    hold(hca,'off')
  end
  hca.CLim = [0 0.5];
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  irf_legend(hca,{sprintf('t = %g',time)},[0.02 0.98],'color',[1 1 1],'fontsize',fontsize)
end

hlinks = linkprop(h,{'XLim'});
h(1).XLim = [65 100];
%compact_panels(h,0.01)
firstrow = 1:nrows:nrows*ncols;
remainingrows = setdiff(1:nrows*ncols,firstrow);

h(nrows+1).Position([2 3 4]) = h(1).Position([2 3 4]);

for ip = 1:npanels
  h(ip).FontSize = fontsize;
%   h(ip).Position(3) = 0.4;
%   h(ip).Position(4) = 0.2;
end
for ip = (nrows+1):npanels
  h(ip).YLabel.String = '';  
  h(ip).YTickLabel = [];  
end

compact_panels(h(firstrow),0.01,0.01)
compact_panels(h(remainingrows),0.01,0.01)
drawnow
c_eval('h(?).Position(4) = h(2).Position(4);',firstrow)
c_eval('h(?).Position(2) = h(?).Position(2)-0.06;',[firstrow remainingrows])
c_eval('h(?).Position(2) = h(?).Position(2)+0.04;',firstrow)

% legends
legends = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','o)'};
legends = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O)'};
legends_left = legends([1 3:7]);
legends_right = legends([2 8:12]);
for ip = 1:6
  irf_legend(h(ip),legends_left{ip},[-0.12 0.98],'color',[0 0 0],'fontsize',17,'verticalalignment','middle','fontweight','bold')
end
for ip = 7:12
  irf_legend(h(ip),legends_right{ip-6},[1.02 0.98],'color',[0 0 0],'fontsize',17,'verticalalignment','middle','fontweight','bold')
end

%
% hpos_y_fields = h(2).Position(2);
% c_eval('h(?).Position(2) = hpos_y_fields-0.1;',[2 nrows+2])
% compact_panels(h([2:nrows (nrows+2):npanels]),0.01,0.01)
% compact_panels(h([1 nrows+1]),0.01,0.01)
% c_eval('h(?).Position(2) = h(?).Position(2)-0.05;',[2:nrows (nrows+2):npanels])
if 0 % extra markings and annotations 
  %%
  % Add DF locations
%   itt = no02m.twpelim(twpe).it;  
%   for ip = 7:9
%     hold(h(ip),'on')
%     hDF = plot(h(ip),xDF(1,itt)*[1 1],h(ip).YLim,'k--');
%     plot(h(ip),xDF(2,itt)*[1 1],h(ip).YLim,'k--');
%     hold(h(ip),'off')
%     if ip == 7
%       hlegdf = legend(hDF,{'DF'},'location','southeast','box','off');
%       hlegdf.Position(1) = 0.32;
%     end    
%   end
  % inflow/outflow
  text(h(1),102,2.5,'inflow','fontsize',13,'horizontalalignment','center')
  text(h(1),128,2.5,'outflow','fontsize',13,'horizontalalignment','center')
  text(h(1),75,2.5,'outflow','fontsize',13,'horizontalalignment','center')
  
  annotation('textarrow',[0.56 0.544],[0.14 0.14],'String','main X line','fontsize',13,'fontweight','light')
  annotation('textarrow',[0.83 0.84],[0.74 0.76],'String',{'thermalization of','cold inbound ions','at separatrices'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.81 0.825],[0.46 0.48],'String',{'inside separatrix,','inbound ions','are thermalized'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.86 0.87],[0.46 0.50],'String',{'deeper within','exhaust, inbound','ions are still cold'},'fontsize',13,'fontweight','light','horizontalalignment','left')
  annotation('textarrow',[0.86 0.87],[0.87 0.84],'String',{'heated exhaust ions'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.78 0.78],[0.82 0.79],'String',{'cold inflow'},'fontsize',13,'fontweight','light')  
  annotation('textarrow',[0.69 0.675],[0.58 0.54],'String',{'cold outbound','beam close to DF'},'fontsize',13,'fontweight','light','horizontalalignment','left')
  annotation('textarrow',[0.46 0.46],[0.20 0.23],'String',{'inbound ions','with v_y<0'},'fontsize',13,'fontweight','light')
  %annotation('textarrow',[0.43 0.43],[0.29 0.27],'String',{'inbound ions','with v_y>0'},'fontsize',13,'fontweight','light')
  %annotation('textarrow',[0.62 0.63],[0.84 0.81],'String',{'inbound ions','with v_y>0'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.61 0.63],[0.86 0.81],'String',{'inbound ions','with v_y>0'},'fontsize',13,'fontweight','light')
  annotation('textarrow',[0.60 0.60],[0.73 0.77],'String',{'inbound ions','with v_y<0'},'fontsize',13,'fontweight','light')
  % islands
  text(h(7),99,2.9,'islands','fontsize',13,'horizontalalignment','center')
  annotation('arrow',[0.260 0.260],[0.33 0.28])
  annotation('arrow',[0.235 0.235],[0.33 0.29]) 
  % secondary X lines
  text(h(8),95,3.5,'secondary X lines','fontsize',13,'horizontalalignment','center')
  annotation('arrow',[0.4795 0.4795],[0.33 0.29])
  annotation('arrow',[0.5117 0.5117],[0.33 0.29])
  % counterstreaming beams
  % alt 1
%   text(h(9),119,-3.25,'counter-streaming beams','fontsize',13,'horizontalalignment','center')
%   annotation('textarrow',0.79*[1 1],[0.17 0.21],'string','cold','fontsize',13) 
%   annotation('textarrow',0.825*[1 1],[0.17 0.21],'string','warmer','fontsize',13)
%   annotation('textarrow',0.875*[1 1],[0.17 0.22],'string','cold','fontsize',13)
  % alt 2
  text(h(9),119,-3.35,'counter-streaming beams','fontsize',13,'horizontalalignment','center')
  annotation('textarrow',[0.78 0.79],[0.17 0.21],'string','cold','fontsize',13) 
  annotation('textarrow',[0.82 0.825],[0.17 0.21],'string','warmer','fontsize',13)
  annotation('textarrow',[0.86 0.870],[0.17 0.19],'string','hot','fontsize',13)
  annotation('textarrow',[0.885 0.875],[0.17 0.22],'string','cold','fontsize',13)  
end

h(1).XTickLabels = [];
h(7).XTickLabels = [];
h(7).XLabel = [];
h(1).XLabel = [];
hcb = findobj(gcf,'type','colorbar');
c_eval('hcb(?).YLabel.FontSize = 14;',1:4)
c_eval('hcb(?).FontSize = 14;',1:4)
%c_eval('h(?).Position(2) = h(?).Position(2)+0.02;',[2:nrows (nrows+2):npanels])
%% Figure 6, ALT1, magnetic curvature
%% Curvature plot, combined, also including f(vx,vz) at z=2 to illustrate inward vs outward beam temperature
clear h;
nrows = 4;
ncols = 5;
h(1) = subplot(nrows,ncols,1:ncols);
h1pos = h(1).Position;
h(2) = subplot(nrows,ncols,ncols+(1:ncols));
h2pos = h(2).Position;
c_eval('h(?+2) = subplot(nrows,ncols,?+2*ncols);',1:2*ncols);
twpe = 24000;
xpos = [72:2:80];
fontsize = 14;
legfontsize = 11;
xlim = [60 110];
xpos2 = [69 71 73 75 77 79];
xpos2 = [71 73 75 77 79];
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','k)','l)','m)','n)','o)'};

ds = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([71 80]).zfind(0).xfind(xpos);
ds2 = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([71 80]).zfind(2).xfind(xpos);
ds2 = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([50 80]).zfind(2).xfind(xpos2);
pic = no02m.twpelim(twpe).xlim(xlim).zlim([-7 7]);

hca = h(1);
%hmap = pic.plot_map(hca,{'log10(curvbrad)'},'A',1,'cmap',pic_colors('blue_red'),'cbarlabels',{'log_{10}r_B'},'smooth',2);
hmap = pic.plot_map(hca,{'n(5)'},'A',1,'cmap',pic_colors('blue_red'),'cbarlabels',{'log_{10}r_B'},'smooth',2);
hmap.CLim = [-1 3];
hca.Position = h1pos;
hold(hca,'on')
ds.plot_boxes(hca,'color',[1 1 1]);
ds2.plot_boxes(hca,'color',[0 0 0]);
hold(hca,'off')

hca = h(2);
comp = 'x';
hline = pic.xlim([63 xlim(2)]).zlim([-0.25 0.25]).plot_line(hca,comp,{{'Babs','curvbrad','curvbrad.*Babs'}},'smooth',10);
hca.Position = h2pos;
hca.YLim = [0 3.99];
%legend(hca,{'|B|','r_B','r_B\omega_{ci}'},'location','north','orientation','horizontal')
legend(hca,{'|B|','r_B','r_B\omega_{ci}'},'location','north','orientation','horizontal','box','off')
hca.Title.String = [];
hca.XLim = xlim;
hca.Position(2) = hca.Position(2) + 0.08;
compact_panels(h(1:2),0)
c_eval('h(?).Position(2) = h(?).Position(2) + 0.03;',1:2)
irf_legend(hca,{'z=0\pm0.25'},[0.98 0.98],'color',[0 0 0],'fontsize',fontsize)

ih0 = 2+ncols;
if 1 % f(v_x,v_y)
  sumdim = 3;
  clim = [-4 1];
  xlim = 3.5*0.99*[-1 1];
  ylim = 3.5*0.99*[-1 1];
  x0 = (ds.xi1{1}+ds.xi2{1})/2;
  hds = ds.plot_map(h(ih0+(1:ncols)),[5],sumdim,'v',no02m,'log','curv',{no02m,1},'nolabel'); % 
  hlinks = linkprop(hds.ax,{'XLim','YLim','CLim','XTick','YTick'});
  c_eval('hds.ax(?).Position(2) = hds.ax(?).Position(2) + 0.10;',1:ncols)
  hds.ax(1).Position(1) = h(1).Position(1);
  compact_panels(hds.ax,0.00,0.00)
  c_eval('irf_legend(hds.ax(?),sprintf(''%s x = %g, z = 0'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',legfontsize);',1:ncols);
  %[hax,hlab] = label_panels(hds.ax);
  hds.ax(1).CLim = clim;
  %h.ax(1).CLim = [0 1];
  hds.ax(1).XLim = xlim;
  hds.ax(1).YLim = ylim;
  %c_eval('hds.ax(?).XLabel.String = []; hds.ax(?).XTickLabels = [];',[1 2 3]);
  c_eval('hds.ax(?).YLabel.String = []; hds.ax(?).YTickLabels = [];',2:ncols);
  pos = hds.ax(end).Position;
  hcb = colorbar('peer',hds.ax(end));
  hcb.YLabel.String = 'f(v_x,v_y)';
  hds.ax(end).Position = pos;
end

ih0 = 2;
if 1 % f(v_x,v_z)
  hds2 = ds2.plot_map(h(ih0+(1:ncols)),[5],2,'v',no02m,'bline',no02m,'log','nolabel','nan','exb',no02m); % 
  hlinks2 = linkprop(hds2.ax,{'XLim','YLim','CLim','XTick','YTick'});
  x0 = (ds2.xi1{1}+ds2.xi2{1})/2;
  c_eval('irf_legend(hds2.ax(?),sprintf(''%s x = %g, z = 2'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',legfontsize);',1:ncols);
  hds2.ax(1).Position(1) = h(1).Position(1);
  %c_eval('hds2.ax(?).Position(2) = hds2.ax(?).Position(2) + 0.10;',1:3)
  compact_panels(hds2.ax,0.00,0.00)
  c_eval('hds2.ax(?).YLabel.String = []; hds2.ax(?).YTickLabels = [];',2:ncols);
  hds2.ax(1).CLim = clim;
  %h.ax(1).CLim = [0 1];
  hds2.ax(1).XLim = xlim;
  hds2.ax(1).YLim = ylim;
  pos = hds2.ax(end).Position;
  hcb = colorbar('peer',hds2.ax(end));
  hcb.YLabel.String = 'f(v_x,v_z)';
  hds2.ax(end).Position = pos;
end

irf_legend(hmap,{'a)'},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)
irf_legend(h(2),{'b)'},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)


hlinks = linkprop(h(3:end),{'XLim','YLim','CLim'});
h(1).Title.String = '';

c_eval('h(?).FontSize = fontsize;',1:numel(h));
c_eval('h(?).Position(1) = h(?).Position(1)-0.04;',1:numel(h));
c_eval('h(?).Position(4) = h(end).Position(4);',3:numel(h))
drawnow
c_eval('h(?).Position(2) = h(?).Position(2)+0.06;',3:(2+ncols));
%compact_panels(h(3:end),0.0,0)

%irf_legend(h(8),{'B'},[0.92 0.27],'color',0.5+[0 0 0],'fontsize',16)

%c_eval('h().Position')
hl = findobj(h(2+ncols),'type','line');
hs = findobj(h(2+ncols),'type','scatter');
legend([hl;hs],{'B','v_{bulk}','v_{ExB}'},'box','off','location','northoutside','orientation','horizontal')
%% Figure 6, ALT2, Curvature plot, combined, also including f(vx,vz) at z=2 to illustrate inward vs outward beam temperature
clear h;
nrows = 4;
ncols = 5;
h(1) = subplot(nrows,ncols,1:ncols);
h1pos = h(1).Position;
h(2) = subplot(nrows,ncols,ncols+(1:ncols));
h2pos = h(2).Position;
c_eval('h(?+2) = subplot(nrows,ncols,?+2*ncols);',1:2*ncols);
twpe = 24000;
xpos = [72:2:80];
fontsize = 14;
legfontsize = 11;
xlim = [60 110];
xpos2 = [69 71 73 75 77 79];
xpos2 = [71 73 75 77 79];
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','k)','l)','m)','n)','o)'};

ds = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([71 80]).zfind(0).xfind(xpos);
ds2 = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([71 80]).zfind(2).xfind(xpos);
ds2 = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([50 80]).zfind(2).xfind(xpos2);
pic = no02m.twpelim(twpe).xlim(xlim).zlim([-7 7]);

hca = h(1);
%hmap = pic.plot_map(hca,{'log10(curvbrad)'},'A',1,'cmap',pic_colors('blue_red'),'cbarlabels',{'log_{10}r_B'},'smooth',2);
hmap = pic.plot_map(hca,{'n(5)'},'A',1,'cmap',pic_colors('blue_red'),'cbarlabels',{'log_{10}r_B'},'smooth',2);
hmap.CLim = [-1 3];
hca.Position = h1pos;
hold(hca,'on')
ds.plot_boxes(hca,'color',[1 1 1]);
ds2.plot_boxes(hca,'color',[0 0 0]);
hold(hca,'off')

hca = h(2);
comp = 'x';
hline = pic.xlim([63 xlim(2)]).zlim([-0.25 0.25]).plot_line(hca,comp,{{'Babs','curvbrad','curvbrad.*Babs'}},'smooth',10);
hca.Position = h2pos;
hca.YLim = [0 3.99];
%legend(hca,{'|B|','r_B','r_B\omega_{ci}'},'location','north','orientation','horizontal')
legend(hca,{'|B|','r_B','r_B\omega_{ci}'},'location','north','orientation','horizontal','box','off')
hca.Title.String = [];
hca.XLim = xlim;
hca.Position(2) = hca.Position(2) + 0.08;
compact_panels(h(1:2),0)
c_eval('h(?).Position(2) = h(?).Position(2) + 0.03;',1:2)
irf_legend(hca,{'z=0\pm0.25'},[0.98 0.98],'color',[0 0 0],'fontsize',fontsize)


xlim = 3.5*0.99*[-1 1];
ylim = 3.5*0.99*[-1 1];
ih0 = 2+ncols;
if 0 % f(v_x,v_y)
  sumdim = 3;
  clim = [-4 1];
  xlim = 3.5*0.99*[-1 1];
  ylim = 3.5*0.99*[-1 1];
  x0 = (ds.xi1{1}+ds.xi2{1})/2;
  hds = ds.plot_map(h(ih0+(1:ncols)),[3],sumdim,'v',no02m,'log','curv',{no02m,1},'nolabel'); % 
  hlinks = linkprop(hds.ax,{'XLim','YLim','CLim','XTick','YTick'});
  c_eval('hds.ax(?).Position(2) = hds.ax(?).Position(2) + 0.10;',1:ncols)
  hds.ax(1).Position(1) = h(1).Position(1);
  compact_panels(hds.ax,0.00,0.00)
  c_eval('irf_legend(hds.ax(?),sprintf(''%s x = %g, z = 0'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',legfontsize);',1:ncols);
  %[hax,hlab] = label_panels(hds.ax);
  hds.ax(1).CLim = clim;
  %h.ax(1).CLim = [0 1];
  hds.ax(1).XLim = xlim;
  hds.ax(1).YLim = ylim;
  %c_eval('hds.ax(?).XLabel.String = []; hds.ax(?).XTickLabels = [];',[1 2 3]);
  c_eval('hds.ax(?).YLabel.String = []; hds.ax(?).YTickLabels = [];',2:ncols);
  pos = hds.ax(end).Position;
  hcb = colorbar('peer',hds.ax(end));
  hcb.YLabel.String = 'f(v_x,v_y)';
  hds.ax(end).Position = pos;
end

ih0 = 2;
if 1 % f(v_x,v_z)
  sumdim = 2;
  hds2 = ds2.plot_map(h(ih0+(1:ncols)),[5],sumdim,'v',no02m,'bline',no02m,'nolabel','log','nan','exb',no02m); % 
  hlinks2 = linkprop(hds2.ax,{'XLim','YLim','CLim','XTick','YTick'});
  x0 = (ds2.xi1{1}+ds2.xi2{1})/2;
  c_eval('irf_legend(hds2.ax(?),sprintf(''%s x = %g, z = 2'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',legfontsize);',1:ncols);
  hds2.ax(1).Position(1) = h(1).Position(1);
  %c_eval('hds2.ax(?).Position(2) = hds2.ax(?).Position(2) + 0.10;',1:3)
  compact_panels(hds2.ax,0.00,0.00)
  c_eval('hds2.ax(?).YLabel.String = []; hds2.ax(?).YTickLabels = [];',2:ncols);
  hds2.ax(1).CLim = clim;
  %hds2.ax(1).CLim = [0 1];
  %h.ax(1).CLim = [0 1];
  hds2.ax(1).XLim = xlim;
  hds2.ax(1).YLim = ylim;
  pos = hds2.ax(end).Position;
  hcb = colorbar('peer',hds2.ax(end));
  hcb.YLabel.String = 'log_{10}f(v_x,v_z)';
  hds2.ax(end).Position = pos;
  %c_eval('colormap(hds2.ax(?),pic_colors(''pasteljet''));',1:numel(hds2.ax))
end
ih0 = 2+ncols;
if 1 % f(v_x,v_z)
  sumdim = 2;
  hds3 = ds2.plot_map(h(ih0+(1:ncols)),[3],sumdim,'ratio',[3 5],'v',no02m,'bline',no02m,'nolabel','nan','exb',no02m); % 
  hlinks2 = linkprop(hds2.ax,{'XLim','YLim','CLim','XTick','YTick'});
  x0 = (ds2.xi1{1}+ds2.xi2{1})/2;
  c_eval('irf_legend(hds3.ax(?),sprintf(''%s x = %g, z = 2'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',legfontsize);',1:ncols);
  hds3.ax(1).Position(1) = h(1).Position(1);
  %c_eval('hds2.ax(?).Position(2) = hds2.ax(?).Position(2) + 0.10;',1:3)
  compact_panels(hds3.ax,0.00,0.00)
  c_eval('hds3.ax(?).YLabel.String = []; hds3.ax(?).YTickLabels = [];',2:ncols);
  %hds3.ax(1).CLim = clim;
  hds3.ax(1).CLim = [0 1];
  %h.ax(1).CLim = [0 1];
  hds3.ax(1).XLim = xlim;
  hds3.ax(1).YLim = ylim;
  pos = hds3.ax(end).Position;
  hcb = colorbar('peer',hds3.ax(end));
  hcb.YLabel.String = 'f(v_x,v_z)';
  hds3.ax(end).Position = pos;
  c_eval('colormap(hds3.ax(?),pic_colors(''pasteljet''));',1:numel(hds3.ax))
end


irf_legend(hmap,{'a)'},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)
irf_legend(h(2),{'b)'},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)


%hlinks = linkprop(h(3:end),{'XLim','YLim','CLim'});
h(1).Title.String = '';

c_eval('h(?).FontSize = fontsize;',1:numel(h));
c_eval('h(?).Position(1) = h(?).Position(1)-0.04;',1:numel(h));
c_eval('h(?).Position(4) = h(end).Position(4);',3:numel(h))
drawnow
c_eval('h(?).Position(2) = h(?).Position(2)+0.06;',3:(2+ncols));
%compact_panels(h(3:end),0.0,0)

%irf_legend(h(8),{'B'},[0.92 0.27],'color',0.5+[0 0 0],'fontsize',16)

%c_eval('h().Position')
hl = findobj(h(2+ncols),'type','line');
hs = findobj(h(2+ncols),'type','scatter');
legend([hl;hs],{'B','v_{bulk}','v_{ExB}'},'box','off','location','northoutside','orientation','horizontal')

%% Figure 6, ALT3, Curvature plot, combined, fred along two selected field lines
% OBS, need to runs econd part to get locations of reduced distributions

doExB = 1;
ExBcol = [0.5 0.5 0.5];
contlev = [-1:0.5:0];
forcelim = [-0.5 0.5];
fclim = [-6 -0.5];

% Prepare plot
clear h hpos;
nrows = 4;
ncols = 5;
nmaps = 3;
h(1) = subplot(nrows,ncols,0*ncols+1:ncols);
hpos{1} = h(1).Position;
h(2) = subplot(nrows,ncols,1*ncols+(1:ncols));
hpos{2} = h(2).Position;
h(3) = subplot(nrows,ncols,2*ncols+(1:ncols));
hpos{3} = h(3).Position;

c_eval('h(?-ncols*nmaps+nmaps) = subplot(nrows,ncols,?)',(ncols*3+1):nrows*ncols);
twpe = 24000;
xpos = [72:2:80];
fontsize = 14;
legfontsize = 12;
xlim = [60 110];
xpos2 = [69 71 73 75 77 79];
xpos2 = [71 73 75 77 79];
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','k)','l)','m)','n)','o)'};
legends = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N'};

%ds = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([71 80]).zfind(0).xfind(xpos);
%ds2 = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([71 80]).zfind(2).xfind(xpos);
%ds2 = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([50 80]).zfind(2).xfind(xpos2);
ds = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([65 80]).zfind(0).xfind(xpos);
%ds2 = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([71 80]).zfind(2).xfind(xpos);
ds2 = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([65 80]).zfind(2).xfind(xpos2);

pic = no02m.twpelim(twpe).xlim(xlim).zlim([-7 7]);

isub = 1;

hca = h(isub); isub = isub + 1;
hmap = pic.plot_map(hca,{'n(3)'},'A',1,'cmap',pic_colors('thermal'),'cbarlabels',{'n_{ic}^{top}'},'smooth',2);
hmap.CLim = [0.0000 0.4999];
hca.Position = hpos{isub-1};
hold(hca,'on')
ds.plot_boxes(hca,'color',[1 1 1]);
%ds2.plot_boxes(hca,'color',[0 0 0]);
hold(hca,'off')

hca = h(isub); isub = isub + 1;
hmap = pic.plot_map(hca,{'log10(curvbrad)'},'A',1,'cmap',pic_colors('blue_red'),'cbarlabels',{'log_{10}r_B'},'smooth',2);
hmap.CLim = [-1 2.99];
hca.Position = hpos{isub-1};
hold(hca,'on')
ds.plot_boxes(hca,'color',[1 1 1]);
%ds2.plot_boxes(hca,'color',[0 0 0]);
hold(hca,'off')

hca = h(isub); isub = isub + 1;
comp = 'x';
hline = pic.xlim([63 xlim(2)]).zlim([-0.25 0.25]).plot_line(hca,comp,{{'Babs','curvbrad','curvbrad.*Babs'}},'smooth',10);
hca.Position = hpos{isub-1};
hca.YLim = [0 3.99];
%legend(hca,{'|B|','r_B','r_B\omega_{ci}'},'location','north','orientation','horizontal')
hlegLlines = legend(hca,{'|B|','r_B','v_{i,\perp}^{\kappa=1}=r_B\omega_{ci}'},'location','north','orientation','horizontal','box','off');
hca.Title.String = [];
hca.XLim = xlim;
ycompr = 0.03;
h(2).Title.String = [];
%compact_panels(h(1:2),0)
%c_eval('h(?).Position(2) = h(?).Position(2) + 0.03;',1:2)
irf_legend(hca,{'z=0\pm0.25'},[0.98 0.98],'color',[0 0 0],'fontsize',fontsize)

compact_panels(h(1:nmaps),0)
% make this panel slightly smaller
h(3).Position(2) = h(3).Position(2) + ycompr;
h(3).Position(4) = h(3).Position(4) - ycompr;
c_eval('h(?).Position(1) = h(?).Position(1) - 0.05;',1:numel(h))

dsslimmin = [65.0 70.8];
for ids = 1:numel(dss)
  ds_tmp = dss{ids}; 
  ds_tmp.xlim([dsslimmin(ids) 100]).plot_boxes(h(1),'color',[1 1 1]);
  ds_tmp.xlim([dsslimmin(ids) 100]).plot_boxes(h(2),'color',[1 1 1]);
  %ds_tmp.plot_map(h(1)'color',[1 1 1]);
end
%
xlim = 3.5*0.99*[-1 1];
ylim = 3.5*0.99*[-1 1];
fclim2d = [-7 0];
ih0 = nmaps;
if 1 % f(v_x,v_y)
  sumdim = 3;
  clim = [-4 1];
  xlim = 3.5*0.99*[-1 1];
  ylim = 3.5*0.99*[-1 1];
  x0 = (ds.xi1{1}+ds.xi2{1})/2;
  hds = ds.plot_map(h(ih0+(1:ncols)),[3],sumdim,'v',no02m,'log','curv',{no02m,1},'nolabel'); % 
  hlinks = linkprop(hds.ax,{'XLim','YLim','CLim','XTick','YTick'});
  %c_eval('hds.ax(?).Position(2) = hds.ax(?).Position(2) + 0.10;',1:ncols)
  hds.ax(1).Position(1) = h(1).Position(1);
  compact_panels(hds.ax,0.00,0.00)
  c_eval('irf_legend(hds.ax(?),sprintf(''%s'',legends{?+ih0}),[0.02 0.99],''color'',[0 0 0],''fontsize'',19,''fontweight'',''bold'');',1:ncols);
  c_eval('irf_legend(hds.ax(?),sprintf(''x = %g, z = 0'',x0(?)),[0.98 0.98],''color'',[0 0 0],''fontsize'',legfontsize);',1:ncols);
  %c_eval('irf_legend(hds.ax(?),sprintf(''%s x = %g, z = 0'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',legfontsize);',1:ncols);
  %[hax,hlab] = label_panels(hds.ax);
  hds.ax(1).CLim = fclim2d;
  %h.ax(1).CLim = [0 1];
  hds.ax(1).XLim = xlim;
  hds.ax(1).YLim = ylim;
  %c_eval('hds.ax(?).XLabel.String = []; hds.ax(?).XTickLabels = [];',[1 2 3]);
  c_eval('hds.ax(?).YLabel.String = []; hds.ax(?).YTickLabels = [];',2:ncols);
  pos = hds.ax(end).Position;
  hcb = colorbar('peer',hds.ax(end));
  hcb.YLabel.String = 'log_{10}f_{ic}^{top}(v_x,v_y)';
  hds.ax(end).Position = pos;
end


%ih0 = nmaps;
if 0 % f(v_x,v_z)
  sumdim = 2;
  hds2 = ds2.plot_map(h(ih0+(1:ncols)),[5],sumdim,'v',no02m,'bline',no02m,'nolabel','log','nan','exb',no02m); % 
  hlinks2 = linkprop(hds2.ax,{'XLim','YLim','CLim','XTick','YTick'});
  x0 = (ds2.xi1{1}+ds2.xi2{1})/2;
  c_eval('irf_legend(hds2.ax(?),sprintf(''%s x = %g, z = 2'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',legfontsize);',1:ncols);
  hds2.ax(1).Position(1) = h(1).Position(1);
  %c_eval('hds2.ax(?).Position(2) = hds2.ax(?).Position(2) + 0.10;',1:3)
  compact_panels(hds2.ax,0.00,0.00)
  c_eval('hds2.ax(?).YLabel.String = []; hds2.ax(?).YTickLabels = [];',2:ncols);
  hds2.ax(1).CLim = fclim2d;
  %hds2.ax(1).CLim = [0 1];
  %h.ax(1).CLim = [0 1];
  hds2.ax(1).XLim = xlim;
  hds2.ax(1).YLim = ylim;
  pos = hds2.ax(end).Position;
  hcb = colorbar('peer',hds2.ax(end));
  hcb.YLabel.String = 'log_{10}f(v_x,v_z)';
  hds2.ax(end).Position = pos;
  %c_eval('colormap(hds2.ax(?),pic_colors(''pasteljet''));',1:numel(hds2.ax))
end
%ih0 = 2+ncols;
if 0 % f(v_x,v_z)
  sumdim = 2;
  hds3 = ds2.plot_map(h(ih0+(1:ncols)),[3],sumdim,'ratio',[3 5],'v',no02m,'bline',no02m,'nolabel','nan','exb',no02m); % 
  hlinks2 = linkprop(hds2.ax,{'XLim','YLim','CLim','XTick','YTick'});
  x0 = (ds2.xi1{1}+ds2.xi2{1})/2;
  c_eval('irf_legend(hds3.ax(?),sprintf(''%s x = %g, z = 2'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',legfontsize);',1:ncols);
  hds3.ax(1).Position(1) = h(1).Position(1);
  %c_eval('hds2.ax(?).Position(2) = hds2.ax(?).Position(2) + 0.10;',1:3)
  compact_panels(hds3.ax,0.00,0.00)
  c_eval('hds3.ax(?).YLabel.String = []; hds3.ax(?).YTickLabels = [];',2:ncols);
  %hds3.ax(1).CLim = clim;
  hds3.ax(1).CLim = [0 1];
  %h.ax(1).CLim = [0 1];
  hds3.ax(1).XLim = xlim;
  hds3.ax(1).YLim = ylim;
  pos = hds3.ax(end).Position;
  hcb = colorbar('peer',hds3.ax(end));
  hcb.YLabel.String = 'f(v_x,v_z)';
  hds3.ax(end).Position = pos;
  c_eval('colormap(hds3.ax(?),pic_colors(''pasteljet''));',1:numel(hds3.ax))
end

for ip = 1:nmaps
  %irf_legend(h(ip),{legends{ip}},[0.01 0.99],'color',[0 0 0],'fontsize',fontsize)
  irf_legend(h(ip),{legends{ip}},[0.01 0.99],'color',[0 0 0],'fontsize',19,'fontweight','bold')
end
%irf_legend(h(2),{'b)'},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)

heigths = [0.12 0.12 0.08];
clear dh
for ih = 1:nmaps
  dh(ih) = heigths(ih) - h(ih).Position(4);
  h(ih).Position(4) = heigths(ih);
  h(ih).Position(2) = h(ih).Position(2) - sum(dh);
end
drawnow
c_eval('h(?).Position(2) = h(?).Position(2)-sum(dh)+0.02;',(nmaps)+[1:ncols]);
%h(1).Position(4)  = heigths(1);
%h(1).Position(2)  = heigths(2);
%h(1).Position(3)  = heigths(3);

%hlinks = linkprop(h(3:end),{'XLim','YLim','CLim'});
h(1).Title.String = '';
%
c_eval('h(?).FontSize = fontsize;',1:numel(h));
c_eval('h(?).FontSize = 12;',1:numel(h));
%c_eval('h(?).Position(1) = h(?).Position(1)-0.04;',1:numel(h));
%c_eval('h(?).Position(4) = h(end).Position(4);',3:numel(h))
%drawnow
%c_eval('h(?).Position(2) = h(?).Position(2)+0.06;',3:(2+ncols));
%compact_panels(h(3:end),0.0,0)

%irf_legend(h(8),{'B'},[0.92 0.27],'color',0.5+[0 0 0],'fontsize',16)
%annotation('textarrow',[0.24 0.22],[0.140 0.15]-sum(dh)+0.02,'String',{'v_{i,\perp}^{\kappa=1}'},'fontsize',13,'fontweight','light')
annotation('textarrow',[0.24 0.22]-0.05,[0.140 0.15]-sum(dh)+0.02,'String',{'v_{i,\perp}^{\kappa=1}'},'fontsize',13,'fontweight','light')

%c_eval('h().Position')
htmp = h(end-1);
tmppos = htmp.Position;
%hl = findobj(htmp,'type','line');
hs = findobj(htmp,'type','scatter');
hleg_dist = legend([hs],{'v_{ExB}','v_{bulk}'},'location','south','orientation','horizontal','edgecolor',[1 1 1]);
htmp.Position = tmppos;

set(gcf,'position',[104         582         616*1.2         610*1.2])

%% SECOND PART, with reduced distributions along field line
% Prepare 1D reduced distributions
twpe = 24000;
pic = no02m.twpelim(twpe).xlim([60 90]).zlim([-8 8]);
tags = {'A=6.0','A=8.0'};
%tags = {'A=6.0','A=7.5','A=8.0'};
tags = tags(end:-1:1);
%tags = {'A=6.0'};

nred = numel(tags);

%tags = {'line vertical'}; xpick = [75 80];
%nred = numel(xpick);

% fpar3 = struct([]);
% fpar5 = struct([]);
dsslimmin = [65.0 70.8];
if 1
clear fx3 fx5 fy3 fy5 fz3 fz5 dss
for itag = 1:nred
  %ds = ds100.twpelim(twpe).dxlim([0.1 0.3]).findtag(tags(itag)).xlim([50 90]);
  ds = ds100.twpelim(twpe).xlim([dsslimmin(itag) 100]).dxlim([0.1 0.3]).findtag(tags(itag)).xlim([50 90]);
  %ds = ds100.twpelim(twpe).dxlim([0.5 1]).findtag(tags).xfind(xpick(itag)).xlim([50 100]);
  iSpecies3 = 3;
  iSpecies5 = 3;
  dss{itag} = ds;
  fx3_tmp = ds.fx(1,:,iSpecies3);
  fx5_tmp = ds.fx(1,:,iSpecies5);
  fx3{itag} = fx3_tmp;
  fx5{itag} = fx5_tmp;
  fy3_tmp = ds.fy(1,:,iSpecies3);
  fy5_tmp = ds.fy(1,:,iSpecies5);
  fy3{itag} = fy3_tmp;
  fy5{itag} = fy5_tmp;
  fz3_tmp = ds.fz(1,:,iSpecies3);
  fz5_tmp = ds.fz(1,:,iSpecies5);
  fz3{itag} = fz3_tmp;
  fz5{itag} = fz5_tmp;
end
end

% 
doExB = 1;
ExBcol = [0.5 0.5 0.5]*0+0;
ExBlinewidth = 0.5;
doCont = 0;
contlev = [-1:0.5:0];
forcelim = [-0.5 0.5];
fclim = [-5 -0.0];
ncols = nred;
nrows = 3;
doPlotSurf = 0;
%nrows = 7;
%ncols = 4;
cmap = pic_colors('pasteljet');
cmap = pic_colors('candy4');
plotxstr = 'arc_z0';
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');
%h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
% Distributions
for itag = 1:nred
  if 0 % arclength on vertical axis
    if 1 % fx5
      hca = h(isub); isub = isub + 1;
      ff1 = fx5{itag};
      pcolor(hca,ff1.v,ff1.(plotxstr),log10(ff1.f)'); 
      shading(hca,'flat'); 
      hca.XLabel.String = 'v (v_A)'; 
      if 1 % contour lines
        hold(hca,'on')      
        contour(hca,ff1.v,ff1.(plotxstr),log10(ff1.f)',contlev,'color',[0 0 0])      
        hold(hca,'off')
      end
      if doExB % ExB
        B2 = ff1.Bx.^2 + ff1.By.^2 + ff1.Bz.^2;
        ExB = (ff1.Ey.*ff1.Bz - ff1.Ez.*ff1.By)./B2; % x-comp
        hold(hca,'on')
        plot(hca,ExB,ff1.(plotxstr),'color',[1 1 1],'linewidth',1)
        hold(hca,'off')
      end
      colormap(hca,pic_colors('candy4')); 
      hca.CLim = fclim;    
      hleg = irf_legend(hca,'f(v_x)',[0.02,0.98],'color',[0 0 0],'backgroundcolor',[1 1 1]); 
    end
    if 1 % fy5
      hca = h(isub); isub = isub + 1;
      ff1 = fy5{itag};
      pcolor(hca,ff1.v,ff1.(plotxstr),log10(ff1.f)'); 
      shading(hca,'flat'); 
      colormap(hca,pic_colors('candy4')); 
      hca.CLim = fclim;
      hca.XLabel.String = 'v (v_A)'; 
      if 1 % contour lines
        hold(hca,'on')      
        contour(hca,ff1.v,ff1.(plotxstr),log10(ff1.f)',contlev,'color',[0 0 0])      
        hold(hca,'off')
      end
      if doExB % ExB
        B2 = ff1.Bx.^2 + ff1.By.^2 + ff1.Bz.^2;
        ExB = (ff1.Ez.*ff1.Bx - ff1.Ex.*ff1.Bz)./B2; % y-comp
        hold(hca,'on')
        plot(hca,ExB,ff1.(plotxstr),'color',[1 1 1],'linewidth',1)
        hold(hca,'off')
      end
      hleg = irf_legend(hca,'f(v_y)',[0.02,0.98],'color',[0 0 0],'backgroundcolor',[1 1 1]); 
    end
    if 1 % fz5
    hca = h(isub); isub = isub + 1;
    ff1 = fz5{itag};
    pcolor(hca,ff1.v,ff1.(plotxstr),log10(ff1.f)'); 
    shading(hca,'flat'); 
    colormap(hca,pic_colors('candy4')); 
    hca.CLim = fclim;
    hca.XLabel.String = 'v (v_A)'; 
    if 1 % contour lines
      hold(hca,'on')      
      contour(hca,ff1.v,ff1.(plotxstr),log10(ff1.f)',contlev,'color',[0 0 0])      
      hold(hca,'off')
    end
    if doExB % ExB
      B2 = ff1.Bx.^2 + ff1.By.^2 + ff1.Bz.^2;
      ExB = (ff1.Ex.*ff1.By - ff1.Ey.*ff1.Bx)./B2; % z-comp
      hold(hca,'on')
      plot(hca,ExB,ff1.(plotxstr),'color',[1 1 1],'linewidth',1)
      hold(hca,'off')
    end    
    hleg = irf_legend(hca,'f(v_z)',[0.02,0.98],'color',[0 0 0],'backgroundcolor',[1 1 1]);     
  end
  end
  if 1 % arclength on horizontal axis
    if 1 % fx5
      hca = h(isub); isub = isub + 1;
      ff1 = fx5{itag};
      plotf = log10(ff1.f);
      plotf(plotf==-Inf) = NaN;
      if doPlotSurf
        xx = ff1.(plotxstr);
        dxx = diff(xx);
        xx_edges = [xx(1)-0.5*dxx(1) xx(1:end-1)+0.5*dxx xx(end)+0.5*dxx(end)];
        yy = ff1.v';
        dyy = diff(yy);
        yy_edges = [yy(1)-0.5*dyy(1) yy(1:end-1)+0.5*dyy yy(end)+0.5*dyy(end)];
        [XX,YY] = meshgrid(xx_edges,yy_edges);            
        surf(hca,XX,YY,YY*0,plotf); 
        view(hca,[0 0 1])
      else
        pcolor(hca,ff1.(plotxstr),ff1.v,plotf); 
      end
      shading(hca,'flat'); 
      hca.YLabel.String = 'v_x'; 
      if doCont % contour lines
        hold(hca,'on')      
        contour(hca,ff1.v,ff1.(plotxstr),log10(ff1.f)',contlev,'color',[0 0 0])      
        hold(hca,'off')
      end
      if doExB % ExB
        B2 = ff1.Bx.^2 + ff1.By.^2 + ff1.Bz.^2;
        ExB = (ff1.Ey.*ff1.Bz - ff1.Ez.*ff1.By)./B2; % x-comp
        hold(hca,'on')
        plot(hca,ff1.(plotxstr),ExB,'color',ExBcol,'linewidth',ExBlinewidth)
        hold(hca,'off')
      end
      colormap(hca,cmap); 
      hca.CLim = fclim;    
      %hleg = irf_legend(hca,'f(v_x)',[0.02,0.98],'color',[0 0 0],'backgroundcolor',[1 1 1]); 
    end
    if 1 % fy5
      hca = h(isub); isub = isub + 1;
      ff1 = fy5{itag};
      plotf = log10(ff1.f);
      plotf(plotf==-Inf) = NaN;
      pcolor(hca,ff1.(plotxstr),ff1.v,plotf); 
      shading(hca,'flat'); 
      colormap(hca,cmap); 
      hca.CLim = fclim;
      hca.YLabel.String = 'v _y'; 
      if doCont % contour lines
        hold(hca,'on')      
        contour(hca,ff1.(plotxstr),ff1.v,log10(ff1.f),contlev,'color',[0 0 0])      
        hold(hca,'off')
      end
      if doExB % ExB
        B2 = ff1.Bx.^2 + ff1.By.^2 + ff1.Bz.^2;
        ExB = (ff1.Ez.*ff1.Bx - ff1.Ex.*ff1.Bz)./B2; % y-comp
        hold(hca,'on')
        plot(hca,ff1.(plotxstr),ExB,'color',ExBcol,'linewidth',ExBlinewidth)
        hold(hca,'off')
      end
      %hleg = irf_legend(hca,'f(v_y)',[0.02,0.98],'color',[0 0 0],'backgroundcolor',[1 1 1]); 
    end
    if 1 % fz5
      hca = h(isub); isub = isub + 1;
      ff1 = fz5{itag};
      plotf = log10(ff1.f);
      plotf(plotf==-Inf) = NaN;
      pcolor(hca,ff1.(plotxstr),ff1.v,plotf); 
      shading(hca,'flat'); 
      colormap(hca,cmap); 
      hca.CLim = fclim;
      hca.YLabel.String = 'v_z'; 
      if doCont % contour lines
        hold(hca,'on')      
        contour(hca,ff1.(plotxstr),ff1.v,log10(ff1.f),contlev,'color',[0 0 0])      
        hold(hca,'off')
      end
      if doExB % ExB
        B2 = ff1.Bx.^2 + ff1.By.^2 + ff1.Bz.^2;
        ExB = (ff1.Ex.*ff1.By - ff1.Ey.*ff1.Bx)./B2; % z-comp
        hold(hca,'on')
        plot(hca,ff1.(plotxstr),ExB,'color',ExBcol,'linewidth',ExBlinewidth)
        hold(hca,'off')
      end    
      %hleg = irf_legend(hca,'f(v_z)',[0.02,0.98],'color',[0 0 0],'backgroundcolor',[1 1 1]);     
    end
  end
end


c_eval('h(?).XLabel.String = ''s'';',1:npanels)
c_eval('h(?).YLabel.String = '''';',4:npanels)
c_eval('h(?).YTickLabels = [];',4:npanels)
c_eval('h(?).Position(1) = h(?).Position(1)-0.05;',1:npanels)
%c_eval('colormap(h(?),pic_colors(''pasteljet''));',1:npanels) 
compact_panels(h,0.0,0.02)
hlinks = linkprop(h,{'XLim','YLim','XTick'});
h(1).YLim = [-2.5 2.5];
h(1).YLim = 1.9*[-1 1];
h(1).XTick = [-10:2:10];
h(1).XLim = 9.5*[-1 1];

hpos = h(end).Position;
hcb = colorbar('peer',h(end));
hcb.Position(4) = hpos(4)*3;
hcb.Position(1) = hpos(1)+hpos(3)+0.01;
hcb.Position(3) = 0.015;
hcb.YLabel.String = 'log_{10}f_{ic}^{top}(s,v)';

c_eval('h(?).XGrid = ''on'';',1:npanels)
c_eval('h(?).YGrid = ''on'';',1:npanels)
c_eval('h(?).Layer = ''top'';',1:npanels)

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
legends = {'i)','j)','k)','l)','m)','n)','o)'};
legends = {'I','J','K','L','M','N','O'};
for ip = 1:npanels
  %irf_legend(h(ip),{legends{ip}},[0.02 0.98],'color',[0 0 0],'fontsize',12)
  irf_legend(h(ip),{legends{ip}},[0.02 0.98],'color',[0 0 0],'fontsize',19,'fontweight','bold')
end

hl_exb = findobj(h(1),'type','line');
hleg_exb = legend(hl_exb,'v_{E\times B}','Box','off','location','southwest');
hleg_exb.Box = 'off';
%
an_fontsize = 12;

hleg = irf_legend(h(1),{'north';'B_x>0'},[0.10,0.98],'color',[0 0 0],'fontsize',an_fontsize);
hleg = irf_legend(h(1),{'south';'B_x<0'},[0.98,0.98],'color',[0 0 0],'fontsize',an_fontsize);
hleg = irf_legend(h(4),{'north';'B_x>0'},[0.10,0.98],'color',[0 0 0],'fontsize',an_fontsize);
hleg = irf_legend(h(4),{'south';'B_x<0'},[0.98,0.98],'color',[0 0 0],'fontsize',an_fontsize);

delete(findall(gcf,'type','annotation'))

hleg = irf_legend(h(2),{'stronger';'magnetization:';'complete gyromotion'},[0.02,0.04],'color',[0 0 0],'fontsize',an_fontsize);
annotation('textarrow',[0.30 0.30],[0.59 0.54],'String',{'v_y>0 \rightarrow -v_yB_x>0'},'fontsize',an_fontsize,'fontweight','light')
annotation('textarrow',[0.34 0.33],[0.46 0.50],'String',{'v_y<0 \rightarrow -v_yB_x<0'},'fontsize',an_fontsize,'fontweight','light','horizontalalignment','center')
annotation('textarrow',[0.29 0.29],[0.31 0.27],'String',{'reflection due to -v_yB_x>0'},'fontsize',an_fontsize,'fontweight','light','horizontalalignment','left')
annotation('textarrow',[0.32 0.32],[0.17 0.21],'String',{'outward force due to -v_yB_x<0'},'fontsize',an_fontsize,'fontweight','light','horizontalalignment','left')

hleg = irf_legend(h(5),{'weaker magnetization:';'incomplete gyromotion'},[0.02,0.04],'color',[0 0 0],'fontsize',an_fontsize);
annotation('textarrow',[0.68 0.68],[0.59 0.55],'String',{'v_y>0 \rightarrow -v_yB_x>0'},'fontsize',an_fontsize,'fontweight','light','horizontalalignment','left')
%annotation('textarrow',[0.70 0.7],[0.17 0.21],'String',{'reflection due to -v_yB_x<0'},'fontsize',an_fontsize,'fontweight','light','horizontalalignment','left')
%annotation('textarrow',[0.72 0.7],[0.30 0.28],'String',{'reflection due to','-v_yB_x<0'},'fontsize',an_fontsize,'fontweight','light','horizontalalignment','center')
%annotation('textarrow',[0.73 0.72],[0.22 0.24],'String',{'simultaneous','scattering'},'fontsize',an_fontsize,'fontweight','light','horizontalalignment','center')
annotation('textarrow',[0.74 0.72],[0.25 0.25],'String',{'reflection','due to','-v_yB_x>0','+','simultaneous','scattering'},'fontsize',an_fontsize,'fontweight','light','horizontalalignment','center')

%annotation('textarrow',[0.30 0.30],[0.15 0.2],'String',{'-v_y+B_x<0'},'fontsize',12,'fontweight','light')
%if strcmp(plotxstr,'arc_z0'), %c_eval('h(?).YDir = ''reverse'';',remainingrows); end
set(gcf,'position',[104         782         616*1.2         331*1.2])

%% Figure 6.5. 2D reduced distributions, showing beams in different perspective
% Pick distributions
twpe = 24000;
ds = ds100.twpelim(twpe).xlim([70 70]).zlim([0 8]).dxlim([0.3 0.7]).findtag('line horizontal');
nrows = 3;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);


%% Figure 7, f(vy), py
%% Reduced distributions, plot
twpe = 24000;
pic = no02m.twpelim(twpe);
fred = fred35_z0;
twpe = 24000;
fontsize = 12;
% What to include
% - overview of whxx  ere boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];
zlim_line = [-0.5 0.5];
xlim_line = [min(fred3_z0.x) max(fred3_z0.x)];

py0 = [8 7 6.5 5.8 5.2 4.6];
py0 = [8 6.9 6.3 5.6 5.1 4.6];
py0 = [8 6.9 6.3 5.6 5.1 4.6];
py0 = [8 7 6.5 5.8 5.4 4.6];
xmin = [70 70 71 73 80 82];
xmax = [80 85 87 93 95 100];


py0 = [8.5 8.1 7.7 7.4 6.9 6.3 5.8 5.2 4.6];
xmin = [71 71 71 71 71 71.5 73 80 83];
xmax = [74 77 79 81 84 87 92 96 97]+00;
npy = numel(py0);

nrows = 1;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 1; colorE = [0 0.8 0.8];
doV = 1; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];
doTraj = 1; colorTraj = [0 0 0];
trs = tr100.find([tr100.z0] == 0,[tr100.vy0] > 0.5,[tr100.x0] > 73,[tr100.xstart] < 100);
trs = tr100.find([tr100.z0] == 0,[tr100.vy0] > 0.5).lim('t',[23000 25000]/200);

cmap_dist = pic_colors('waterfall');

%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};

if 0 % A(x),
  hca = h(isub); isub = isub + 1; 
  pic_tmp = pic.xlim(xlim_line).zlim(zlim_line);
  plot(hca,pic_tmp.xi,mean(pic_tmp.A,2))  
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'A_{y}';    
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % fi(v_y)
  hca = h(isub); isub = isub + 1;
  fred = fred35_z0;    
  pcolor(hca,fred.x,fred.v,log10(fred.fvy)')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('candy4'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'f_{ic}(x,v_{y})';  
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  hca.CLim = 0.99*[-4 2];
  hca.CLim = 0.99*[-4 2];
  hca.YLim = 0.99*4*[-1 1];
  
  for ipy = 1:npy
    hold(hca,'on')
    pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
    A = squeeze(mean(pic_tmp.A,2));    
    xx = pic_tmp.xi;
    yy = f_vy_A(py0(ipy),A);
    plot(hca,xx,yy,':k','linewidth',0.5)
    hold(hca,'off')
    ht = text(hca,xx(end),yy(end),sprintf(' %.1f',py0(ipy)),'fontsize',fontsize);
    ht.HorizontalAlignment = 'left';
    dy = (yy(end)-yy(end-1))./diff(hca.YLim);
    dx = (xx(end)-xx(end-1))./diff(hca.XLim);
    ht.Rotation = atand(dy/dx);
    
  end
  if 0*doE
    hold(hca,'on')
    plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
    hold(hca,'off')
  end
  if 0*doV
    hold(hca,'on')
    plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if 0 % doExB
    hold(hca,'on')
    xx = eval(['x_z' num2str(unique(fred.z))]);
    vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
    plot(hca,xx,vv,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi(v_y)
  hca = h(isub); isub = isub + 1;
  fred = fred3_z0;    
  pcolor(hca,fred.x,fred.v,log10(fred.fvy)')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('candy4'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'f_{i,cold}(x,v_{y})';  
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  hca.CLim = 0.99*[-4 2];
  hca.CLim = 0.99*[-2 1];
  hca.YLim = 0.99*4*[-1 1];
  
%   for ipy = 1:npy
%     hold(hca,'on')
%     pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
%     A = squeeze(mean(pic_tmp.A,2));
%     plot(hca,pic_tmp.xi,f_vy_A(py0(ipy),A),'--k')
%     hold(hca,'off')
%   end

   for itr = 1:trs.ntr
      hold(hca,'on')
      plot(hca,trs(itr).x,trs(itr).vy,'k')
      hold(hca,'off')
   end
end
if 0 % fi_top(v_y)/fi_tot(v_y)
  hca = h(isub); isub = isub + 1;
  fred = fred3_z0;    
  pcolor(hca,fred.x,fred.v,fred3_z0.fvy'./fred35_z0.fvy')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('waterfall'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'f_{i,top}/f_{i,tot}(x,v_{y})';  
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  hca.CLim = 1*[0 1];  
  hca.YLim = 0.99*4*[-1 1];
  
%   for ipy = 1:npy
%     hold(hca,'on')
%     pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
%     A = squeeze(mean(pic_tmp.A,2));
%     plot(hca,pic_tmp.xi,f_vy_A(py0(ipy),A),'--k')
%     hold(hca,'off')
%   end

   for itr = 1:trs.ntr
      hold(hca,'on')
      plot(hca,trs(itr).x,trs(itr).vy,'k')
      hold(hca,'off')
   end
end
if 0 % A0(v_y,A)
  hca = h(isub); isub = isub + 1; 
  A_ = no02m.interp(fred.x,fred.z,no02m.twpelim(twpe).twci,'A');
  [VY,A] = meshgrid(fred.v,A_);  
  A0map = f_A0_vy_py0(A,VY);
  A0map(fred.fvy==0) = NaN;
  %[Ccont,hcont] = contourf(hca,fred.x,fred.v,A0map',-30:1:30);
  [Ccont,hcont] = contourf(hca,fred.x,fred.v,A0map',[0 sort(py0)]);
  clabel(Ccont,hcont,'LabelSpacing',72,'Color','k','FontWeight','bold');
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('waterfall'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'A_{y0}(A_{y,loc},v_y,v_{y0}=0)';  
  irf_legend(hca,{'A_{y0} = A_{y,loc} - v_{y,loc}';'assuming v_{y0}=0'},[0.98 0.1],'color',[0 0 0])
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  %hca.CLim = 0.99*[-4 2];
  %hca.CLim = 0.99*[-2 1];
  hca.YLim = 0.99*4*[-1 1];
  
  for ipy = 1:npy
    hold(hca,'on')
    pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
    A = squeeze(mean(pic_tmp.A,2));
    xx = pic_tmp.xi;
    yy = f_vy_A(py0(ipy),A);
    plot(hca,xx,yy,'--k')
    hold(hca,'off')
    ht = text(hca,xx(end),yy(end),sprintf('p_y = %.1f',py0(ipy)),'fontsize',fontsize);
  end
  if 0*doE
    hold(hca,'on')
    plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
    hold(hca,'off')
  end
  if 0*doV
    hold(hca,'on')
    plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if 0 % doExB
    hold(hca,'on')
    xx = eval(['x_z' num2str(unique(fred.z))]);
    vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
    plot(hca,xx,vv,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % ExB_x
  hca = h(isub); isub = isub + 1; 
  pic_tmp = pic.xlim(xlim_line).zlim(zlim_line);
  plot(hca,pic_tmp.xi,mean(pic_tmp.vExBx,2),pic_tmp.xi,mean(pic_tmp.vix,2))
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{x}';    
  legend(hca,{'v_{ExB}','v_i'},'location','best')
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 0 % A(x), for different times
  hca = h(isub); isub = isub + 1; 
  pic_tmp = no02m.twpelim(19000:1000:24000,'exact').xlim(xlim_line).zlim(zlim_line);
  plot(hca,pic_tmp.xi,squeeze(mean(pic_tmp.A,2)))  
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'A_{y}';    
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim'});
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
for ip = 1:nrows*ncols
  axwidth(ip) = h(ip).Position(3);
  h(ip).FontSize = 14;
end
for ip = 1:nrows*ncols
  h(ip).Position(3) = 0.7;%min(axwidth);
end
% for ip = 1:nrows*ncols
%   h(ip).Position(2) = h(ip).Position(2)-0.05;
% end

%c_eval('h(?).YTickLabel = []; h(?).YLabel = [];',[2 3 5 6 8 9])
%c_eval('h(?).YTick = -10:1:10;',1)


%% Figure 8, Energy content
nrows = 1;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

Utot = no02m.UB+no02m.Uke+no02m.Uki+no02m.Ute+no02m.Uti;
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,no02m.twci,no02m.UB,...
           no02m.twci,no02m.Uti,...
           no02m.twci,no02m.Ute,...
           no02m.twci,no02m.Uki,...
           no02m.twci,no02m.Uke,...
           no02m.twci,Utot)
legend(hca,{'U_B','U_{ti}','U_{te}','U_{ki}','U_{ke}','U_{tot}'},'location','best')
end

%% Figure 9, streaming instability
%% Figure 9, prepare data
twpe = 24000; xlim = [130 110]; zlim = [-8 8]; % later time
twpe = 23000; xlim = [64 80]; zlim = [-8 8]; % earlier time

ds = ds100.twpelim(twpe).findtag({'A=7.5'});
%fpar3 = ds.fpar()

xdist = (ds.xi1{1}+ds.xi2{1})/2;
zdist = (ds.zi1{1}+ds.zi2{1})/2;

pic_lim = no02m.xlim(xlim).zlim(zlim).twpelim(twpe);
pic = no02m.twpelim(twpe);
Bx_ = pic.Bx;
By_ = pic.By;
Bz_ = pic.Bz;

Bx = interpfield(pic.xi,pic.zi,Bx_,xdist,zdist); 
By = interpfield(pic.xi,pic.zi,By_,xdist,zdist); 
Bz = interpfield(pic.xi,pic.zi,Bz_,xdist,zdist);
%get_points(obj,x,z,t,range,field,varargin)

% reduced distributions
if 0 % saved reduced distributions
  % save('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/matlab/fred_twpe23000_A75.mat','fred35_A75','fred46_A75','fred3_A75')
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/matlab//fred_twpe23000_A75.mat')
elseif 0 % make reduced distributions
  fred35_A75 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  fred3_A75 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  fred46_A75 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});  
elseif 1 
  fpar3 = ds.fpar(1,:,3);
  fpar5 = ds.fpar(1,:,5);
  fpar35 = ds.fpar(1,:,[3 5]);
end

% Scalar quantities
fredi_str = '3'; iSpecies = [3];
frede_str = '3'; eSpecies = [3];
fredi = eval(['fred' fredi_str '_A75']);
frede = eval(['fred' frede_str '_A75']);

arclength = [0; cumsum(sqrt(diff(fredi.x).^2 + diff(fredi.z).^2))];
arclength = fpar3.arc;
%arclength = arclength(end:-1:1);
if 1; arclength = arclength - arclength(find(abs(fredi.z)==min(abs(fredi.z)))); end

arclength = fpar3.arc_z0;

darc = arclength(2)-arclength(1);
arcedges = [arclength(1)-0.5*darc; arclength+0.5*darc];
narc = 1800; % 900
arclength_interp = linspace(arclength(1),arclength(end),narc);
darc_interp = (arclength_interp(end)-arclength_interp(1))/narc;
xinterp = interp1(arclength,xdist,arclength_interp);
zinterp = interp1(arclength,zdist,arclength_interp);

ni = interpfield(pic.xi,pic.zi,pic.ni,xinterp,zinterp); 
ne = interpfield(pic.xi,pic.zi,pic.ne,xinterp,zinterp); 
vepar = interpfield(pic.xi,pic.zi,pic.vpar([3 5]),xinterp,zinterp); 
Epar = interpfield(pic.xi,pic.zi,pic.Epar,xinterp,zinterp);
Babs = interpfield(pic.xi,pic.zi,pic.Babs,xinterp,zinterp); 
intEpar = [-cumsum(Epar*darc_interp)];

% % Set up for instability analysis
arcval = 7; % twpe=23000, at neutral plane
%arccenter = (arclength(end)-arclength(1))/2;
%idist = find(abs(arclength-arcval-arccenter)==min(abs(arclength-arcval-arccenter)));
idist = find(abs(arclength-arcval)==min(abs(arclength-arcval)));
xzarcdist = [xdist' zdist' arclength];
xval = xdist(idist);
zval = zdist(idist);
ds_pick = ds.twpelim(twpe).xfind(xval).zfind(zval);
%
% Fit to distribution
qe = 1;
n = [0.05 0.08 0.15 0.11 0.13 0.04]; % n0
vd = [0.5 1.15 0.05 -2.1 2.7 0.5]; % vA0
vt = [1.5 0.3 0.07 2.8 3.0 4.0]; % vA0
m = [1 1 1 1/100 1/100 1/100]*1; % m0
q = [1 1 1 -1 -1 -1]*qe;
iIncl = [1 2 3];
eIncl = [4 5 6];
% n = sqrt(n*e^2/m*eps0)
wpewce = 2;
mime = 200;
eps0 = 1/(wpewce^2*mime);
wp = sqrt(n./m)*sqrt(1/eps0);

% Prepare simulation distributions
clear f
for iSpecies = 1:nSpecies
  f(iSpecies) = ds_pick.fxyz(1,1,iSpecies);
end
for iSpecies = 1:nSpecies
  f(iSpecies).n_map = mean(mean(pic.xlim(f(iSpecies).x).zlim(f(iSpecies).z).n(iSpecies)));
  f(iSpecies).n_dist = sum(f(iSpecies).f(:))*f(iSpecies).dv^3;
end
if 1 % takes some time, so not necessary to do many times if the dist is the same
%%
  clear f_rot
for iSpecies = 1:nSpecies
  x_tmp = f(iSpecies).x;
  z_tmp = f(iSpecies).z;
  pic_tmp = pic.xlim(x_tmp).zlim(z_tmp);
  Bx_tmp = mean(mean(pic_tmp.Bx));
  By_tmp = mean(mean(pic_tmp.By));
  Bz_tmp = mean(mean(pic_tmp.Bz));
  r1 = [Bx_tmp,By_tmp,Bz_tmp]; r1 = r1/sqrt(sum(r1.^2));
  r2 = cross(r1,cross([0 1 0],r1)); r2 = r2/sqrt(sum(r2.^2));
  r3 = cross(r1,r2); r3 = r3/sqrt(sum(r3.^2));
  f_rot(iSpecies) = rotate_dist(f(iSpecies),r1,r2,r3,1); % last entry after r3 is doInterp;
end
end
ve = f_rot(iEle(1)).v;    
iSpecies = 2; fehot = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
iSpecies = 4; fecoldtop = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
iSpecies = 6; fecoldbot = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
fetot = fehot + fecoldtop + fecoldbot;
vi = f_rot(iIon(1)).v;   
iSpecies = 1; fihot = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
iSpecies = 3; ficoldtop = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
iSpecies = 5; ficoldbot = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
fitot = fehot + ficoldtop + ficoldbot;

fefit = fmax1D(ve,n(eIncl),vd(eIncl),vt(eIncl));
fifit = fmax1D(vi,n(iIncl),vd(iIncl),vt(iIncl));
    
%% Plot 
fi_clim = [0 0.0499];
fe_clim = [0 1.3e-2];
colors = pic_colors('matlab');

nrows = 4;
ncols = 2;
h = setup_subplots(nrows,ncols,'vertical');
%[h,h2] = initialize_combined_plot(nrows,2,2,0.4,'vertical')
isub = 1;
doE = 0; colorE = [0.0 0.0 0.0];
doV = 0; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
isMap = [];


if 1 % Epar
  hca = h(isub); isub = isub + 1;  
  %plot(hca,arclength_interp,Epar_,arclength_interp,ni_,arclength_interp,ne_,arclength_interp,(ni_-ne_)*10)  
  %legend(hca,{'E_{||}','n_{i,cold}','n_{e,cold}','ni-ne'},'location','eastoutside')   
  plot(hca,arclength_interp,Epar,'k')  
  hold(hca,'on')
  %plot(hca,smooth(arclength_interp,50),smooth(Epar,50)*1,'linewidth',2,'color',colors(6,:))
  hold(hca,'off')
  %plot(hca,arclength_interp,Epar,'k',smooth(arclength_interp,100),smooth(Epar,100),'linewidth',1)  
  %plot(hca,arclength_interp,Epar,'k',arclength_interp,intEpar*3,arclength,Babs/2)  
  %legend(hca,{'E_{||}'},'location','eastoutside')   
  hca.YLabel.String = 'E_{||}';   
  hca.XLabel.String = 's_{||} (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLim = 0.399*[-1 1];
end

if 0 % dn
  hca = h(isub); isub = isub + 1;  
  %plot(hca,arclength_interp,Epar_,arclength_interp,ni_,arclength_interp,ne_,arclength_interp,(ni_-ne_)*10)  
  %legend(hca,{'E_{||}','n_{i,cold}','n_{e,cold}','ni-ne'},'location','eastoutside')   
  plot(hca,arclength_interp,(ni_-ne_))  
  legend(hca,{'n_{i,cold}-n_{e,cold}'},'location','eastoutside')   
  hca.YLabel.String = 'n';   
  hca.XLabel.String = 'arclength (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Epar, n
  hca = h(isub); isub = isub + 1;  
  %plot(hca,arclength_interp,Epar_,arclength_interp,ni_,arclength_interp,ne_,arclength_interp,(ni_-ne_)*10)  
  %legend(hca,{'E_{||}','n_{i,cold}','n_{e,cold}','ni-ne'},'location','eastoutside')   
  plot(hca,arclength_interp,Epar_,arclength_interp,(ni_-ne_)*10)  
  legend(hca,{'E_{||}','n_{i,cold}-n_{e,cold}'},'location','eastoutside')   
  hca.YLabel.String = 'E_{||}, n';   
  hca.XLabel.String = 'arclength (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred35_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fred.fvpar')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,cold}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(fred.fvpar(:),99);
  hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi3(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fred.fvpar')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,cold}^{top}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(fred.fvpar(:),99);
  hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred35_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  
  if 0 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.vpar_center,log10(fred.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{i,cold}(s_{||},v_{||})'];
  hcb.YLabel.String = ['log_{10}f_{ic}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = [-5 -0];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % log 10 fi3(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{i,cold}^{top}(l_{||},v_{||})'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = [-6 -1];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % fi3/fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fredall = fred35_A75;
  fred = fred3_A75;
  fplot = (fred.fvpar./fredall.fvpar);
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fplot')
  view(hca,[0 0 1]);   
  shading(hca,'flat')
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.vpar_center,log10(fredall.fvpar)',[-1 -1],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 's_{||} (d_i)';
  hca.YLabel.String = 'v_{||}';
  %colormap(hca,pic_colors('candy')) 
  %colormap(hca,pic_colors('blue_red'))
  colormap(hca,pic_colors('pasteljet'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,cold}^{top}/f_{i,cold}^{tot}(s_{||},v_{||})'];
  hcb.YLabel.String = ['f_{ic}^{top}/f_{ic}^{tot}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.CLim = [0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fe(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred46_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')  
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fred.fvpar')
  view(hca,[0 0 1]);  
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{e,cold}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(fred.fvpar(:),99);
  hca.CLim = fe_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-7 7];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vepar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fe(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred46_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')  
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]);  
  shading(hca,'flat')
  hca.XLabel.String = 's_{||} (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{e,cold}(s_{||},v_{||})'];
  hcb.YLabel.String = ['log_{10}f_{ec}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fe_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-7 7];
  hca.CLim = [-3 -0.8];
  irf_legend(hca,{'v_{ec,||}'},[0.02 0.4],'color',[0 0 0])
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if 1%doV
    hold(hca,'on')
    plot(hca,arclength_interp,vepar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end

if 1 % line position on map, Epar
  isMap(end+1) = isub; 
  hca = h(isub); isub = isub + 1;
  %imagesc(hca,pic_lim.xi,pic_lim.zi,smooth2(pic_lim.Epar,3)');
  %imagesc(hca,pic_lim.xi,pic_lim.zi,pic_lim.Epar');
  imagesc(hca,pic_lim.xi,pic_lim.zi,pic_lim.n(5)');
  colormap(hca,pic_colors('blue_red'));
  hcb = colorbar('peer',hca);
  hca.CLim = max(max(get(findobj(hca.Children,'type','Image'),'CData')))*[-1 1];
  hcb.YLabel.String = 'E_{||}';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.CLim = 0.5*[-1 1];
  if 1 % plot_boxes
    %%
    hold(hca,'on')
    plot(hca,xdist,zdist,'k','linewidth',1)
    %ds.plot_boxes(hca)
    %ds_pick.plot_boxes(hca,'color',colors(3,:),'linewidth',1)
    %ds_pick.plot_boxes(hca)
    hl = plot(hca,xval,zval,'Marker','s','MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
    hold(hca,'off')
  end
  if 1 % plot A
    hold(hca,'on')
    A = pic_lim.A;
    clim = hca.CLim;
    contour(hca,pic_lim.xi,pic_lim.zi,A',[0:1:25],'color',0.5*[1 1 1])
    hca.CLim = clim;
    hold(hca,'off')
  end
end

% Distributions
if 1 % f(vpar), ions
  hca = h(isub); isub = isub + 1;
  
  plot(hca,vi,fihot,vi,ficoldtop,vi,ficoldbot,vi,fitot,vi,fifit,'--','linewidth',1)  
  hca.XLabel.String = 'v_{||} (v_{A0})';
  hca.XLabel.String = 'v_{||}';
  hca.YLabel.String = 'f_i(v_{||}) (d_{i0}^3v_{A0})';
  hca.YLabel.String = 'f_i(v_{||})';
  %legend(hca,{'f_{hot}','f_{cold}^{top}','f_{cold}^{top}','f_{tot}','f_{fit}'},'location','northeast','box','off')
  irf_legend(hca,{'hot','cold, top','cold, bot','tot','fit'}',[0.95 0.95])
  hca.XLim = 0.99*[-1 2];
  irf_legend(hca,{'ions'},[0.02 0.98],'color',[0 0 0])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLim = [0 1.2];
  hca.XLim = [-0.5 1.9];
end
if 1 % f(vpar), electrons
  hca = h(isub); isub = isub + 1;
  plot(hca,ve,fehot,ve,fecoldtop,ve,fecoldbot,ve,fetot,ve,fefit,'--','linewidth',1) 
  hca.XLabel.String = 'v_{||} (v_{A0})';
  hca.XLabel.String = 'v_{||}';
  %hca.YLabel.String = 'f_e(v_{||}) (d_{i0}^3v_{A0})';
  hca.YLabel.String = 'f_e(v_{||})';
  %legend(hca,{'f_{hot}','f_{cold}^{top}','f_{cold}^{top}','f_{tot}','f_{fit}'},'location','east','box','off')  
  hca.XLim = 0.99*[-10 10];
  irf_legend(hca,{'electrons'},[0.02 0.98],'color',[0 0 0])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end

% Solution to dispersion relation
if 1 % 
  hca = h(isub); isub = isub + 1;
  vph_extrap = interp1(kvec, vph_store, [0 kvec], 'spline', 'extrap');

  plot(hca,[0 kvec],[0 wr_store/10],[0 kvec],[0 wi_store],[0 kvec],vph_extrap,'linewidth',1)
  hca.XLabel.String = 'k (d_i^{-1})';
  hca.YLabel.String = '\omega, v_{ph}';
  hleg_disp = legend(hca,{'\omega_r/10','\gamma','v_{ph}'},'location','east','box','off');
  %hca.XLim = [-1 2];
  irf_legend(hca,{'dispersion relation'},[0.02 0.98],'color',[0 0 0]);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)'};
for ip = [1 2 3 4 5 6 8]
  irf_legend(h(ip),legends{ip},[-0.15 0.99],'fontsize',14,'color','k')
end
for ip = 7
  irf_legend(h(7),legends{ip},[-0.2 0.99],'fontsize',14,'color','k')
end
compact_panels(h(1:4),0.005,0)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
%fig = gcf; h_ = findobj(fig.Children,'type','axes');
%hlinks_dist = linkprop(h(setdiff(1:nrows*ncols,isMap)),{'XLim'});
h_all = findobj(get(gcf,'Children'));
c_eval('try; h_all(?).FontSize = 12; end',1:numel(h_all))
hleg_disp.FontSize = 10;

hlinks_arc = linkprop(h(1:4),{'XLim'});
hlinks_arc.Targets(1).XLim = 14.0*[-1 1];

%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
for ip = 1:npanels
  h(ip).FontSize = 13;
end
for ip = 1:4
  h(ip).Position(1) = 0.08;
  h(ip).Position(3) = 0.4;
end
for ip = 6:8
  h(ip).Position(1) = 0.70;
  h(ip).Position(3) = 0.25;
end
h(5).Position(1) = 0.58;
h(5).Position(3) = 0.25;

for ip = 1:4 % mark single distribution location
  hold(h(ip),'on');
  plot(h(ip),arcval*[1 1],h(ip).YLim,'k--')
  hold(h(ip),'off');
end
if 1 % plot vphmax  and lambda max
  %%
%   hold(h(2),'on')
%   hvph2 = plot(h(2),arcval,vphmax,'k*');
%   hold(h(2),'off')
%   hold(h(3),'on')
%   hvph3 = plot(h(3),arcval,vphmax,'k*');
%   hold(h(3),'off')
  lambdamax = 2*pi/kmax;
  hold(h(2),'on')
  hlam2 = errorbar(h(2),arcval,vphmax,lambdamax/2,...
    'horizontal','color',[0 0 0],'linewidth',1.5);  
  hold(h(2),'off')
  hold(h(3),'on')
%   hlam_dummy = errorbar(h(3),arcval,vphmax,lambdamax/2,...
%     'horizontal','color',[0 0 0]);
%   hvph3 = plot(h(3),arcval,vphmax,'k*');
  hold(h(3),'off')
  
  %htext = text(h(2),0.6,0.2,'v_{ph} and \lambda at \gamma_{max}','Units','normalized','fontsize',13);
  %annotation('arrow',0.38*[1 1],[0.58 0.61])
  %htext = annotation('textarrow',0.38*[1 1],[0.58 0.61],'String','v_{ph} and \lambda at \gamma_{max}','fontsize',13);
  htext = annotation('textarrow',0.38*[1 1],[0.68 0.65],'String','v_{ph} and \lambda at \gamma_{max}','fontsize',13);
%   hold(h(2),'on')
%   quiver(h(2),arcval,-1,0,0.8,10)
%   hold(h(2),'off')
  %hlegvph = legend([hvph2 hlam2],{'v_{ph} at \gamma_{max}','\lambda at \gamma_{max}'},'edgecolor','white','location','southeast');
end

if 1 % annotations
  %%
  hann_trap = annotation('textarrow',0.40*[1 1],[0.38 0.4],'String',{'only cold ions','from the bottom'},'fontsize',13);
  hann_cbot = annotation('textarrow',0.40*[1 1],[0.38 0.4],'String',{'only cold ions','from the bottom'},'fontsize',13);
  hann_ctop = annotation('textarrow',0.152*[1 1],[0.45 0.43],'String',{'only cold ions','from the top'},'fontsize',13);
  hann_trap = annotation('textarrow',0.365*[1 1],[0.58 0.61],'String',{'ion-wave','interaction'},'fontsize',13);
  hann_esw = annotation('textarrow',[0.33 0.36],[0.826 0.81],'String',{'ESWs'},'fontsize',13);
  hann_esw2 = annotation('arrow',[0.27 0.25],[0.826 0.81]);
  
end
if 0 % potential velocity in fipar plot
  %%
  hold(h(2),'on')
  plot(h(2),arclength_interp,sqrt(2*abs(intEpar)))
  hold(h(2),'off')
end
% for ip = 1:nrows*ncols
%   axwidth(ip) = h(ip).Position(3); 
% end
% for ip = 1:nrows*ncols
%   h(ip).Position(3) = min(axwidth);
% end
if 1 % compact figure to ake Epar arcpanel shorter, waves show off better like that
  %%
  h(1).Position(4) = 0.12;
  
  h(5).Position(4) = 0.14;
  h(6).Position(4) = 0.13;
  h(7).Position(4) = 0.13;
  h(8).Position(4) = 0.13;
  
  h(5).Position(2) = 0.71;
  h(6).Position(2) = 0.51;
  h(7).Position(2) = 0.31;
  h(8).Position(2) = 0.11;
  
end
%% Figure 9, ALT2, prepare data + plot, with premade fpar
twpe = 24000; xlim = [130 110]; zlim = [-8 8]; % later time
twpe = 23000; xlim = [64 80]; zlim = [-8 8]; % earlier time

ds = ds100.twpelim(twpe).findtag({'A=7.5'});
%fpar3 = ds.fpar()

xdist = (ds.xi1{1}+ds.xi2{1})/2;
zdist = (ds.zi1{1}+ds.zi2{1})/2;

pic_lim = no02m.xlim(xlim).zlim(zlim).twpelim(twpe);
pic = no02m.twpelim(twpe);
Bx_ = pic.Bx;
By_ = pic.By;
Bz_ = pic.Bz;

Bx = interpfield(pic.xi,pic.zi,Bx_,xdist,zdist); 
By = interpfield(pic.xi,pic.zi,By_,xdist,zdist); 
Bz = interpfield(pic.xi,pic.zi,Bz_,xdist,zdist);

% reduced distributions
fpar3 = ds.fpar(1,:,3);
fpar5 = ds.fpar(1,:,5);
fepar3 = ds.fpar(1,:,4);
fepar5 = ds.fpar(1,:,6);


% Scalar quantities
arclength = fpar3.arc_z0;

darc = arclength(2)-arclength(1);
arcedges = [arclength(1)-0.5*darc, arclength+0.5*darc];
narc = 1800; % 900
arclength_interp = linspace(arclength(1),arclength(end),narc);
darc_interp = (arclength_interp(end)-arclength_interp(1))/narc;
xinterp = interp1(arclength,xdist,arclength_interp);
zinterp = interp1(arclength,zdist,arclength_interp);

ni = interpfield(pic.xi,pic.zi,pic.ni,xinterp,zinterp); 
ne = interpfield(pic.xi,pic.zi,pic.ne,xinterp,zinterp); 
vepar = interpfield(pic.xi,pic.zi,pic.vpar([3 5]),xinterp,zinterp); 
Epar = interpfield(pic.xi,pic.zi,pic.Epar,xinterp,zinterp);
Babs = interpfield(pic.xi,pic.zi,pic.Babs,xinterp,zinterp); 
intEpar = [-cumsum(Epar*darc_interp)];

% % Set up for instability analysis
arcval = 7; % twpe=23000, at neutral plane
%arccenter = (arclength(end)-arclength(1))/2;
%idist = find(abs(arclength-arcval-arccenter)==min(abs(arclength-arcval-arccenter)));
idist = find(abs(arclength-arcval)==min(abs(arclength-arcval)));
xzarcdist = [xdist' zdist' arclength'];
xval = xdist(idist);
zval = zdist(idist);
ds_pick = ds.twpelim(twpe).xfind(xval).zfind(zval);
%
% Fit to distribution
qe = 1;
n = [0.05 0.08 0.15 0.11 0.13 0.04]; % n0
vd = [0.5 1.15 0.05 -2.1 2.7 0.5]; % vA0
vt = [1.5 0.3 0.07 2.8 3.0 4.0]; % vA0
m = [1 1 1 1/100 1/100 1/100]*1; % m0
q = [1 1 1 -1 -1 -1]*qe;
iIncl = [1 2 3];
eIncl = [4 5 6];
% n = sqrt(n*e^2/m*eps0)
wpewce = 2;
mime = 200;
eps0 = 1/(wpewce^2*mime);
wp = sqrt(n./m)*sqrt(1/eps0);

% Prepare simulation distributions
clear fpar_pick
c_eval('fpar_pick{?} = ds_pick.fpar(1,:,?);',1:6)

ve_hot = fpar_pick{2}.v;
ve_cold = fpar_pick{4}.v;
iSpecies = 2; fehot = fpar_pick{iSpecies}.f;
iSpecies = 4; fecoldtop = fpar_pick{iSpecies}.f;
iSpecies = 6; fecoldbot = fpar_pick{iSpecies}.f;
fetot = fehot + fecoldtop + fecoldbot;

vi_hot = fpar_pick{1}.v;
vi_cold = fpar_pick{3}.v;
iSpecies = 1; fihot = fpar_pick{iSpecies}.f;
iSpecies = 3; ficoldtop = fpar_pick{iSpecies}.f;
iSpecies = 5; ficoldbot = fpar_pick{iSpecies}.f;
fitot = fehot + ficoldtop + ficoldbot;

fefit = fmax1D(ve_cold',n(eIncl),vd(eIncl),vt(eIncl));
fifit = fmax1D(vi_cold',n(iIncl),vd(iIncl),vt(iIncl));
    
%% Plot 
fi_clim = [0 0.0499];
fe_clim = [0 1.3e-2];
fcontlevel = [-1 1];
colors = pic_colors('matlab');

nrows = 4; 
ncols = 2;
npanels = nrows*ncols;
clear h;
h = setup_subplots(nrows,ncols,'vertical');
%[h,h2] = initialize_combined_plot(nrows,2,2,0.4,'vertical')
isub = 1;
doE = 0; colorE = [0.0 0.0 0.0];
doV = 0; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
isMap = [];


if 1 % Epar
  hca = h(isub); isub = isub + 1;  
  %plot(hca,arclength_interp,Epar_,arclength_interp,ni_,arclength_interp,ne_,arclength_interp,(ni_-ne_)*10)  
  %legend(hca,{'E_{||}','n_{i,cold}','n_{e,cold}','ni-ne'},'location','eastoutside')   
  plot(hca,arclength_interp,Epar,'k')  
  hold(hca,'on')
  %plot(hca,smooth(arclength_interp,50),smooth(Epar,50)*1,'linewidth',2,'color',colors(6,:))
  hold(hca,'off')
  %plot(hca,arclength_interp,Epar,'k',smooth(arclength_interp,100),smooth(Epar,100),'linewidth',1)  
  %plot(hca,arclength_interp,Epar,'k',arclength_interp,intEpar*3,arclength,Babs/2)  
  %legend(hca,{'E_{||}'},'location','eastoutside')   
  hca.YLabel.String = 'E_{||}';   
  hca.XLabel.String = 's (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLim = 0.399*[-1 1];
end

if 0 % dn
  hca = h(isub); isub = isub + 1;  
  %plot(hca,arclength_interp,Epar_,arclength_interp,ni_,arclength_interp,ne_,arclength_interp,(ni_-ne_)*10)  
  %legend(hca,{'E_{||}','n_{i,cold}','n_{e,cold}','ni-ne'},'location','eastoutside')   
  plot(hca,arclength_interp,(ni_-ne_))  
  legend(hca,{'n_{i,cold}-n_{e,cold}'},'location','eastoutside')   
  hca.YLabel.String = 'n';   
  hca.XLabel.String = 'arclength (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Epar, n
  hca = h(isub); isub = isub + 1;  
  %plot(hca,arclength_interp,Epar_,arclength_interp,ni_,arclength_interp,ne_,arclength_interp,(ni_-ne_)*10)  
  %legend(hca,{'E_{||}','n_{i,cold}','n_{e,cold}','ni-ne'},'location','eastoutside')   
  plot(hca,arclength_interp,Epar_,arclength_interp,(ni_-ne_)*10)  
  legend(hca,{'E_{||}','n_{i,cold}-n_{e,cold}'},'location','eastoutside')   
  hca.YLabel.String = 'E_{||}, n';   
  hca.XLabel.String = 'arclength (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred35_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fred.fvpar')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,cold}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(fred.fvpar(:),99);
  hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi3(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fred.fvpar')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,cold}^{top}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(fred.fvpar(:),99);
  hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vpar)
  hca = h(isub); isub = isub + 1;
  ff1 = fpar3;
  ff2 = fpar5;  
  fplot = ff1.f + 1*ff2.f;
  fplot = log10(fplot);
  pcolor(hca,ff1.arc_z0,ff1.v,fplot)
  %surf(hca,arcedges,ff.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',log10(fred.fvpar)')
  %view(hca,[0 0 1]); 
  shading(hca,'flat')
  
  if 0 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.vpar_center,log10(fred.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{i,cold}(s_{||},v_{||})'];
  hcb.YLabel.String = ['log_{10}f_{ic}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = [-4.99 -0];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % log 10 fi3(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{i,cold}^{top}(l_{||},v_{||})'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = [-6 -1];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % fi3/fi35(vpar)
  hca = h(isub); isub = isub + 1;
  ff1 = fpar3;
  ff2 = fpar5;  
  fplot = ff1.f./(ff1.f + ff2.f);  
  pcolor(hca,ff1.arc_z0,ff1.v,fplot)  
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  %surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fplot')
  %view(hca,[0 0 1]);   
  shading(hca,'flat')
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,ff1.arc_z0,ff1.v,log10(ff1.f + ff2.f),fcontlevel,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 's_{||} (d_i)';
  hca.YLabel.String = 'v_{||}';
  %colormap(hca,pic_colors('candy')) 
  %colormap(hca,pic_colors('blue_red'))
  colormap(hca,pic_colors('pasteljet'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,cold}^{top}/f_{i,cold}^{tot}(s_{||},v_{||})'];
  hcb.YLabel.String = ['f_{ic}^{top}/f_{ic}^{tot}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.CLim = [0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fe(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred46_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')  
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fred.fvpar')
  view(hca,[0 0 1]);  
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{e,cold}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(fred.fvpar(:),99);
  hca.CLim = fe_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-7 7];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vepar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fe(vpar)
  hca = h(isub); isub = isub + 1;  
  ff1 = fepar3;
  ff2 = fepar5;  
  fplot = ff1.f + 1*ff2.f;
  fplot = log10(fplot);
  pcolor(hca,ff1.arc_z0,ff1.v,fplot)
  shading(hca,'flat')
  hca.XLabel.String = 's';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{e,cold}(s_{||},v_{||})'];
  hcb.YLabel.String = ['log_{10}f_{ec}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fe_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-7 7];
  hca.CLim = [-6 -0.8];
  hca.CLim = [-2 -1.4]; % high saturation to show off structures at low v
  
  irf_legend(hca,{'v_{ec,||}'},[0.02 0.4],'color',[0 0 0])
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if 1%doV
    hold(hca,'on')
    plot(hca,arclength_interp,vepar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end

if 1 % line position on map, Epar
  isMap(end+1) = isub; 
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic_lim.xi,pic_lim.zi,smooth2(pic_lim.Epar,3)');
  %imagesc(hca,pic_lim.xi,pic_lim.zi,pic_lim.Epar');
  %imagesc(hca,pic_lim.xi,pic_lim.zi,pic_lim.n(5)');
  colormap(hca,pic_colors('blue_red'));
  hcb = colorbar('peer',hca);
  hca.CLim = max(max(get(findobj(hca.Children,'type','Image'),'CData')))*[-1 1];
  hcb.YLabel.String = 'E_{||}';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.CLim = 0.5*[-1 1];
  if 1 % plot_boxes
    %%
    hold(hca,'on')
    plot(hca,xdist,zdist,'k','linewidth',1)
    %ds.plot_boxes(hca)
    %ds_pick.plot_boxes(hca,'color',colors(3,:),'linewidth',1)
    %ds_pick.plot_boxes(hca)
    hl = plot(hca,xval,zval,'Marker','s','MarkerFaceColor',[0 0 0],'Color',[0 0 0]);
    hold(hca,'off')
  end
  if 1 % plot A
    hold(hca,'on')
    A = pic_lim.A;
    clim = hca.CLim;
    contour(hca,pic_lim.xi,pic_lim.zi,A',[0:1:25],'color',0.5*[1 1 1])
    hca.CLim = clim;
    hold(hca,'off')
  end
end

% Distributions
if 1 % f(vpar), ions
  hca = h(isub); isub = isub + 1;
  
  hlines = plot(hca,vi_hot,fihot,vi_cold,ficoldtop,vi_cold,ficoldbot,vi_cold,fitot,vi_cold,fifit,'--','linewidth',1);
  hca.XLabel.String = 'v_{||} (v_{A0})';
  hca.XLabel.String = 'v_{||}';
  hca.YLabel.String = 'f_i(v_{||}) (d_{i0}^3v_{A0})';
  hca.YLabel.String = 'f_i(v_{||})';
  %legend(hca,{'f_{hot}','f_{cold}^{top}','f_{cold}^{top}','f_{tot}','f_{fit}'},'location','northeast','box','off')
  irf_legend(hca,{'hot','cold, top','cold, bot','tot','fit'}',[0.95 0.95])
  hca.XLim = 0.99*[-1 2];
  irf_legend(hca,{'ions'},[0.02 0.98],'color',[0 0 0])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLim = [0 1.2];
  hca.XLim = [-0.5 1.9];
end
if 1 % f(vpar), electrons
  hca = h(isub); isub = isub + 1;
  plot(hca,ve_hot,fehot,ve_cold,fecoldtop,ve_cold,fecoldbot,ve_cold,fetot,ve_cold,fefit,'--','linewidth',1) 
  hca.XLabel.String = 'v_{||} (v_{A0})';
  hca.XLabel.String = 'v_{||}';
  %hca.YLabel.String = 'f_e(v_{||}) (d_{i0}^3v_{A0})';
  hca.YLabel.String = 'f_e(v_{||})';
  %legend(hca,{'f_{hot}','f_{cold}^{top}','f_{cold}^{top}','f_{tot}','f_{fit}'},'location','east','box','off')  
  hca.XLim = 0.99*[-10 10];
  irf_legend(hca,{'electrons'},[0.02 0.98],'color',[0 0 0])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end

% Solution to dispersion relation
if 1 % 
  hca = h(isub); isub = isub + 1;
  vph_extrap = interp1(kvec, vph_store, [0 kvec], 'spline', 'extrap');

  plot(hca,[0 kvec],[0 wr_store/10],[0 kvec],[0 wi_store],[0 kvec],vph_extrap,'linewidth',1)
  hca.XLabel.String = 'k (d_i^{-1})';
  hca.YLabel.String = '\omega, v_{ph}';
  hleg_disp = legend(hca,{'\omega_r/10','\gamma','v_{ph}'},'location','east','box','off');
  %hca.XLim = [-1 2];
  irf_legend(hca,{'dispersion relation'},[0.02 0.98],'color',[0 0 0]);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)'};
for ip = [1 2 3 4 5 6 8]
  irf_legend(h(ip),legends{ip},[-0.15 0.99],'fontsize',14,'color','k')
end
for ip = 7
  irf_legend(h(7),legends{ip},[-0.2 0.99],'fontsize',14,'color','k')
end
compact_panels(h(1:4),0.005,0)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
%fig = gcf; h_ = findobj(fig.Children,'type','axes');
%hlinks_dist = linkprop(h(setdiff(1:nrows*ncols,isMap)),{'XLim'});
h_all = findobj(get(gcf,'Children'));
c_eval('try; h_all(?).FontSize = 12; end',1:numel(h_all))
hleg_disp.FontSize = 10;

hlinks_arc = linkprop(h(1:4),{'XLim'});
hlinks_arc.Targets(1).XLim = 14.0*[-1 1];

%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
c_eval('h(ip).FontSize = 12;',1:npanels)
for ip = 1:4
  h(ip).Position(1) = 0.08;
  h(ip).Position(3) = 0.4;
end
for ip = 6:8
  h(ip).Position(1) = 0.70;
  h(ip).Position(3) = 0.25;
end
h(5).Position(1) = 0.58;
h(5).Position(3) = 0.25;

for ip = 1:4 % mark single distribution location
  hold(h(ip),'on');
  plot(h(ip),arcval*[1 1],h(ip).YLim,'k--')
  hold(h(ip),'off');
end
if 1 % plot vphmax  and lambda max
  %%
%   hold(h(2),'on')
%   hvph2 = plot(h(2),arcval,vphmax,'k*');
%   hold(h(2),'off')
%   hold(h(3),'on')
%   hvph3 = plot(h(3),arcval,vphmax,'k*');
%   hold(h(3),'off')
  lambdamax = 2*pi/kmax;
  hold(h(2),'on')
  hlam2 = errorbar(h(2),arcval,vphmax,lambdamax/2,...
    'horizontal','color',[0 0 0],'linewidth',1.5);  
  hold(h(2),'off')
  hold(h(3),'on')
%   hlam_dummy = errorbar(h(3),arcval,vphmax,lambdamax/2,...
%     'horizontal','color',[0 0 0]);
%   hvph3 = plot(h(3),arcval,vphmax,'k*');
  hold(h(3),'off')
  
  %htext = text(h(2),0.6,0.2,'v_{ph} and \lambda at \gamma_{max}','Units','normalized','fontsize',13);
  %annotation('arrow',0.38*[1 1],[0.58 0.61])
  %htext = annotation('textarrow',0.38*[1 1],[0.58 0.61],'String','v_{ph} and \lambda at \gamma_{max}','fontsize',13);
  htext = annotation('textarrow',0.38*[1 1],[0.68 0.65],'String','v_{ph} and \lambda at \gamma_{max}','fontsize',13);
%   hold(h(2),'on')
%   quiver(h(2),arcval,-1,0,0.8,10)
%   hold(h(2),'off')
  %hlegvph = legend([hvph2 hlam2],{'v_{ph} at \gamma_{max}','\lambda at \gamma_{max}'},'edgecolor','white','location','southeast');
end
drawnow

if 1 % annotations
  %%
  annotation('textarrow',[0.67 0.71],[0.765 0.74]+0.005,'string',{'input to','dispersion','analysis'},'fontsize',12);
  annotation('textarrow',0.38*[1 1],[0.87 0.85],'string',{'input to dispersion analysis'},'fontsize',13);
  irf_legend(h(1),{'north'},[0.02 0.98],'color',[0 0 0],'fontsize',13)
  irf_legend(h(1),{'south'},[0.98 0.98],'color',[0 0 0],'fontsize',13)  
  hann_cbot = annotation('textarrow',0.40*[1 1],[0.38 0.4],'String',{'only cold ions','from the bottom'},'fontsize',13);
  hann_ctop = annotation('textarrow',0.152*[1 1],[0.45 0.43],'String',{'only cold ions','from the top'},'fontsize',13);
  hann_trap = annotation('textarrow',0.365*[1 1],[0.58 0.61],'String',{'ion-wave','interaction'},'fontsize',13);
  hann_etrap = annotation('textarrow',0.365*[1 1],[0.14 0.16],'String',{'electron-wave interaction'},'fontsize',13);
  hann_esw = annotation('textarrow',[0.33 0.36],[0.826 0.81],'String',{'ESWs'},'fontsize',13);
  hann_esw2 = annotation('arrow',[0.27 0.25],[0.826 0.81]);
  
end
if 0 % potential velocity in fipar plot
  %%
  hold(h(2),'on')
  plot(h(2),arclength_interp,sqrt(2*abs(intEpar)))
  hold(h(2),'off')
end
% for ip = 1:nrows*ncols
%   axwidth(ip) = h(ip).Position(3); 
% end
% for ip = 1:nrows*ncols
%   h(ip).Position(3) = min(axwidth);
% end
if 1 % compact figure to ake Epar arcpanel shorter, waves show off better like that
  %%
  h(1).Position(4) = 0.12;
  
  h(5).Position(4) = 0.14;
  h(6).Position(4) = 0.13;
  h(7).Position(4) = 0.13;
  h(8).Position(4) = 0.13;
  
  h(5).Position(2) = 0.71;
  h(6).Position(2) = 0.51;
  h(7).Position(2) = 0.31;
  h(8).Position(2) = 0.11;
  
end

hcball=findall(gcf,'type','colorbar');
c_eval('hcball(?).Position(3) = 0.012;',1:numel(hcball))
h(5).Position(3) = 0.28;
hcball(end).Position(1) = h(5).Position(1) + h(5).Position(3)+0.004;
hcball(end).Position(2) = h(5).Position(2);
hcball(end).Position(4) = h(5).Position(4);

h(5).CLim = 0.399*[-1 1];

%% Figure 9.2 reduced parallel ion distributiona at two different times
%% Figure 9.2 prepare data
% Distribution at 1st time
twpe1 = 23000;
pic = no02m.twpelim(twpe1);
ds1 = ds100.twpelim(twpe1).findtag({'A=7.5'});
xdist1 = (ds1.xi1{1}+ds1.xi2{1})/2;
zdist1 = (ds1.zi1{1}+ds1.zi2{1})/2;
Bx1_ = pic.Bx;
By1_ = pic.By;
Bz1_ = pic.Bz;
Bx1 = interpfield(pic.xi,pic.zi,Bx1_,xdist1,zdist1); 
By1 = interpfield(pic.xi,pic.zi,By1_,xdist1,zdist1); 
Bz1 = interpfield(pic.xi,pic.zi,Bz1_,xdist1,zdist1);
fred35_A75_1 = ds1.reduce_1d_new('x',[3 5],[],'vpar',{Bx1,By1,Bz1},'pitch',{Bx1,By1,Bz1});
fred3_A75_1 = ds1.reduce_1d_new('x',[3],[],'vpar',{Bx1,By1,Bz1},'pitch',{Bx1,By1,Bz1});
arclength1 = [0 cumsum(sqrt(diff(xdist1).^2 + diff(zdist1).^2))];
arc01 = arclength1(find(abs(zdist1)==min(abs(zdist1))));
arclength1 = arclength1 - arc01;
darcs1 = diff(arclength1);
arcedges1 = [arclength1(1)-0.5*darcs1(1) arclength1(1:end-1)+0.5*darcs1 arclength1(end)+0.5*darcs1(end)];

% Distribution at 2nd time
twpe2 = 24000;
pic = no02m.twpelim(twpe2);
ds2 = ds100.twpelim(twpe2).findtag({'A=7.5'});
xdist2 = (ds2.xi1{1}+ds2.xi2{1})/2;
zdist2 = (ds2.zi1{1}+ds2.zi2{1})/2;
Bx2_ = pic.Bx;
By2_ = pic.By;
Bz2_ = pic.Bz;
Bx2 = interpfield(pic.xi,pic.zi,Bx2_,xdist2,zdist2); 
By2 = interpfield(pic.xi,pic.zi,By2_,xdist2,zdist2); 
Bz2 = interpfield(pic.xi,pic.zi,Bz2_,xdist2,zdist2);
fred35_A75_2 = ds2.reduce_1d_new('x',[3 5],[],'vpar',{Bx2,By2,Bz2},'pitch',{Bx2,By2,Bz2});
fred3_A75_2 = ds2.reduce_1d_new('x',[3],[],'vpar',{Bx2,By2,Bz2},'pitch',{Bx2,By2,Bz2});
arclength2 = [0 cumsum(sqrt(diff(xdist2).^2 + diff(zdist2).^2))];
arc02 = arclength2(find(abs(zdist2)==min(abs(zdist2))));
arclength2 = arclength2 - arc02;
darcs2 = diff(arclength2);
arcedges2 = [arclength2(1)-0.5*darcs2(1) arclength2(1:end-1)+0.5*darcs2 arclength2(end)+0.5*darcs2(end)];

%% Parallel electric field between the two times
% obtain point along field line
twpe = [23000 24000];
xlim = [60 80];
zlim = [-6 6];
varnames = {'Epar','nic_top','nic_bot'};
varstrs = {'Epar','n(3)','n(5)'};
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
Aval = 7.5;
clear S
for it = 1:pic.nt 
  it
  pic_tmp = pic(it);
  A = squeeze(pic_tmp.A);
  S_tmp = contourcs(pic_tmp.xi,pic_tmp.zi,A',Aval*[1 1]);
  % Interpolate equidistant
  d_arc_x = diff(S_tmp.X);
  d_arc_y = diff(S_tmp.Y);
  d_arc_distance = sqrt(d_arc_x.^2 + d_arc_y.^2);
  arc_distance = [0 cumsum(d_arc_distance)];
  d_arcdist = 0.03;
  new_arc_distance = arc_distance(1):d_arcdist:arc_distance(end); % try to get atleast one box at start
  x_new = interp1(arc_distance,S_tmp.X,new_arc_distance);
  z_new = interp1(arc_distance,S_tmp.Y,new_arc_distance);
  % Collect data into table array
  S(it).twpe = pic_tmp.twpe;
  S(it).twci = pic_tmp.twci;
  S(it).A = S_tmp.Level;
  S(it).np_orig = S_tmp.Length;
  S(it).x_orig = S_tmp.X;
  S(it).z_orig = S_tmp.Y;  
  S(it).s_orig = arc_distance;  
  S(it).np = numel(x_new);
  S(it).x = x_new;
  S(it).z = z_new;
  S(it).s = new_arc_distance;
  for ivar = 1:numel(varstrs)
    xlim_tmp = [min(S(it).x) max(S(it).x)] + [-1 1];
    zlim_tmp = [min(S(it).z) max(S(it).z)] + [-1 1];
    pic_tmptmp = no02m.twpelim(S(it).twpe).xlim(xlim_tmp).zlim(zlim_tmp);
    %varstr = varstrs{ivar};
    %var = pic_tmptmp.(varstr);
    var_line = pic_tmptmp.interp(S(it).x,S(it).z,S(it).twci,varstrs{ivar});
    S(it).(varnames{ivar}) = var_line';
  end
  
end

xmin = min([S.x]);
%
for iS = 1:numel(S) % add some stuffs
  % center arclength where s=0 is at z=0
  S(iS).i0 = find(abs(S(iS).z)==min(abs(S(iS).z)));
  S(iS).s0 = S(iS).s(S(iS).i0);
  S(iS).s_centered = S(iS).s-S(iS).s0;
  % interpolate to common s where s=0 is at z=0
  s_common = S(1).s_centered; % first time has longest arclength
  S(iS).s_common = s_common;
  for ivar = 1:numel(varnames) 
    tmp_var = S(iS).(varnames{ivar});
    new_var = interp1(S(iS).s_centered,tmp_var,S(iS).s_common);
    S(iS).([varnames{ivar} '_common']) = new_var;
  end
end
%% 
times = [S.twci];
ss = S(1).s_common;
stE = cat(1,S.Epar_common);
stNitop = cat(1,S.nic_top_common);
stNibot = cat(1,S.nic_bot_common);

nrows = 4;
ncols = 1;
npaels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % Epar
  hca = h(isub); isub = isub + 1;
  pcolor(hca,ss,times,smooth2(stE,2))
  shading(hca,'flat')
  hca.XLabel.String = 's';
  hca.YLabel.String = 'twpe';
  hcb = colorbar('peer',hca);
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = [-0.5 0.5];
end
if 1 % ntop
  hca = h(isub); isub = isub + 1;
  pcolor(hca,ss,times,stNitop)
  shading(hca,'flat')
  hca.XLabel.String = 's';
  hca.YLabel.String = 'twpe';
  hcb = colorbar('peer',hca);
  colormap(hca,pic_colors('pasteljet'))
  hca.CLim = [0 0.5];
end
if 1 % nbot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,ss,times,stNibot)
  shading(hca,'flat')
  hca.XLabel.String = 's';
  hca.YLabel.String = 'twpe';
  hcb = colorbar('peer',hca);
  colormap(hca,pic_colors('pasteljet'))
  hca.CLim = [0 0.5];
end
if 1 % ntop/ntot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,ss,times,stNitop./(stNibot+stNitop))
  shading(hca,'flat')
  hca.XLabel.String = 's';
  hca.YLabel.String = 'twpe';
  hcb = colorbar('peer',hca);
  colormap(hca,pic_colors('pasteljet'))
  hca.CLim = [0 1];
end
compact_panels(0.01,0.0)
hlinks = linkprop(h,{'XLim','YLim'});
% for iS = 1:numel(S)
%   plot(S(iS).s_centered,S(iS).Epar)
%   if iS == 1
%     hold on
%   end
% end
% hold off

% % reduced distributions
% if 0 % saved reduced distributions
%   % save('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/matlab/fred_twpe23000_A75.mat','fred35_A75','fred46_A75','fred3_A75')
%   load('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/matlab//fred_twpe23000_A75.mat')
% elseif 0 % make reduced distributions
%   fred35_A75 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
%   fred3_A75 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
%   fred46_A75 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});  
% end
% 
% % Scalar quantities
% 
% arclength = [0; cumsum(sqrt(diff(fredi.x).^2 + diff(fredi.z).^2))];
% %arclength = arclength(end:-1:1);
% if 1; arclength = arclength - arclength(find(abs(fredi.z)==min(abs(fredi.z)))); end
% darc = arclength(2)-arclength(1);
% arcedges = [arclength(1)-0.5*darc; arclength+0.5*darc];
% narc = 1800; % 900
% arclength_interp = linspace(arclength(1),arclength(end),narc);
% darc_interp = (arclength_interp(end)-arclength_interp(1))/narc;
% xinterp = interp1(arclength,xdist,arclength_interp);
% zinterp = interp1(arclength,zdist,arclength_interp);
% 
% ni = interpfield(pic.xi,pic.zi,pic.ni,xinterp,zinterp); 
% ne = interpfield(pic.xi,pic.zi,pic.ne,xinterp,zinterp); 
% vepar = interpfield(pic.xi,pic.zi,pic.vpar([3 5]),xinterp,zinterp); 
% Epar = interpfield(pic.xi,pic.zi,pic.Epar,xinterp,zinterp);
% Babs = interpfield(pic.xi,pic.zi,pic.Babs,xinterp,zinterp); 
% intEpar = [-cumsum(Epar*darc_interp)];
%% Figure 9.2 Plot
colors = pic_colors('matlab');

nrows = 5;
ncols = 1;
%h = setup_subplots(nrows,ncols,'vertical');
h = setup_subplots(nrows,ncols,'horizontal');
%[h,h2] = initialize_combined_plot(nrows,2,2,0.4,'vertical')
isub = 1;
doE = 0; colorE = [0.0 0.0 0.0];
doV = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 0; colorExB = 0*[1 1 1]+0.5;
isMap = [];


if 1 % log 10 fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred35_A75_1;  
  surf(hca,arcedges1,fred.vpar_edges,zeros(numel(arcedges1),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  
  if 0 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.vpar_center,log10(fred.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{ic}'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = [-6 -1];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end
if 1 % fi3/fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fredall = fred35_A75_1;
  fred = fred3_A75_1;
  fplot = (fred.fvpar./fredall.fvpar);
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges1,fred.vpar_edges,zeros(numel(arcedges1),numel(fred.vpar_edges))',fplot')
  view(hca,[0 0 1]);   
  shading(hca,'flat')
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength1,fred.vpar_center,log10(fredall.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 's_{||} (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('pasteljet'))
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['f_{ic}^{top}/f_{ic}^{tot}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.CLim = [0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end

if 1 % Epar
  hca = h(isub); isub = isub + 1;
  pcolor(hca,ss,times,smooth2(stE,2))
  shading(hca,'flat')
  hca.XLabel.String = 's';
  hca.YLabel.String = 'twpe';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_{||}';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = [-0.5 0.5];  
  hca.XGrid = 'on';
  %hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YDir = 'reverse';
end

if 1 % log 10 fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred35_A75_2;  
  surf(hca,arcedges2,fred.vpar_edges,zeros(numel(arcedges2),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  
  if 0 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength2,fred.vpar_center,log10(fred.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{ic}'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = [-6 -1];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end
if 1 % fi3/fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fredall = fred35_A75_2;
  fred = fred3_A75_2;
  fplot = (fred.fvpar./fredall.fvpar);
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges2,fred.vpar_edges,zeros(numel(arcedges2),numel(fred.vpar_edges))',fplot')
  view(hca,[0 0 1]);   
  shading(hca,'flat')
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength2,fred.vpar_center,log10(fredall.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 's_{||} (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('pasteljet'))
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['f_{ic}^{top}/f_{ic}^{tot}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.CLim = [0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end
compact_panels(0.01)

hlinksx = linkprop(h,{'XLim'});
hlinksy = linkprop(h([1 2 4 5]),{'YLim'});

%% Figure 9.0 top, map of Epar att different times
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)'};

twpe = 23000:200:24000;
xlim = [67 76];
zlim = [-6 -1];
xlim = [67 76]+[-3 3];
zlim = [-6 -1]+[-2 2];
h = setup_subplots(2,3);

for it = 1:numel(twpe)
  hca = h(it);
  pic_tmp = no02m.twpelim(twpe(it)).xlim(xlim).zlim(zlim);
  %imagesc(hca,pic_tmp.xi,pic_tmp.zi,smooth2(pic_tmp.Epar,2)')
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,smooth2(pic_tmp.Ez,2)')
  %imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.n([3 5])')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';  
  irf_legend(hca,['t\omega_{ci}=' sprintf('%g',pic_tmp.twci)],[0.02 0.98],'color',[0 0 0],'fontsize',13)%,'fontweight','bold')
  colormap(hca,pic_colors('blue_red'));
  %colormap(hca,pic_colors('thermal'));
  hca.CLim = 0.99*[-0.5 0.5];
  %hca.CLim = [0 0.5];
  hca.YDir = 'normal';
  hold(hca,'on')
  contour(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.A',[0:1:25])
  contour(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.A',[7.5 7.5],'linewidth',1)
  %axis(hca,'equal')
  hold(hca,'off')
end
%c_eval('h(?).Position(1) = h(?).Position(1)-0.01;',1:6)
c_eval('h(?).Position(2) = h(4).Position(2)+h(4).Position(4);',[1 2 3])
c_eval('h(?).Position(2) = h(?).Position(2)+0.1;',1:6)
c_eval('h(?).Position(1) = h(?).Position(1)-0.06;',1:6)

compact_panels(h,0,0);
drawnow
c_eval('h(?).FontSize = 13;',1:6)

hold(h(1),'on')
[~,hbox] = ds_pick.plot_boxes(h(1));
hbox.LineWidth = 2;
%plot(h(1),xval,zval,'sk')
hold(h(1),'off')
drawnow
hb = colorbar('peer',h(end));
hb.YLabel.String = 'E_{||}';
hb.Position(1) = h(end).Position(1)+h(end).Position(3)+0.005;
hb.Position(4) = h(end).Position(4)*2;
c_eval('h(?).YLabel.String = [];',[2 3 5 6])
c_eval('h(?).YTickLabels = [];',[2 3 5 6])

drawnow
hb.Position(1) = h(end).Position(1)+h(end).Position(3)+0.005;
if 1 % annotate
  %%
  hann_box = annotation('textarrow',[0.20 0.24]-0.035,[0.75 0.71]-0.05,...
    'string',{'input for','instability','analysis'},...
    'fontsize',12,'TextBackgroundColor',[1 1 1],...
    'TextEdgeColor',[0 0 0]);
  
  hann(1) = annotation('arrow',[0.22 0.24]+0.00,[0.79 0.71]+0.01);
  hann(2) = annotation('arrow',[0.22 0.24]+0.239,[0.79 0.71]+0.04);
  hann(3) = annotation('arrow',[0.22 0.24]+0.479,[0.79 0.71]+0.06);
  hann(4) = annotation('arrow',[0.22 0.24]-0.06,[0.79 0.71]-0.26);
  hann(5) = annotation('arrow',[0.22 0.24]+0.185,[0.79 0.71]-0.26);
  hann(6) = annotation('arrow',[0.22 0.24]+0.43,[0.79 0.71]-0.26);  
end  

%% Figure 10, forces on fermi ions
if 0 % prepare data
twpe = 23000; xlim = [64 80]; zlim = [-8 8]; 
pic_lim = no02m.xlim(xlim).zlim(zlim).twpelim(twpe);

Bx_ = pic_lim.Bx;
By_ = pic_lim.By;
Bz_ = pic_lim.Bz;
Bx_interp = interpfield(pic_lim.xi,pic_lim.zi,Bx_,xinterp,zinterp); 
By_interp = interpfield(pic_lim.xi,pic_lim.zi,By_,xinterp,zinterp); 
Bz_interp = interpfield(pic_lim.xi,pic_lim.zi,Bz_,xinterp,zinterp);
Ex_ = pic_lim.Ex;
Ey_ = pic_lim.Ey;
Ez_ = pic_lim.Ez;
Ex = interpfield(pic_lim.xi,pic_lim.zi,Ex_,xdist,zdist); 
Ey = interpfield(pic_lim.xi,pic_lim.zi,Ey_,xdist,zdist); 
Ez = interpfield(pic_lim.xi,pic_lim.zi,Ez_,xdist,zdist);
end 
colors = pic_colors('matlab');
dv = diff(fred.v(1:2));
v_edges = [fred.v(1)-0.5*dv fred.v+0.5*dv];

doB = 1;
Bcolors = [0 0 0; colors([1 4],:)];
fcont = [-1:1:2];
fredcont = fred35_A75;


nrows = 3;
ncols = 4;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

if 0 % log 10 fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred35_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  
  if 0 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.vpar_center,log10(fred.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{i,cold}(s_{||},v_{||})'];
  hcb.YLabel.String = ['log_{10}f_{ic}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = [-6 -1];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end

if 1 % log 10 fi35(vx)
  hca = h(isub); isub = isub + 1;
  fred = fred35_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,v_edges,zeros(numel(arcedges),numel(v_edges))',log10(fred.fvx)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{x}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['log_{10}f_{ic}'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.v,log10(fred.fvx)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end  
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doB
    hold(hca,'on')
    hBy = plot(hca,arclength_interp,By_interp,'color',Bcolors(2,:));
    hBz = plot(hca,arclength_interp,Bz_interp,'color',Bcolors(3,:));
    hold(hca,'off')
    legend([hBy hBz],{'B_y','B_z'},'location','northeast','box','off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vy)
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,v_edges,zeros(numel(arcedges),numel(v_edges))',log10(fred.fvy)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  
  if 0 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.vpar_center,log10(fred.fvy)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{i,cold}(s_{||},v_{||})'];
  hcb.YLabel.String = ['log_{10}f_{ic}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  %hca.CLim = [-6 -1];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.v,log10(fred.fvy)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doB
    hold(hca,'on')
    hBx = plot(hca,arclength_interp,Bx_interp,'color',Bcolors(1,:));
    hBz = plot(hca,arclength_interp,Bz_interp,'color',Bcolors(3,:));
    hold(hca,'off')
    legend([hBx hBz],{'B_x','B_z'},'location','southeast','box','off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vz)
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,v_edges,zeros(numel(arcedges),numel(v_edges))',log10(fred.fvz)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  
  if 0 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.vpar_center,log10(fred.fvz)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{z}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{i,cold}(s_{||},v_{||})'];
  hcb.YLabel.String = ['log_{10}f_{ic}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  %hca.CLim = [-6 -1];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.v,log10(fred.fvz)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
    
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doB
    hold(hca,'on')
    plot(hca,arclength_interp,Bx_interp,'color',Bcolors(1,:))
    plot(hca,arclength_interp,By_interp,'color',Bcolors(2,:))
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end

if 1 % log 10 fi35(vx) vx*By
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  [S,V] = meshgrid(arclength,fred.v);
  [BY,V] = meshgrid(By,fred.v);  
  F = V.*BY;
  F(fred.fvx'<1e-2) = NaN;
  pcolor(hca,S,V,F)
  %view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{x}';
  colormap(hca,pic_colors('blue_red'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['F_z=v_xB_y'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = 0.5*[-1 1];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fredcont.v,log10(fredcont.fvx)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doB
    hold(hca,'on')
    plot(hca,arclength_interp,By_interp,'color',Bcolors(2,:))
    plot(hca,arclength_interp,Bz_interp,'color',Bcolors(3,:))
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vx) -vy*Bx
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  [S,V] = meshgrid(arclength,fred.v);
  [BX,V] = meshgrid(Bx,fred.v);  
  F = -V.*BX;
  F(fred.fvy'<1e-2) = NaN;
  pcolor(hca,S,V,F)
  %view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('blue_red'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['F_z=-v_yB_x'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = 0.5*[-1 1];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fredcont.v,log10(fredcont.fvy)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doB
    hold(hca,'on')
    plot(hca,arclength_interp,By_interp,'color',Bcolors(2,:))
    plot(hca,arclength_interp,Bz_interp,'color',Bcolors(3,:))
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vx) vz*Bx
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  [S,V] = meshgrid(arclength,fred.v);
  [BX,V] = meshgrid(Bx,fred.v);  
  F = V.*BX;
  F(fred.fvz'<1e-2) = NaN;
  pcolor(hca,S,V,F)
  %view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{z}';
  colormap(hca,pic_colors('blue_red'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['F_y=v_zB_x'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = 0.5*[-1 1];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fredcont.v,log10(fredcont.fvz)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doB
    hold(hca,'on')
    plot(hca,arclength_interp,By_interp,'color',Bcolors(2,:))
    plot(hca,arclength_interp,Bz_interp,'color',Bcolors(3,:))
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vx) -vx*Bz
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  [S,V] = meshgrid(arclength,fred.v);
  [BZ,V] = meshgrid(Bz,fred.v);  
  F = -V.*BZ;
  F(fred.fvx'<1e-2) = NaN;
  %F(smooth2(fred.fvx,2)<1e-2) = NaN;
  pcolor(hca,S,V,F)
  %view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{x}';
  colormap(hca,pic_colors('blue_red'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['F_y=-v_xB_z'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = 0.5*[-1 1];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fredcont.v,log10(fredcont.fvx)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doB
    hold(hca,'on')
    plot(hca,arclength_interp,By_interp,'color',Bcolors(2,:))
    plot(hca,arclength_interp,Bz_interp,'color',Bcolors(3,:))
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vx) vy*Bz
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  [S,V] = meshgrid(arclength,fred.v);
  [BZ,V] = meshgrid(Bz,fred.v);  
  F = V.*BZ;
  F(fred.fvy'<1e-2) = NaN;
  %F(smooth2(fred.fvx,2)<1e-2) = NaN;
  pcolor(hca,S,V,F)
  %view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('blue_red'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['F_x=v_yB_z'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = 0.5*[-1 1];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fredcont.v,log10(fredcont.fvy)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doB
    hold(hca,'on')
    plot(hca,arclength_interp,By_interp,'color',Bcolors(2,:))
    plot(hca,arclength_interp,Bz_interp,'color',Bcolors(3,:))
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vx) -vz*By
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  [S,V] = meshgrid(arclength,fred.v);
  [BY,V] = meshgrid(By,fred.v);  
  F = -V.*BY;
  F(fred.fvz'<1e-2) = NaN;
  %F(smooth2(fred.fvx,2)<1e-2) = NaN;
  pcolor(hca,S,V,F)
  %view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{z}';
  colormap(hca,pic_colors('blue_red'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['F_x=-v_zB_y'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = 0.5*[-1 1];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fredcont.v,log10(fredcont.fvz)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doB
    hold(hca,'on')
    plot(hca,arclength_interp,By_interp,'color',Bcolors(2,:))
    plot(hca,arclength_interp,Bz_interp,'color',Bcolors(3,:))
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vx) Ex
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;  
  [S,V] = meshgrid(arclength,fred.v);
  [EX,V] = meshgrid(Ex,fred.v);    
  F = EX;
  F(fred.fvz'<1e-2) = NaN;
  pcolor(hca,S,V,smooth2(F,2))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{x}';
  colormap(hca,pic_colors('blue_red'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['F_x=E_x'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = 0.5*[-1 1];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fredcont.v,log10(fredcont.fvz)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vx) Ey
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;  
  [S,V] = meshgrid(arclength,fred.v);
  [EY,V] = meshgrid(Ey,fred.v);    
  F = EY;
  F(fred.fvz'<1e-2) = NaN;
  pcolor(hca,S,V,smooth2(F,2))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('blue_red'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['F_y=E_y'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = 0.5*[-1 1];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fredcont.v,log10(fredcont.fvz)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vx) Ex
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;  
  [S,V] = meshgrid(arclength,fred.v);
  [EZ,V] = meshgrid(Ez,fred.v);    
  F = EZ;
  F(fred.fvz'<1e-2) = NaN;
  pcolor(hca,S,V,smooth2(F,2))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{z}';
  colormap(hca,pic_colors('blue_red'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['F_z=E_z'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = 0.5*[-1 1];
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fredcont.v,log10(fredcont.fvz)',fcont,'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end

%hlinks_dist = linkprop(h(1:3),{'XLim','YLim','CLim'});
hlinks_forces = linkprop(h(4:end),{'XLim','YLim','CLim'});

hlinks_forces.Targets(1).CLim = [-0.4 0.4];

%% Figure, last, illustration/sketch of thermalization process
nv = 1000;
v = linspace(-3,4,nv);
nf = 5;
n = [1 1 1 1 1];
vdstep = 0.5;
vd = 0.0:vdstep:nf*vdstep;
vt = [0.2 0.3 0.4 0.5 0.6];

fcold = fmax1D(v,n(1),0,vt(1));
fhot = fmax1D(v,n(1),0,vt(3));
clear f
for ii = 1:nf
  f(ii,:) = fmax1D(v,n(ii),vd(ii),vt(ii));
end

% Plot
h(1) = subplot(1,3,1);
h(2) = subplot(1,3,[2 3]);

% Illustration of single population thermalization
hca = h(1);
plot(hca,v,fcold,v,fhot,'linewidth',1.5)
hca.Box = 'off';
hca.Visible = 'off';
hca.XLim = 1.1*[-1 1];
hca.YLim = [0 3];
hca.Position(3) = 0.25;

hca = h(2);
plot(hca,v,f,'linewidth',1.5)
hold(hca,'on')
plot(hca,v,sum(f(1:end,:),1),'linestyle','-','linewidth',1.5,'color',[0 0 0]+0.5)
hold(hca,'off')
hca.Box = 'off';
hca.XLim = [-1 3.5];
hca.YLim = [0 3];
hca.Visible = 'off';

