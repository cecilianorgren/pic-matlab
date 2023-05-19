%sus = PIC('/Volumes/DataRaid/susanne-stuff/fields.h5');

%% Reconnection electric field for varying By

twci = sus.twci;
RE = sus.RE;

xz_xline = sus.xline_position;

for it = 1:size(xz_xline,1)
  xl = xz_xline(it,1) + 0.2*[-1 1];
  zl = xz_xline(it,2) + 0.2*[-1 1];

  bytmp = sus(it).xlim(xl).zlim(zl).By;
  BY_xline(it) = mean(bytmp(:));

  dz = 0.5;
  bxtmp_top = sus(it).xlim(xl).zlim(zl+dz).Bx; bxtmp_top = mean(bxtmp_top(:));
  bxtmp_bot = sus(it).xlim(xl).zlim(zl-dz).Bx; bxtmp_bot = mean(bxtmp_bot(:));
  bytmp_top = sus(it).xlim(xl).zlim(zl+dz).By; bytmp_top = mean(bytmp_top(:));
  bytmp_bot = sus(it).xlim(xl).zlim(zl-dz).By; bytmp_bot = mean(bytmp_bot(:));
  btmp_top = sqrt(bxtmp_top^2+bytmp_top^2);
  bytmp_top = bytmp_top/btmp_top;
  bxtmp_top = bxtmp_top/btmp_top;
  btmp_bot = sqrt(bxtmp_bot^2+bytmp_bot^2);
  bytmp_bot = bytmp_bot/btmp_bot;
  bxtmp_bot = bxtmp_bot/btmp_bot;

  B_shear(it) = acosd(bytmp_top*bytmp_bot+bxtmp_top*bxtmp_bot);
end

% Plot
fontsize = 14;

ncols = 1;
nrows = 3;
ipanel = 0;
h = gobjects([nrows,ncols]);
for icol = 1:ncols, for irow = 1:nrows, ipanel = ipanel + 1; h(irow,icol) = subplot(nrows,ncols,ipanel); end, end
isub = 1;

if 1 % initial conditions
  hca = h(isub); isub = isub + 1;
  By = sus(1).xlim([1 2]).By;
  By = mean(By,1);
  Bx = sus(1).xlim([1 2]).Bx;
  Bx = mean(Bx,1);
  zi = sus.zi;

  plot(hca,zi,Bx,sus.zi,By)
  hca.XLim = [-30 30];
  hca.XLabel.String = 'z/d_i';
  hca.YLim = 1.1*[-1 1];
  hca.YLabel.String = 'B';

  legend(hca,{'B_x','B_y'},'Box','off','location','northwest')
  hca.Title.String = 'Simulation initial conditions';
end
if 1
  its = 5:numel(twci);
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,twci(its),RE(its),twci,BY_xline);
  AX(1).YLim = [0 0.25]*0.99;
  AX(2).YLim = [0 0.5]*0.99;
  AX(1).YTick = 0:0.05:1;
  AX(2).YTick = 0:0.1:1;
  AX(1).YLabel.String = 'Reconnection rate';
  AX(2).YLabel.String = 'B_y at X line';
  
  c_eval('AX(?).XLim = [0 twci(end)];',1:numel(AX))
  c_eval('AX(?).FontSize = 14;',1:numel(AX))
  %hca.Position(2) = 0.2;
  %hca.Position(4) = 0.7;
  
  H1.LineWidth = 1;
  H2.LineWidth = 1;
  hca.LineWidth = 1;
  
  AX(1).XLabel.String = 't\omega_{ce}';
  hca.Title.String = 'Magnetic shear effect on reconnection rate';
end
if 0 % rec rate and magnetic shear
  its = 5:numel(twci);
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,twci(its),RE(its),twci,B_shear);
  AX(1).YLim = [0 0.3];
  AX(2).YLim = [0 180];
  AX(1).YTick = 0:0.05:1;
  AX(2).YTick = 0:30:180;
  AX(1).YLabel.String = 'Reconnection rate';
  AX(2).YLabel.String = 'Magnetic shear across current sheet';
  %irf_legend(hca,'cos^{-1}(B\cdot B)',[0.02 0.9])
  
  c_eval('AX(?).XLim = [0 twci(end)];',1:numel(AX))
  c_eval('AX(?).FontSize = 14;',1:numel(AX))
  %hca.Position(2) = 0.2;
  %hca.Position(4) = 0.7;
  
  H1.LineWidth = 1;
  H2.LineWidth = 1;
  hca.LineWidth = 1;
  
  AX(1).XLabel.String = 't\omega_{ce}';
  hca.Title.String = 'Magnetic shear effect on reconnection rate';
end
if 0 % outflow at later times
  hca = h(isub); isub = isub + 1;
  twci = 100;
  By = sus.twcilim(twci).zlim(3*[-1 1]).By;
  By = mean(By,2);
  Bz = sus.twcilim(twci).zlim(3*[-1 1]).Bz;
  Bz = mean(Bz,2);
  ni = sus.twcilim(twci).zlim(3*[-1 1]).ni;
  ni = mean(ni,2);
  
  xi = sus.xi;

  plot(hca,xi,ni,xi,By,xi,Bz)
  hca.XLim = [20 180];
  hca.XLabel.String = 'x/d_i';
  %hca.YLim = 1.1*[-1 1];
  hca.YLabel.String = 'B';

  legend(hca,{'n','B_y','B_z'},'Box','off','location','northwest')
  hca.Title.String = 'Simulation initial conditions';
end

if 1 % timex map outflow at later times
  hca = h(isub); isub = isub + 1;
  
  [hax,hb] = sus.zlim([-1 1]).plot_timemap(hca,'tx',{'By'},'A',1);
  hold(hca,'on')
  plot(hca,twci,xz_xline(:,1),'k--','linewidth',1)
  hold(hca,'off')
  hb.YLabel.String = 'B_y at z = 0';
  hb.FontSize = fontsize;
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = [-1 1];
  hca.YLim = [20 180];
  hca.YLabel.String = 'x/d_i';
  %hca.YLim = 1.1*[-1 1];
  hca.XLabel.String = 't\omega_{ci}';


  %legend(hca,{'n','B_y','B_z'},'Box','off','location','northwest')
  hca.Title.String = 'Outflow structure';
  hca.Position(3) = h(1).Position(3);
end

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).Position(3) = 0.75;',1:numel(h))
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))
hca.LineWidth = 2;






%% Reconnection electric field for varying By

pic = no02m;
twci = pic.twci;
RE = pic.RE;

xz_xline = pic.xline_position;

for it = 1:size(xz_xline,1)
  xl = xz_xline(it,1) + 0.2*[-1 1];
  zl = xz_xline(it,2) + 0.2*[-1 1];


  ntmp_h = pic(it).xlim(xl).zlim(zl).n(1);
  ntmp_c = pic(it).xlim(xl).zlim(zl).n([3 5]);
  nc_xline(it) = mean(ntmp_c(:));
  nh_xline(it) = mean(ntmp_h(:));
end

% Plot
fontsize = 14;
colors = pic_colors('matlab');

ncols = 1;
nrows = 3;
ipanel = 0;
h = gobjects([nrows,ncols]);
for icol = 1:ncols, for irow = 1:nrows, ipanel = ipanel + 1; h(irow,icol) = subplot(nrows,ncols,ipanel); end, end
isub = 1;

if 1 % initial conditions
  hca = h(isub); isub = isub + 1;
  ni_hot = pic(1).xlim([1 6]).n(1);
  ni_hot = mean(ni_hot,1);
  ni_cold = pic(1).xlim([1 6]).n([3 5]);
  ni_cold = mean(ni_cold,1);
  Bx = pic(1).xlim([1 4]).Bx;
  Bx = mean(Bx,1);
  zi = pic.zi;

  plot(hca,zi,smooth(ni_cold,10),zi,smooth(ni_hot,10))
  hca.XLim = [-12.5 12.5];
  hca.XLabel.String = 'z/d_i';
  %hca.YLim = 1.1*[-1 1];
  hca.YLabel.String = 'n';

  legend(hca,{'n_{cold}','n_{hot}'},'Box','off','location','northwest')
  hca.Title.String = 'Simulation initial conditions';
end
if 1
  its = 1:numel(twci);
  hca = h(isub); isub = isub + 1;
  
  [AX,H1,H2] = plotyy(hca,twci(its),RE(its),twci,[nc_xline;nh_xline;nc_xline+nh_xline]);
  H1.Color = [0 0 0];
  AX(1).YColor = [0 0 0];
  H2(1).Color = colors(1,:);
  H2(2).Color = colors(2,:);
  H2(3).Color = [0.5 0.5 0.5];
  AX(1).YLabel.String = 'Reconnection rate';
  AX(2).YLabel.String = 'Density at X line';
  
 
  AX(1).XLim = [0 twci(end)];
  AX(2).XLim = [0 twci(end)];
  hca.FontSize = fontsize;
  %hca.Position(2) = 0.2;
  %hca.Position(4) = 0.7;
  
  hca.LineWidth = 1;
  
  hca.XLabel.String = 't\omega_{ce}';
  %hca.Title.String = 'Reconnection r';
  legend(hca,{'R_E','n_{cold}','n_{hot}','n_{tot}'},'Box','off','location','northwest')
end

if 1 % initial conditions
  hca = h(isub); isub = isub + 1;
  twpe = 23000;
  ni_hot = pic.twpelim(twpe).zlim([-1 1]).n(1);
  ni_hot = mean(ni_hot,2);
  ni_cold = pic.twpelim(twpe).zlim([-1 1]).n([3 5]);
  ni_cold = mean(ni_cold,2);
  Bz = pic.twpelim(twpe).zlim([-1 1]).Bz;
  Bz = mean(Bz,2);
  xi = pic.xi;

  hl = plot(hca,xi,smooth(ni_cold,1),xi,smooth(ni_hot,1),xi,Bz);
  hca.XLim = [50 150];
  hca.XLabel.String = 'x/d_i';
  %hca.YLim = 1.1*[-1 1];
  hca.YLabel.String = 'n, B';

  legend(hca,{'n_{cold}','n_{hot}','B_z'},'Box','off','location','northwest')
  hca.Title.String = ['Outflow structure at t\omega_{pe}' sprintf(' = %g',pic.twpelim(twpe).twci)];
end
if 0 % rec rate and magnetic shear
  its = 5:numel(twci);
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,twci(its),RE(its),twci,B_shear);
  AX(1).YLim = [0 0.3];
  AX(2).YLim = [0 180];
  AX(1).YTick = 0:0.05:1;
  AX(2).YTick = 0:30:180;
  AX(1).YLabel.String = 'Reconnection rate';
  AX(2).YLabel.String = 'Magnetic shear across current sheet';
  %irf_legend(hca,'cos^{-1}(B\cdot B)',[0.02 0.9])
  
  c_eval('AX(?).XLim = [0 twci(end)];',1:numel(AX))
  c_eval('AX(?).FontSize = 14;',1:numel(AX))
  %hca.Position(2) = 0.2;
  %hca.Position(4) = 0.7;
  
  H1.LineWidth = 1;
  H2.LineWidth = 1;
  hca.LineWidth = 1;
  
  AX(1).XLabel.String = 't\omega_{ce}';
  hca.Title.String = 'Magnetic shear effect on reconnection rate';
end
if 0 % outflow at later times
  hca = h(isub); isub = isub + 1;
  twci = 100;
  By = sus.twcilim(twci).zlim(3*[-1 1]).By;
  By = mean(By,2);
  Bz = sus.twcilim(twci).zlim(3*[-1 1]).Bz;
  Bz = mean(Bz,2);
  ni = sus.twcilim(twci).zlim(3*[-1 1]).ni;
  ni = mean(ni,2);
  
  xi = sus.xi;

  plot(hca,xi,ni,xi,By,xi,Bz)
  hca.XLim = [20 180];
  hca.XLabel.String = 'x/d_i';
  %hca.YLim = 1.1*[-1 1];
  hca.YLabel.String = 'B';

  legend(hca,{'n','B_y','B_z'},'Box','off','location','northwest')
  hca.Title.String = 'Simulation initial conditions';
end

if 0 % timex map outflow at later times
  hca = h(isub); isub = isub + 1;
  
  [hax,hb] = sus.zlim([-1 1]).plot_timemap(hca,'tx',{'By'},'A',1);
  hold(hca,'on')
  plot(hca,twci,xz_xline(:,1),'k--','linewidth',1)
  hold(hca,'off')
  hb.YLabel.String = 'B_y at z = 0';
  hb.FontSize = fontsize;
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = [-1 1];
  hca.YLim = [20 180];
  hca.YLabel.String = 'x/d_i';
  %hca.YLim = 1.1*[-1 1];
  hca.XLabel.String = 't\omega_{ci}';


  %legend(hca,{'n','B_y','B_z'},'Box','off','location','northwest')
  hca.Title.String = 'Outflow structure';
  hca.Position(3) = h(1).Position(3);
end

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).Position(3) = 0.75;',1:numel(h))
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))
%hca.LineWidth = 2;




