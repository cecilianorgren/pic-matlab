%sus = PIC('/Volumes/DataRaid/susanne-stuff/fields.h5');
tint_cavity = irf.tint('2017-07-25T22:09:46.00Z/2017-07-25T22:09:56.00Z');
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

%% Illustration of difficulty of estimating current sheet thickness
colors = pic_colors('matlab');

B = @(z,z0,L,B0) B0*tanh((z-z0)/L);
J = @(z,z0,L,B0) (B0/L)*(1-(B(z,z0,L,B0)/B0).^2);

z = linspace(-3.5,2.5,100);

hca = subplot(2,1,1);


plot(hca,z,B(z,0,1,1),'k',...
         z,B(z,0.5,1,1),'--',...
         z,B(z,0,2,1.0),'-.')
irf_legend(hca,{'Magnetic field'},[0.1 0.98],'color','k','fontsize',16)
set(hca,'ColorOrder',[0 0 0; colors])
irf_legend(hca,{'Original','Flapping','Thickening'}',[0.1 0.8],'fontsize',16)


hca = subplot(2,1,2);
plot(hca,z,J(z,0,1,1),'k',...
         z,J(z,0.5,1,1),'--',...
         z,J(z,0,2,1.0),'-.')
irf_legend(hca,{'Current density'},[0.1 0.98],'color','k','fontsize',16)       
set(hca,'ColorOrder',[0 0 0; colors])
irf_legend(hca,{'Original','Flapping','Thickening'}',[0.1 0.8],'fontsize',16)

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 2;',1:numel(hl))

h = findobj(gcf,'type','axes');
c_eval('h(?).LineWidth = 2;',1:numel(h))
c_eval('h(?).Visible = ''off'';',1:numel(h))
c_eval('h(?).XLim = z([1 end]);',1:numel(h))






%% Oppoaite Hall fields, 1
ic = 1;

time_df1 = irf_time('2017-07-25T22:10:06.00Z','utc>EpochTT');
time_df2 = irf_time('2017-07-25T22:10:28.00Z','utc>EpochTT');

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
feeps_ion_omni1234 = feeps_ion_omni1;
%feeps_ion_omni1234.data = (feeps_ion_omni1.data + feeps_ion_omni2.resample(feeps_ion_omni1).data + feeps_ion_omni3.resample(feeps_ion_omni1).data + feeps_ion_omni4.resample(feeps_ion_omni1).data)/4;
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

fontsize = 14;
fontsize_leg = 13;

npanels = 4;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

isub = 0;
zoomy = [];
cs = 'gse';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'B_x','B_y','B_z'}',[1.02 0.9],'fontsize',fontsize_leg);
end 
if 0 % E gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gse lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gseE?.filt(0,1,[],3).x,gseE?.filt(0,1,[],3).y,gseE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 1 % E gse downsampled  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gse downsamp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gseE?.resample(gseVi1).x,gseE?.resample(gseVi1).y,gseE?.resample(gseVi1).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
    irf_legend(hca,{'E_x','E_y','E_z'}',[1.02 0.9],'fontsize',fontsize_leg);
  %irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i{\perp}x}','v_{i{\perp}y}','v_{i{\perp}z}','v_{i||}'}',[1.02 0.9],'fontsize',fontsize_leg);
end
if 0 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?.x,%sVi?.y,%sVi?.z},''comp'');',cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'}',[1.02 0.9],'fontsize',fontsize_leg);
end
if 0 % Vi xyz, fpi, hpca
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi fpi hpca');
  colors_plot = [mms_colors('xyza'); mms_colors('x').^0.5];
  set(hca,'ColorOrder',colors_plot)
  c_eval(sprintf('irf_plot(hca,{%sVi?.x,%sVi?.y,%sVi?.z},''comp'');',cs,cs,cs),ic)  
  c_eval('irf_plot(hca,{gseVHp?_srvy.x,gseVHp?_srvy.y,gseVHp?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',colors_plot)
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{x}','v_{y}','v_{z}'}',[1.02 0.9],'fontsize',fontsize_leg);

  irf_legend(hca,{'FPI'},[.17 0.77],'fontsize',fontsize_leg,'color','k');
  irf_legend(hca,{'HPCA'},[.03 0.8],'fontsize',fontsize_leg,'color','k');
end

if 1 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',[mms_colors('z'); 0.1000    0.4000    1.0000; 0 0 0])
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    c_eval(sprintf('ne = ne?.resample(%sVi?);',cs),ic);    
    %ve.data(abs(ve.data)>5000) = NaN;
    ve.data(abs(ne.data)<0.04,:) = NaN;
    
    %c_eval(sprintf('irf_plot(hca,{ve.resample(%sVi?),%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp),gseVHp?_srvy.(comp)},''comp'');',cs,cs,cs,cs),ic)  
    c_eval(sprintf('irf_plot(hca,{ve.resample(%sVi?),%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp%s}',comp),'(km/s)'};
    set(hca,'ColorOrder',[mms_colors('z'); mms_colors('x'); 0 0 0])
    irf_legend(hca,{'v_e','v_i','ExB'}',[1.02 0.9],'fontsize',fontsize_leg);    
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end


if 0 % JotE
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('JdotE');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{JdotE?.resample(gseVi?)},''comp'');',ic)  
  hca.YLabel.String = {'J\cdot E','(nW/m^3)'}; % 
  set(hca,'ColorOrder',mms_colors('xyza'))
end


%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0],'fontsize',fontsize_leg)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize_leg;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);

c_eval('h(?).FontSize = fontsize;',1:numel(h))
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

c_eval('h(?).LineWidth = 1;',1:numel(h))


c_eval('h(?).XGrid = ''on'';',1:numel(h))
c_eval('h(?).YGrid = ''on'';',1:numel(h))
c_eval('h(?).GridLineStyle = '':'';',1:numel(h))
c_eval('h(?).GridAlpha = 0.3;',1:numel(h))

h(end).XTickLabelRotation = 0;

irf_zoom(h,'x',tint_action)
%irf_zoom(h,'y')
%c_eval('hmark_df!(?) = irf_pl_mark(h(?),time_df!.epochUnix,''k''); hmark_df!(?).LineStyle = ''--'';',1:numel(h),1:2)
hca = irf_panel('Vi fpi hpca');
irf_zoom(hca,'y')

%hca = irf_panel('feeps omni mms 1234'); hca.CLim = [-1 100];
%hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
hca = irf_panel('B gsm'); hca.YLim = 0.99*[-17 22];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [1e-2 1e3];
%hca = irf_panel('Tepar/Teperp'); hca.YLim = [0.5 1.5];

%% Oppoaite Hall fields, 2
tint_cavity = irf.tint('2017-07-25T22:09:46.00Z/2017-07-25T22:09:56.00Z');
tint_cavity_longer = tint_cavity + [-3 3];
ic = 1;

time_df1 = irf_time('2017-07-25T22:10:06.00Z','utc>EpochTT');
time_df2 = irf_time('2017-07-25T22:10:28.00Z','utc>EpochTT');

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
feeps_ion_omni1234 = feeps_ion_omni1;
%feeps_ion_omni1234.data = (feeps_ion_omni1.data + feeps_ion_omni2.resample(feeps_ion_omni1).data + feeps_ion_omni3.resample(feeps_ion_omni1).data + feeps_ion_omni4.resample(feeps_ion_omni1).data)/4;
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

fontsize = 14;
fontsize_leg = 13;

npanels = 5;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab')

plot_colors = [0.6000    0.6000    0.6000;...
               0         0         0;...
               1.0000    0.2000         0];
matlab_colors = pic_colors('matlab');
j_colors = [mms_colors('123'); matlab_colors(6,:)];
j_colors = [mms_colors('123'); matlab_colors(3,:)];
j_colors = [mms_colors('xyz'); 0 0 0];

plot_colors = [0.8 0.8 0.8; j_colors([4 3],:)];
leg_colors = [0.6 0.6 0.6; j_colors([4 3],:)];

isub = 0;
zoomy = [];
cs = 'gse';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'B_x','B_y','B_z'}',[1.02 0.9],'fontsize',fontsize_leg);
end 
if 0 % E gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gse lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gseE?.filt(0,1,[],3).x,gseE?.filt(0,1,[],3).y,gseE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 1 % E gse downsampled  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gse downsamp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gseE?.resample(gseVi1).x,gseE?.resample(gseVi1).y,gseE?.resample(gseVi1).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
    irf_legend(hca,{'E_x','E_y','E_z'}',[1.02 0.9],'fontsize',fontsize_leg);
  %irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % JotE
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('JdotE');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{JdotE?.resample(gseVi?)},''comp'');',ic)  
  hca.YLabel.String = {'J\cdot E','(nW/m^3)'}; % 
  set(hca,'ColorOrder',mms_colors('xyza'))
end
if 1 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',[plot_colors; 0 0 1])
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    c_eval(sprintf('ne = ne?.resample(%sVi?);',cs),ic);    
    %ve.data(abs(ve.data)>5000) = NaN;
    ve.data(abs(ne.data)<0.04,:) = NaN;
    
    %c_eval(sprintf('irf_plot(hca,{ve.resample(%sVi?),%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp),gseVHp?_srvy.(comp)},''comp'');',cs,cs,cs,cs),ic)  
    c_eval(sprintf('irf_plot(hca,{ve.resample(%sVi?),%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp%s}',comp),'(km/s)'};
    set(hca,'ColorOrder',[leg_colors; 0 0 1])
    irf_legend(hca,{'v_e','v_i','ExB'}',[1.02 0.9],'fontsize',fontsize_leg);    
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end

if 1 % JxBne gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('JxB gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseJxBne?_mVm.x,gseJxBne?_mVm.y,gseJxBne?_mVm.z},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJxBne?_mVm.x.resample(gseVi?),gseJxBne?_mVm.y.resample(gseVi?),gseJxBne?_mVm.z.resample(gseVi?)},''comp'');',ic)
  hca.YLabel.String = {'J\times B/ne','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'J\times{B_x}','J\times{B_y}','J\times{B_z}'}',[1.02 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % vepar
  isub = isub + 1;
  hca = irf_panel('vepar');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gseVe?par},''comp'')',ic)
  hold(hca,'on')
  c_eval('irf_patch(hca,{gseVe?par,0})',ic)
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('12'))  
  hca.YLabel.String = {'v_{e||}','(km/s)'};
  hca.YLabel.Interpreter = 'tex';
end

%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0],'fontsize',fontsize_leg)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize_leg;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
%irf_zoom(h(zoomy),'y')
drawnow
irf_zoom(h,'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);

c_eval('h(?).FontSize = fontsize;',1:numel(h))
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

c_eval('h(?).LineWidth = 1;',1:numel(h))

c_eval('h(?).XGrid = ''on'';',1:numel(h))
c_eval('h(?).YGrid = ''on'';',1:numel(h))
c_eval('h(?).GridLineStyle = '':'';',1:numel(h))
c_eval('h(?).GridAlpha = 0.3;',1:numel(h))




h(end).XTickLabelRotation = 0;

irf_zoom(h,'x',tint_cavity_longer)
irf_zoom(h,'y')
%c_eval('hmark_df!(?) = irf_pl_mark(h(?),time_df!.epochUnix,''k''); hmark_df!(?).LineStyle = ''--'';',1:numel(h),1:2)
%hca = irf_panel('Vi fpi hpca');
%irf_zoom(hca,'y')

%hca = irf_panel('feeps omni mms 1234'); hca.CLim = [-1 100];
%hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
%hca = irf_panel('B gsm'); hca.YLim = 0.99*[-17 22];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [1e-2 1e3];
%hca = irf_panel('Tepar/Teperp'); hca.YLim = [0.5 1.5];

%%
%recon1 = 1;
twpe = 17000;
varstrs = {'By','Ez','ne','JxBz','vex','vepar'}';
clims = {0.2*[-1 1],1*[-1 1],[0 1],[-1 1],5*[-1 1],5*[-1 1]};
h = no02m.twpelim(twpe).xlim([80 120]).zlim([-7 7]).plot_map(varstrs,'A',0.5,'clim',clims);
%%
pic = rec2;
twpe = 10000;
varstrs = {'By','Ez','ne','JxBz','vex','vepar'}';
clims = {0.2*[-1 1],1*[-1 1],[0 1],[-1 1],5*[-1 1],5*[-1 1]};
xlim = mean(pic.xi) + [-20 20];
zlim = [-7 7];

varstrs = {'By','Ez'}';
clims = {0.5*[-1 1],3*[-1 1],[0 1],[-1 1],5*[-1 1],5*[-1 1]};
xlim = mean(pic.xi) + [-20 20];
zlim = [-4 4];

pic = no02m;
twpe = 21700;
varstrs = {'By','Ez','vepar'}';
clims = {0.5*[-1 1],1*[-1 1],3*[-1 1],[-1 1],5*[-1 1],5*[-1 1]};
xlim = mean(pic.xi) + [-20 20];
xlim = [100 140];
zlim = [-7 7];

h = pic.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map(varstrs,'A',0.5,'clim',clims,'smooth',2);

%% timemap of Bz, showing islands
pic = no02m;
varstrs = {'Bz'}';
clims = {0.2*[-1 1],1*[-1 1],3*[-1 1],[-1 1],5*[-1 1],5*[-1 1]};
xlim = mean(pic.xi) + [-20 20];
xlim = [100 140];

h(1) = subplot(3,1,1);
h(2) = subplot(3,1,2);
h(3) = subplot(3,1,3);

%ax1 = no02m.xlim(mean(no02m.xi)+[-45 45]).zlim([-0.5 0.5]).plot_timemap(h(1),'xt',varstrs,'clim',clims);
%ax2 =  rec2.zlim([-0.5 0.5]).plot_timemap(h(2),'xt',varstrs,'clim',clims);
%ax3 =  rec1.zlim([-0.5 0.5]).plot_timemap(h(3),'xt',varstrs,'clim',clims);

ax1 = no02m.xlim(mean(no02m.xi)+[-45 45]).zlim([-0.5 0.5]).plot_timemap(h(1),'xt',varstrs,'clim',{[-0.99 0.99]});
ax2 =  rec2.zlim([-0.5 0.5]).plot_timemap(h(2),'xt',varstrs,'clim',{[-0.7 0.7]});
ax3 =  sus.zlim([-0.5 0.5]).plot_timemap(h(3),'xt',varstrs,'clim',{[-0.099 0.099]});

h(1).YLim(1) = 61;
h(2).YLim(1) = 11;
h(3).YLim(1) = 21;

c_eval('h(?).Title = [];',1:numel(h))
c_eval('h(?).Box = ''on'';',1:numel(h))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).FontSize = 14;',1:numel(h))
compact_panels(0.01,0.05)

%% timemap of Bz, showing density
pic = no02m;
varstrs = {'ni'}';
clims = {[0 1.3],[0 1.3],[0 1.3]};
xlim = mean(pic.xi) + [-20 20];
xlim = [100 140];

h(1) = subplot(3,1,1);
h(2) = subplot(3,1,2);
h(3) = subplot(3,1,3);

%ax1 = no02m.xlim(mean(no02m.xi)+[-45 45]).zlim([-0.5 0.5]).plot_timemap(h(1),'xt',varstrs,'clim',clims);
%ax2 =  rec2.zlim([-0.5 0.5]).plot_timemap(h(2),'xt',varstrs,'clim',clims);
%ax3 =  rec1.zlim([-0.5 0.5]).plot_timemap(h(3),'xt',varstrs,'clim',clims);

ax1 = no02m.xlim(mean(no02m.xi)+[-45 45]).zlim([-0.5 0.5]).plot_timemap(h(1),'xt',varstrs,'clim',{[0 1.4]});
ax2 =  rec2.zlim([-0.5 0.5]).plot_timemap(h(2),'xt',varstrs,'clim',{[0 1.4]});
ax3 =  rec1.zlim([-0.5 0.5]).plot_timemap(h(3),'xt',varstrs,'clim',{[0 1.4]});

h(1).YLim(1) = 61;
h(2).YLim(1) = 11;
h(3).YLim(1) = 21;

c_eval('h(?).Title = [];',1:numel(h))
c_eval('h(?).Box = ''on'';',1:numel(h))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).FontSize = 14;',1:numel(h))
compact_panels(0.01,0.05)

%% Denstiy fraction
pic = no02m;
twci = 85;

varstrs = {'n(1)./n([1 3 5])'}';
clims = {[0 1]};
xlim = mean(pic.xi) + [-80 80];
zlim = [-8 8]*0.99;

h(1) = subplot(2,1,1);
h(2) = subplot(2,1,2);

ax1 = pic.twcilim(60).xlim(xlim).zlim(zlim).plot_map(h(1),varstrs,'A',1,'clim',clims,'cbarlabels',{'n_{hot}/n_{tot}'},'sep');
ax2 = pic.twcilim(110).xlim(xlim).zlim(zlim).plot_map(h(2),varstrs,'A',1,'clim',clims,'cbarlabels',{'n_{hot}/n_{tot}'},'sep');



%%
pic = rec2;
varstrs = {'By','Ez','ne','JxBz','vex','vix'}';
clims = {0.2*[-1 1],1*[-1 1],[0 1],[-1 1],5*[-1 1],5*[-1 1]};
xlim = mean(pic.xi) + [-20 20];
h = pic.xlim(xlim).zlim([-7 7]).movie(varstrs,'A',0.5,'clim',clims,'filename',[printpath 'E05_By-Ez_ne_JxBz_vec_vix']);

%%
pic = rec2;
varstrs = {'Bz','vex','vix','Jx','Ey.*Bz'}';
clims = {[-1 1],2*[-1 1],2*[-1 1],2*[-1 1],2*[-1 1]};
xlim = mean(pic.xi) + [-20 20];
h = pic.xlim(xlim).zlim([-1 1]).plot_timemap('xt',varstrs,'A',0.5,'clim',clims,'filename',[printpath 'E05_By-Ez_ne_JxBz_vec_vix']);

%% Tsyganenko illustration of tail
% Set up input parameters
%clear all
%clear C;
ic = 0;


season = 'summer solstice';
activity = 'active';
plane = 'xy';
switch [season,activity,plane]
  case ['winter solstice','active','xz']
    ic = ic + 1; % Winter solstice, this one looksgood now, using mlat/mlon
    C(ic).date = [1990,12,22,00,00,00];
    C(ic).doy = day(datetime(C(ic).date(1),C(ic).date(2),C(ic).date(3)),'dayofyear');
    C(ic).parmod = [8.7,-128,0.3,-10.9,33.1,28.5]; % Pdyn, Dst, ByGSM, BzGSM, G1, G2
    C(ic).geolat = [25:5:165];
    C(ic).geolon = 180; % in the middle of the night, geolon = 180 is at midnight
    C(ic).mlat = [-142 35:5:170 -140:5:-105]; %-145 140 
    C(ic).mlon = [180 0*ones(size(35:5:170)) 180*ones(size(-140:5:-105))];% 180 180     
  case ['winter solstice','active','xy']
    ic = ic + 1; 
    C(ic).date = [1990,12,22,06,00,00];
    C(ic).doy = day(datetime(C(ic).date(1),C(ic).date(2),C(ic).date(3)),'dayofyear');
    C(ic).parmod = [8.7,-128,0.3,-10.9,33.1,28.5]; % Pdyn, Dst, ByGSM, BzGSM, G1, G2
    C(ic).geolat = [25:5:165];
    C(ic).geolon = 180; % in the middle of the night, geolon = 180 is at midnight
    C(ic).mlat = [80 -90];
    C(ic).mlon = [0 180];
    
    C(ic).mlat = [45*ones(1,18), 70*ones(1,18)]; %-145 140 
    C(ic).mlon = [0:20:350, 0:20:350];% 180 180 
    
    C(ic).mlat = [50*ones(1,36)]; %-145 140 
    C(ic).mlon = [0:10:350];% 180 180 
    
    C(ic).mlat = [55*ones(1,36)]; %-145 140 
    C(ic).mlon = [0:10:350];% 180 180   
  case ['winter solstice','inactive','xy']
    ic = ic + 1; 
    C(ic).date = [1990,12,22,06,00,00];
    C(ic).doy = day(datetime(C(ic).date(1),C(ic).date(2),C(ic).date(3)),'dayofyear');
    C(ic).parmod = [2,-70,0.3,-5,33.1,28.5]; % Pdyn, Dst, ByGSM, BzGSM, G1, G2
    C(ic).geolat = [25:5:165];
    C(ic).geolon = 180; % in the middle of the night, geolon = 180 is at midnight
    C(ic).mlat = [80 -90]; 
    C(ic).mlon = [0 180];
    
    C(ic).mlat = [45*ones(1,18), 70*ones(1,18)]; %-145 140 
    C(ic).mlon = [0:20:350, 0:20:350];% 180 180  
    
    
    C(ic).mlat = [50*ones(1,36)]; %-145 140 
    C(ic).mlon = [0:10:350];% 180 180 
    
    C(ic).mlat = [55*ones(1,36)]; %-145 140 
    C(ic).mlon = [0:10:350];% 180 180     
  case ['summer solstice','active','xy']
    ic = ic + 1; % Summer solstice, looks good, using mlat/mlon
    C(ic).date = [1990,06,21,00,00,00];
    C(ic).doy = day(datetime(C(ic).date(1),C(ic).date(2),C(ic).date(3)),'dayofyear');
    C(ic).parmod = [8.7,-128,0.3,-10.9,33.1,28.5]; % Pdyn, Dst, ByGSM, BzGSM, G1, G2
    C(ic).geolat = [25:2:165];
    C(ic).geolon = 180; % in the middle of the night, geolon = 180 is at midnight
    C(ic).mlat = [45*ones(1,18), 70*ones(1,18)]; %-145 140 
    C(ic).mlon = [0:20:350, 0:20:350];% 180 180 
    
    C(ic).mlat = [50*ones(1,36)]; %-145 140 
    C(ic).mlon = [0:10:350];% 180 180 
    
    C(ic).mlat = [60*ones(1,36)]; %-145 140 
    C(ic).mlon = [0:10:350];% 180 180 
end

    ic = ic + 1; % Summer solstice, looks good, using mlat/mlon
    C(ic).date = [1990,07,25,22,08,00];
    C(ic).doy = day(datetime(C(ic).date(1),C(ic).date(2),C(ic).date(3)),'dayofyear');
    C(ic).parmod = [8.7,-20,0.3,-10.9,33.1,28.5]; % Pdyn, Dst, ByGSM, BzGSM, G1, G2
    C(ic).geolat = [25:2:165];
    C(ic).geolon = 180; % in the middle of the night, geolon = 180 is at midnight
    C(ic).mlat = [45*ones(1,18), 70*ones(1,18)]; %-145 140 
    C(ic).mlon = [0:20:350, 0:20:350];% 180 180 
    
    C(ic).mlat = [50*ones(1,36)]; %-145 140 
    C(ic).mlon = [0:10:350];% 180 180 
    
    C(ic).mlat = [60*ones(1,36)]; %-145 140 
    C(ic).mlon = [0:10:350];% 180 180 
    

colors = [obs_colors('matlab'); 0 0 0; obs_colors('matlab')];
colors = repmat([0 0 0],10,1);


% Run Tsyganenko model
units = irf_units;
method_cs = 'new_mlat'; % 'geo;
doDownload = 0;

for ic = 1%[5 6 9 10];%[5 6 7]%1:2%:numel(C)
  C(ic).C = [];

  % Collect input data
  GEOPACK_RECALC(C(ic).date(1),C(ic).doy,C(ic).date(4),C(ic).date(5),C(ic).date(6));
  global GEOPACK1
  
  if doDownload % Base parmod on downloaded omnidata based on time interval
    tint = irf_time(C(ic).date)+[-60*60 0];    
    disp(['Loading omni data: ' irf_time(tint(1),'epoch>utc_yyyy-mm-ddTHH:MM:SS') ' to ' irf_time(tint(2),'epoch>utc_yyyy-mm-ddTHH:MM:SS')])
    omni_data = irf_get_data_omni(tint,'n,v,By,Bz','omni_min');
    omni_hour = irf_get_data_omni(tint,'dst','omni_hour');

    t   = omni_data(:,1); % s
    n   = omni_data(:,2); % cc
    v   = omni_data(:,3); v(v>1500) = NaN; % km/s
    By  = omni_data(:,4); % nT
    Bz  = omni_data(:,5); % nT
    dst = omni_hour(:,2); % 
    %Ma  = omni_data(:,7); % Mach number 

    mp = units.mp; % proton mass
    m = 0.98*mp+0.02*2*mp; % 98 percent hydrogen, 2 percent helium
    Dp = n.*1e6*m.*(v*1e3).^2*1e9; % Dynamic pressure, nPa
    Dp(isnan(Dp)) = [];
    [G1,G2] = G1G2par(By,Bz,v);
    
    PARMOD = [Dp(end), dst(end), By(end), Bz(end), (G1), G2];
    C(ic).parmod = PARMOD; 
  end
  
  
  PARMOD = C(ic).parmod;
  
  GEOLAT = C(ic).geolat;
  lon = C(ic).geolon;
  %k=1.25;

  % Set up coordinate system
  switch method_cs
    case 'geo'
      col=(90-GEOLAT);
      f=pi/180;
      lon=lon*f;
      col=col*f;

      l=length(col);
      s=length(lon);

      for j=1:length(lon)
          long=lon(j);
          for i=1:length(col)
          [xgeo(i,j), ygeo(i,j), zgeo(i,j)] = GEOPACK_SPHCAR(1,col(i),long,1);
          end
      end

      xgeof=reshape(xgeo,s*l,1);
      ygeof=reshape(ygeo,s*l,1);
      zgeof=reshape(zgeo,s*l,1);

      for i=1:length(xgeof)
         [xgsm(i), ygsm(i), zgsm(i)] = GEOPACK_GEOGSM(xgeo(i), ygeo(i), zgeo(i),1);
      end
    case 'ml'
      %mlat = [35:2:145]; XLAT = mlat;
      mlat = C(ic).mlat; XLAT = mlat;
      mlon = C(ic).mlon; XLON = mlon;
      T    = (90.-XLAT)*(pi/180); % from elevation angle to polar angle and deg to rad
      XL   = XLON*(pi/180); % deg to rad
      xgsm = sin(T).*cos(XL);
      ygsm = sin(T).*sin(XL);
      zgsm = cos(T);
      
    case 'new_mlat'
      mlat = C(ic).mlat; XLAT = mlat;
      mlon = C(ic).mlon; XLON = mlon;
      
      COLAT    = (90.-XLAT)*(pi/180); % from elevation angle to polar angle and deg to rad
      %disp([COLAT, XLON])
      [XGEO,YGEO,ZGEO] = GEOPACK_SPHCAR(1.,COLAT,XLON*(pi/180),1);
      %disp([XGEO,YGEO,ZGEO])
      % C
      % C   TRANSFORM GEOGRAPHICAL GEOCENTRIC COORDS INTO SOLAR MAGNETOSPHERIC ONES:
      % C
      [XGSM,YGSM,ZGSM] = GEOPACK_GEOGSM (XGEO,YGEO,ZGEO,1);
      xgsm = XGSM;
      ygsm = YGSM;
      zgsm = ZGSM;
  end
  
  % SPECIFY TRACING PARAMETERS:  
  lonix = length(xgsm); % CN: not used?
  
  %DIR   = 1; % (TRACE THE LINE WITH A FOOTPOINT IN THE NORTHERN HEMISPHERE, THAT IS,  ANTIPARALLEL TO THE MAGNETIC FIELD)  
  RLIM  = 40; % (LIMIT THE TRACING REGION WITHIN R=80 Re)
  R0    = 1; %(LANDING POINT WILL BE CALCULATED ON THE SPHERE R=1,I.E. ON THE EARTH'S SURFACE)   
  IOPT  = 0; % (IN THIS EXAMPLE IOPT IS JUST A DUMMY PARAMETER,  WHOSE VALUE DOES NOT MATTER)

  % TRACE THE FIELD LINE:  
  figure(100)
  for i=1:length(xgsm)
    if C(ic).mlat(i) && C(ic).mlon(i) == 0, DIR = -1;
    else DIR = -1;
    end
    if C(ic).mlon(i) == 0, DIR = 1;
    else DIR = -1;
    end
    if C(ic).mlat(i) > 0, DIR = 1; else, DIR = -1; end
    
      [XF,YF,ZF,XXf,YYf,ZZf,M] = GEOPACK_TRACE(xgsm(i),ygsm(i),zgsm(i),DIR,RLIM,R0,IOPT,PARMOD,'T01','GEOPACK_IGRF_GSM');
      %[Xi,Yi,Zi,XXi,YYi,ZZi,Mi] = GEOPACK_TRACE(xgsm(i),ygsm(i),zgsm(i),DIR,RLIM,R0,IOPT,PARMOD1,'T01','GEOPACK_IGRF_GSM');

      [XSEf,YSEf,ZSEf] = GEOPACK_GSMGSE(XXf,YYf,ZZf,1);
      hold on
      %plot3(XSEf,YSEf,ZSEf,'color',colors(ic,:))
      %view([0 1 0])
      
          %XSEf = XXf;
          %YSEf = YYf;
          %ZSEf = ZZf;
      
      switch plane
        case 'xz'
          XSEf = XF;
          YSEf = YF;
          ZSEf = ZF;
          
          hca = subplot(3,1,1);
          plot(hca,XSEf,ZSEf,'color',colors(ic,:))        
          hold(hca,'on')
          plot(hca,XSEf(1),ZSEf(1),'o','color',[0 1 0])
          xlabel(hca,'X_{GSE} [R_E]')
          ylabel(hca,'Z_{GSE} [R_E]')
          axis(hca,'equal')  
          hca.Box = 'on';
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          
          hca = subplot(3,1,2);
          plot(hca,XSEf,YSEf,'color',colors(ic,:))        
          hold(hca,'on')
          plot(hca,XSEf(1),YSEf(1),'o','color',[0 1 0])
          xlabel(hca,'X_{GSE} [R_E]')
          ylabel(hca,'Y_{GSE} [R_E]')
          axis(hca,'equal')          
          hold(hca,'on')
          hca.Box = 'on';
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          
          hca = subplot(3,1,2);
          plot3(hca,XSEf,YSEf,ZSEf,'color',colors(ic,:))        
          hold(hca,'on')
          plot3(hca,XSEf(1),YSEf(1),ZSEf(1),'o','color',[0 1 0])
          xlabel(hca,'X_{GSE} [R_E]')
          ylabel(hca,'Y_{GSE} [R_E]')
          zlabel(hca,'Z_{GSE} [R_E]')
          axis(hca,'equal')          
          hold(hca,'on')
          hca.Box = 'on';
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          hca.ZGrid = 'on';
        case 'xy'
          
          hca = subplot(1,2,1);
          plot(hca,XSEf,ZSEf,'color',colors(ic,:))     
          hold(hca,'on')
          plot(hca,XSEf(1),ZSEf(1),'.','color',[0 0 0])
          plot(hca,XSEf(end),ZSEf(end),'.','color',[0 0 0])
          xlabel(hca,'X_{GSE} [R_E]')
          ylabel(hca,'Z_{GSE} [R_E]')
          %axis(hca,'equal')          
          hold(hca,'on')
          hca.Box = 'on';
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          
          hca = subplot(1,2,2);
          plot(hca,XSEf,YSEf,'color',colors(ic,:))     
          hold(hca,'on')
          plot(hca,XSEf(1),YSEf(1),'.','color',[0 0 0])
          plot(hca,XSEf(end),YSEf(end),'.','color',[0 0 0])
          xlabel(hca,'X_{GSE} [R_E]')
          ylabel(hca,'Y_{GSE} [R_E]')
          axis(hca,'equal')          
          hold(hca,'on')
          hca.Box = 'on';
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          
          if 0
          hca = subplot(3,1,3);
          plot3(hca,XSEf,YSEf,ZSEf,'color',colors(ic,:))        
          hold(hca,'on')
          plot3(hca,XSEf(1),YSEf(1),ZSEf(1),'.','color',[0 0 0])
          plot3(hca,XSEf(end),YSEf(end),ZSEf(end),'.','color',[0 0 0])
          xlabel(hca,'X_{GSE} [R_E]')
          ylabel(hca,'Y_{GSE} [R_E]')
          zlabel(hca,'Z_{GSE} [R_E]')
          %axis(hca,'equal')          
          hold(hca,'on')
          hca.Box = 'on';
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          hca.ZGrid = 'on';
          end
      end
          
      % 
      drawnow
      
      C(ic).C(i).X_GSM = XXf;
      C(ic).C(i).Y_GSM = YYf;
      C(ic).C(i).Z_GSM = ZZf;
      C(ic).C(i).X_GSE = XSEf;
      C(ic).C(i).Y_GSE = YSEf;
      C(ic).C(i).Z_GSE = ZSEf;

  end

  %axis([0 2 -1 1])
  % set(gca,'xtick',0:0.4:2);
  % set(gca,'ytick',-1:0.4:1);
  %zlabel('Z_{GSE} [R_E]')

end

% Add bowshock and magnetopause
h = findobj(gcf,'type','axes'); h = h(end:-1:1);

Bz = C(1).parmod(4);
Dp = C(1).parmod(1);

%
[x_mp,y_mp,omni_mp] = irf_magnetosphere('mp_shue1998',Dp,Bz);
[x_bs,y_bs,omni_bs] = irf_magnetosphere('bs',Dp,Bz);


hca = h(1);
hold(hca,'on')
%xx = [x_bs(end:-1:1),x_bs];
xx = [x_bs(end:-1:1),x_bs];
yy = [y_bs(end:-1:1),-y_bs];
incl = intersect(find(yy>-2.1),find(yy<4.7));
plot(hca,xx(incl),yy(incl),'k')


xx = [x_mp(end:-1:1),x_mp];
yy = [y_mp(end:-1:1),-y_mp];
incl = intersect(find(yy>-2.1),find(yy<4.7));
plot(hca,xx(incl),yy(incl),'k')
%plot(hca,[x_mp(end:-1:1),x_mp],[y_mp(end:-1:1),-y_mp],'k')

plot(hca,[-8 19.6 19.6 -8 -8],[-2.1 -2.1 4.7 4.7 -2.1],'k');
hold(hca,'off')


hca = h(2);
hold(hca,'on')

xx = [x_bs(end:-1:1),x_bs];
yy = [y_bs(end:-1:1),-y_bs];
incl = intersect(find(yy>-11.1),find(yy<11.1));
plot(hca,xx(incl),yy(incl),'k')

xx = [x_mp(end:-1:1),x_mp];
yy = [y_mp(end:-1:1),-y_mp];
incl = intersect(find(yy>-11.1),find(yy<11.1));
plot(hca,xx(incl),yy(incl),'k')

plot(hca,[-8 19.6 19.6 -8 -8],[-11.1 -11.1 11.1 11.1 -11.1],'k');

hold(hca,'off')

c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:2)
c_eval('h(?).Visible = ''off'';',1:2)


%% The two DFs or DF/island
ic = 1;

time_df1 = irf_time('2017-07-25T22:10:06.10Z','utc>EpochTT');
time_df2 = irf_time('2017-07-25T22:10:28.00Z','utc>EpochTT');

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
feeps_ion_omni1234 = feeps_ion_omni1;
%feeps_ion_omni1234.data = (feeps_ion_omni1.data + feeps_ion_omni2.resample(feeps_ion_omni1).data + feeps_ion_omni3.resample(feeps_ion_omni1).data + feeps_ion_omni4.resample(feeps_ion_omni1).data)/4;
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

fontsize = 14;
fontsize_leg = 13;

npanels = 6;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

isub = 0;
zoomy = [];
cs = 'gse';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'B_x','B_y','B_z'}',[1.02 0.9],'fontsize',fontsize_leg);
end 
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',12);  
end
if 0 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i{\perp}x}','v_{i{\perp}y}','v_{i{\perp}z}','v_{i||}'}',[1.02 0.9],'fontsize',fontsize_leg);
end
if 0 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?.x,%sVi?.y,%sVi?.z},''comp'');',cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize_leg);
end
if 0 % beta
  hca = irf_panel('beta');
  betalow = 1/0.3;
  betalow = 3;
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,betalow},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-5:1:5];
  hca.YLabel.Interpreter = 'tex';
  hca.YLim = [1.01e-2 1e3];
  irf_legend(hca,['\beta < ' num2str(betalow,'%.2f')],[0.02 0.62],'color',hp.FaceColor,'fontsize',fontsize_leg)
end
if 0 % JotE
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('JdotE');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{JdotE?.resample(gseVi?)},''comp'');',ic)  
  hca.YLabel.String = {'J\cdot E','(nW/m^3)'}; % 
  set(hca,'ColorOrder',mms_colors('xyza'))
end

elim_feeps = [8e4 Inf];
if 0 % FEEPS Pitch all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps pitch mms 1234');  
  %elim = [0 32000];
  specrec = feeps_ion_pa1234.elim(elim_feeps).specrec('pitchangle');
  %specrec.f = f*1e-3;
  irf_spectrogram(hca,specrec);
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta^{FEEPS}_i','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  hca.YTick = 0:60:180;
end
if 1 % i DEF feeps omni all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps omni mms 1234'); 
  specrec = feeps_ion_omni1.elim(elim_feeps).specrec('energy');
  specrec.f = specrec.f*1e-3;
  specrec.f_label =  {'E_i (keV)'};
  [hout,hcb] = irf_spectrogram(hca,specrec,'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(keV)'};   
end
if 0 % i DEF omni, FPI
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 0 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','M','N'];
  comps = ['x','y','z'];
  for icomp = 1:3 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end


if 1 % i psd x,y,z, 3 panels
  for comp = ['x','y','z']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?_700;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 0 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
    %hca.XGrid = 'on';
    %hca.YGrid = 'on';
    %hca.Layer = 'top';
    %irf_legend(hca,sprintf('E > %.0f eV',if1D.ancillary.energy(1,1)-if1D.ancillary.delta_energy_minus(1,1)),[0.02 0.08],'color','k','fontsize',fontsize_leg)
  end
end

if 0 % e DEF feeps omni all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps e omni mms 1234'); 
  specrec = feeps_ele_omni1.elim(elim_feeps).specrec('energy');
  specrec.f = specrec.f*1e-3;
  specrec.f_label =  {'E_e (keV)'};
  [hout,hcb] = irf_spectrogram(hca,specrec,'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_e^{FEEPS}','(keV)'};   
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
  hca.CLim = [-5 -1];
end

if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{Te?par/Te?perp},''comp'');',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 5];
end
if 0 % Tepar/Teperp log
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp patch');
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_patch(hca,{Te?par/Te?perp,1});',ic) 
  c_eval('tsdata = Te?par.resample(gseVi?)/Te?perp.resample(gseVi?);',ic)
  c_eval('irf_patch(hca,{tsdata,1});',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 2];
  hca.YLabel.Interpreter = 'tex';
end
if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp patch');
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_patch(hca,{Te?par/Te?perp,1});',ic) 
  c_eval('tsdata = Te?par.resample(gseVi?)/Te?perp.resample(gseVi?);',ic)
  tsdata = tsdata.tlim(irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:11:00.00Z'));
  if 1
    tsdata.data(ne1.resample(tsdata).data<0.01) = NaN; % doesn't work well with patch
    tsdata.data = log10(tsdata.data);
    c_eval('hp = irf_patch(hca,{tsdata,0});',ic)  
    hp.FaceColor = [0 0 0];
    hp.EdgeColor = [0 0 0];
  else
    tsdata.data(ne1.resample(tsdata).data<0.01) = NaN;
    tsdata.data = log10(tsdata.data);
    irf_plot(hca,tsdata,'k')
  end
  hca.YLabel.String = 'log_{10}(T_{e||}/T_{e\perp})';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [-0.201 0.201];
  hca.YLabel.Interpreter = 'tex';
end
%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0],'fontsize',fontsize_leg)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize_leg;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);

c_eval('h(?).FontSize = fontsize;',1:numel(h))


h(end).XTickLabelRotation = 0;

irf_zoom(h,'x',tint_action)
%irf_zoom(h,'y')
c_eval('hmark_df!(?) = irf_pl_mark(h(?),time_df!.epochUnix,''k''); hmark_df!(?).LineStyle = ''--'';',1:numel(h),1:2)


%hca = irf_panel('feeps omni mms 1234'); hca.CLim = [-1 100];
hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
hca = irf_panel('B gsm'); hca.YLim = 0.99*[-16 25];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [1e-2 1e3];
%hca = irf_panel('Tepar/Teperp'); hca.YLim = [0.5 1.5];


%hca1 = irf_panel('feeps omni mms 1234');  
%hca2 = irf_panel('i DEF omni'); 
%hca2.YLim(2) = hca1.YLim(1);
%h(end).YLim  = [-350 1150];


% text: 'Moderate beta...'
if 0 % pressure panel included
  %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]-dy);
elseif 0
  %%
  annotation('textarrow',[0.3 0.3],[.7 .72],'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]);
  annotation('textarrow',[0.63 0.63],[.7 .72],'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]);  
elseif 0
  %%
  dy = 0.08;
  dx = 0.07; 
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63]-dx,[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63]-dx,[.68 .64]-dy);
end

%% Physical scales
npanels = 5;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');
fontsize_leg = 14;
fontsize = 14;
isub = 1;

cs = 'gse';
if 1 % B cs
  isub = isub + 1;
  
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'B_x','B_y','B_z'}',[1.02 0.9],'fontsize',fontsize_leg);
end 
if 1 % beta
  hca = irf_panel('beta');
  betalow = 1/0.3;
  betalow = 3;
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  %hp = irf_patch(hca,{beta1,betalow},'smaller');
  %hp.EdgeColor = 'none';
  %hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-5:1:5];
  hca.YLabel.Interpreter = 'tex';
  hca.YLim = [1.01e-2 1e3];
  %irf_legend(hca,['\beta < ' num2str(betalow,'%.2f')],[0.02 0.62],'color',hp.FaceColor,'fontsize',fontsize_leg)
end
if 1 % frequencies
  hca = irf_panel('freq');
  fpe1_ = fpe1;
  fpe1_.data(ne1.data<0.02) = NaN;
  set(hca,'ColorOrder',pic_colors('matlab'))  
  irf_plot(hca,{fce1,fpe1_,fcp1.resample(ni1),fpp1.resample(ni1)},'comp')
  hca.YLabel.String = {'Frequency','(Hz)'};
  hca.YScale = 'log';
  set(hca,'ColorOrder',pic_colors('matlab'))  
  irf_legend(hca,{'\omega_{ce}','\omega_{pe}','\omega_{cp}','\omega_{pp}'}',[1.02,0.98])
end
if 1 % lengths
  hca = irf_panel('length');
  Lp1_ = Lp1;
  Lp1_.data(ne1.data<0.1) = NaN;
  Le1_ = Le1;
  Le1_.data(ne1.data<0.02) = NaN;
  re1_ = re1;
  re1_.data(ne1.data<0.02) = NaN;
  set(hca,'ColorOrder',pic_colors('matlab'))  
  irf_plot(hca,{re1_,Le1_,rp1.resample(ni1),Lp1_.resample(ni1)},'comp')
  hca.YLabel.String = {'Length','(km)'};
  %hca.YScale = 'log';
  set(hca,'ColorOrder',pic_colors('matlab'))  
  irf_legend(hca,{'\rho_e','d_e','\rho_p','d_p'}',[1.02,0.98])
  hca.YScale = 'log';
end
if 0 % lengths
  hca = irf_panel('length');
  Lp1_ = Lp1;
  Lp1_.data(ne1.data<0.1) = NaN;
  Le1_ = Le1;
  Le1_.data(ne1.data<0.02) = NaN;
  re1_ = re1;
  re1_.data(ne1.data<0.02) = NaN;
  set(hca,'ColorOrder',pic_colors('matlab'))  
  irf_plot(hca,{re1_,Le1_,rp1.resample(ni1),Lp1_.resample(ni1)},'comp')
  hca.YLabel.String = {'Length','(km)'};
  %hca.YScale = 'log';
  set(hca,'ColorOrder',pic_colors('matlab'))  
  irf_legend(hca,{'\rho_e','d_e','\rho_p','d_p'}',[1.02,0.98])
end
if 1 % speeds
  hca = irf_panel('speed');
  vA1_ = vA1;
  vA1_.data(ne1.data<0.1) = NaN;
  set(hca,'ColorOrder',pic_colors('matlab'))  
  irf_plot(hca,{vA1_.resample(ni1)},'comp')
  hca.YLabel.String = {'Speed','(km/s)'};
  %hca.YScale = 'log';
  set(hca,'ColorOrder',pic_colors('matlab'))  
  irf_legend(hca,{'v_{A}'}',[1.02,0.98])
end

irf_zoom(h,'x',tint)

hca = irf_panel('freq');
hca.YTick = 10.^[-2:5];
hca.YLim = [1e-2 1e4]*0.99;

hca = irf_panel('speed');
hca.YLim = [0 4000]*0.99;
%hca.YLim = [3 1000]*0.99;


hca = irf_panel('length');
hca.YLim = [4 4000]*0.99; 
hca.YTick = 10.^[-2:5];

hca = irf_panel('length');
hca.YLim = [0 3000]*0.99;  

irf_plot_axis_align