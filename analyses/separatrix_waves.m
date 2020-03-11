df04n = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/data_h5/fields.h5');

%% Plot timeline of ni, ne and Ex at given z.
% turb = PIC('/Volumes/Fountain/Data/PIC/turbulencerun/data_h5/fields.h5');
% df04n = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/data_h5/fields.h5');
% df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
pic = df04;
zilim = 3 + [-0.1 0.1];
xilim = mean(pic.xi) + 50*[-1 1];
pic = pic.xlim(xilim).zlim(zilim).twcilim([00 500]);
meandim = 2; % z
ni = squeeze(mean(pic.ni,meandim));
ne = squeeze(mean(pic.ne,meandim));
vix = squeeze(mean(pic.vix,meandim));
vex = squeeze(mean(pic.vex,meandim));
%vixc = squeeze(mean(pic.vx([3 5]),meandim));
%vexc = squeeze(mean(pic.vx([4 6]),meandim));
%vixh = squeeze(mean(pic.vx([1]),meandim));
%vexh = squeeze(mean(pic.vx([2]),meandim));
Ex = squeeze(mean(pic.Ex,meandim));

nrows = 5;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % ni
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.twci,ni');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_i (n_0)';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
  
end
if 1 % ne
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.twci,ne');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e (n_0)';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
end
if 1 % vix
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.twci,vix');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_{ix,tot} (v_{A0})';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';  
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.Colormap = pic_colors('blue_red');
end
if 1 % vex
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.twci,vex');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_{ex,tot} (v_{A0})';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.Colormap = pic_colors('blue_red');
end
if 0 % vix cold
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.twci,vixc');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_{ix,cold} (v_{A0})';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';  
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.Colormap = pic_colors('blue_red');
end
if 0 % vex cold
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.twci,vexc');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_{ex,cold} (v_{A0})';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.Colormap = pic_colors('blue_red');
end
if 1 % Ex
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.twci,Ex');
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_x (B_0v_{A0})';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.Colormap = pic_colors('blue_red');
end

h(1).Title.String = sprintf('z = [%g, %g]',zilim(1),zilim(2));
hlinks = linkprop(h,{'XLim','YLim'});
compact_panels(0.01)

for ip = 1:npanels
  h(ip).XTick = 0:10:max(pic.xi);
  h(ip).YTick = 0:10:500;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).Layer = 'top';
end

%% Plot the distribution boxes locations on top of density maps
ds = ds04.dxlim([0 0.21]);
pic = df04;
twpe = ds.twpe;
nt = ds.nt;
h = setup_subplots(nt,1);

for itime = 1:nt
  hca = h(itime);
  imagesc(hca,pic.xi,pic.zi,pic.twpelim(twpe(itime)).n(5)')
  hca.YDir = 'normal';
  ds.twpelim(twpe(itime)).plot_boxes(hca);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_i/n_0';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'z/d_{i0}';
end

compact_panels(0.01)
hlinks = linkprop(h,{'CLim','XLim','YLim'});
h(1).XLim = [80 240];
h(1).YLim = [-7 7];
h(1).CLim = [0 2.2];

%% Chose one time and z to plot distribution (f, def) along x
twpe = 6000;
zvals = 0;
zi  = zvals;

ds = ds08.dxlim([0 0.21]).zfind(zi).twpelim(twpe);
x0 = (ds.xi2{:}+ds.xi1{:})/2;

xvals = sort(x0);%170:0.2:235;

vmax = [5 15 2 10 2 10];
vmax = [5 15 2 10 2 10];
nv = 101;

c_eval('f? = ds.reduce_1d(''x'',xvals,zvals,linspace(-vmax(?),vmax(?),nv),?,''vabs'');',1:1:4)
%%
f35 = ds.reduce_1d('x',xvals,zvals,linspace(-vmax(3),vmax(3),nv),[3 5],'vabs');
f46 = ds.reduce_1d('x',xvals,zvals,linspace(-vmax(4),vmax(4),nv),[4 6],'vabs');
% f1 = ds.reduce_1d('x',xvals,zvals,linspace(-vmax(1),vmax(1),nv),1,'vabs');
% f2 = ds.reduce_1d('x',xvals,zvals,linspace(-vmax(2),vmax(2),nv),2,'vabs');
% f3 = ds.reduce_1d('x',xvals,zvals,linspace(-vmax(3),vmax(3),nv),3,'vabs');
% f4 = ds.reduce_1d('x',xvals,zvals,linspace(-vmax(4),vmax(4),nv),4,'vabs');
% f5 = ds.reduce_1d('x',xvals,zvals,linspace(-vmax(5),vmax(5),nv),5,'vabs');
% f6 = ds.reduce_1d('x',xvals,zvals,linspace(-vmax(6),vmax(6),nv),6,'vabs');
disp('Done.')

%% Plot single species
xlim = xvals([1 end]) + [-10 10];%[170 240];
zlim = [-8 8];
pic = df08.twpelim(twpe).xlim(xlim).zlim(zlim);

fdist = 3;
fdiststr = {'hot ions','hot electrons',...
  'cold ions from the north','cold electrons from the north',...
  'cold ions from the south','cold electrons from the south'};

tlims = [0.7 0.2 0.2 0.2 0.2 0.2];

doLogRED = 0;
doLogDEF = 1;
isDEF = [];
isRED = [];
isMAP = [];

doV = 1;
doVExB = 1;
doE = 1;

nrows = 7;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;


fftmp = eval(sprintf('f%g(1,1)',fdist));
dx = fftmp.x(2)-fftmp.x(1);
% just plotting solution
%fftmp.x = [fftmp.x - (fftmp.x(2)-fftmp.x(1)) fftmp.x + (fftmp.x(2)-fftmp.x(1))];

z = unique(fftmp.z);

v.x = pic.xi;
v.vx = squeeze(mean(pic.zlim(z+[-0.1 0.1]).vx(fdist),2));
v.vy = squeeze(mean(pic.zlim(z+[-0.1 0.1]).vy(fdist),2));
v.vz = squeeze(mean(pic.zlim(z+[-0.1 0.1]).vz(fdist),2));
Bx = pic.zlim(z+[-0.1 0.1]).Bx; 
By = pic.zlim(z+[-0.1 0.1]).By;
Bz = pic.zlim(z+[-0.1 0.1]).Bz;
Babs = sqrt(Bx.^2+By.^2+Bz.^2);
Ex = pic.zlim(z+[-0.1 0.1]).Ex;
Ey = pic.zlim(z+[-0.1 0.1]).Ey;
Ez = pic.zlim(z+[-0.1 0.1]).Ez;
ExB = cross_product(Ex,Ey,Ez,Bx,By,Bz);
v.ExBx = mean(ExB.x./Babs.^2,2);
v.ExBy = mean(ExB.y./Babs.^2,2);
v.ExBz = mean(ExB.z./Babs.^2,2);
vlinestyle = '-';
vlinecolor = [1 1 1];
ExBlinestyle = '-';
ExBlinecolor = [0 0 0];


if 1 % map of where the distribution boxes are, n
  hca = h(isub); isub = isub + 1;
  isMAP(end+1) = isub - 1;
  imagesc(hca,pic.xi,pic.zi,pic.twpelim(twpe).n(fdist)')
  hca.YDir = 'normal';
  ds.zfind(zvals).twpelim(twpe).plot_boxes(hca);
  
  if 1 % A    
    hold(hca,'on')
    A = pic.twpelim(twpe).A;
    iplx = 1:5:pic.nx;
    iplz = 1:5:pic.nz;
    contour(hca,pic.xi(iplx),pic.zi(iplz),A(iplx,iplz)',-25:1:0,'k')
    hold(hca,'off')
  end
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('n_%g/n_0',fdist);
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'z/d_{i0}';
end
if 1 % map of where the distribution boxes are, Ex
  hca = h(isub); isub = isub + 1;
  isMAP(end+1) = isub - 1;
  imagesc(hca,pic.xi,pic.zi,pic.twpelim(twpe).Ex')
  hca.YDir = 'normal';
  ds.zfind(zvals).twpelim(twpe).plot_boxes(hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.Colormap = pic_colors('blue_red');
  
  if 1 % A    
    hold(hca,'on')
    A = pic.twpelim(twpe).A;
    iplx = 1:5:pic.nx;
    iplz = 1:5:pic.nz;
    contour(hca,pic.xi(iplx),pic.zi(iplz),A(iplx,iplz)',-25:1:0,'k')
    hold(hca,'off')
  end
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_x/B_0v_{A0}';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'z/d_{i0}';
end
if 1 % map of where the distribution boxes are, T
  hca = h(isub); isub = isub + 1; isMAP(end+1) = isub - 1;  
  toplot = pic.twpelim(twpe).p(fdist)'./pic.twpelim(twpe).n(fdist)';
  imagesc(hca,pic.xi,pic.zi,toplot)
  hca.YDir = 'normal';
  ds.zfind(zvals).twpelim(twpe).plot_boxes(hca);
  hca.CLim = [0 prctile(toplot(:),99.9)];
  hca.CLim = [0 tlims(fdist)];
  
  if 1 % A    
    clim = hca.CLim;
    hold(hca,'on')
    A = pic.twpelim(twpe).A;
    iplx = 1:5:pic.nx;
    iplz = 1:5:pic.nz;
    contour(hca,pic.xi(iplx),pic.zi(iplz),A(iplx,iplz)',-25:1:0,'k')
    hold(hca,'off')
    hca.CLim = clim;
  end
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'T/m_pv_{A0}^2';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'z/d_{i0}';
end

if 1 %f(vx)
  hca = h(isub); isub = isub + 1;  
  if doLogRED
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,log10(fftmp.fvx')); shading(hca,'flat'); 
  else
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,fftmp.fvx'); shading(hca,'flat'); 
  end
  isRED(end+1) = isub - 1;
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('f(x,v_x)',twpe,zi);
  hca.YLabel.String = 'v/v_{A0}';
  if doE % Ex
    hold(hca,'on')
    plot(hca,pic.xi,mean(Ex,2),'linewidth',1,'color',[0 1 1],'linestyle','-')
    hold(hca,'off')
  end
  if doV % vx
    hold(hca,'on')
    plot(hca,v.x,v.vx,'linewidth',1,'color',vlinecolor,'linestyle',vlinestyle)
    hold(hca,'off')
  end
  if doVExB % vExBx
    hold(hca,'on')
    plot(hca,v.x,v.ExBx,'linewidth',1,'color',ExBlinecolor,'linestyle',ExBlinestyle)
    hold(hca,'off')
  end
end 
if 1 % f(vy)
  hca = h(isub); isub = isub + 1;  
  if doLogRED
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,log10(fftmp.fvy')); shading(hca,'flat'); 
  else
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,fftmp.fvy'); shading(hca,'flat'); 
  end
  isRED(end+1) = isub - 1;
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('f(x,v_y)',twpe,zi);
  hca.YLabel.String = 'v/v_{A0}';
  
  if doE % Ey
    hold(hca,'on')
    plot(hca,pic.xi,mean(Ey,2),'linewidth',1,'color',[0 1 1],'linestyle','-')
    hold(hca,'off')
  end
  if doV % vy
    hold(hca,'on')
    plot(hca,v.x,v.vy,'linewidth',1,'color',vlinecolor,'linestyle',vlinestyle)
    hold(hca,'off')
  end
  if doVExB % vExB
    hold(hca,'on')
    plot(hca,v.x,v.ExBy,'linewidth',1,'color',ExBlinecolor,'linestyle',ExBlinestyle)
    hold(hca,'off')
  end
end
if 1 % f(vz)
  hca = h(isub); isub = isub + 1;  
  if doLogRED
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,log10(fftmp.fvz')); shading(hca,'flat'); 
  else
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,fftmp.fvz'); shading(hca,'flat'); 
  end
  isRED(end+1) = isub - 1;
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('f(x,v_z)',twpe,zi);
  hca.YLabel.String = 'v/v_{A0}';
  
  if doE % Ez
    hold(hca,'on')
    plot(hca,pic.xi,mean(Ez,2),'linewidth',1,'color',[0 1 1],'linestyle','-')
    hold(hca,'off')
  end
  if doV % vz
    hold(hca,'on')
    plot(hca,v.x,v.vz,'linewidth',1,'color',vlinecolor,'linestyle',vlinestyle)
    hold(hca,'off')
  end
  if doVExB % vExBx
    hold(hca,'on')
    plot(hca,v.x,v.ExBz,'linewidth',1,'color',ExBlinecolor,'linestyle',ExBlinestyle)
    hold(hca,'off')
  end
end
if 0 % def, on vscale
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,fftmp.x,fftmp.v(fftmp.v>=0),log10(fftmp.def')); shading(hca,'flat'); 
  isDEF(end+1) = isub - 1;
end
if 1 % def on log10(v^2) scale
  hca = h(isub); isub = isub + 1;
  pcolor(hca,fftmp.x,fftmp.v(fftmp.v>=0).^2/2,log10(fftmp.def')); shading(hca,'flat'); 
  isDEF(end+1) = isub - 1;
  hca.YScale = 'log';
  hca.YLabel.String = 'log_{10}(v^2/2)';  
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('DEF',twpe,zi);  
  %hca.YLim(1) = min(fftmp.v(fftmp.v>=0).^2/2);
  
  if 1 % vx
    hold(hca,'on')
    plot(hca,v.x,(v.vx.^2+v.vy.^2+v.vz.^2)/2,'linewidth',1,'color',vlinecolor,'linestyle',vlinestyle)
    hold(hca,'off')
  end
  if 1 % vExBabs
    ylim = hca.YLim;
    hold(hca,'on')
    plot(hca,v.x,(v.ExBx.^2+v.ExBy.^2+v.ExBz.^2)/2,'linewidth',1,'color',ExBlinecolor,'linestyle',ExBlinestyle)
    hold(hca,'off')
    hca.YLim = ylim;
  end
end

h(1).Title.String = sprintf('twpe = %g, iSpecies = %g (%s)',twpe,fdist,fdiststr{fdist});

h(1).CLim = [0 2.5];

h(isDEF).YLim(1) = 0.0005;

hlinksRED = linkprop(h(isRED),{'CLim','YLim'});
%hlinksRED.Targets(1).YLim = 3*[-1 1];
hlinksRED.Targets(1).CLim = [0 0.15];
hlinksRED.Targets(1).CLim = [0 0.1];

hlinksALL = linkprop(h,{'XLim'});

hlinksMAP = linkprop(h(isMAP),{'YLim'});
hlinksMAP.Targets(1).YLim = zlim;

for ipanel = 1:npanels
  h(ipanel).XTick = 0:10:500;
  h(ipanel).YTick = -20:1:20;
  h(ipanel).XGrid = 'on';
  h(ipanel).YGrid = 'on';
  h(ipanel).Layer = 'top';
  h(ipanel).FontSize = 14;
end

drawnow
compact_panels(0.01)

%% Plot combined species

pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);

if 1
  fdist = 46;
  fftmp = f46;
  iSpecies = [4 6];
  sp_str = '[4 6]';
else
  fdist = 35;
  fftmp = f35;
  iSpecies = [3 5];
  sp_str = '[3 5]';
end

tlims = [];
fdiststrall = {'hot ions','hot electrons',...
  'cold ions from the north','cold electrons from the north',...
  'cold ions from the south','cold electrons from the south'};

tlims = [0.7 0.2 0.2 0.2 0.2 0.2];


fdiststr = [fdiststrall{iSpecies(1)} ', ' fdiststrall{iSpecies(2)}];


doLogRED = 0;
doLogDEF = 1;
isDEF = [];
isRED = [];


z = unique(fftmp.z);
v.x = pic.xi;
v.vx = squeeze(mean(pic.zlim(z+[-0.1 0.1]).vx(iSpecies),2));
v.vy = squeeze(mean(pic.zlim(z+[-0.1 0.1]).vy(iSpecies),2));
v.vz = squeeze(mean(pic.zlim(z+[-0.1 0.1]).vz(iSpecies),2));
Bx = pic.zlim(z+[-0.1 0.1]).Bx; 
By = pic.zlim(z+[-0.1 0.1]).By;
Bz = pic.zlim(z+[-0.1 0.1]).Bz;
Babs = sqrt(Bx.^2+By.^2+Bz.^2);
Ex = pic.zlim(z+[-0.1 0.1]).Ex;
Ey = pic.zlim(z+[-0.1 0.1]).Ey;
Ez = pic.zlim(z+[-0.1 0.1]).Ez;
ExB = cross_product(Ex,Ey,Ez,Bx,By,Bz);
v.ExBx = mean(ExB.x./Babs.^2,2);
v.ExBy = mean(ExB.y./Babs.^2,2);
v.ExBz = mean(ExB.z./Babs.^2,2);
vlinestyle = '-';
vlinecolor = [1 1 1];
ExBlinestyle = '-';
ExBlinecolor = [0 0 0];

% Figure
nrows = 7;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % map of where the distribution boxes are, n
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.twpelim(twpe).n(iSpecies)')
  hca.YDir = 'normal';
  ds.zfind(zvals).twpelim(twpe).plot_boxes(hca);
  
  if 1 % A    
    hold(hca,'on')
    A = pic.twpelim(twpe).A;
    iplx = 1:5:pic.nx;
    iplz = 1:5:pic.nz;
    contour(hca,pic.xi(iplx),pic.zi(iplz),A(iplx,iplz)',-25:1:0,'k')
    hold(hca,'off')
  end
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('n_%g/n_0',fdist);
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'z/d_{i0}';
end
if 1 % map of where the distribution boxes are, Ex
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.twpelim(twpe).Ex')
  hca.YDir = 'normal';
  ds.zfind(zvals).twpelim(twpe).plot_boxes(hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.Colormap = pic_colors('blue_red');
  
  if 1 % A    
    hold(hca,'on')
    A = pic.twpelim(twpe).A;
    iplx = 1:5:pic.nx;
    iplz = 1:5:pic.nz;
    contour(hca,pic.xi(iplx),pic.zi(iplz),A(iplx,iplz)',-25:1:0,'k')
    hold(hca,'off')
  end
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('E_x/v_{A0}B_0',fdist);
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'z/d_{i0}';
end
if 1 % map of where the distribution boxes are, T
  hca = h(isub); isub = isub + 1; isMAP(end+1) = isub - 1;  
  toplot = pic.twpelim(twpe).p(iSpecies)'./pic.twpelim(twpe).n(iSpecies)';
  imagesc(hca,pic.xi,pic.zi,toplot)
  hca.YDir = 'normal';
  ds.zfind(zvals).twpelim(twpe).plot_boxes(hca);
  hca.CLim = [0 prctile(toplot(:),99.9)];
  hca.CLim = [0 tlims(iSpecies(1))];
  
  if 1 % A    
    clim = hca.CLim;
    hold(hca,'on')
    A = pic.twpelim(twpe).A;
    iplx = 1:5:pic.nx;
    iplz = 1:5:pic.nz;
    contour(hca,pic.xi(iplx),pic.zi(iplz),A(iplx,iplz)',-25:1:0,'k')
    hold(hca,'off')
    hca.CLim = clim;
  end
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'T/m_pv_{A0}^2';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'z/d_{i0}';
end

if 1 %f(vx)
  hca = h(isub); isub = isub + 1;  
  if doLogRED
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,log10(fftmp.fvx')); shading(hca,'flat'); 
  else
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,fftmp.fvx'); shading(hca,'flat'); 
  end
  isRED(end+1) = isub - 1;
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('f(x,v_x)',twpe,zi);
  hca.YLabel.String = 'v/v_{A0}';
  
  
  if doE % Ex
    hold(hca,'on')
    plot(hca,pic.xi,mean(Ex,2),'linewidth',1,'color',[0 1 1],'linestyle','-')
    hold(hca,'off')
  end
  if doV % vx
    hold(hca,'on')
    plot(hca,v.x,v.vx,'linewidth',1,'color',vlinecolor,'linestyle',vlinestyle)
    hold(hca,'off')
  end
  if doVExB % vExBx
    hold(hca,'on')
    plot(hca,v.x,v.ExBx,'linewidth',1,'color',ExBlinecolor,'linestyle',ExBlinestyle)
    hold(hca,'off')
  end
  
end 
if 1 % f(vy)
  hca = h(isub); isub = isub + 1;  
  if doLogRED
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,log10(fftmp.fvy')); shading(hca,'flat'); 
  else
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,fftmp.fvy'); shading(hca,'flat'); 
  end
  isRED(end+1) = isub - 1;
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('f(x,v_y)',twpe,zi);
  hca.YLabel.String = 'v/v_{A0}';
  
  
  if doE % Ey
    hold(hca,'on')
    plot(hca,pic.xi,mean(Ey,2),'linewidth',1,'color',[0 1 1],'linestyle','-')
    hold(hca,'off')
  end
  if doV % vy
    hold(hca,'on')
    plot(hca,v.x,v.vy,'linewidth',1,'color',vlinecolor,'linestyle',vlinestyle)
    hold(hca,'off')
  end
  if doVExB % vExBy
    hold(hca,'on')
    plot(hca,v.x,v.ExBy,'linewidth',1,'color',ExBlinecolor,'linestyle',ExBlinestyle)
    hold(hca,'off')
  end
end
if 1 % f(vz)
  hca = h(isub); isub = isub + 1;  
  if doLogRED
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,log10(fftmp.fvz')); shading(hca,'flat'); 
  else
    pcolor(hca,fftmp.x-0.5*dx,fftmp.v,fftmp.fvz'); shading(hca,'flat'); 
  end
  isRED(end+1) = isub - 1;
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('f(x,v_z)',twpe,zi);
  hca.YLabel.String = 'v/v_{A0}';
  
  
  if doE % Ez
    hold(hca,'on')
    plot(hca,pic.xi,mean(Ez,2),'linewidth',1,'color',[0 1 1],'linestyle','-')
    hold(hca,'off')
  end
  if doV % vz
    hold(hca,'on')
    plot(hca,v.x,v.vz,'linewidth',1,'color',vlinecolor,'linestyle',vlinestyle)
    hold(hca,'off')
  end
  if doVExB % vExBx
    hold(hca,'on')
    plot(hca,v.x,v.ExBz,'linewidth',1,'color',ExBlinecolor,'linestyle',ExBlinestyle)
    hold(hca,'off')
  end
end
if 0 % def, on vscale
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,fftmp.x,fftmp.v(fftmp.v>=0),log10(fftmp.def')); shading(hca,'flat'); 
  isDEF(end+1) = isub - 1;
end
if 1 % def on log10(v^2) scale
  hca = h(isub); isub = isub + 1;
  pcolor(hca,fftmp.x,fftmp.v(fftmp.v>=0).^2/2,log10(fftmp.def')); shading(hca,'flat'); 
  isDEF(end+1) = isub - 1;
  hca.YScale = 'log';
  hca.YLabel.String = 'log_{10}(v^2/2)';  
  hb = colorbar('peer',hca);
  hb.YLabel.String = sprintf('DEF',twpe,zi);  
  
  if doV % vx
    hold(hca,'on')
    plot(hca,v.x,(v.vx.^2+v.vy.^2+v.vz.^2)/2,'linewidth',1,'color',vlinecolor,'linestyle',vlinestyle)
    hold(hca,'off')
  end
  if doVExB % vExBabs
    ylim = hca.YLim;
    hold(hca,'on')
    plot(hca,v.x,(v.ExBx.^2+v.ExBy.^2+v.ExBz.^2)/2,'linewidth',1,'color',ExBlinecolor,'linestyle',ExBlinestyle)
    hold(hca,'off')
    hca.YLim = ylim;
  end
  
end

h(1).Title.String = sprintf('twpe = %g, iSpecies = %s (%s)',twpe,sp_str,fdiststr);

h(1).CLim = [0 2.5];

hlinksRED = linkprop(h(isRED),{'CLim','YLim'});
%hlinksRED.Targets(1).YLim = 3*[-1 1];
hlinksRED.Targets(1).CLim = [0 0.15];

hlinksALL = linkprop(h,{'XLim'});

hlinksMAP = linkprop(h(isMAP),{'YLim'});
hlinksMAP.Targets(1).YLim = zlim;

for ipanel = 1:npanels
  h(ipanel).XTick = 0:1:500;
  h(ipanel).YTick = -20:1:20;
  h(ipanel).XGrid = 'on';
  h(ipanel).YGrid = 'on';
  h(ipanel).Layer = 'top';
end

drawnow
compact_panels(0.01)

%% Instability analysis

