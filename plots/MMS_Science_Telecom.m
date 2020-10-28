%% Load PIC and PICDist objects
no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
tr100 = PICTraj('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories.h5');
ds100 = PICDist('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');
no02 = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02/data_h5/fields.h5');
bs = PIC('/Volumes/Fountain/Data/PIC/baselinev4/data_h5/fields.h5');

%% Initial conditions
pic = no02m;
comp = 'z';
twpe = 1000;
xlim = 53+0.5*[-1 1];
zlim = pic.zi([1 end]);
pic = pic.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');
varstrs = {{'Bx'};{'Jy'};{'n(1)','n([3 5])','n([1 3 5])'};{'Ey','Ez'};{'txx(1)','txx(3)','txx(5)','txx([1 3 5])'}};
varstrs = {{'Bx'};{'Jy'};{'n([1 3 5])'};{'Ey'};{'Ez'};{'viz'};{'vez'}};
varstrs = {{'Ey'};{'Ez'};{'vez'};{'tzz(4)'};{'txx(4)'}};
varstrs = {{'tzz([1])'};{'tzz(2)'};{'ni'}};
%varstrs = {{'Ey'};{'n(4)'};{'tzz(4)'}};
varstrs = {{'n([1])','n([3 5])'}};

h = pic.plot_line(comp,varstrs,'smooth',10);

%% Overview plot, before going to distributions
% Plotmap
twpe = 24000;
xlim = [50 150];
zlim = [-10 10];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);



varstrs = {'By';'Ey';'Jy';'n([1])';'n([3 5])';'A'};
varstrs = {'By';'Jy';'Jx';'Ey';'ni'};
varstrs = {'Epar';'log10(ne)';'ti';'te';'n([3 5])'};
varstrs = {'log10(ne)';'Ez';'viy';'vex';'n([3 5])'};
clims = {[-2 0.2],[-2 2],[-2 2],[-6 6],[0 0.4],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1]};
clims = {[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1]};
varstrs = {'Epar';'log10(ne)';'ti';'te';'n([3 5])'};
varstrs = {'By';'Jy';'Jx';'Ey';'ni'};
varstrs = {'Ez','By','vx([4 6])'}';
varstrs = {'Bz';'ti';'te';'log10(ni)'};
clims = {1*[-1 1],[0 0.6],[0 0.3],[-4 1.5]};
%varstrs = {'Jx';'Epar';'log10(abs(n([1])))';'log10(abs(n([3 5])))'};
%varstrs = {'Ey';'Ey+vexBy';'Ey+vixBy';'vex';'hvey';'vez';'vix'};
%varstrs = {'Babs','By','log10(ni)','log10(magmom([3 5]))','log10(magmom([4 6]))'}';
%varstrs = {'Ez';'pxy([2 4 6])';'pyz([2 4 6])'};
%varstrs = {'Epar';'ne';'ni';'log10(ni)'};
%clims = {0.5*[-1 1],[0 1],[0 1],[-1.5 1]};
cmapbr = pic_colors('blue_red');
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr};
h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps);

%% Map plot, magnetic curvature
% Plotmap
twpe = 24000;
xlim = [60 110];
zlim = 0.99*[-8 8];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);

varstrs = {'Babs';'curvbabs';'log10(curvbabs)';'log10(curvbrad)'};
clims = {[0 1.2];[0 1];[-2 1];[-2 3]};
varstrs = {'Babs';'log10(curvbabs)'};
clims = {[0 1.2];[-2 1]};

varstrs = {'Ez';'vepar';'log10(ne)';'ti';'log10(curvbabs)';'divpz([3 5])'};
varstrs = {'n([3 5])';'ni';'log10(ni)'};
varstrs = {'log10(ni)';'log10(curvbabs)'};
clims = {[-2 0.2];[-2 2]};


cmapbr = pic_colors('blue_red');
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr};
h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps,'smooth',5);

%% Reduced distributions, make
twpe = 24000; xlim = [50 155]; zlim = [-15 15];

for zpick = [4]
  ds = ds100.twpelim(twpe).zfind(zpick).xlim(xlim).findtag({'line horizontal'});
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
  %fred3_tmp = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred3_z%g = fred3_tmp;',zpick))
  fred5_tmp = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred5_z%g = fred5_tmp;',zpick))
  fred3_tmp = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred3_z%g = fred3_tmp;',zpick))
  fred35_tmp = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred35_z%g = fred35_tmp;',zpick))
  %fred46_tmp = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred46_z%g = fred46_tmp;',zpick))    
end

%% Reduced distributions, plot 1
xlim_fred = [min(fred35_z0.x) max(fred35_z0.x)];
zlim_fred = min(fred35_z0.z) + 0.25*[-1 1];
x_z0 = pic.xlim(xlim_fred).zlim(zlim_fred).xi;
vExBx_z0 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBx,2)); 
vExBy_z0 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBy,2)); 
vExBz_z0 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBz,2)); 

xlim_fred = [min(fred35_z2.x) max(fred35_z2.x)];
zlim_fred = min(fred35_z2.z) + 0.25*[-1 1];
x_z2 = pic.xlim(xlim_fred).zlim(zlim_fred).xi;
vExBx_z2 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBx,2)); 
vExBy_z2 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBy,2)); 
vExBz_z2 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBz,2)); 

xlim_fred = [min(fred35_z4.x) max(fred35_z4.x)];
zlim_fred = min(fred35_z4.z) + 0.25*[-1 1];
x_z4 = pic.xlim(xlim_fred).zlim(zlim_fred).xi;
vExBx_z4 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBx,2)); 
vExBy_z4 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBy,2)); 
vExBz_z4 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBz,2)); 
%%
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 3;
ncols = 1;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 1; colorE = [0 0.8 0.8];
doV = 1; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];

cmap_dist = pic_colors('waterfall');

%freds = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
freds = {fred35_z0,fred35_z0,fred35_z0};
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};

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
    irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
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
    if 0 % doExB
      hold(hca,'on')
      xx = eval(['x_z' num2str(unique(fred.z))]);
      vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
      plot(hca,xx,vv,'color',colorExB,'linewidth',1.5)
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end
  end
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
h(1).CLim = 0.99*[-4 2];
h(1).YLim = 0.99*4*[-1 1];
%%
for ip = 1:nrows*ncols
  axwidth(ip) = h(ip).Position(3);
end
for ip = 1:nrows*ncols
  h(ip).Position(3) = min(axwidth);
end
for ip = 1:nrows*ncols
  h(ip).Position(2) = h(ip).Position(2)-0.05;
end

if 1
  hpos = h(1).Position;
  hcbvx = colorbar('peer',h(1),'location','northoutside','fontsize',14);  
  hcbvx.YLabel.String = 'f_{i,cold}(x,v_{x})';
  h(1).Position = hpos;
  hcbvx.Position(2) = hpos(2) + hpos(4)+0.01;

  hpos = h(2).Position;
  hcbvx = colorbar('peer',h(2),'location','northoutside','fontsize',14);  
  hcbvx.YLabel.String = 'f_{i,cold}(x,v_{y})';
  h(2).Position = hpos;
  hcbvx.Position(2) = hpos(2) + hpos(4)+0.01;

  hpos = h(3).Position;
  hcbvx = colorbar('peer',h(3),'location','northoutside','fontsize',14);  
  hcbvx.YLabel.String = 'f_{i,cold}(x,v_{z})';
  h(3).Position = hpos;
  hcbvx.Position(2) = hpos(2) + hpos(4)+0.01;
end
%c_eval('h(?).YTickLabel = []; h(?).YLabel = [];',[2 3 5 6 8 9])
c_eval('h(?).YTick = -10:1:10;',1:9)

%% Reduced distributions, plot 2
pic = no02m;
xlim_fred = [min(fred35_z0.x) max(fred35_z0.x)];
zlim_fred_z0 = min(fred35_z0.z) + 0.25*[-1 1];
zlim_fred_z2 = min(fred35_z2.z) + 0.25*[-1 1];
x_z0 = pic.xlim(xlim_fred).zlim(zlim_fred).xi;
twpe = 24000;
vExBx_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z0).vExBx,2)); 
vExBy_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z0).vExBy,2)); 
vExBz_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z0).vExBz,2)); 
vicx_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z0).vx([3 5]),2)); 
vicy_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z0).vy([3 5]),2)); 
vicz_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z0).vz([3 5]),2)); 
Bz_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z0).Bz,2)); 
Ey_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z0).Ey,2)); 

vExBx_z2 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z2).vExBx,2)); 
vExBy_z2 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z2).vExBy,2)); 
vExBz_z2 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z2).vExBz,2)); 
vicx_z2 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z2).vx([3 5]),2)); 
vicy_z2 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z2).vy([3 5]),2)); 
vicz_z2 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z2).vz([3 5]),2)); 
Bz_z2 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z2).Bz,2)); 
Ey_z2 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred_z2).Ey,2)); 


% xlim_fred = [min(fred35_z2.x) max(fred35_z2.x)];
% zlim_fred = min(fred35_z2.z) + 0.25*[-1 1];
% x_z2 = pic.xlim(xlim_fred).zlim(zlim_fred).xi;
% vExBx_z2 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBx,2)); 
% vExBy_z2 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBy,2)); 
% vExBz_z2 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBz,2)); 
% 
% xlim_fred = [min(fred35_z4.x) max(fred35_z4.x)];
% zlim_fred = min(fred35_z4.z) + 0.25*[-1 1];
% x_z4 = pic.xlim(xlim_fred).zlim(zlim_fred).xi;
% vExBx_z4 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBx,2)); 
% vExBy_z4 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBy,2)); 
% vExBz_z4 = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBz,2)); 

% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 3;
ncols = 1;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 1; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 0; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];

cmap_dist = pic_colors('waterfall');

%freds = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
freds = {fred35_z2,fred35_z2,fred35_z2};
freds = {fred35_z0,fred35_z0,fred35_z0};
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};

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
    irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 12;
    if 0 %*doE
      hold(hca,'on')
      plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
      hold(hca,'off')
    end
    if 0 %*doV
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
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
h(1).CLim = 0.99*[-4 2];
h(1).YLim = 0.99*4*[-1 1];

%
axwidth = [];
for ip = 1:nrows*ncols
  axwidth(ip) = h(ip).Position(3);
end
for ip = 1:nrows*ncols
  h(ip).Position(3) = min(axwidth);
end
for ip = 1:nrows*ncols
  %h(ip).Position(2) = h(ip).Position(2)-0.05;
  h(ip).FontSize = 14;
end

if 1
  hpos = h(1).Position;
  hcbvx = colorbar('peer',h(1),'fontsize',14);  
  hcbvx.YLabel.String = 'log_{10}f_{i,cold}(x,v_{x})';
  
  hpos = h(2).Position;
  hcbvx = colorbar('peer',h(2),'fontsize',14);  
  hcbvx.YLabel.String = 'log_{10}f_{i,cold}(x,v_{y})';
  
  hpos = h(3).Position;
  hcbvx = colorbar('peer',h(3),'fontsize',14);  
  hcbvx.YLabel.String = 'log_{10}f_{i,cold}(x,v_{z})';
  
end
%c_eval('h(?).YTickLabel = []; h(?).YLabel = [];',[2 3 5 6 8 9])
c_eval('h(?).YTick = -10:1:10;',1:3)

% Add vExB and v velocities
if 0
  np_smooth = 50;
  vi_col = [1 1 1]*0.0;
  hold(h(1),'on')
  h1 = plot(h(1),x_z0,vicx_z0,'color',vi_col,'linestyle','-.');
  h2 = plot(h(1),x_z0,smooth(vExBx_z0,np_smooth),'k');
  hold(h(1),'off')
  
  hold(h(2),'on')
  plot(h(2),x_z0,vicy_z0,'color',vi_col,'linestyle','-.')
  plot(h(2),x_z0,smooth(vExBy_z0,np_smooth),'k')
  hold(h(2),'off')
  
  hold(h(3),'on')
  plot(h(3),x_z0,vicz_z0,'color',vi_col,'linestyle','-.')
  plot(h(3),x_z0,smooth(vExBz_z0,np_smooth),'k')
  hold(h(3),'off')
  
  hold(h(1),'on')
  %plot(h(1),x_z0,Bz_z0,'color',[1 1 0])  
  %plot(h(1),x_z0,Ey_z0,'color',[0 1 1])  
  hold(h(1),'off')
  legend(findobj(h(1),'type','line'),{'v_{ix}','v_{ExB,x}'},'location','best')
end

%% Reduced distributions, fraction
xlim_fred = [min(fred35_z0.x) max(fred35_z0.x)];
zlim_fred = min(fred35_z0.z) + 0.25*[-1 1];
x_z0 = pic.xlim(xlim_fred).zlim(zlim_fred).xi;
twpe = 24000;
vExBx_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred).vExBx,2)); 
vExBy_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred).vExBy,2)); 
vExBz_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred).vExBz,2)); 

vicx_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred).vx([3 5]),2)); 
vicy_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred).vy([3 5]),2)); 
vicz_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred).vz([3 5]),2)); 

Bz_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred).Bz,2)); 
Ey_z0 = squeeze(mean(pic.twpelim(twpe,'exact').xlim(xlim_fred).zlim(zlim_fred).Ey,2)); 

% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 3;
ncols = 1;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 1; colorE = [0 0.8 0.8];
doV = 1; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];

cmap_dist = pic_colors('waterfall');

freds_den = {fred35_z0,fred35_z0,fred35_z0};
freds_num = {fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;
    fred_den = freds_den{ifred};
    fred_num = freds_num{ifred};
    labstr = labstrs{ifred};
    fredplot_den = eval(['fred_den.fv' labstr]);
    fredplot_num = eval(['fred_num.fv' labstr]);
    fredplot = fredplot_num./fredplot_den;
    pcolor(hca,fred.x,fred.v,fredplot')
    shading(hca,'flat')
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
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
    if 0 % doExB
      hold(hca,'on')
      xx = eval(['x_z' num2str(unique(fred.z))]);
      vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
      plot(hca,xx,vv,'color',colorExB,'linewidth',1.5)
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end
  end
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
h(1).CLim = [0 1];
colormap(pic_colors('waterfall'))
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
%
axwidth = [];
for ip = 1:nrows*ncols
  axwidth(ip) = h(ip).Position(3);
end
for ip = 1:nrows*ncols
  h(ip).Position(3) = min(axwidth);
end
for ip = 1:nrows*ncols
  %h(ip).Position(2) = h(ip).Position(2)-0.05;
  h(ip).FontSize = 14;
end

h(1).CLim = [0 1];
h(1).YLim = 0.99*4*[-1 1];

if 1
  hpos = h(1).Position;
  hcbvx = colorbar('peer',h(1),'fontsize',14);  
  hcbvx.YLabel.String = 'f_{i,cold,top}/f_{i,cold}(x,v_{x})';
  
  hpos = h(2).Position;
  hcbvx = colorbar('peer',h(2),'fontsize',14);  
  hcbvx.YLabel.String = 'f_{i,cold,top}/f_{i,cold}(x,v_{y})';
  
  hpos = h(3).Position;
  hcbvx = colorbar('peer',h(3),'fontsize',14);  
  hcbvx.YLabel.String = 'f_{i,cold,top}/f_{i,cold}(x,v_{z})';  
end
%c_eval('h(?).YTickLabel = []; h(?).YLabel = [];',[2 3 5 6 8 9])
c_eval('h(?).YTick = -10:1:10;',1:9)

%% Reduced distributions, locations of distribution boxes
varstrs = {'log10(ni)'};
clims = {0.99*[-1.5 0.5],0.99*0.5*[-1 1]};
cmapbr = pic_colors('blue_red');
cmaps = {cmapbr,cmapbr};
xlim = [50 150];
zlim = [-12 12];
h = setup_subplots(1,1,'vertical');
h1 = no02m.twpelim(24000).xlim(xlim).zlim(zlim).plot_map(h(1),varstrs,'A',1,'clim',clims,'cmap',cmaps);
hold(h,'on')
ds100.twpelim(24000).xlim([50 105]).zfind([0 2 4]).findtag({'line horizontal'}).plot_boxes(h);
hold(h,'off')

%%
varstrs = {'Babs';'curvbabs';'log10(curvbabs)';'log10(curvbrad)'};
clims = {[0 1.2];[0 1];[-2 1];[-2 3]};
varstrs = {'log10(curvbabs)'};
clims = {0.99*[-2 1]};

cmapbr = pic_colors('blue_red');
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr};
xlim = [60 110];
zlim = [-7 7];
h = setup_subplots(1,1,'vertical');
h1 = no02m.twpelim(24000).xlim(xlim).zlim(zlim).plot_map(h(1),varstrs,'A',1,'clim',clims,'cmap',cmaps);
hold(h,'on')
ds100.twpelim(24000).xlim([50 105]).zfind([0 2 4]).findtag({'line horizontal'}).plot_boxes(h);
hold(h,'off')

hcbar = findobj(gcf,'type','Colorbar');
c_eval('hcbar(?).YLabel.FontSize = 14;',1)
c_eval('h(?).YLabel.FontSize = 14;',1)
c_eval('h(?).XLabel.FontSize = 14;',1)
c_eval('h(?).FontSize = 14;',1)

h(1).Position(2) = 0.1;
h(1).Position(4) = 0.8;

hlines = findobj(gcf,'type','line');
c_eval('hlines(?).LineWidth = 1.5;',1:numel(hlines))

%% Reduced distribution, horizontal, with some overview plots
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 5;
ncols = 1;
npanels  = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 1; colorE = [0 0.8 0.8];
doV = 1; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];

cmap_dist = pic_colors('waterfall');

%freds = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
freds = {fred35_z0,fred35_z0,fred35_z0};
freds = {fred3_z4,fred3_z4,fred3_z4};
freds = {fred5_z4,fred5_z4,fred5_z4};
%freds = {fred5_z2,fred5_z2,fred5_z2};
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};

xlim = [64 140];
zlim = [-10 10]*0.99;
pic = no02m.twpelim(24000).xlim(xlim).zlim(zlim);
A = pic.A;
[Ainds,Avals] = saddle(A,'sort');
if 1 % log10(ni)
  hca = h(isub); isub = isub + 1;
  ni = smooth2(pic.ni,3);
  imagesc(hca,pic.xi,pic.zi,log10(ni)')
  colormap(hca,pic_colors('blue_red'))
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.CLim = [-2 0.2];
  hca.YDir = 'normal';
  hcb = colorbar(hca,'peer',hca);
  hcb.YLabel.String = 'log_{10} n_i';
  if 1
    hold(hca,'on')
    contour(hca,pic.xi,pic.zi,A',0:0.5:25,'k')    
    %hsep = contour(hca,pic.xi,pic.zi,A',Avals(1)*[1 1],'color',[0 0 0],'linewidth',1,'linestyle','-');
    hold(hca,'off')
  end
end
if 1 % |K_B|
  hca = h(isub); isub = isub + 1;
  KB = smooth2(pic.curvbabs,3);
  imagesc(hca,pic.xi,pic.zi,log10(KB)')
  colormap(hca,pic_colors('blue_red'))
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.CLim = [-2 1];
  hca.YDir = 'normal';
  hcb = colorbar(hca,'peer',hca);
  hcb.YLabel.String = 'log_{10}|K_B|';
  if 1
    hold(hca,'on')
    contour(hca,pic.xi,pic.zi,A',0:0.5:25,'k')
    hold(hca,'off')
  end
end

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
    irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    hcb = colorbar('peer',hca);  
    hcb.YLabel.String = sprintf('f_{i,cold}(x,v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
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
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks_fred = linkprop(h(3:end),{'YLim','CLim'});
hlinks = linkprop(h,{'XLim'});
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
h(3).CLim = 0.99*[-4 2];
h(3).YLim = 0.99*3.5*[-1 1];

hb = findobj(gcf,'type','colorbar');

for ip = 1:npanels
  h(ip).FontSize = 13;
  h(ip).Position(3) = 0.7;
  hb(ip).FontSize = 13;
end

if 1 % add boxes
  %%
  ds = ds100.twpelim(twpe).zfind(4).xlim(xlim).findtag({'line horizontal'});
  ds.plot_boxes(h(1))
  ds.plot_boxes(h(2))
%   ds = ds100.twpelim(twpe).xfind(75).findtag({'line vertical'});
%   ds.plot_boxes(h(1))
%   ds.plot_boxes(h(2))
%   ds = ds100.twpelim(twpe).xfind(85).findtag({'line vertical'});
%   ds.plot_boxes(h(1))
%   ds.plot_boxes(h(2))
  
end
h(1).XLim = [65 138];

%% Reduced distributions, make vertical
twpe = 24000; xlim = [50 155]; zlim = [-15 15];

for xpick = 75:10:95
  ds = ds100.twpelim(twpe).xfind(xpick).findtag({'line vertical'});
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
  %fred1_tmp = ds.reduce_1d_new('x',[1],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred1_x%g = fred1_tmp;',xpick))  
%   fred3_tmp = ds.reduce_1d_new('x',[3],[]); eval(sprintf('fred3_x%g = fred3_tmp;',xpick))
%   fred5_tmp = ds.reduce_1d_new('x',[5],[]); eval(sprintf('fred5_x%g = fred5_tmp;',xpick))
%   fred35_tmp = ds.reduce_1d_new('x',[3 5],[]); eval(sprintf('fred35_x%g = fred35_tmp;',xpick))
%   fred3_tmp = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred3_x%g = fred3_tmp;',xpick))
%   fred5_tmp = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred5_x%g = fred5_tmp;',xpick))
%   fred35_tmp = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred35_x%g = fred35_tmp;',xpick))
  
  fred46_tmp = ds.reduce_1d_new('x',[4 6],[]); eval(sprintf('fred46_x%g = fred46_tmp;',xpick))    
  %fred46_tmp = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred46_x%g = fred46_tmp;',xpick))    
end

%% Reduced distributions, vertical plot 1
%xpick = ;
for xpick = 75:10:95
  eval(sprintf('zlim_fred = [min(fred35_x%g.z) max(fred35_x%g.z)];',xpick,xpick))
  eval(sprintf('xlim_fred = min(fred35_x%g.x) + 0.25*[-1 1];',xpick,xpick))
  eval(sprintf('z_x%g = pic.xlim(xlim_fred).zlim(zlim_fred).zi;',xpick))
%   eval(sprintf('vExBx_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBx,1));',xpick))
%   eval(sprintf('vExBy_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBy,1));',xpick))
%   eval(sprintf('vExBz_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBz,1));',xpick))
%   eval(sprintf('vix_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vx([3 5]),1));',xpick))
%   eval(sprintf('viy_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vy([3 5]),1));',xpick))
%   eval(sprintf('viz_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vz([3 5]),1));',xpick))
%   eval(sprintf('vex_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vx([4 6]),1));',xpick))
%   eval(sprintf('vey_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vy([4 6]),1));',xpick))
%   eval(sprintf('vez_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vz([4 6]),1));',xpick))
  eval(sprintf('Ex_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Ex,1));',xpick))
  eval(sprintf('Ey_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Ey,1));',xpick))
  eval(sprintf('Ez_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Ez,1));',xpick))
  eval(sprintf('phiz_x%g = -cumtrapz(z_x%g,Ez_x%g);',xpick,xpick,xpick))
end
sep = no02m.twpelim(twpe).separatrix_location;
disp('Done.')
%%
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 3;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1; 
doE = 0; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doVe = 0; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
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
freds = {fred35_x75,fred35_x85,fred35_x95,fred35_x75,fred35_x85,fred35_x95,fred35_x75,fred35_x85,fred35_x95};
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
    irf_legend(hca,{sprintf('x = %g',unique(fred.x))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
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
      hExB = plot(hca,smooth(vv,50),zz,'color',0*colorExB,'linewidth',1,'linestyle','--');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end  
    if doPhi
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
      hphi = plot(hca,vv_phi,zz,'color',0*[1 1 1],'linewidth',1,'linestyle','-');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      end
      hold(hca,'off')
    end   
  end
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
position_hcb_peer = h(end).Position;
hcb = colorbar('peer',h(end));
h(end).Position = position_hcb_peer;
hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})',labstr);
hcb.YLabel.FontSize = 14;
if all([doV doExB])
  legend([hv,hExB],{sprintf('v_{i,%s}',labstr),sprintf('v_{ExB,%s}',labstr)})
elseif doV
  legend([hv],{sprintf('v_{i,%s}',labstr)})
elseif doExB
  legend([hExB],{sprintf('v_{ExB,%s}',labstr)})
end
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
%compact_panels(0.005,0.01)
h(1).CLim = 0.99*[-4 2];
h(1).CLim = [-3.5 1.8];
%h(1).CLim = 0.99*[-6 0]; % electrons
%h(1).XLim = 0.99*3*[-1 1];
for ip = 2:npanels
  h(ip).YTickLabels = '';
  h(ip).YLabel.String = '';  
end
for ip = 1:npanels
  h(ip).FontSize = 14;
  %h(ip).Position(2) = 0.17;
end

%% Comparison of density structure
varstrs = {'log10(ni)','Epar'};
clims = {0.99*[-1.5 0.5],0.99*0.5*[-1 1]};
cmapbr = pic_colors('blue_red');
cmaps = {cmapbr,cmapbr};
xlim = [100 200];
zlim = [-12 12];
h = setup_subplots(2,2,'vertical');
h1 = bs.twpelim(7300).xlim(xlim).zlim(zlim).plot_map(h(1:2),varstrs,'A',1,'clim',clims,'cmap',cmaps);
h2 = no02.twpelim(3000).xlim(xlim).zlim(zlim).plot_map(h(3:4),varstrs,'A',1,'clim',clims,'cmap',cmaps);

hlinks = linkprop([h1 h2],{'XLim','YLim'});
compact_panels(0.01,0.01)
h2(1).YTick = [];
h2(2).YTick = [];
h2(1).YLabel = [];
h2(2).YLabel = [];
c_eval('h(?).FontSize = 12;',1:4)

hcbar = findobj(gcf,'type','Colorbar');
c_eval('hcbar(?).FontSize = 12;',1:4)
delete(hcbar([2 4]))
hcbar(1).Position(1) = hcbar(3).Position(1);

%% Integrated density as a function of x, comparison
varstrs = {{'ti'}}';

xlim = [100 200];
zlim = [-5 5];
h = setup_subplots(1,1,'vertical');
h1 = bs.twpelim(7300).xlim(xlim).zlim(zlim).plot_line(h(1),'x',varstrs);
hold(h,'on')
h2 = no02.twpelim(3000).xlim(xlim).zlim(zlim).plot_line(h(1),'x',varstrs);
hold(h,'off')

%% Fermi acceleration (slingshot) comparison
%input_varstr = 'B.z';
%input_ivar = find(cellfun(@(x)strcmp(x,input_varstr),varstrs_ts_line_x));
%input_data_for_velocity = squeeze(cell_ts_line_x{input_ivar}(:,11,:));
xlim = [0 200];
zlim = [-2 2];
zpicklim = [3 4];
pic = no02.xlim(xlim);
input_data_for_velocity = squeeze(mean(pic.zlim(zlim).Bz,2));
x = pic.xi;
times = pic.twci;
ntimes = pic.nt;

%%
[front_velocity,front_location,value_at_location] = expansion_velocity(x,times,input_data_for_velocity,'pos',1);
n_vel = numel(front_velocity);

% If a particle was reflected at the front at time = t1, where would it be
% at a alter time t2? First we get the velocity as a function of the front
% speed: v2 = -v1+2*v_front, we assume v1 = 0
velocity_before_reflection = 0.5;
velocity_after_reflection = cellfun(@(x)2*x+velocity_before_reflection,front_velocity,'UniformOutput',false);

% n_particles = ntimes;
% particle_velocity = zeros(ntimes,n_particles);
% for itime = 1:ntimes
%   particle_velocity(itime,itime:end) = velocity_after_reflection{1}(itime);
%   particle_location_{iparticle}(1) = front_location{1}(iparticle);
% end

% reflection location is front location
% particle velocity = front_location+(t-t_reflection)*velocity_after_reflection
% calculate this for each time step

% set location prior to/at the time of the reflection = location of front
particle_location = cell(ntimes,1);
n_particles = ntimes;
for iparticle = 1:n_particles
  particle_location{iparticle}(1) = NaN;%front_location{1}(iparticle);
end
% Integrate particle position 
dt = times(2)-times(1);
for itime = 1:ntimes
  %fprintf('itime = %g, time: %g \n',itime,times(itime))
  for iparticle = 1:n_particles
   % fprintf('Particle: %g \n',iparticle)
    %fprintf('Particle location: %g, Front location: %g \n',particle_location{iparticle}(itime),front_location{1}(itime))
    if iparticle > itime % particle has not yet encountered the front
      particle_location{iparticle}(itime) = NaN;
    elseif iparticle == itime
      particle_location{iparticle}(itime) = front_location{1}(itime);   
    else % particle has encountered the front, advance position            
      particle_location{iparticle}(itime) = particle_location{iparticle}(itime-1) + dt*velocity_after_reflection{1}(iparticle);
    end
  end
end

% Plot
if exist('h','var'); delete(h); end
nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);

isub = 1;
if 0 % Data used to get the expansion velocity
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,times,x,input_data_for_velocity);
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = input_varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  for i_vel = 1:n_vel
    hold(hca,'on')
    plot(hca,times,front_location{i_vel},'k')
    hold(hca,'off')
  end  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Location of expansion front
  hca = h(isub); isub = isub + 1;
  prefix = {'1','1'};
  operator = cellfun(@(x)eval(x),prefix);
  for i_vel = 1:n_vel
    %velocity_tmp = eval(sprintf('%s,%s;',velocity{i_vel}));
    if i_vel == 1, hold(hca,'on'); end
    plot(hca,times,operator(i_vel)*front_location{i_vel});
    if i_vel == n_vel, hold(hca,'off'); end
  end
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'front position';
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
end
if 0 % Velocity of expansion
  hca = h(isub); isub = isub + 1;
  prefix = {'1','1'};
  operator = cellfun(@(x)eval(x),prefix);
  for i_vel = 1:n_vel
    %velocity_tmp = eval(sprintf('%s,%s;',velocity{i_vel}));
    if i_vel == 1, hold(hca,'on'); end
    plot(hca,times,operator(i_vel)*front_velocity{i_vel});
    if i_vel == n_vel, hold(hca,'off'); end
  end
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'velocity';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Velocity of particle after having encountered the front at the given time
  hca = h(isub); isub = isub + 1;
  prefix = {'1','1'};
  operator = cellfun(@(x)eval(x),prefix);
  for i_vel = 1:n_vel
    %velocity_tmp = eval(sprintf('%s,%s;',velocity{i_vel}));
    if i_vel == 1, hold(hca,'on'); end
    plot(hca,times,operator(i_vel)*velocity_after_reflection{i_vel});
    if i_vel == n_vel, hold(hca,'off'); end
  end
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = {'particle velocity','after reflection at front'};
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  irf_legend(hca,{'velocity after reflection = - velocity before reflection + 2*vfront';...
    sprintf('velocity before reflection = %g',velocity_before_reflection)},[0.02 0.98],'k')
end
if 1 % Velocity of front and particle after having encountered the front at a given time
  hca = h(isub); isub = isub + 1;
  prefix = {'1','1'};
  operator = cellfun(@(x)eval(x),prefix);
  for i_vel = 1:n_vel
    %velocity_tmp = eval(sprintf('%s,%s;',velocity{i_vel}));
    if i_vel == 1, hold(hca,'on'); end
    plot(hca,times,operator(i_vel)*front_velocity{i_vel});
    plot(hca,times,operator(i_vel)*velocity_after_reflection{i_vel});
    if i_vel == n_vel, hold(hca,'off'); end
  end
  
%   for i_vel = 1:n_vel
%     %velocity_tmp = eval(sprintf('%s,%s;',velocity{i_vel}));
%     if i_vel == 1, hold(hca,'on'); end
%     plot(hca,times,operator(i_vel)*front_velocity{i_vel});
%     if i_vel == n_vel, hold(hca,'off'); end
%   end
  legend(hca,{'velocity of front','particle velocity after reflection at front'},'box','off','location','northwest')
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = {'particle velocity','after reflection at front'};
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hleg = irf_legend(hca,{'velocity after reflection = - velocity before reflection + 2*vfront';...
      sprintf('velocity before reflection = %g',velocity_before_reflection)},[0.1 0.5]);
  arrayfun(@(x)eval(sprintf('x.Color = ''k'';'),x),hleg)
end

if 0 % Location of particle
  hca = h(isub); isub = isub + 1;
  hold(hca,'on');
  for i_particle = 1:n_particles    
    plot(hca,times,particle_location{i_particle});
  end
  hold(hca,'off');
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'front position';
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
end
  zzz = 2;
if 1 % Plot particle lines starting from front
  hca = h(isub); isub = isub + 1;
  %zzind = find_closest_ind(zpicks,zzz);
  varstr = 'E.z';
  plot_data = squeeze(mean(pic.zlim(zpicklim).Ez,2));
  himag = imagesc(hca,times,x,plot_data');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  for i_vel = 1:n_vel
    hold(hca,'on')
    hfront = plot(hca,times,front_location{i_vel},'color','k','linewidth',1);
    hold(hca,'off')
  end  
%  irf_legend(hca,{sprintf('z = %g',zpicks(zzind))},[0.02 0.98],'k')
  
  hold(hca,'on');  
  for i_particle = 1:n_particles
    plot(hca,times,particle_location{i_particle},'k','linewidth',0.5);
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.YLim(1) = 0;
  %hca.YLim(2) = 100;
end
if 1 % Plot particle lines starting from front
  hca = h(isub); isub = isub + 1;
  zzind = find_closest_ind(zpicks,zzz);
  varstr = 'vi2.x';
  plot_data = squeeze(cell_ts_line_x{find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_line_x))}(:,zzind,:));
  himag = imagesc(hca,times,x,plot_data);
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  for i_vel = 1:n_vel
    hold(hca,'on')
    hfront = plot(hca,times,front_location{i_vel},'color','k','linewidth',1);
    hold(hca,'off')
  end  
  irf_legend(hca,{sprintf('z = %g',zpicks(zzind))},[0.02 0.98],'k')
  
  hold(hca,'on');  
  for i_particle = 21:2:n_particles
    plot(hca,times,particle_location{i_particle},'k','linewidth',0.5);
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.YLim(1) = 0;
  %hca.YLim(2) = 100;
end

arrayfun(@(x)eval(sprintf('x.Position(3) = 0.7;'),x),h)
arrayfun(@(x)eval(sprintf('x.XLabel.String = [];'),x),h(1:end-1))
arrayfun(@(x)eval(sprintf('x.Box = ''on'';'),x),h)
hlink = linkprop(h,{'XLim'});
hlink.Targets(1).XLim(1) = 80;
hlink.Targets(1).XLim(2) = times(end);
arrayfun(@(x)eval(sprintf('x.YLim = [0 200];'),x),h(2:3))
compact_panels

%% Particle trajectories
%% Get r0,v0 from peaks of distributions

it = 2;
if 1
  iSpecies = [3];
  %ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([-0.1 0.6]);
%  ds = ds04(it).xlim([190 203]).zlim([-0.1 0.1]); % top row.
else % electrons
  iSpecies = [4];
  %ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([-0.1 0.6]);
 % ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([0.4 0.6]); % top row.

end
% ds = ds04.twcilim(120).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([180:1:210]);
% ds = ds04.twcilim(140).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([177:1:210]);
% ds = ds04.twcilim(160).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([166:1:205]);
% ds = ds04.twcilim(160).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([167:0.2:175]);
% ds = ds04.twcilim(160).zlim(2+[-0.2 0.2]).dxlim([0 0.25]).xfind([165:0.2:175]);
ds = ds100.twpelim(24000).zfind(0).findtag({'line horizontal'}).xfind([82:1:92]);

nPeaks = 1;
spacingPeaks = 6; % for ions its 2 equals 0.2 vA
fpeaks = ds.get_peaks(nPeaks,spacingPeaks,iSpecies,'vz',[-1.0 0.0]); % ,'vz',[-0.19 0.19]
%fpeaks = ds.get_peaks(nPeaks,spacingPeaks,iSpecies); % ,'vz',[-0.19 0.19]

nDists = ds.nd;
doPlot = 1;
doPrint = 0;
if doPlot
  % plot results
  for id = ds.nd{1}:-1:1
    f = ds.f(1,id,iSpecies);
    
    figure(27)
    h = setup_subplots(1,3);
    
    
    hca = h(1);
    imagesc(hca,f.v,f.v,log10(f.fxy)')
    hca.YDir = 'normal';
    colormap(pic_colors('candy'))
    hold(hca,'on')
    plot(hca,[fpeaks(:,id).vx],[fpeaks(:,id).vy],'k.')
    for iPeak = 1:nPeaks
      text(hca,[fpeaks(iPeak,id).vx],[fpeaks(iPeak,id).vy],sprintf('%g',iPeak))
    end
    hold(hca,'off')
    hca.XGrid = 'on'; hca.YGrid = 'on';
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_y';
    hca.Title.String = sprintf('x = [%.1f,%.1f], z = [%.1f,%.1f]',ds.xi1{1}(id),ds.xi2{1}(id),ds.zi1{1}(id),ds.zi2{1}(id));
    

    hca = h(2);
    imagesc(hca,f.v,f.v,log10(f.fxz)')
    hca.YDir = 'normal';
    colormap(pic_colors('candy'))
    hold(hca,'on')
    plot(hca,[fpeaks(:,id).vx],[fpeaks(:,id).vz],'k.')
    for iPeak = 1:nPeaks
      text(hca,[fpeaks(iPeak,id).vx],[fpeaks(iPeak,id).vz],sprintf('%g',iPeak))
    end
    hold(hca,'off')
    hca.XGrid = 'on'; hca.YGrid = 'on';
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_z';
    

    hca = h(3);
    imagesc(hca,f.v,f.v,log10(f.fyz)')
    hca.YDir = 'normal';
    colormap(pic_colors('candy'))
    hold(hca,'on')
    plot(hca,[fpeaks(:,id).vy],[fpeaks(:,id).vz],'k.')
    for iPeak = 1:nPeaks
      text(hca,[fpeaks(iPeak,id).vy],[fpeaks(iPeak,id).vz],sprintf('%g',iPeak))
    end
    hold(hca,'off')
    hca.XGrid = 'on'; hca.YGrid = 'on';
    hca.XLabel.String = 'v_y';
    hca.YLabel.String = 'v_z';
    
    
    for ip = 1:3
      axis(h(ip),'square')
      h(ip).XTick = h(ip).YTick;
      h(ip).FontSize = 14;
       h(ip).Position(2) = 0.2;
      h(ip).Position(4) = 0.7;
      h(ip).XLim = [-2.5 2];
      h(ip).YLim = [-2 2];
    end
    hb = colorbar('peer',h(3));
    hb.Position(1) = h(3).Position(1) + h(3).Position(3) + 0.01;
    if doPrint
      cn.print(sprintf('fpeaks_z=0_id=%04.0f',id))
    end
    pause
  end
end

%% Integrate trajectories based on fpeaks (which may be based or moments)
pic = no02m.twpelim(0:200:25000,'exact');
t0 = 120-0.01; 
tspan = [t0,70];
t0 = 80;
tspan = [80,90];
%t0 = 60;
%tspan = [t0,239];
m = 1; 
q = 1;
istart = 1;
%ntr_pre = trs.ntr-(istart-1); % the +(istart-1) is there becuse i already did the 1 and added it to the datafile
%ntr_pre = 80; % the +(istart-1) is there becuse i already did the 1 and added it to the datafile
ntr_pre = 1;

for iTr = istart%:numel(fpeaks)
  tic
  fprintf('iTr = %g/%g\n',iTr,numel(fpeaks))
  r0 = [fpeaks(iTr).x, fpeaks(iTr).y, fpeaks(iTr).z];
  v0 = [fpeaks(iTr).vx, fpeaks(iTr).vy, fpeaks(iTr).vz];
  %tr_tmp = no02m.twpelim(0:200:25000,'exact').integrate_trajectory(r0,v0,tspan,m,q);
  tr_tmp = no02m.integrate_trajectory(r0,v0,tspan,m,q);
  
  hca = subplot(3,1,1);
  plot(hca,tr_tmp.x,tr_tmp.z,tr_tmp.x0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca = subplot(3,1,2);
  plot(hca,tr_tmp.y,tr_tmp.z,tr_tmp.y0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  hca = subplot(3,1,3);
  plot(hca,tr_tmp.x,tr_tmp.y,tr_tmp.x0,tr_tmp.y0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  drawnow
%   [Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(tr_tmp.x,tr_tmp.z,tr_tmp.t);  % interpolate
% 
%   tr_tmp.t0 = t0;
%   tr_tmp.Ex = Ex;
%   tr_tmp.Ey = Ey;
%   tr_tmp.Ez = Ez;
%   tr_tmp.Bx = Bx;
%   tr_tmp.By = By;
%   tr_tmp.Bz = Bz;
  %0
  %h5write_trajs('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories.h5',tr_tmp,'fpeaks',fpeaks(iTr),'id',ntr_pre+iTr)
  h5write_trajs('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories_test.h5',tr_tmp,'fpeaks',fpeaks(iTr),'id',ntr_pre+iTr)
  %h5write_trajs('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr_tmp)
  %tr(iPeak,id) = tr_tmp;
  toc
  %catch
  %  continue
  %end
end

%% Movie with trajectories
tr = tr100;
trs = tr.find([tr.z0]==4);
trs = tr100(2:10);
trs = tr100.find([tr100.ncross]>5,[tr100.Ustart]<1,tr100.xstart>95);

twpe = [18000 24000];
xlim = [70 110];
zlim = [-4 4];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'A'}';
clims = {[-1 1];[-1 1]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmapc2 = pic_colors('candy2');
cmapma = pic_colors('matlab');
cmaps = {cmapbr,cmapbr};
cmaps = {cmapwa,cmapbr};

clims = {[2 8]};
filename = [printpath 'no02m_Ay_3traj_pm4'];
%pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename,'tr',trs,'trajcolordot',[trs.x0]);
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename,'tr',trs);

%% Movie, overview, introductory
twpe = 7000:1000:12000;
%twpe = 24000;
twpe = 17000:100:24000;
%twpe = [18000 20000];
xlim = no02m.xi([1 end])+[40 -40]';
zlim = [-8 8];
pic = no02m.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');

cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
%cmapjet = colormap('jet');


varstrs = {'ni','log10(ni)'}';
varstrs = {'jepar','vepar','log10(ne)','Eparx','Epary','Eparz'}';
varstrs = {'ni','ti'}';
clims = {[0 1],[0 0.6],[-2 0.5],[-1 1],[-1 1],[-1 1]};
varstrs = {'n([3 5])./n([1 3 5])'}';
clims = {[0 1]};
% varstrs = {'n([3 5])','t([3 5])'}';
% clims = {[0 0.5],[0 0.5]};
cmaps = {flipdim(cmapbr,1),cmapbr,cmapbr,cmapbr,cmapbr,cmapbr};

filename = [printpath 'no02m_ncold_blue'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename);

%% Reconnection rate vs alfven speed
xlim = no02m.xi([1 end])+[40 -40]';
zlim = 0.5*[-1 1];
pic = no02m.xlim(xlim).zlim(zlim);
Bxline_z1 = interp(no02m,no02m.x_xline,no02m.z_xline+1,no02m.twci,'Bx');
nxline_z1 = interp(no02m,no02m.x_xline,no02m.z_xline+1,no02m.twci,'ne');
nxline_z0 = interp(no02m,no02m.x_xline,no02m.z_xline+0,no02m.twci,'ni');
ncold_z0_ = pic.n([3 5]); ncold_z0 = squeeze(mean(ncold_z0_,2));
A_z0_ = pic.A; A_z0 = squeeze(mean(A_z0_,2));
Bz_z0_ = pic.Bz; Bz_z0 = squeeze(mean(Bz_z0_,2));

%% Figure, nr 1
ylim_R = [0 0.15];
ytick_R = [0 0.05 0.1 0.15];
pic = no02m.xlim(xlim).zlim(zlim);
xlim_plot = no02m.xi([1 end])+[65 -65]';
ylim_n = [0 1];

nrows = 4;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % ER, B/sqrt(n)
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,pic.twci,Bxline_z1./sqrt(nxline_z1),pic.twci,pic.RA);
  
  AX(1).YLabel.String = 'v_A(z=1) (v_{A0})';
  AX(1).XLabel.String = 't (\omega_{ci}^{-1})';
  AX(1).YLim = [0.0 1.5];
  AX(1).YTick = [0 0.5 1];
  AX(2).YLabel.String = 'E_R (v_{A0}B_0)';
  AX(2).YLim = ylim_R;
  AX(2).YTick = ytick_R;
  legend(hca,{'v_A(x_x,z_x+1)','E_R'},'location','northwest')
  %legend(hca,'n(x_x,z_x)','n(x_x,z_x+1)')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % Bx, n
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.twci,Bxline_z1,pic.twci,sqrt(nxline_z1));
  hca.YLabel.String = 'n^{1/2} (n_0), B_x (B_0)';
  hca.XLabel.String = 't (\omega_{ci}^{-1})';  
  legend(hca,'B_x(x_x,z_x+1)','n^{1/2}(x_x,z_x+1)')
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
end
if 0 % ER, n
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,pic.twci,([nxline_z1])',pic.twci,pic.RA);
  AX(1).YLabel.String = 'n (n_0)';
  AX(1).XLabel.String = 't (\omega_{ci}^{-1})';
  AX(1).YLim = ylim_n;
  AX(2).YLabel.String = 'E_R (v_{A0}B_0)';
  AX(2).YLim = ylim_R;
  AX(2).YTick = ytick_R;
  legend(hca,'n(x_x,z_x)','n(x_x,z_x+1)')
end
if 1 % Bz(x,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.twci,pic.xi,Bz_z0);
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'B_z';
  hca.YLabel.String = 'x (d_i)';
  hca.XLabel.String = 't (\omega_{ci}^{-1})';  
  hold(hca,'on')
  clim = hca.CLim;
  contour(hca,pic.twci,pic.xi,A_z0,[0:1:25],'k');
  hca.CLim = clim;
  hold(hca,'on')
  hca.YLim = xlim_plot;
  hca.XGrid = 'on';
  hca.YGrid = 'on';    
  hca.Layer = 'top';  
  hca.CLim = [-1 1];
end
if 1 % ncold(x,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.twci,pic.xi,ncold_z0);
  shading(hca,'flat')
  colormap(hca,pic_colors('candy4'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'n_{i,cold}';
  hca.YLabel.String = 'x (d_i)';
  hca.XLabel.String = 't (\omega_{ci}^{-1})';  
  hold(hca,'on')
  clim = hca.CLim;
  contour(hca,pic.twci,pic.xi,A_z0,[0:1:25],'k');
  hca.CLim = clim;
  hold(hca,'on')
  hca.YLim = xlim_plot;
  hca.XGrid = 'on';
  hca.YGrid = 'on';    
  hca.Layer = 'top';  
  hca.CLim = [0 0.8];
end

if 0 % Er, n
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,df04.twci,[df04_nxline_z0,df04_nxline_z1]',df04.twci,df04.RA);
  AX(1).YLabel.String = 'n (n_0)';
  AX(1).XLabel.String = 't (\omega_{ci}^{-1})';
  AX(1).YLim = ylim_n;
  AX(2).YLabel.String = 'E_R (v_{A0}B_0)';
  AX(2).YLim = ylim_R;
  AX(2).YTick = ytick_R;
  legend(hca,'n(x_x,z_x)','n(x_x,z_x+1)')
end

fig = gcf;
hlinks = linkprop(findobj(gcf,'type','axes'),{'XLim'});
compact_panels(0.01)
fontsize = 12;

for ip = 1:npanels
  h(ip).Position(3) = 0.7;  
  h(ip).FontSize = fontsize;
end
hb = findobj(gcf,'type','colorbar');
for ip = 1:numel(hb)
  hb(ip).Position(1) = 0.85;  
  hb(ip).FontSize = fontsize;
end
%
if 0 % Rotate figure so that the text is ok for it to be upright
  h(1).XTick = [0:20:120];
  h(1).XTickLabel = {'0','20','40','60','80','100','120'};
  h(1).XLabel = h(4).XLabel;
  h(1).XAxisLocation = 'top';
  h(1).XTickLabelRotation = 270;
  h(1).XLabel.String = 't (\omega_{ci}^{-1})';
  h(1).XLabel.Rotation = 270;
  h(1).XLabel.HorizontalAlignment = 'right';
  
  h(4).XTickLabel = [];
  
  
  %h(ip).YLabel.Rotation = 270;
  %h(ip).YAxisLocation = 'top';
  
  AX(2).YLabel.VerticalAlignment = 'bottom';
  AX(2).YLabel.Rotation = 270;
    AX(2).YTickLabelRotation= 270;
    AX(2).YDir = 'reverse';
  
  for ip = 1:npanels
    h(ip).YLabel.VerticalAlignment = 'top';
    h(ip).YLabel.Rotation = 270;
    h(ip).YTickLabelRotation= 270;
    h(ip).YDir = 'reverse';
  end
  for ip = 1:numel(hb)
    hb(ip).YLabel.VerticalAlignment = 'bottom';
    hb(ip).YLabel.Rotation = 270;
    %hb(ip).YTickLabelRotation = 270;
  end
  
  %hleg = findobj(gcf,'type','legend');
  
  for ip = 1:npanels
    h(ip).Position(2) = h(ip).Position(2)-0.05;    
  end
  for ip = 1:numel(hb)
    hb(ip).Position(2) = hb(ip).Position(2)-0.05;    
  end
end

%irf_plot_axis_align

%% Figure, nr 2
ylim_R = [0 0.15];
ytick_R = [0 0.05 0.1 0.15];
pic = no02m.xlim(xlim).zlim(zlim);
xlim_plot = no02m.xi([1 end])+[65 -65]';
ylim_n = [0 1];

nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % ER
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.twci,pic.RA,'linewidth',1);
  
  hca.XLabel.String = 't (\omega_{ci}^{-1})';
  %AX(1).YLim = [0.0 1.5];
  %AX(1).YTick = [0 0.5 1];
  hca.YLabel.String = 'E_R (v_{A0}B_0)';
  %AX(2).YLim = ylim_R;
  %AX(2).YTick = ytick_R;
  %legend(hca,{'v_A(x_x,z_x+1)','E_R'},'location','northwest')
  %legend(hca,'n(x_x,z_x)','n(x_x,z_x+1)')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % ER, B/sqrt(n)
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,pic.twci,Bxline_z1./sqrt(nxline_z1),pic.twci,pic.RA);
  
  AX(1).YLabel.String = 'v_A(z=1) (v_{A0})';
  AX(1).XLabel.String = 't (\omega_{ci}^{-1})';
  AX(1).YLim = [0.0 1.5];
  AX(1).YTick = [0 0.5 1];
  AX(2).YLabel.String = 'E_R (v_{A0}B_0)';
  AX(2).YLim = ylim_R;
  AX(2).YTick = ytick_R;
  legend(hca,{'v_A(x_x,z_x+1)','E_R'},'location','northwest')
  %legend(hca,'n(x_x,z_x)','n(x_x,z_x+1)')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Bx, n
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.twci,Bxline_z1,pic.twci,sqrt(nxline_z1));
  hca.YLabel.String = 'n^{1/2} (n_0), B_x (B_0)';
  hca.XLabel.String = 't (\omega_{ci}^{-1})';  
  legend(hca,'B_x(x_x,z_x+1)','n^{1/2}(x_x,z_x+1)')
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
end
if 0 % ER, n
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,pic.twci,([nxline_z1])',pic.twci,pic.RA);
  AX(1).YLabel.String = 'n (n_0)';
  AX(1).XLabel.String = 't (\omega_{ci}^{-1})';
  AX(1).YLim = ylim_n;
  AX(2).YLabel.String = 'E_R (v_{A0}B_0)';
  AX(2).YLim = ylim_R;
  AX(2).YTick = ytick_R;
  legend(hca,'n(x_x,z_x)','n(x_x,z_x+1)')
end
if 1 % Bz(x,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.twci,pic.xi,Bz_z0);
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'B_z';
  hca.YLabel.String = 'x (d_i)';
  hca.XLabel.String = 't (\omega_{ci}^{-1})';  
  hold(hca,'on')
  clim = hca.CLim;
  contour(hca,pic.twci,pic.xi,A_z0,[0:1:25],'k');
  hca.CLim = clim;
  hold(hca,'on')
  hca.YLim = xlim_plot;
  hca.XGrid = 'on';
  hca.YGrid = 'on';    
  hca.Layer = 'top';  
  hca.CLim = [-1 1];
  
  if 1
    hold(hca,'on')
    plot(hca,pic.twci,pic.x_xline,'k','linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % ncold(x,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.twci,pic.xi,ncold_z0);
  shading(hca,'flat')
  colormap(hca,pic_colors('candy4'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'n_{i,cold}';
  hca.YLabel.String = 'x (d_i)';
  hca.XLabel.String = 't (\omega_{ci}^{-1})';  
  hold(hca,'on')
  clim = hca.CLim;
  contour(hca,pic.twci,pic.xi,A_z0,[0:1:25],'k');
  hca.CLim = clim;
  hold(hca,'on')
  hca.YLim = xlim_plot;
  hca.XGrid = 'on';
  hca.YGrid = 'on';    
  hca.Layer = 'top';  
  hca.CLim = [0 0.8];
  
  if 1
    hold(hca,'on')
    plot(hca,pic.twci,pic.x_xline,'k','linewidth',1.5)
    hold(hca,'off')
  end
end

if 0 % Er, n
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,df04.twci,[df04_nxline_z0,df04_nxline_z1]',df04.twci,df04.RA);
  AX(1).YLabel.String = 'n (n_0)';
  AX(1).XLabel.String = 't (\omega_{ci}^{-1})';
  AX(1).YLim = ylim_n;
  AX(2).YLabel.String = 'E_R (v_{A0}B_0)';
  AX(2).YLim = ylim_R;
  AX(2).YTick = ytick_R;
  legend(hca,'n(x_x,z_x)','n(x_x,z_x+1)')
end

fig = gcf;
hlinks = linkprop(findobj(gcf,'type','axes'),{'XLim'});
compact_panels(0.01)
fontsize = 12;

for ip = 1:npanels
  h(ip).Position(3) = 0.7;  
  h(ip).FontSize = fontsize;
end
hb = findobj(gcf,'type','colorbar');
for ip = 1:numel(hb)
  hb(ip).Position(1) = 0.85;  
  hb(ip).FontSize = fontsize;
end
%
if 0 % Rotate figure so that the text is ok for it to be upright
  h(1).XTick = [0:20:120];
  h(1).XTickLabel = {'0','20','40','60','80','100','120'};
  h(1).XLabel = h(4).XLabel;
  h(1).XAxisLocation = 'top';
  h(1).XTickLabelRotation = 270;
  h(1).XLabel.String = 't (\omega_{ci}^{-1})';
  h(1).XLabel.Rotation = 270;
  h(1).XLabel.HorizontalAlignment = 'right';
  
  h(4).XTickLabel = [];
  
  
  %h(ip).YLabel.Rotation = 270;
  %h(ip).YAxisLocation = 'top';
  
  AX(2).YLabel.VerticalAlignment = 'bottom';
  AX(2).YLabel.Rotation = 270;
    AX(2).YTickLabelRotation= 270;
    AX(2).YDir = 'reverse';
  
  for ip = 1:npanels
    h(ip).YLabel.VerticalAlignment = 'top';
    h(ip).YLabel.Rotation = 270;
    h(ip).YTickLabelRotation= 270;
    h(ip).YDir = 'reverse';
  end
  for ip = 1:numel(hb)
    hb(ip).YLabel.VerticalAlignment = 'bottom';
    hb(ip).YLabel.Rotation = 270;
    %hb(ip).YTickLabelRotation = 270;
  end
  
  %hleg = findobj(gcf,'type','legend');
  
  for ip = 1:npanels
    h(ip).Position(2) = h(ip).Position(2)-0.05;    
  end
  for ip = 1:numel(hb)
    hb(ip).Position(2) = hb(ip).Position(2)-0.05;    
  end
end

%irf_plot_axis_align
%% Figure, nr 2, vertical

ylim_R = [0 0.15];
ytick_R = [0 0.05 0.1 0.15];
pic = no02m.xlim(xlim).zlim(zlim);
xlim_plot = no02m.xi([1 end])+[65 -65]';
ylim_n = [0 1];

nrows = 1;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % ER
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.RA,pic.twci,'linewidth',1);  
  hca.YLabel.String = 't (\omega_{ci}^{-1})';
  hca.XLabel.String = 'E_R (v_{A0}B_0)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % Bz(x,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.twci,Bz_z0');
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hcb = colorbar('peer',hca,'Location','northoutside');
  hcb.YLabel.String = 'B_z';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';  
  hold(hca,'on')
  clim = hca.CLim;
  contour(hca,pic.xi,pic.twci,A_z0',[0:1:25],'k');
  hca.CLim = clim;
  hold(hca,'on')
  hca.XLim = xlim_plot;
  hca.XGrid = 'on';
  hca.YGrid = 'on';    
  hca.Layer = 'top';  
  hca.CLim = 0.99*[-1 1];
  
  if 1
    hold(hca,'on')
    plot(hca,pic.x_xline,pic.twci,'k','linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % ncold(x,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.twci,ncold_z0');
  shading(hca,'flat')
  colormap(hca,pic_colors('candy4'))
  hcb = colorbar('peer',hca,'Location','northoutside');
  hcb.YLabel.String = 'n_{i,cold}';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';  
  hold(hca,'on')
  clim = hca.CLim;
  contour(hca,pic.xi,pic.twci,A_z0',[0:1:25],'k');
  hca.CLim = clim;
  hold(hca,'on')
  hca.XLim = xlim_plot;
  hca.XGrid = 'on';
  hca.YGrid = 'on';    
  hca.Layer = 'top';  
  hca.CLim = [0.001 0.8]*0.99;
  
  if 1
    hold(hca,'on')
    plot(hca,pic.x_xline,pic.twci,'k','linewidth',1.5)
    hold(hca,'off')
  end
end


fig = gcf;
hlinks = linkprop(findobj(gcf,'type','axes'),{'YLim'});
%compact_panels(0.01)
fontsize = 14;

for ip = 1:npanels
  h(ip).Position(2) = 0.15;
  h(ip).Position(3) = 0.275;
  h(ip).Position(4) = 0.63;
  h(ip).FontSize = fontsize;
end
for ip = 2:npanels
  h(ip).YTickLabels = [];
  h(ip).YLabel.String = '';
end

hb = findobj(gcf,'type','colorbar');
for ip = 1:numel(hb)
  %hb(ip).Position(1) = 0.85;  
  hb(ip).Location = 'northoutside';  
  hb(ip).FontSize = fontsize;
end

%% Bz(x,t), nicol(x,t)
xlim = no02m.xi([1 end])+[40 -40]';
zlim = 0.5*[-1 1];
pic = no02m.xlim(xlim).zlim(zlim);
ncold_z0_ = pic.n([3 5]); ncold_z0 = squeeze(mean(ncold_z0_,2));
tic; A_z0_ = pic.A; A_z0 = squeeze(mean(A_z0_,2)); toc
Bz_z0_ = pic.Bz; Bz_z0 = squeeze(mean(Bz_z0_,2));

%% Figure
ylim_R = [0 0.15];
ytick_R = [0 0.05 0.1 0.15];
pic = no02m.xlim(xlim).zlim(zlim);
xlim_plot = no02m.xi([1 end])+[65 -65]';
ylim_n = [0 1];

nrows = 1;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;


if 1 % Bz(x,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.twci,Bz_z0');
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'B_z';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';  
  hold(hca,'on')
  clim = hca.CLim;
  contour(hca,pic.xi,pic.twci,A_z0',[0:1:25],'k');
  hca.CLim = clim;
  hold(hca,'on')
  hca.XLim = xlim_plot;
  hca.XGrid = 'on';
  hca.YGrid = 'on';    
  hca.Layer = 'top';  
  hca.CLim = [-1 1];
  if 1
    hold(hca,'on')
    plot(hca,pic.x_xline,pic.twci,'k','linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % ncold(x,t)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,pic.twci,ncold_z0');
  shading(hca,'flat')
  colormap(hca,pic_colors('candy4'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'n_{i,cold}';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 't (\omega_{ci}^{-1})';  
  hold(hca,'on')
  clim = hca.CLim;
  contour(hca,pic.xi,pic.twci,A_z0',[0:1:25],'k');
  hca.CLim = clim;
  hold(hca,'on')
  hca.XLim = xlim_plot;
  hca.XGrid = 'on';
  hca.YGrid = 'on';    
  hca.Layer = 'top';  
  hca.CLim = [0 0.8];
  if 1
    hold(hca,'on')
    plot(hca,pic.x_xline,pic.twci,'k','linewidth',1.5)
    hold(hca,'off')
  end
end

fig = gcf;
hlinks = linkprop(findobj(gcf,'type','axes'),{'XLim'});
compact_panels(0.01,0.01)
fontsize = 16;
h(2).YTickLabels = [];
h(2).YLabel.String = '';

for ip = 1:npanels
  h(ip).Position(2) = 0.15;  
  h(ip).Position(4) = 0.65;  
  h(ip).FontSize = fontsize;
end
hb = findobj(gcf,'type','colorbar');
for ip = 1:numel(hb)  
  hb(ip).FontSize = fontsize;
  hb(ip).Location = 'northoutside';
  hb(ip).Position(2) = 0.81;  
end

%% plotline, vertical, Ez balance
comp = 'z';
twpe = [24000];
xlim = 75+1*[-1 1];
zlim = [-1 8];
% 
% twpe = [21500];
% xlim = 80+1*[-1 1];
% %xlim = 90+1*[-1 1];


pic = no02m.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);

% twpe = 10000;
% xlim = 160+1*[-1 1];

% pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);

varstrs = {{'Bx','Bz','Ey'};{'Ez','-vxBz([1 3 5])','divpz([1 3 5])','dvzdt([1 3 5])'};{'Ez','-vxBz([2 4 6])','-divpz([2 4 6])'};{'Ey','-vxBy([2 4 6])','divpy([2 4 6])','divpy([1 3 5])'};{'pxx([1 3 5])','pyy([1 3 5])','pzz([1 3 5])','pxy([1 3 5])','pxz([1 3 5])','pyz([1 3 5])'}};
varstrs = {{'n(1)','n([3 5])','n(2)','n([4 6])'};{'Ez','-vxBz([1 3 5])','-vxBz([3 5])','-vxBz([2])','-vxBz([4 6])'};{'Ez','-vxBz([1 3 5])','divpz([1 3 5])','dvzdt([1 3 5])','vdvz([1 3 5])','-vxBz([1 3 5])+divpz([1 3 5])+dvzdt([1 3 5])+vdvz([1 3 5])'};{'Ez','-vxBz([2 4 6])','-divpz([2 4 6])'}};
varstrs = {{'PB','pi','pe','PB+pi+pe'};{'n(1)','n([3 5])','n(2)','n([4 6])'};{'Ez','-vxBz([1 3 5])','-vxBz([3 5])','-vxBz([2])','-vxBz([4 6])'};{'Ez','-vxBz([1])','divpz([1])','dvzdt([1])','vdvz([1])','-vxBz([1])+divpz([1])+dvzdt([1])+vdvz([1])'};{'Ez','-vxBz([3 5])','divpz([3 5])','dvzdt([3 5])','vdvz([3 5])','-vxBz([3 5])+divpz([3 5])+dvzdt([3 5])+vdvz([3 5])'};{'Ez','-vxBz([2 4 6])','-divpz([2 4 6])'}};
varstrs = {...%{'Ez+vxBz([3 5])','divpz([3 5])./n([3 5])'};...
  {'Ez','-vxBz([3 5])'};...
  {'Ez+vxBz([3 5])','divpz([3 5])./n([3 5])+dvzdt([3 5])+vdvz([3 5])'};...
  ...%{'Ez+vxBz([3 5])','divpz([3 5])./n([3 5])','dvzdt([3 5])'};...
  ...%{'Ez','-vxBz([3 5])','divpz([3 5])./n([3 5])','dvzdt([3 5])','vdvz([3 5])','-vxBz([3 5])+divpz([3 5])./n([3 5])+dvzdt([3 5])+vdvz([3 5])'};...
  {'divpz([3 5])./n([3 5])','dvzdt([3 5])','vdvz([3 5])'};...
  };
% varstrs = {{'Ez+vxBz([3 5])','divpz([3 5])./n([3 5])'};...
%   {'Ez'};...
%   };

h = pic.plot_line(comp,varstrs','smooth',50,'vertical');
%
colors = pic_colors('matlab');
h_lines = h(2).Children(end:-1:1);

for ip = 1:numel(h_lines)
  h_lines(ip).Color = colors(2+ip,:);  
end

for ip = 1:numel(h)
  h(ip).FontSize = 14;
  h(ip).XLim = 0.37*[-1 1];
  h(ip).XLabel.String = 'E_z (v_{A}B_0)';
end

%% Find some nice trajectories
tr = tr100;
trs = tr.find([tr.z0]==4);
trs = tr100(2:10);
trs = tr100.find([tr100.ncross]<2,[tr100.Ustart]<1,tr100.xstart>95);
trs = tr100.find([tr100.ncross]<2,[tr100.ncross]>0,[tr100.Ustart]<1,[tr100.z0]<1);
%trs.plot_all_xz

nrows = 2;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
holdon = 0;

for itr = 1:trs.ntr
  isub = 1;
  tr = trs(itr);
  if 1
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t,tr.Ex)
  end
  if 1
    hca = h(isub); isub = isub + 1;
    Erms = sqrt(sum(tr.Ex.^2));
    Erms = rms(tr.Ex.^2);
    plot(hca,abs(tr.zstart),Erms,'o')
  end
  if not(holdon)
    for ip = 1:npanels
      hold(h(ip),'on')
    end
    holdon = 1; 
  end
end
if holdon
  for ip = 1:npanels
    hold(h(ip),'off')
  end
  holdon = 0; 
end
  


