% ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists_new.h5');

%% Define boxes
twpe = 24000;
xlim = [105 112];
zlim = [-8 8];
xlim_ds = [105 112];
zlim_ds = [-1 2.5];


twpe = 16000;
xlim = [98 107];
zlim = [-2 2];
xlim_ds = [99 106];
zlim_ds = [-1 1];


%zpick = ;


ds = ds100.twpelim(twpe).xlim(xlim_ds).zlim(zlim_ds); 

x0 = (ds.xi2{1}+ds.xi1{1})/2;
z0 = (ds.zi2{1}+ds.zi1{1})/2;
zpick = -0.4:0.1:0.4;
xpick = 99:0.4:102.4;

zpick = -0.2:0.05:0.2;
xpick = 100.8:0.2:102.4;

ds = ds.xfind(xpick).zfind(zpick);
% ds.plot_map(4,3,'off-diag')
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim); 

%% Plot box locations
figure(501)

varstrs = {'pexy','pexz','peyz','log10(abs(tepar)./abs(teperp))'}';
cmapbr = pic_colors('blue_red');
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr};
clims = {[-1 1]*10e-3,[-1 1]*10e-3,[-1 1]*10e-3,[-1 1]};

h = pic.plot_map(varstrs,'cmap',cmaps,'clim',clims,'A',1,'sep');

for ip = 1:numel(h)
  hd = ds.plot_boxes(h(ip));
end

%% Plot off-diagonal components
fig = figure(502); delete(fig.Children);
h = ds.plot_map([2 4 6],3,'off-diag','bline',pic);
colormap(pic_colors('blue_red'));
compact_panels(0.00,0.00)

for ip = 1:numel(h.ax)
  axis(h.ax(ip),'square')
  h.ax(ip).XGrid = 'on'; 
  h.ax(ip).YGrid = 'on';
  h.ax(ip).YTick = [-20:5:20];
  h.ax(ip).XTick = [-20:5:20];
end
%% Plot distributions
fig = figure(503); delete(fig.Children);
h = ds.plot_map([2 4 6],2);
compact_panels(0.00,0.00)

for ip = 1:numel(h.ax)
  axis(h.ax(ip),'square')
  h.ax(ip).XGrid = 'on'; 
  h.ax(ip).YGrid = 'on';
  h.ax(ip).YTick = [-20:5:20];
  h.ax(ip).XTick = [-20:5:20];
end
%% Plot forces
fig = figure(504); delete(fig.Children);
%h = ds.plot_map([2 4 6],3,'bline',pic,'force','vx*Bz-Ey',pic);
h = ds.plot_map([2 4 6],3,'bline',pic,'force','-vy*Bz-Ex',pic);

colormap(pic_colors('blue_red'));
compact_panels(0.00,0.00)

for ip = 1:numel(h.ax)
  axis(h.ax(ip),'square')
  h.ax(ip).XGrid = 'on'; 
  h.ax(ip).YGrid = 'on';
  h.ax(ip).YTick = [-20:5:20];
  h.ax(ip).XTick = [-20:5:20];
end
hlinks = linkprop(h.ax,{'CLim'});
hlinks.Targets(1).CLim = max(abs(hlinks.Targets(1).CLim))*[-1 1];