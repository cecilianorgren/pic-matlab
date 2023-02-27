%% Map plot of electorn Ohm's law
twpe = 16000;
varstrs = {'Ey','vxBy([2 4 6])','Ey+vxBy([2 4 6])','-dxPexy./ne','-dzPezy./ne','pexy','peyz'}';
climE = 0.3*[-1 1];
climP = [-0.03 0.03];
clims = {climE,climE,climE,climE,climE,climP,climP};
no02m.twpelim(twpe).xlim(mean(no02m.xi)+[-5 5]).zlim(2*[-1 1]).plot_map(varstrs,'A',0.2,'smooth',3)
colormap(pic_colors('blue_red'))

%% Line plot
varstrs = {{'Bz','vex'},{'Ey','vxBy([2 4 6])','Ey+vxBy([2 4 6])','-dxPexy./ne','-dzPezy./ne'}}';
no02m.twpelim(twpe).xlim(mean(no02m.xi)+[-8 8]).zlim(0.05*[-1 1]).plot_line('x',varstrs,'smooth',20)

%% Distribution plot
%ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists.h5');
sumdim = 1;
species = 1;
twpe = 24000;
%ds = ds100.twpelim(twpe).xlim([66 80]);
ds = ds.update_inds({1}); % just get one
xlim = [ds.xi1{1} ds.xi2{1}];
zlim = [ds.zi1{1} ds.zi2{1}];
n = mean(mean(no02m.twpelim(twpe).xlim(xlim).zlim(zlim).n(species)));
vx = mean(mean(no02m.twpelim(twpe).xlim(xlim).zlim(zlim).vx(species)));
vy = mean(mean(no02m.twpelim(twpe).xlim(xlim).zlim(zlim).vy(species)));
vz = mean(mean(no02m.twpelim(twpe).xlim(xlim).zlim(zlim).vz(species)));
disp(sprintf('n = %g, vx = %g, vy = %g, vz = %g',n,vx,vy,vz))

nrows = 2;
ncols = 1;
ipanel = 0;
for irow = 1:nrows
  for icol = 1:ncols
    ipanel = ipanel + 1;
    h(ipanel) = subplot(nrows,ncols,ipanel);
  end
end

isub = 1;

hca = h(isub); isub = isub + 1;
hout = ds.plot_map(hca,species,sumdim);
linkprop(hout.ax,{'CLim'})
hout.ax(1).CLim = 0.2*[-1 1];
colormap(pic_colors('blue_red'))

hca = h(isub); isub = isub + 1;
hout = ds.plot_map(hca,species,sumdim,'off-diag');
linkprop(hout.ax,{'CLim'})
hout.ax(1).CLim = 0.2*[-1 1];
colormap(pic_colors('blue_red'))


%% Moments from distribution
%ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists.h5');
sumdim = 1;
species = 1;
twpe = 23000;
ds = ds100.twpelim(twpe).xlim([66 80]);
ds = ds.update_inds({1}); % just get one

f = ds.f(1,1,species);
f2 = ds.fxyz(1,1,species,1);

xlim = f.x;
zlim = f.z;
n = mean(mean(no02m.twpelim(twpe).xlim(xlim).zlim(zlim).n(species)));
vx = mean(mean(no02m.twpelim(twpe).xlim(xlim).zlim(zlim).vx(species)));
vy = mean(mean(no02m.twpelim(twpe).xlim(xlim).zlim(zlim).vy(species)));
vz = mean(mean(no02m.twpelim(twpe).xlim(xlim).zlim(zlim).vz(species)));
disp(sprintf('n = %g, vx = %g, vy = %g, vz = %g',n,vx,vy,vz))


[V1,V2,V3] = ndgrid(f.v,f.v,f.v);  
n_d = sum(sum(sum(f.fyz)))*f.dv; % the f.dv multi, is for the already compressed direction
vbulk1 = sum(sum(f.f.*V1))*(f.dv^2)/n;
vbulk2 = sum(sum(f.f.*V2))*(f.dv^2)/n;
f.f = f.f.*(V1-vbulk1).*(V2-vbulk2);

%%
ds = ds100.findx(95).zlim([-1 1]);
inds = ds.indices{1};
ds = ds.update_inds({3:(numel(inds)-2)});

ds.plot_map([4 6],1,'off-diag')