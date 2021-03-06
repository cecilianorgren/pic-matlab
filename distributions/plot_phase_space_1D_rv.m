% Plot cuts of 1D phase space f(x,vx), f(x,vy), f(x,vz), etc
df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
df08 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data_h5/fields.h5');
ds04 = PICDist('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/dists.h5');
tr04 = PICTraj('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5');
savedir  = ['/Users/' localuser '/GoogleDrive/DF_Cold_ions_figures_for_AGU/'];

%%
zlim = 0+[-0.2 0.2];
xlim = [0 205];
its = 2;
twci = 160;
iss = [2];
icold = [3 5];

rdim = 1; % x
vdim = 2; % y
f1_x_vx = make_space_time_1d(ds04,its,xlim,zlim,rdim,1,1);
f1_x_vy = make_space_time_1d(ds04,its,xlim,zlim,rdim,2,1);
f1_x_vz = make_space_time_1d(ds04,its,xlim,zlim,rdim,3,1);
f35_x_vx = make_space_time_1d(ds04,its,xlim,zlim,rdim,1,icold);
f35_x_vy = make_space_time_1d(ds04,its,xlim,zlim,rdim,2,icold);
f35_x_vz = make_space_time_1d(ds04,its,xlim,zlim,rdim,3,icold);

% Get X line location
A = df04.twcilim(twci).A;
[saddle_locations,saddle_values] = saddle(A);

nrows = 3;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

if 1 % f(x,vx)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,f1_x_vx.x,f1_x_vx.v(1,:),f1_x_vx.f'); 
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'v_x/v_A';
end
if 1 % f(x,vy)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,f1_x_vy.x,f1_x_vy.v(1,:),f1_x_vy.f'); 
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'v_y/v_A';
end
if 1 % f(x,vz)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,f1_x_vz.x,f1_x_vz.v(1,:),f1_x_vz.f'); 
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'v_z/v_A';
end

if 1 % f(x,vx)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,f35_x_vx.x,f35_x_vx.v(1,:),f35_x_vx.f'); 
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'v_x/v_A';
end
if 1 % f(x,vy)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,f35_x_vy.x,f35_x_vy.v(1,:),f35_x_vy.f'); 
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'v_y/v_A';
end
if 1 % f(x,vz)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,f35_x_vz.x,f35_x_vz.v(1,:),f35_x_vz.f'); 
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'v_z/v_A';
end

maxc = 0;
for ip = 1:npanels
  hca = h(ip);
  hb(ip) = colorbar('peer',h(ip));
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  maxc = max([maxc hca.CLim(2)]);
  hca.YLim = [-3 3];
  hold(hca,'on')
  for isaddle = 1:size(saddle_locations,1)
    plot(hca,df04.xi(3181)*[1 1],hca.YLim,'--k')
  end
  hold(hca,'off')
end

h(1).Title.String = sprintf('z = [%.1f %.1f] d_i\nhot ions',zlim(1),zlim(2));
h(4).Title.String = sprintf('z = [%.1f %.1f] d_i\ncold ions',zlim(1),zlim(2));
colormap(pic_colors('candy'))
compact_panels(0.01)
%%
hlinks_1 = linkprop(h(1:3),{'CLim'});
hlinks_1.Targets(1).CLim = [0 0.2*maxc];

hlinks_35 = linkprop(h(4:6),{'CLim'});
hlinks_35.Targets(1).CLim = [0 0.5*maxc];


%%
it = 2; % 1600 (8000)
ds = ds04(it).zlim(zlim);

f = ds.f(1,1,3);