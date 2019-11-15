%% load data
timesteps = 00200:200:06000;
data_dir_resave = '/Volumes/pic/finished_runs/turbulencerun/data_separated/';
[sim_info,E_ts] = fun_load_resaved_data(data_dir_resave,{'E'},timesteps);
[sim_info,B_ts] = fun_load_resaved_data(data_dir_resave,{'B'},timesteps);
[sim_info,vi12_ts] = fun_load_resaved_data(data_dir_resave,{'vi12'},timesteps);
[sim_info,ve12_ts] = fun_load_resaved_data(data_dir_resave,{'ve12'},timesteps);
[sim_info,ni12_ts] = fun_load_resaved_data(data_dir_resave,{'ni12'},timesteps);
[sim_info,ne12_ts] = fun_load_resaved_data(data_dir_resave,{'ne12'},timesteps);

x = sim_info.x-mean(sim_info.x);
z = sim_info.z;
%%
A_ts = vector_potential(x,z,B_ts(:,:,:,1),B_ts(:,:,:,3)); % vector potential
Babs_ts = sqrt(B_ts(:,:,:,1).^2 + B_ts(:,:,:,2).^2 + B_ts(:,:,:,3).^2);
ve12par_ts = (ve12_ts(:,:,:,1).*B_ts(:,:,:,1) + ve12_ts(:,:,:,2).*B_ts(:,:,:,2) + ve12_ts(:,:,:,3).*B_ts(:,:,:,3))./Babs_ts;
Epar_ts = (E_ts(:,:,:,1).*B_ts(:,:,:,1) + E_ts(:,:,:,2).*B_ts(:,:,:,2) + E_ts(:,:,:,3).*B_ts(:,:,:,3))./Babs_ts;

[sep_x,sep_z,sep_vepar] = interpolate_to_separatrix(ve12par_ts,x,z,x,[],A_ts);
[sep_x,sep_z,sep_ne] = interpolate_to_separatrix(ne12_ts,x,z,x,[],A_ts);
[sep_x,sep_z,sep_ni] = interpolate_to_separatrix(ni12_ts,x,z,x,[],A_ts);
[sep_x,sep_z,sep_Epar] = interpolate_to_separatrix(Epar_ts,x,z,x,[],A_ts);

%% Example plot, non-adaptive
npanels = 8;
h = setup_subplots(npanels,1);

isub = 1;

if 1 % R
  hca = h(isub); isub = isub + 1;
  plot(hca,timesteps/200,R.E,(timesteps(2:end)-0.5*timesteps(2)+0.5*timesteps(1))/200,R.A)
  hca.YLim = [0 0.5]; hca.XGrid = 'on'; hca.YGrid = 'on';
  legend(hca,{'E_y','A'})
  hca.YLabel.String = 'R';
end
if 1 % vepar
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,timesteps/200,x,sep_vepar')
  hca.CLim = [-11 11];
  hb = colorbar('peer',hca);
  hca.YLabel.String = 'x/d_i';
  hb.YLabel.String = 'v_{e,||}';
end
if 1 % ne
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,timesteps/200,x,sep_ne')
  hca.CLim = 0.11+1*0.11*[-1 1];
  hb = colorbar('peer',hca);
  hca.YLabel.String = 'x/d_i';
  hb.YLabel.String = 'n_{e}';
end
if 1 % ni
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,timesteps/200,x,sep_ni')
  hb = colorbar('peer',hca);
  hca.CLim = 0.11+1*0.11*[-1 1];
  hca.YLabel.String = 'x/d_i';
  hb.YLabel.String = 'n_{i}';
end
if 1 % ni - ne
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,timesteps/200,x,sep_ni'-sep_ne')
  hb = colorbar('peer',hca);
  hca.CLim = 1*0.2*0.11*[-1 1];
  hca.YLabel.String = 'x/d_i';
  hb.YLabel.String = 'n_{i}-n_e';
end
if 1 % Epar
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,timesteps/200,x,sep_Epar')
  hca.CLim = 0.3*[-1 1];
  hb = colorbar('peer',hca);
  hca.YLabel.String = 'x/d_i';
	hb.YLabel.String = 'E_{||}';
end
if 1 % ne*vepar
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,timesteps/200,x,sep_ne'.*sep_vepar')
  hb = colorbar('peer',hca);
  hca.CLim = 0.5*[-1 1];
  hca.YLabel.String = 'x/d_i';
  hb.YLabel.String = 'n_{e}v_{e,||}';
end
if 1 % z_sep
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,timesteps/200,x,sep_z')
  hb = colorbar('peer',hca);
  hca.CLim = 10*[-1 1];
  hca.YLabel.String = 'x/d_i';
  hb.YLabel.String = 'z';
end


colormap(pic_colors('blue_red'))
h(1).Title.String = 'Northern separatrix';
h(end).XLabel.String = 't\omega_{ci}';


c_eval('h(?).YTick = [-40:20:40];',2:npanels)
c_eval('h(?).YDir = ''normal'';',2:npanels)
c_eval('h(?).XLim = timesteps([1 end])/200;',1:npanels)

compact_panels(0.010)
irf_plot_axis_align(h)
%cn.print(sprintf('northern_serparatrix_R_vepar_ne_Epar_nevepar'))


%%

varstrs_map = {'ve12par_ts(it,:,:)'};
varstrs_time = {'sep_vepar','sep_ne'};

[X,Z] = meshgrid(x,z);
hca = subplot(1,1,1);

ntimes = size(B_ts,1);
xlim = [-50 50];
ylim = [-10 10];



for it = 1:ntimes  
  A = squeeze(A_ts(it,:,:));
  [saddle_locations,saddle_values] = saddle(A,'sort');
  AX = saddle_values(1);
  S = contourcs(x,z,A',AX*[1 1]*0.999);

  variable = squeeze(ve12par_ts(it,:,:,1));
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLim = xlim;
  hca.YLim = ylim;
  %hca.Title.String = varstrs{ivar};
  hca.YLabel.String = 'v_{e,||}/v_{te,0}';
  hca.XLabel.String = 'x/d_{di0}';
  pause(0.1)
end

%%
% plot for presentation
[saddle_locations,saddle_values] = saddle(A,'sort');
AX = saddle_values(1);
S = contourcs(x,z,A',AX*[1 1]*0.999);
[X,Z] = meshgrid(x,z);

xlim = [x(1) x(end)]+[3 -3]; 
zlim = [z(1) z(end)]; zlim = [-10 10];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:1:ipx2;
ipz = ipz1:1:ipz2;

% plot
nrows = 6;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1
  hca = h(isub); isub = isub + 1;  
  variable = ve12.x;
  imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)')
  %hca.Title.String = varstrs{ivar};
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_{e,x}/v_{A0}';
  hca.XLabel.String = 'x/d_{di0}';
  hca.YLabel.String = 'z/d_{di0}';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = 6*[-1 1];
  hold(hca,'on')
  contour(hca,x(1:5:end),z(1:5:end),A(1:5:end,1:5:end)',-25:2:0,'k')
  plot(hca,S(1).X,S(1).Y,'linewidth',3,'color',0*pic_colors('3')+[0 0 0])
  hold(hca,'off')
  hca.YDir = 'normal';
end

if 1
  hca = h(isub); isub = isub + 1;
  variable = ve12.par/vte12.par(100,100);
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);      
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.Title.String = varstrs{ivar};
  hca.YLabel.String = 'v_{e,||}/v_{te,0}';
  hca.XLabel.String = 'x/d_{di0}';
end
if 1
  hca = h(isub); isub = isub + 1;
  variable = ne12/ne12(100,100);
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);      
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.Title.String = varstrs{ivar};
  hca.YLabel.String = 'n_{e}/n_{0}';
  hca.XLabel.String = 'x/d_{di0}';
end
if 1 % flux/flux0
  hca = h(isub); isub = isub + 1;
  variable = ne12.*ve12.par/ne12(100,100)/vte12.par(100,100);
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);      
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.Title.String = varstrs{ivar};
  hca.YLabel.String = 'n_ev_{e,||}/n_0v_{te,0}';
  hca.XLabel.String = 'x/d_{di0}';
end
if 1 % vtepar/vt0
  hca = h(isub); isub = isub + 1;
  variable = vte12.par/vte12.par(100,100);
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);      
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.Title.String = varstrs{ivar};
  hca.YLabel.String = 'v_{te,||}/v_{te,0}';
  hca.XLabel.String = 'x/d_{di0}';
end
if 1 % Epar
  hca = h(isub); isub = isub + 1;
  variable = E.par;
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);      
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.Title.String = varstrs{ivar};
  hca.YLabel.String = 'E_{||}/B_0v_{A0}';
  hca.XLabel.String = 'x/d_{di0}';
end

irf_plot_axis_align
compact_panels(0.01)
hlink = linkprop(h,{'XLim'});
%%
varstrs = {'smooth2(vte12.par,np_smooth)/vte12.par(100,100)','ve12.par/vte12.par(100,100)','ne12/ne12(100,100)','ne12.*ve12.par/ne12(100,100)/vte12.par(100,100)'};




nvars = numel(varstrs);
vars = cell(nvars,1);

h = setup_subplots(2,nvars,'vertical');
isub = 1;
for ivar = 1:nvars
  variable = eval(varstrs{ivar});
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);
  vars{ivar} = Vq;
  hca = h(isub); isub = isub + 1;
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Title.String = varstrs{ivar};
  
  hca = h(isub); isub = isub + 1;
  imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)')
  hold(hca,'on')
  plot(hca,S(1).X,S(1).Y,'linewidth',2,'color',[0 0 0])
  hca.YDir = 'normal';
  hold(hca,'off')
  colormap(hca,pic_colors('blue_red'))
  hca.Title.String = varstrs{ivar};
  
end

hlink = linkprop(h,{'XLim'});