%% Define times
timesteps = 00200:200:10800;
ntimes = numel(timesteps);
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/';
data_dir = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data/';
screensize = get( groot, 'Screensize' );

%% Define xlim zlim, common for all figures, do this below, after having loaded x and z
% xlim = [x(1) x(end)] + [100 -100];
% zlim = [-25 25];
% xlim = [x(1) x(end)];
% zlim = [z(1) z(end)];

%% Quantities to plot or save
clear subdirs_all varstrs_all clims_all cylims_all nrows_all plot_structure      
iplot = 0;
if 0 % vex, vepar
  iplot = iplot + 1;

  subdirs_all{iplot} = 'vex_vepar';
  varstrs_all{iplot} = {'ve1.x','ve2.x','ve1.par','ve2.par'};
  clims_all{iplot} = [-3 3];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 2; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % vix, vipar
  iplot = iplot + 1;

  subdirs_all{iplot} = 'vix_vipar';
  varstrs_all{iplot} = {'vi1.x','vi2.x','vi1.par','vi2.par'};
  clims_all{iplot} = 0.5*[-1 1];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 2; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % vex vey vez
  iplot = iplot + 1;

  subdirs_all{iplot} = 've_xyz';
  varstrs_all{iplot} = {'ve1.x','ve2.x','ve1.y','ve2.y','ve1.z','ve2.z'};
  clims_all{iplot} = [-3 3];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 3; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % vix viy viz
  iplot = iplot + 1;

  subdirs_all{iplot} = 'vi_xyz';
  varstrs_all{iplot} = {'vi1.x','vi2.x','vi1.y','vi2.y','vi1.z','vi2.z'};
  clims_all{iplot} = 0.5*[-1 1];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 3; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end  
if 0 % ni, ne
  iplot = iplot + 1;

  subdirs_all{iplot} = 'n';
  varstrs_all{iplot} = {'ne1','ne2','ni1','ni2'};
  clims_all{iplot} = [-3 3];        
  cylims_all{iplot} = [0 clims_all{iplot}(2)];
  nrows_all{iplot} = 2; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % Pe1 tensor
  iplot = iplot + 1;

  subdirs_all{iplot} = 'pe1_tensor';
  varstrs_all{iplot} = {'pe1.xx','pe1.xy','pe1.yy','pe1.xz','pe1.zz','pe1.yz'};
  clims_all{iplot} = 0.25*[-1 1];        
  cylims_all{iplot} = [0 clims_all{iplot}(2)];
  nrows_all{iplot} = 3; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % Pi1 tensor
  iplot = iplot + 1;

  subdirs_all{iplot} = 'pi1_tensor';
  varstrs_all{iplot} = {'pi1.xx','pi1.xy','pi1.yy','pi1.xz','pi1.zz','pi1.yz'};
  clims_all{iplot} = 1*[-1 1];        
  cylims_all{iplot} = [0 clims_all{iplot}(2)];
  nrows_all{iplot} = 3; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % Pe2 tensor
  iplot = iplot + 1;

  subdirs_all{iplot} = 'pe2_tensor';
  varstrs_all{iplot} = {'pe2.xx','pe2.xy','pe2.yy','pe2.xz','pe2.zz','pe2.yz'};
  clims_all{iplot} = 0.25*[-1 1];        
  cylims_all{iplot} = [0 clims_all{iplot}(2)];
  nrows_all{iplot} = 3; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % Pi2 tensor
  iplot = iplot + 1;

  subdirs_all{iplot} = 'pi2_tensor';
  varstrs_all{iplot} = {'pi2.xx','pi2.xy','pi2.yy','pi2.xz','pi2.zz','pi2.yz'};
  clims_all{iplot} = 1*[-1 1];        
  cylims_all{iplot} = [0 clims_all{iplot}(2)];
  nrows_all{iplot} = 3; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % E,B forces on hot electrons
  iplot = iplot + 1;

  subdirs_all{iplot} = '-ve1xB_E_sum_xyz';
  varstrs_all{iplot} = {'-ve1xB.x','-ve1xB.y','-ve1xB.z','-E.x','-E.y','-E.z','-ve1xB.x-E.x','-ve1xB.y-E.y','-ve1xB.z-E.z'};
  clims_all{iplot} = [-1 1];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 3; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % E,B forces on cold electrons
  iplot = iplot + 1;

  subdirs_all{iplot} = '-ve2xB_E_sum_xyz';
  varstrs_all{iplot} = {'-ve2xB.x','-ve2xB.y','-ve2xB.z','-E.x','-E.y','-E.z','-ve2xB.x-E.x','-ve2xB.y-E.y','-ve2xB.z-E.z'};
  clims_all{iplot} = [-1 1];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 3; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % E,B forces on hot electrons
  iplot = iplot + 1;

  subdirs_all{iplot} = 'vi1xB_E_sum_xyz';
  varstrs_all{iplot} = {'vi1xB.x','vi1xB.y','vi1xB.z','E.x','E.y','E.z','vi1xB.x+E.x','vi1xB.y+E.y','vi1xB.z+E.z'};
  clims_all{iplot} = [-1 1];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 3; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % E,B forces on cold electrons
  iplot = iplot + 1;

  subdirs_all{iplot} = 'vi2xB_E_sum_xyz';
  varstrs_all{iplot} = {'vi2xB.x','vi2xB.y','vi2xB.z','E.x','E.y','E.z','vi2xB.x+E.x','vi2xB.y+E.y','vi2xB.z+E.z'};
  clims_all{iplot} = [-1 1];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 3; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % E,B,divP forces on hot electrons
  iplot = iplot + 1;

  subdirs_all{iplot} = 'e1_force_terms';
  varstrs_all{iplot} ={'-ne1.*ve1xB.x','-ne1.*ve1xB.y','-ne1.*ve1xB.z','-ne1.*E.x','-ne1.*E.y','-ne1.*E.z',...
           '-ne1.*(ve1xB.x+E.x)','-ne1.*(ve1xB.y+E.y)','-ne1.*(ve1xB.z+E.z)',...
           '-gradpe1_smooth.x','-gradpe1_smooth.y','-gradpe1_smooth.z',...
           '-ne1.*(ve1xB.x+E.x)-gradpe1_smooth.x','-ne1.*(ve1xB.y+E.y)-gradpe1_smooth.y','-ne1.*(ve1xB.z+E.z)-gradpe1_smooth.z'...
           }; 
  clims_all{iplot} = 0.5*[-1 1];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 5; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % E,B forces on cold electrons
  iplot = iplot + 1;

  subdirs_all{iplot} = 'e2_force_terms';
  varstrs_all{iplot} ={'-ne2.*ve2xB.x','-ne2.*ve2xB.y','-ne2.*ve2xB.z','-ne2.*E.x','-ne2.*E.y','-ne2.*E.z',...
           '-ne2.*(ve2xB.x+E.x)','-ne2.*(ve2xB.y+E.y)','-ne2.*(ve2xB.z+E.z)',...
           '-gradpe2_smooth.x','-gradpe2_smooth.y','-gradpe2_smooth.z',...
           '-ne2.*(ve2xB.x+E.x)-gradpe2_smooth.x','-ne2.*(ve2xB.y+E.y)-gradpe2_smooth.y','-ne2.*(ve2xB.z+E.z)-gradpe2_smooth.z'...
           }; 
  clims_all{iplot} = 0.5*[-1 1];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 5; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % E,B forces on hot ions
  iplot = iplot + 1;

  subdirs_all{iplot} = 'i1_force_terms';
  varstrs_all{iplot} ={'ni1.*vi1xB.x','ni1.*vi1xB.y','ni1.*vi1xB.z','ni1.*E.x','ni1.*E.y','ni1.*E.z',...
           'ni1.*(vi1xB.x+E.x)','ni1.*(vi1xB.y+E.y)','ni1.*(vi1xB.z+E.z)',...
           '-gradpi1_smooth.x','-gradpi1_smooth.y','-gradpi1_smooth.z',...
           'ni1.*(vi1xB.x+E.x)-gradpi1_smooth.x','ni1.*(vi1xB.y+E.y)-gradpi1_smooth.y','ni1.*(vi1xB.z+E.z)-gradpi1_smooth.z'...
           }; 
  clims_all{iplot} = 0.5*[-1 1];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 5; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % E,B forces on cold ions
  iplot = iplot + 1;

  subdirs_all{iplot} = 'i2_force_terms';
  varstrs_all{iplot} ={'ni2.*vi2xB.x','ni2.*vi2xB.y','ni2.*vi2xB.z','ni2.*E.x','ni2.*E.y','ni2.*E.z',...
           'ni2.*(vi2xB.x+E.x)','ni2.*(vi2xB.y+E.y)','ni2.*(vi2xB.z+E.z)',...
           '-gradpi2_smooth.x','-gradpi2_smooth.y','-gradpi2_smooth.z',...
           'ni2.*(vi2xB.x+E.x)-gradpi2_smooth.x','ni2.*(vi2xB.y+E.y)-gradpi2_smooth.y','ni2.*(vi2xB.z+E.z)-gradpi2_smooth.z'...
           }; 
  clims_all{iplot} = 0.5*[-1 1];        
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 5; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end
if 0 % vi1xB.y, vi2xB.y, -vi1xB.y+vi2xB.y (y-forces on cold ions)
  iplot = iplot + 1;

  subdirs_all{iplot} = 'i2_force_terms';
  varstrs_all{iplot} = {'vi1xB.y','vi2xB.y','-vi1xB.y+vi2xB.y','vi1xB.y_zx','vi2xB.y_zx','-vi1xB.y_zx+vi2xB.y_zx','vi1xB.y_xz','vi2xB.y_xz','-vi1xB.y_xz+vi2xB.y_xz'};
  clims_all{iplot} = 0.2*[-1 1];          
  cylims_all{iplot} = clims_all{iplot};
  nrows_all{iplot} = 3; % ncols is calculated from nrows and nvars

  plot_structure.subdir = subdirs_all{iplot};
  plot_structure.varstrs = varstrs_all{iplot};
  plot_structure.clim = clims_all{iplot};
  plot_structure.cylim = cylims_all{iplot};
  plot_structure.nrows = nrows_all{iplot};
  plot_structures_all{iplot} = plot_structure;
end


  
% Time series of quantities
clear varstrs_ts_adapted
varstrs_ts = {'Uke1','Uke1','Uki1','Uki2',...
              'Ute1','Ute2','Uti1','Uti2',...
              'UB.tot','UB.x','UB.y','UB.z',...              
              'pe1_mean','pe2_mean','pi1_mean','pi2_mean',...              
              'pe1_std','pe2_std','pi1_std','pi2_std',...
              'E_mean','E_std'};
nvars_ts = numel(varstrs_ts);
for ivar_ts = 1:nvars_ts
  varstr_ts = varstrs_ts{ivar_ts};
  varstr_ts_adapted = [varstr_ts '_ts'];
  varstr_ts_adapted(strfind(varstr_ts_adapted,'.')) = '_';
  varstrs_ts_adapted{ivar_ts,:} = varstr_ts_adapted;
  eval([varstr_ts_adapted ' = nan(1,ntimes);']);
end

% Time stacked line plots (for example to see how Bz or vx spread out)
clear varstrs_ts_stacked_z0
zval_collect = [-10:1:10];
nz_collect = numel(zval_collect);
varstrs_ts_stacked = {'B.z','E.y',...
              've1.x','ve2.x','vi1.x','vi2.x',...
              've1.y','ve2.y','vi1.y','vi2.y',...
              'ne1','ne2','ni1','ni2',...
              'jtot.y',...
              'pe1.scalar','pe2.scalar','pi1.scalar','pi2.scalar',...
              'te1.scalar','te2.scalar','ti1.scalar','ti2.scalar'...
              };
nx = 6400;
nvars_ts_stacked = numel(varstrs_ts_stacked);
nvars_ts = numel(varstrs_ts);
for ivar_ts = 1:nvars_ts_stacked
  varstr_ts_stacked = varstrs_ts_stacked{ivar_ts};
  varstr_ts_stacked_adapted = [varstr_ts_stacked '_ts_stacked'];
  varstr_ts_stacked_adapted(strfind(varstr_ts_stacked_adapted,'.')) = '_';
  varstrs_ts_stacked_adapted{ivar_ts,:} = varstr_ts_stacked_adapted;
  eval([varstr_ts_stacked_adapted ' = nan(nx,nz_collect,ntimes);']);
end

%% Time loop
doTs = 1;
doPatch = 0;

for itime = 1:ntimes
  %% Load data
  timestep = timesteps(itime);
  txtfile = sprintf('%s/fields-%05.0f.dat',data_dir,timestep); % michael's perturbation
  tic; [x,z,E,B,...
        ni1,ne1,ni2,ne2,...
        vi1,ve1,vi2,ve2,...
        ji1,je1,ji2,je2,...
        pi1,pe1,pi2,pe2,...
        ti1,te1,ti2,te2,...
        dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q]... 
        = read_fields(txtfile); toc

  ind_z0 = find_closest_ind(z,zval_collect);
  
  xlim = [x(1) x(end)];
  zlim = [z(1) z(end)];
  
  %% Calculate auxillary quantities
  A = vector_potential(x,z,B.x,B.z); % vector potential
  [saddle_locations,saddle_values] = saddle(A);
  if 0
    levels_edges = -30:0.5:1;
    levels_centers = levels_edges(1:end-1) + 0.5*diff(levels_edges(1:end));
    tic; [A_volume,A_map,A_levels] = fluxtube_volume(A,levels_edges); toc
    xline_ind = find(saddle_values == max(saddle_values));
    xline_A = saddle_values(xline_ind);
    ind_inner = find(A<xline_A);
    ind_outer = find(A>xline_A);
    A_outer = A; A_outer(ind_inner) = NaN;
    A_inner = A; A_inner(ind_outer) = NaN;
    A_map_outer = A_map; A_map_outer(ind_inner) = NaN;
    A_map_inner = A_map; A_map_inner(ind_outer) = NaN;
    A_levels_outer = A_levels; A_levels_outer(ind_inner) = NaN;
    A_levels_inner = A_levels; A_levels_inner(ind_outer) = NaN;
  end
  pb = B.abs.^2/2; % magnetic pressure
  bcurv = magnetic_field_curvature(x,z,B.x,B.y,B.z); % magnetic curvature
  c_eval('ve?xB = cross_product(ve?.x,ve?.y,ve?.z,B.x,B.y,B.z);',1:2) % electron motional electric field
  c_eval('vi?xB = cross_product(vi?.x,vi?.y,vi?.z,B.x,B.y,B.z);',1:2) % ion motional electric field
  ExB = cross_product(E.x,E.y,E.z,B.x,B.y,B.z); % Poynting flux
  c_eval('E_ve?xB.x = E.x + ve?xB.x; E_ve?xB.y = E.y + ve?xB.y; E_ve?xB.z = E.z + ve?xB.z;',1:2) % electron motional electric field
  c_eval('E_vi?xB.x = E.x + vi?xB.x; E_vi?xB.y = E.y + vi?xB.y; E_vi?xB.z = E.z + vi?xB.z;',1:2) % electron motional electric field
  c_eval('je?E = je?.x.*E.x + je?.y.*E.y + je?.y.*E.z;',1:2)
  c_eval('ji?E = ji?.x.*E.x + ji?.y.*E.y + ji?.y.*E.z;',1:2)
  UB.x = 0.5*B.x.^2;
  UB.y = 0.5*B.y.^2;
  UB.z = 0.5*B.z.^2;
  UB.tot = 0.5*B.abs.^2; 
  c_eval('gradpe? = grad_scalar(x,z,pe?.scalar);',1:2) %
  c_eval('gradpi? = grad_scalar(x,z,pi?.scalar);',1:2) %
  c_eval('gradpe?_smooth = grad_scalar(x,z,smooth2(pe?.scalar,1));',1:2) %
  c_eval('gradpi?_smooth = grad_scalar(x,z,smooth2(pi?.scalar,1));',1:2) %
  c_eval('Uke? = mass(2)/mass(1)*0.5*ne?.*(ve?.x.^2 + ve?.y.^2 + ve?.z.^2);',1:2)
  c_eval('Uki? = mass(1)/mass(1)*0.5*ni?.*(vi?.x.^2 + vi?.y.^2 + vi?.z.^2);',1:2)
  c_eval('Ute? = pe?.scalar;',1:2)
  c_eval('Uti? = pi?.scalar;',1:2)
  Uke = Uke1 + Uke2;
  Uki = Uki1 + Uki2;
  Uktot = Uke + Uki;
  Ute = Ute1 + Ute2;
  Uti = Uti1 + Uti2;
  Uttot = Ute + Uti;
  jtot.x = ji1.x + ji2.x - je1.x - je2.x;
  jtot.y = ji1.y + ji2.y - je1.y - je2.y;
  jtot.z = ji1.z + ji2.z - je1.z - je2.z;
  
  % noise levels
  xbox = [20 25];
  zbox = [20 25];

  pe1_box = pe1.scalar(lim2ind(x,xbox),lim2ind(z,zbox));
  pe2_box = pe2.scalar(lim2ind(x,xbox),lim2ind(z,zbox));
  pi1_box = pi1.scalar(lim2ind(x,xbox),lim2ind(z,zbox));
  pi2_box = pi2.scalar(lim2ind(x,xbox),lim2ind(z,zbox));
  E_box = E.abs(lim2ind(x,xbox),lim2ind(z,zbox));
  
  pe1_mean = mean(pe1_box(:));
  pe2_mean = mean(pe2_box(:));
  pi1_mean = mean(pi1_box(:));
  pi2_mean = mean(pi2_box(:));
  E_mean = mean(E_box(:));

  pe1_std = std(pe1_box(:));
  pe2_std = std(pe2_box(:));
  pi1_std = std(pi1_box(:));
  pi2_std = std(pi2_box(:));
  E_std = std(E_box(:));
  
  %% Collect time series
  for ivar_ts = 1:nvars_ts
    disp([varstrs_ts_adapted{ivar_ts} '(1,itime) = sum(' varstrs_ts{ivar_ts} '(:));']);    
    eval([varstrs_ts_adapted{ivar_ts} '(1,itime) = sum(' varstrs_ts{ivar_ts} '(:));']);    
  end
  
  %% Collect stacked time series
  for ivar_ts = 1:nvars_ts_stacked
    disp([varstrs_ts_stacked_adapted{ivar_ts} '(:,:,itime) = ' varstrs_ts_stacked{ivar_ts} '(:,ind_z0);']);    
    eval([varstrs_ts_stacked_adapted{ivar_ts} '(:,:,itime) = ' varstrs_ts_stacked{ivar_ts} '(:,ind_z0);']);    
    if ivar_ts == 15
      subdir = 'jtot_at_x=0';
      savedir = [savedir_root,subdir];
      mkdir(savedir)
      savestr = sprintf('%s_t%05.0f',subdir,timestep);    
      figure(33)
      hca = subplot(2,1,1);
      plot(hca,x,eval([varstrs_ts_stacked{ivar_ts} '(:,11);']))      
      hca.XLabel.String = 'x (di)';
      hca.YLabel.String = varstrs_ts_stacked_adapted{ivar_ts};
      hca.Title.String = [varstrs_ts_stacked_adapted{ivar_ts} 'at z = 0'];
      hca = subplot(2,1,2);
      imagesc(hca,timesteps/wpewce/mass(1),x,squeeze(eval([varstrs_ts_stacked_adapted{ivar_ts} '(:,11,:)'])))
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = varstrs_ts_stacked_adapted{ivar_ts};
      hca.XLabel.String = 'time (1/wci)';
      hca.YLabel.String = 'x (di)';
      print('-dpng','-r200',[savedir '/' savestr '.png']);
      
    end
  end
  
  %% Plots
  if 0 % Plot, energy densities
    %% Save and print info
    subdir = 'energy_density_1';
    savedir = [savedir_root,subdir];
    mkdir(savedir)
    savestr = sprintf('%s_t%05.0f',subdir,timestep);
    % Define what variables to plot
    %varstrs = {'ve1.x','ve2.x','ve1.z','ve2.z','ve1.par','ve2.par','-ve1xB.x','-ve2xB.x','-ve1xB.z','-ve2xB.z','E.x','E.z'};
    varstrs = {'UB.tot','Uke1','Uke2','Uki1','Uki2','Ute1','Ute2','Uti1','Uti2'};
    clim = [-1 1];    
    nvars = numel(varstrs);

    % Initialize figure
    fig = figure(101);
    fig.Position = [screensize(1) screensize(2) screensize(3)*0.4 screensize(4)*0.7];
    npanels = nvars + doTs;
    nrows = 5;
    ncols = ceil(npanels/nrows);
    npanels = nrows*ncols;
    isub = 1; 
    for ipanel = 1:npanels  
      h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
    end
    clear hb;

    doA = 0;
    if doA    
      cA = [0.8 0.8 0.8];
      nA = 20;
      nA = [0:-2:min(A(:))];
    end
    
    % Plot part of data
    xlim = [x(1) x(end)] + [100 -100];
    zlim = [-10 10];
    ix1 = find(x>xlim(1),1,'first');
    ix2 = find(x<xlim(2),1,'last');
    iz1 = find(z>zlim(1),1,'first');
    iz2 = find(z<zlim(2),1,'last');
    ipx = ix1:2:ix2;
    ipz = iz1:2:iz2;
    
    % Panels
    isub = 1;
    if doTs % ts plot of energy
      hca = h(isub); isub = isub + 1;
      ts_varstrs = {'UB','Uke1','Uke2','Uki1','Uki2','Ute1','Ute2','Uti1','Uti2'};
      variables = nan(numel(ts_varstrs),ntimes);
      for ivar = 1:numel(ts_varstrs)
        variables(ivar,:) = eval([ts_varstrs{ivar} '_ts']);
      end         
      hlines = plot(hca,timesteps/wpewce/mass(1),[UB_ts; Uke1_ts; Uke2_ts; Uki1_ts; Uki2_ts; Ute1_ts; Ute2_ts; Uti1_ts; Uti2_ts]);      
      for iline = 1:numel(hlines)
        hlines(iline).Marker = '.';
      end
      legend(hca,ts_varstrs,'location','eastoutside')
%      labels = arrayfun(@(x,y) {[num2str(x) ' > Q_{||} > ' num2str(y)]}, edgesQ(end:-1:2),edgesQ(end-1:-1:1));
      hca.XLim = [0 (timesteps(end)+200)/wpewce/mass(1)];
      hca.XLabel.String = 'time (omega_{ci})';
      hca.YLabel.String = 'Energy density (...)';
    end
    for ivar = 1:nvars  
      hca = h(isub); isub = isub + 1;
      varstr = varstrs{ivar};
      variable = eval(varstr);  
      himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
      hca.Title.String = sprintf('%s',varstr); 
      hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      hcb = colorbar('peer',hca);
      hb(ivar) = hcb;
      %hcb.YLim = hca.CLim(2)*[-1 1];
      colormap(hca,pic_colors('blue_red'));

      if doA
        hold(hca,'on')
        hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    for ipanel = 2:npanels
      h(ipanel).YDir = 'normal';
      h(ipanel).XLim = xlim;
      h(ipanel).YLim = zlim;
      h(ipanel).CLim = clim;
      hb(ipanel-1).Limits(1) = 0;
    end
    
    print('-dpng','-r200',[savedir '/' savestr '.png']);
  end
  if 0 % Plot, vector potential
    %% Save and print info
    subdir = 'vector_potential';
    savedir = [savedir_root,subdir];
    mkdir(savedir)
    savestr = sprintf('%s_t%05.0f',subdir,timestep);                  

    % Initialize figure
    fig = figure(103);
    fig.Position = [screensize(1) screensize(2) screensize(3)*0.4 screensize(4)*0.7];
    npanels = nvars + doTs;
    nrows = 3;
    ncols = ceil(npanels/nrows);
    npanels = nrows*ncols;
    isub = 1; 
    for ipanel = 1:npanels  
      h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
    end
    clear hb;

    doA = 0;
    if doA    
      cA = [0.8 0.8 0.8];
      nA = 20;
      nA = [0:-2:min(A(:))];
    end
    
    % Plot part of data
    xlim = [x(1) x(end)];
    zlim = [-10 10];
    ix1 = find(x>xlim(1),1,'first');
    ix2 = find(x<xlim(2),1,'last');
    iz1 = find(z>zlim(1),1,'first');
    iz2 = find(z<zlim(2),1,'last');
    ipx = ix1:2:ix2;
    ipz = iz1:2:iz2;
    
    % Panels
    isub = 1;
    if 1 % plot of fluxtube volume vs A
      hca = h(isub); isub = isub + 1;
      hline = plot(hca,levels_centers,A_volume);
      hline.Marker = '.';
      hca.XLabel.String = 'A_level';
      hca.XLabel.Interpreter = 'none';  
      hca.YLabel.String = 'A_volume';
      hca.YLabel.Interpreter = 'none';  
    end
    if 1 % plot of fluxtube volume vs A, hold on between time steps
      hca = h(isub); isub = isub + 1;
      hold(hca,'on')
      hline = plot(hca,levels_centers,A_volume);
      hline.Marker = '.';
      hca.XLabel.String = 'A_level';
      hca.XLabel.Interpreter = 'none';  
      hca.YLabel.String = 'A_volume';
      hca.YLabel.Interpreter = 'none';  
      hca.Box = 'on';
    end
    if 0 % A
      hca = h(isub); isub = isub + 1;
      imagesc(hca,x,z,A')
      hca.Title.String = 'A';
      hca.Title.Interpreter = 'none';
      hca.XLabel.String = 'x (di)';
      hca.YLabel.String = 'z (di)';
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = 'rel. vol.';
    end
    if 1 % A_levels
      hca = h(isub); isub = isub + 1;
      imagesc(hca,x,z,A_levels')
      hcb = colorbar('peer',hca);
      hca.Title.String = 'A_levels';
      hca.Title.Interpreter = 'none';
      hca.XLabel.String = 'x (di)';
      hca.YLabel.String = 'z (di)';
    end
    if 1 % A_volume, automatic caxis
      hca = h(isub); isub = isub + 1;
      xlim = [x(2) x(end-1)];
      zlim = [z(2) z(end-1)];
      himag = imagesc(hca,x(lim2ind(x,xlim)),z(lim2ind(z,zlim)),A_map(lim2ind(x,xlim),lim2ind(z,zlim))');
      hca.CLim = [0.019 0.04];
      hca.Title.String = 'A_map';
      hca.Title.Interpreter = 'none';
      hca.XLabel.String = 'x (di)';
      hca.YLabel.String = 'z (di)';
      hcb = colorbar('peer',hca);  
      hcb.YLabel.String = 'rel. fluxtube vol.';
    end
    if 1 % A_volume, automatic caxis
      hca = h(isub); isub = isub + 1;
      xlim = mean(x) + [-50 50];
      zlim = mean(z) + [-10 10];
      himag = imagesc(hca,x(lim2ind(x,xlim)),z(lim2ind(z,zlim)),A_map(lim2ind(x,xlim),lim2ind(z,zlim))');
      hca.CLim = [0.019 0.04];
      hca.Title.String = 'A_map';
      hca.Title.Interpreter = 'none';
      hca.XLabel.String = 'x (di)';
      hca.YLabel.String = 'z (di)';
      hcb = colorbar('peer',hca);  
      hcb.YLabel.String = 'rel. fluxtube vol.';
      hca.XLim = xlim;
      hca.YLim = zlim;
    end
    if 1 % A_volume, automatic caxis
      hca = h(isub); isub = isub + 1;
      xlim = mean(x) + [-25 25];
      zlim = mean(z) + [-5 5];
      himag = imagesc(hca,x(lim2ind(x,xlim)),z(lim2ind(z,zlim)),A_map(lim2ind(x,xlim),lim2ind(z,zlim))');
      hca.CLim = [0.019 0.04];
      hca.Title.String = 'A_map';
      hca.Title.Interpreter = 'none';
      hca.XLabel.String = 'x (di)';
      hca.YLabel.String = 'z (di)';
      hcb = colorbar('peer',hca);  
      hcb.YLabel.String = 'rel. fluxtube vol.';
      hca.XLim = xlim;
      hca.YLim = zlim;
    end
    if 0 % A_volume, outside outermost saddle point (main x line)
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,A_map_outer');
      hca.CLim(1) = 0.018;
      hcb = colorbar('peer',hca);
      hca.XLabel.String = 'x (di)';
      hca.YLabel.String = 'z (di)';
    end
    if 0 % A_volume, inside outermost saddle point (main x line)
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,A_map_inner');
      hcb = colorbar('peer',hca);
      hca.XLabel.String = 'x (di)';
      hca.YLabel.String = 'z (di)';
    end
    if 0 % A_levels, outside outermost saddle point (main x line)
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,A_levels_outer');
      hcb = colorbar('peer',hca);
      hca.XLabel.String = 'x (di)';
      hca.YLabel.String = 'z (di)';
    end
    if 0 % A_levels, inside outermost saddle point (main x line)
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,A_levels_inner');
      hcb = colorbar('peer',hca);
      hca.XLabel.String = 'x (di)';
      hca.YLabel.String = 'z (di)';
    end


    for ipanel = 1:npanels
      h(ipanel).YDir = 'normal';
    end
    
    print('-dpng','-r200',[savedir '/' savestr '.png']);
  end  
    
  nplots = iplot;
  for iplot = 1:nplots % Plot, adaptive
    %% Save and print info
    % subdirs_all varstrs_all clims_all cylims_all nrows_all plot_structure     
    subdir = subdirs_all{iplot};
    varstrs = varstrs_all{iplot};
    clim = clims_all{iplot};
    cylim = cylims_all{iplot};
    nrows = nrows_all{iplot};
    
    nvars = numel(varstrs);
    npanels = nvars;        
    ncols = ceil(npanels/nrows);
    % set figure position
    screensize = get(groot,'Screensize');
    figure_position(1) = 1;
    figure_position(2) = 1;
    figure_position(3) = screensize(3)/4*ncols;
    figure_position(4) = figure_position(3)*nrows/ncols*0.5;
    
    savedir = [savedir_root,subdir];
    mkdir(savedir)
    savestr = sprintf('%s_t%05.0f',subdir,timestep);            

    % Initialize figure
    fig = figure(102);
    fig.Position = figure_position;

    isub = 1; 
    for ipanel = 1:npanels  
      h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
    end
    clear hb;

    doA = 1;
    if doA    
      cA = [0.8 0.8 0.8];
      cA = [0.7 0.7 0.7];
      nA = 40;
      nA = [0:-1:min(A(:))];
      ipxA = ipx1:20:ipx2;
      ipzA = ipz1:20:ipz2;
    end
    
    % Plot part of data    
    ix1 = find(x>xlim(1),1,'first');
    ix2 = find(x<xlim(2),1,'last');
    iz1 = find(z>zlim(1),1,'first');
    iz2 = find(z<zlim(2),1,'last');
    ipx = ix1:2:ix2;
    ipz = iz1:2:iz2;
    
    % Panels
    isub = 1;
    for ivar = 1:nvars
      hca = h(isub); isub = isub + 1;
      varstr = varstrs{ivar};
      variable = eval(varstr);  
      himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
      hca.Title.String = sprintf('%s',varstr); 
      hca.Title.Interpreter = 'none';
      if abs(himag.CData(:)) % dont do if is zero
        hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      end
      hcb = colorbar('peer',hca);
      hb(ivar) = hcb;
      %hcb.YLim = hca.CLim(2)*[-1 1];
      colormap(hca,pic_colors('blue_red'));

      if doA
        hold(hca,'on')
        hcont = contour(hca,x(ipxA),z(ipzA),A(ipxA,ipzA)',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    for ipanel = 1:npanels
      h(ipanel).YDir = 'normal';
      h(ipanel).XLim = xlim;
      h(ipanel).YLim = zlim;
      h(ipanel).CLim = clim;      
    end
    
    print('-dpng','-r200',[savedir '/' savestr '.png']);
  end
end