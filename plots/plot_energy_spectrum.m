%% Read and plot distributions
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/distributions/';
timestep = 08000;
str_timestep = sprintf('%05.0f',timestep);
txttime = sprintf('timestep = %05.0f',timestep); 
   
clear f_dist_all
idist = 0;
tic
for distnumber = 1:281%281%:281%:281%30:40%39%180:200%:250%:100%:100%:10%40%:40%:4%:40
  read_sub_dir = '/1/';
  txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber); % michael's perturbation
  if not(exist(txtfile,'file'))
    warning(sprintf('File not found: %s', txtfile))
    continue
  end  
    
  idist = idist + 1;
  %f_dist_all = [];
  %idist = distnumber;
  
  % Load data  
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
      = read_distributions(txtfile);
    
  if x(1)<-150
    xlo = xlo-x0;
    xhi = xhi-x0;
    zlo = zlo;
    zhi = zhi;    
  end
  xc = (xlo+xhi)/2;
  zc = (zlo+zhi)/2;
  
  vx = axes;
  vy = axes;
  vz = axes;
  Bloc.x = B.x(find_closest_ind(x,0.5*(xlo+xhi)),find_closest_ind(z,0.5*(zlo+zhi)));
  Bloc.y = B.y(find_closest_ind(x,0.5*(xlo+xhi)),find_closest_ind(z,0.5*(zlo+zhi)));
  Bloc.z = B.z(find_closest_ind(x,0.5*(xlo+xhi)),find_closest_ind(z,0.5*(zlo+zhi)));
  disp(sprintf('distnumber: %g, %7.2f %7.2f %7.2f %7.2f',distnumber,xlo,xhi,zlo,zhi))  
  
  vabs = sqrt(vx.^2+vy.^2+vz.^2); %vabs = vabs(51:end,:);
  f_energy_edges = cell(4,1);
  f_energy_centers = cell(4,1);
  f_dist_mean = cell(4,1);
  f_dist_sum = cell(4,1);  
  for ispecies = 1:4
    energy = mass(ispecies)/mass(1)*vabs.^2/2;
    ftmp = fxyz(:,:,:,ispecies);
    ftot = sum(ftmp(:));
    [VX,VY,VZ] = ndgrid(vx(:,ispecies),vy(:,ispecies),vz(:,ispecies));
    VABS = sqrt(VX.^2+VY.^2+VZ.^2); %vabs = vabs(51:end,:);
    ENERGY = mass(ispecies)/mass(1)*VABS.^2/2;
    energy_edges = logspace(-2,log10(1.0*max(energy(:,ispecies))),20);
    energy_edges = linspace(0,0.2*max(energy(:,ispecies)),100);
    %ienergy = hist(energy(:,ispecies))
    [N,EDGES,BIN] = histcounts(ENERGY(:),energy_edges);
    f_energy_edges{ispecies} = tocolumn(EDGES);
    f_energy_centers{ispecies} = tocolumn((EDGES(2:end)+EDGES(1:end-1))*0.5);
    nbins = (numel(energy_edges)-1);
    for ibin = 1:nbins
      ind_bin = find(BIN==ibin);      
      f_dist_tmp = ftmp(ind_bin);
      f_dist_mean{ispecies}(ibin,1) = mean(f_dist_tmp);
      f_dist_sum{ispecies}(ibin,1) = sum(f_dist_tmp);
    end
    f_dist_all.distnumber(idist,1) = distnumber;
    f_dist_all.x(idist,1) = xc;
    f_dist_all.z(idist,1) = zc;
    f_dist_all.f_energy_centers(idist,ispecies,:) = f_energy_centers{ispecies};
    f_dist_all.f_energy_edges(idist,ispecies,:) = f_energy_edges{ispecies};
    f_dist_all.f_dist_sum(idist,ispecies,:) = f_dist_sum{ispecies};
    f_dist_all.f_dist_mean(idist,ispecies,:) = f_dist_mean{ispecies};
  end
end
toc
 
%% Plot spectrogram
ivar = find(cellfun(@(x)strcmp(x,'B.z'),varstrs_ts_line_x));
zz = 1;
izpick = find_closest_ind(zpicks,zz);
Bz = squeeze(cell_ts_line_x{ivar}(:,izpick,find_closest_ind(timesteps,08000))); % Bz(x,z=0,t)

ndists = size(f_dist_all.x,1);
save_ind = find(f_dist_all.z == zz);
f_dist_fields = fields(f_dist_all);

for ifields = 1:numel(f_dist_fields)
  eval(['f_dist_z0.' f_dist_fields{ifields} ' = f_dist_all.' f_dist_fields{ifields} '(save_ind,:,:,:,:);'])
end

ispecies = [1 3];
xplot = f_dist_z0.x;
yplot = squeeze(f_dist_z0.f_energy_centers(1,ispecies(1),:));
cplot = squeeze(sum(f_dist_z0.f_dist_mean(:,ispecies,:),2));
cplot = squeeze(sum(f_dist_z0.f_dist_mean(:,ispecies,:),2).*f_dist_z0.f_energy_centers(:,ispecies(1),:).^2);

h = setup_subplots(2,1); isub = 1;
hca = h(isub); isub = isub + 1;
plot(hca,x,Bz)
hca = h(isub); isub = isub + 1;
pcolor(hca,xplot,yplot,log10(cplot)'); shading flat;
hcb = colorbar('peer',hca);
hca.CLim = [-6 -4];

hlink = linkprop(h,{'XLim'});
hlink.Targets(1).XLim = [-55 0];

%% Read and plot distributions, df04
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_n04/distributions/';
timestep = 05000;
str_timestep = sprintf('%05.0f',timestep);
txttime = sprintf('timestep = %05.0f',timestep); 
nss = 6;
mass = df04.mass;
z = df04.zi;
x = df04.xi;
%Bx = df04.twpelim(5000).Bx; By = df04.twpelim(5000).By; Bz = df04.twpelim(5000).Bz;
%B.x = Bx; B.y = By; B.z = Bz;

%clear f_dist_all
idist = 0;
tic
for distnumber = 13:259%281%:281%:281%30:40%39%180:200%:250%:100%:100%:10%40%:40%:4%:40
  read_sub_dir = '/1/';
  txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber); % michael's perturbation
  txtfile = sprintf('/Users/cno062/tesla/cno062/df_cold_protons_n04/distributions/%05.0f/%.0f.dat',timestep,distnumber); % df04
  if not(exist(txtfile,'file'))
    warning(sprintf('File not found: %s', txtfile))
    continue
  end  
    
  idist = idist + 1;
  %f_dist_all = [];
  %idist = distnumber;
  
  % Load data  
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
      = read_distributions(txtfile,nss);
    
  if x(1)<-150
    xlo = xlo-x0;
    xhi = xhi-x0;
    zlo = zlo;
    zhi = zhi;    
  end
  xc = (xlo+xhi)/2;
  zc = (zlo+zhi)/2;
  
  vx = axes;
  vy = axes;
  vz = axes;
  Bloc.x = B.x(find_closest_ind(x,0.5*(xlo+xhi)),find_closest_ind(z,0.5*(zlo+zhi)));
  Bloc.y = B.y(find_closest_ind(x,0.5*(xlo+xhi)),find_closest_ind(z,0.5*(zlo+zhi)));
  Bloc.z = B.z(find_closest_ind(x,0.5*(xlo+xhi)),find_closest_ind(z,0.5*(zlo+zhi)));
  disp(sprintf('distnumber: %g, %7.2f %7.2f %7.2f %7.2f',distnumber,xlo,xhi,zlo,zhi))  
  
  vabs = sqrt(vx.^2+vy.^2+vz.^2); %vabs = vabs(51:end,:);
  f_energy_edges = cell(nss,1);
  f_energy_centers = cell(nss,1);
  f_dist_mean = cell(nss,1);
  f_dist_sum = cell(nss,1);  
  for ispecies = 1:nss
    energy = mass(ispecies)/mass(1)*vabs.^2/2;
    ftmp = fxyz(:,:,:,ispecies);
    ftot = sum(ftmp(:));
    [VX,VY,VZ] = ndgrid(vx(:,ispecies),vy(:,ispecies),vz(:,ispecies));
    VABS = sqrt(VX.^2+VY.^2+VZ.^2); %vabs = vabs(51:end,:);
    ENERGY = mass(ispecies)/mass(1)*VABS.^2/2;
    energy_edges = logspace(-2,log10(1.0*max(energy(:,ispecies))),20);
    energy_edges = linspace(0,0.2*max(energy(:,ispecies)),100);
    %ienergy = hist(energy(:,ispecies))
    [N,EDGES,BIN] = histcounts(ENERGY(:),energy_edges);
    f_energy_edges{ispecies} = tocolumn(EDGES);
    f_energy_centers{ispecies} = tocolumn((EDGES(2:end)+EDGES(1:end-1))*0.5);
    nbins = (numel(energy_edges)-1);
    for ibin = 1:nbins
      ind_bin = find(BIN==ibin);      
      f_dist_tmp = ftmp(ind_bin);
      f_dist_mean{ispecies}(ibin,1) = mean(f_dist_tmp);
      f_dist_sum{ispecies}(ibin,1) = sum(f_dist_tmp);
    end
    f_dist_all.distnumber(idist,1) = distnumber;
    f_dist_all.x(idist,1) = xc;
    f_dist_all.z(idist,1) = zc;
    f_dist_all.f_energy_centers(idist,ispecies,:) = f_energy_centers{ispecies};
    f_dist_all.f_energy_edges(idist,ispecies,:) = f_energy_edges{ispecies};
    f_dist_all.f_dist_sum(idist,ispecies,:) = f_dist_sum{ispecies};
    f_dist_all.f_dist_mean(idist,ispecies,:) = f_dist_mean{ispecies};
  end
end
toc
 
%% Plot spectrogram
%ivar = find(cellfun(@(x)strcmp(x,'B.z'),varstrs_ts_line_x));
zz = 1;
%izpick = find_closest_ind(zpicks,zz);
%Bz = squeeze(cell_ts_line_x{ivar}(:,izpick,find_closest_ind(timesteps,08000))); % Bz(x,z=0,t)
Bz = mean(df04.twpelim(timestep).zlim([-1 1]).Bz,2);

ndists = size(f_dist_all.x,1);
save_ind = find(f_dist_all.z == zz);
f_dist_fields = fields(f_dist_all);

for ifields = 1:numel(f_dist_fields)
  eval(['f_dist_z0.' f_dist_fields{ifields} ' = f_dist_all.' f_dist_fields{ifields} '(save_ind,:,:,:,:);'])
end

ispecies = [1 3];
xplot = f_dist_z0.x;
yplot = squeeze(f_dist_z0.f_energy_centers(1,ispecies(1),:));
cplot = squeeze(sum(f_dist_z0.f_dist_mean(:,ispecies,:),2));
cplot = squeeze(sum(f_dist_z0.f_dist_mean(:,ispecies,:),2).*f_dist_z0.f_energy_centers(:,ispecies(1),:).^2);

h = setup_subplots(2,1); isub = 1;
hca = h(isub); isub = isub + 1;
plot(hca,x,Bz)
hca = h(isub); isub = isub + 1;
pcolor(hca,xplot,yplot,log10(cplot)'); shading flat;
hcb = colorbar('peer',hca);
hca.CLim = [-6 -4];

hlink = linkprop(h,{'XLim'});
hlink.Targets(1).XLim = [-55 0];



