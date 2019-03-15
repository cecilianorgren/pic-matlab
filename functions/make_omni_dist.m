function varargout = make_omni_dist(dir,distnumbers,nspecies,mass,energy_edges)
tic;
ndists = numel(distnumbers);

% First check how many files exist and remove distnumbers who does not
% exist
ndists = 0;
for idist = distnumbers
  txtfile = sprintf('%s/%.0f.dat',dir,idist); 
  if not(exist(txtfile,'file'))
    warning(sprintf('File not found: %s', txtfile))
    distnumbers(idist) = []; % remove from array
    continue
  else
    ndists = ndists + 1;
  end
end

% If energy is not given, define it from velocity axes
energy_set = 0;
nbins = 100; % default
if not(isempty(energy_edges))
  energy_centers = tocolumn((energy_edges(2:end)+energy_edges(1:end-1))*0.5);
  nbins = numel(energy_edges)-1;
  energy_set = 1;
end

% Initialize arrays
f = cell(1,nspecies);
for ispecies = 1:nspecies
  f{ispecies}.distnumber = distnumbers;
  f{ispecies}.x = zeros(ndists,1);
  f{ispecies}.z = zeros(ndists,1);
  f{ispecies}.dist_sum = zeros(ndists,nbins);
  f{ispecies}.dist_mean = zeros(ndists,nbins);
end

for idist = 1:ndists
  distnumber = distnumbers(idist);   
  
  % Load data  
  txtfile = sprintf('%s/%.0f.dat',dir,distnumber);    
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
      = read_distributions(txtfile);
    
  % Center coordinate of bin
  xc = (xlo+xhi)/2;
  zc = (zlo+zhi)/2;
  
  vx = axes;
  vy = axes;
  vz = axes;
  disp(sprintf('distnumber: %4.f, [xlo xhi zlo zhi] = [%7.2f %7.2f %7.2f %7.2f]',distnumber,xlo,xhi,zlo,zhi))  
  
  vabs = sqrt(vx.^2+vy.^2+vz.^2); %vabs = vabs(51:end,:);
  for ispecies = 1:4
    if idist == 1 
      if not(energy_set) % Base energy on velocity axes
        energy = mass(ispecies)/mass(1)*vabs.^2/2;    
        f{ispecies}.energy_edges = linspace(0,0.2*max(energy(:,ispecies)),nbins);
        f{ispecies}.energy_centers = torow((f{ispecies}.energy_edges(2:end)+f{ispecies}.energy_edges(1:end-1))*0.5);      
      else
        f{ispecies}.energy_edges = energy_edges;
        f{ispecies}.energy_centers = torow(energy_centers);  
      end      
    end
  
    f{ispecies}.x(idist,1) = xc;
    f{ispecies}.z(idist,1) = zc;
    %energy = mass(ispecies)/mass(1)*vabs.^2/2;
    ftmp = fxyz(:,:,:,ispecies);
    ftot = sum(ftmp(:));
    [VX,VY,VZ] = ndgrid(vx(:,ispecies),vy(:,ispecies),vz(:,ispecies));
    VABS = sqrt(VX.^2+VY.^2+VZ.^2); %vabs = vabs(51:end,:);
    ENERGY = mass(ispecies)/mass(1)*VABS.^2/2;
    %energy_edges = logspace(-2,log10(1.0*max(energy(:,ispecies))),20);
    %energy_edges = linspace(0,0.2*max(energy(:,ispecies)),100);
    %ienergy = hist(energy(:,ispecies))
    
    [N,EDGES,BIN] = histcounts(ENERGY(:),f{ispecies}.energy_edges);            
    for ibin = 1:nbins
      ind_bin = find(BIN==ibin);      
      f_dist_tmp = ftmp(ind_bin);
      f{ispecies}.dist_mean(idist,ibin) = mean(f_dist_tmp);
      f{ispecies}.dist_sum(idist,ibin) = sum(f_dist_tmp);
    end    
  end
end
toc

for ispecies = 1:nspecies
  varargout{ispecies} = f{ispecies};
end