sim = df04;
for it = 1:sim.length
  iter = sim.iteration(it);
  str_iteration = sprintf('%010.0f',iter); % same group format as SMILEI
  %it
  %tic
  Bx = sim(it).Bx; %toc
  Bz = sim(it).Bz; %toc
  A = vector_potential(sim.xi,sim.zi,squeeze(Bx),squeeze(Bz)); %toc
  %[saddle_locations,saddle_values] = saddle(A,'sort'); %toc
  %Ax(it) = saddle_values(1);
  %Eyx(it) = 
  %contour(sim.xi(1:10:end),sim.zi(1:10:end),A(1:10:end,1:10:end)',[-25:0.5:8],'k');  %toc
  %title(sprintf('t=%g',sim.twci(it))); drawnow;  
  dataset_name = ['/data/' str_iteration '/A'];
  disp(dataset_name)
  %h5create(sim.file, dataset_name, size(A));
  h5write( sim.file, dataset_name, A);
end

%% Thermal and kinetic energy n08
sim = df08;
tic;
clear UT UK
for it = 1:sim.length  
  sim_tmp = sim(it); 
  for iSpecies = 1:numel(sim.mass)    
    disp(sprintf('it = %g/%g, sp = %g ',it,sim.length,iSpecies))
    
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp(iSpecies);
    %toc;
    pdyn = sim.mass(iSpecies)/sim.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
    p = (pxx+pyy+pzz)/3; % scalar pressure
    UT(it,iSpecies) = 3/2*nansum(p(:));
    UK(it,iSpecies) = nansum(pdyn(:));
    %imagesc(sim.xi,sim.zi,squeeze(p)')
    %drawnow
  end
end
toc

%% Thermal and kinetic energy of combined groups n08
sim = df08;
tic;
species_group = {[1 3],[2 4]};
clear UT UK
for it = 1:sim.length  
  sim_tmp = sim(it); 
  for iSpecies = 1:numel(species_group)
    disp(sprintf('it = %g/%g, sp = %g ',it,sim.length,iSpecies))
    
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp(species_group{iSpecies});
    %toc;
    pdyn = sim.mass(iSpecies)/sim.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
    p = (pxx+pyy+pzz)/3; % scalar pressure
    UT(it,iSpecies) = 3/2*nansum(p(:));
    UK(it,iSpecies) = nansum(pdyn(:));
    %imagesc(sim.xi,sim.zi,squeeze(p)')
    %drawnow
  end
end
toc

%% Thermal end kinetic energy of separate species n04
sim = df04;
tic;

for it = 1:sim.length  
  sim_tmp = sim(it); 
  for iSpecies = 1:numel(sim.mass)
    disp(sprintf('it = %g/%g, sp = %g ',it,sim.length,iSpecies))
    
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp(iSpecies);
    %toc;
    pdyn = sim.mass(iSpecies)/sim.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
    p = (pxx+pyy+pzz)/3; % scalar pressure
    UT(it,iSpecies) = 3/2*nansum(p(:));
    UK(it,iSpecies) = nansum(pdyn(:));
    %imagesc(sim.xi,sim.zi,squeeze(p)')
    %drawnow
  end
end
toc

%% Thermal and kinetic energy of combined groups n04
sim = df04;
tic;
clear UT UK
species_group = {[1 3 5],[3 5],[2 4 6],[4 6]};
for it = 1:sim.length  
  sim_tmp = sim(it); 
  for iGroup = 1:numel(species_group)
    disp(sprintf('it = %g/%g, grp = %g ',it,sim.length,iGroup))
    
    group = species_group{iGroup};
    mass_group = sim.mass(species_group{iGroup}(1));
    
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp(species_group{iGroup});
    %toc;
    
    
    pdyn = mass_group/sim.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
    p = (pxx+pyy+pzz)/3; % scalar pressure
    UT(it,iGroup) = 3/2*nansum(p(:));
    UK(it,iGroup) = nansum(pdyn(:));
    %imagesc(sim.xi,sim.zi,squeeze(p)')
    %drawnow
  end
end
toc

%%
sim = df04;
species_str = {'135','35','246','46'};
for iSpecies = 1:numel(species_str)
  %h5create(sim.file,sprintf('/scalar_timeseries/U/T/%s',species_str{iSpecies}), [1,sim.length]);
  h5write( sim.file,sprintf('/scalar_timeseries/U/T/%s',species_str{iSpecies}), UT(:,iSpecies)');
  %h5create(sim.file,sprintf('/scalar_timeseries/U/K/%s',species_str{iSpecies}), [1,sim.length]);
  h5write( sim.file,sprintf('/scalar_timeseries/U/K/%s',species_str{iSpecies}), UK(:,iSpecies)');
end
%%
sim = df08;
species_str = {'13','24'};
for iSpecies = 1:numel(species_str)
  %h5create(sim.file,sprintf('/scalar_timeseries/U/T/%s',species_str{iSpecies}), [1,sim.length]);
  h5write( sim.file,sprintf('/scalar_timeseries/U/T/%s',species_str{iSpecies}), UT(:,iSpecies)');
  %h5create(sim.file,sprintf('/scalar_timeseries/U/K/%s',species_str{iSpecies}), [1,sim.length]);
  h5write( sim.file,sprintf('/scalar_timeseries/U/K/%s',species_str{iSpecies}), UK(:,iSpecies)');
end

%%

for iSpecies = 1:numel(sim.mass)      
  %h5create(sim.file,sprintf('/scalar_timeseries/U/T/%g',iSpecies), [1,sim.length]);
  h5write( sim.file,sprintf('/scalar_timeseries/U/T/%g',iSpecies), UT(:,iSpecies)');
  %h5create(sim.file,sprintf('/scalar_timeseries/U/K/%g',iSpecies), [1,sim.length]);
  h5write( sim.file,sprintf('/scalar_timeseries/U/K/%g',iSpecies), UK(:,iSpecies)');
end
