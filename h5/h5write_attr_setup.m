%df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
pic = no02m;
missingAttr = pic.get_missing_attributes('RE');
%pic = pic.twcilim(pic.twci(missingAttr(1:end-2)),'exact');
%pic = pic.twcilim(pic.twci(missingAttr),'exact');
pic = pic.twcilim(pic.twci([missingAttr(1)-2:missingAttr(1)+2]),'exact');
%pic = no02m;
%pic = pic.twpelim(200:200:6000);
%pic = turb.twcilim(1:1:30,'exact');

%% Energy partitioning, UB, UK, UT
times = pic.twci;
for it = 1:pic.nt
  pic_tmp = pic.twcilim(times(it));
  Bx = pic_tmp.Bx;
  By = pic_tmp.By;
  Bz = pic_tmp.Bz;
  Babs = sqrt(Bx.^2+By.^2+Bz.^2);
  UB(it) = sum(Babs(:).^2)/2;  
  disp(sprintf('%g/%g',it,pic.nt))
end

%% Energy partitioning, UB, UK, UT
times = pic.twci;
%% Thermal and kinetic energy n08
sim = no02m;
tic;
clear UT UK
for it = 1 %sim.twpelim(2000:1000:3000,'exact').indices%:sim.length  
  sim_tmp = sim(it); 
  for iSpecies = 1:numel(sim.mass)    
    disp(sprintf('it = %g/%g, sp = %g ',it,sim.length,iSpecies))
    
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp(iSpecies);
    %toc;
    pdyn = sim.mass(iSpecies)/sim.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
    p = (pxx+pyy+pzz)/3; % scalar pressure
    UT(it,iSpecies) = 3/2*nansum(p(:));
    UK(it,iSpecies) = nansum(pdyn(:));
    imagesc(sim.xi,sim.zi,squeeze(p)')
    drawnow
  end
  h5write_attr(sim_tmp,sim_tmp.twci,'UK',UK(it,:))
  h5write_attr(sim_tmp,sim_tmp.twci,'UT',UT(it,:))
  toc
end

%% Thermal and kinetic energy n02m, cold
sim = no02m;
tic;
clear UT UK
species_groups = {[3 5],[4 6]};
species_groups_str = {'cold_ions','cold_electrons'};
for it = sim.twpelim(7000:1000:14000,'exact').indices%:sim.length  
  sim_tmp = sim(it); 
  for iSpecies_group = 1:numel(species_groups)
    iSpecies = species_groups{iSpecies_group};
    species_str = species_groups_str{iSpecies_group};
    disp(sprintf('it = %g/%g, sp = %g ',it,sim.length,iSpecies))
    
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp(iSpecies);
    %toc;
    pdyn = sim.mass(iSpecies(1))/sim.mass(1)*0.5*(jx.^2 + jy.^2 + jz.^2)./n;
    p = (pxx+pyy+pzz)/3; % scalar pressure
    UT(it,iSpecies_group) = 3/2*nansum(p(:));
    UK(it,iSpecies_group) = nansum(pdyn(:));
    imagesc(sim.xi,sim.zi,squeeze(p)')
    drawnow
    h5write_attr(sim_tmp,sim_tmp.twci,['UK_' species_str],UK(it,iSpecies_group))
    h5write_attr(sim_tmp,sim_tmp.twci,['UT_' species_str],UT(it,iSpecies_group))
  end
  toc
end

%% Fix att
for it = sim.twpelim(15000:1000:25000,'exact').indices
  pic_tmp = sim(it);
  str_iteration = sprintf('%010.0f',pic_tmp.iteration); % same group format as SMILEI      
  dataset_name = ['/data/' str_iteration '/'];
  
  % Load data
  UK_cold_ions = h5readatt(no02m.file,dataset_name,'UK_cold_ions');
  UK_cold_electrons = h5readatt(no02m.file,dataset_name,'UK_cold_electrons');
  UT_cold_ions = h5readatt(no02m.file,dataset_name,'UT_cold_ions');
  UT_cold_electrons = h5readatt(no02m.file,dataset_name,'UT_cold_electrons');
  
  attrstr = {'UK_cold_ions','UK_cold_electrons','UT_cold_ions','UT_cold_electrons'};
  
  for iatt = 1:numel(attrstr)
    % Remove attribute
    h5file        =  pic_tmp.file;
    location      =  dataset_name;
    attributeName = attrstr{iatt};

    % Open the file (ensure to close it automatically when done)
    fileID = H5F.open(h5file,'H5F_ACC_RDWR','H5P_DEFAULT');
    fileIDCleanUp = onCleanup(@()H5F.close(fileID));
    % Open the dataset/group
    locID  = H5O.open(fileID, location,'H5P_DEFAULT');
    locIDCleanUp = onCleanup(@()H5O.close(locID));
    try %to open the attribute.
       attID = H5A.open(locID, attributeName, 'H5P_DEFAULT');
       H5A.close(attID);
       H5A.delete(locID, attributeName);
    catch ALL
        % do nothing if the attribute does not exist.
    end
    
  end
  
  h5writeatt(no02m.file,dataset_name,'UK_cold_ions',UK_cold_ions(1));  
  h5writeatt(no02m.file,dataset_name,'UK_cold_electrons',UK_cold_electrons(1));
  h5writeatt(no02m.file,dataset_name,'UT_cold_ions',UT_cold_ions(1));  
  h5writeatt(no02m.file,dataset_name,'UT_cold_electrons',UT_cold_electrons(1));
end

%% Rem att
for it = 1:14
  pic_tmp = sim(it);
  str_iteration = sprintf('%010.0f',pic_tmp.iteration); % same group format as SMILEI      
  dataset_name = ['/data/' str_iteration '/'];
  
    
  attrstr = {'UK','UT'};
  
  for iatt = 1:numel(attrstr)
    % Remove attribute
    h5file        =  pic_tmp.file;
    location      =  dataset_name;
    attributeName = attrstr{iatt};

    % Open the file (ensure to close it automatically when done)
    fileID = H5F.open(h5file,'H5F_ACC_RDWR','H5P_DEFAULT');
    fileIDCleanUp = onCleanup(@()H5F.close(fileID));
    % Open the dataset/group
    locID  = H5O.open(fileID, location,'H5P_DEFAULT');
    locIDCleanUp = onCleanup(@()H5O.close(locID));
    try %to open the attribute.
       attID = H5A.open(locID, attributeName, 'H5P_DEFAULT');
       H5A.close(attID);
       H5A.delete(locID, attributeName);
    catch ALL
        % do nothing if the attribute does not exist.
    end
    
  end
    
end

%% X line position, Ey at X line, A at X line

xXlineAll = [];
zXlineAll = [];
EyXlineAll = [];
xXline = [];
zXline = [];
Aval = [];
Aall = nan(pic.nx,pic.nz,pic.nt);
plot(nan,nan); hold on
times = pic.twci;

for it = 1:pic.nt
  pic_tmp = pic.twcilim(times(it));
  if 0
    Bx = pic_tmp.Bx;
    Bz = pic_tmp.Bz;
    if 1 % also magnetic pressure
      By = pic_tmp.By;
      Babs = sqrt(Bx.^2+By.^2+Bz.^2);
      UB(it) = sum(Babs(:).^2)/2;  
    end
    A = vector_potential(pic_tmp.xi,pic_tmp.zi,Bx,Bz);
    try
    h5write_fields_ancillary(pic_tmp,pic_tmp.twpe,'A',A)  
    catch
      disp('Could not write A.')
    end
  else
    A = pic_tmp.A;
  end
  Aall(:,:,it) = A;
  
  [Ainds,Avals] = saddle(A,'sort');
  xXline(it) = pic_tmp.xi(Ainds(1,1));
  zXline(it) = pic_tmp.zi(Ainds(1,2));  
  Aval(it) = Avals(1);
  EyXline(it) = mean(mean(pic_tmp.xlim(xXline(it)+[-0.1 0.1]).zlim(zXline(it)+[-0.1 0.1]).Ey));
  for iX = 1:numel(Avals)
    xXlineAll = [xXlineAll pic_tmp.xi(Ainds(iX,1))];
    zXlineAll = [zXlineAll pic_tmp.zi(Ainds(iX,2))];
    EyXlineAll = [EyXlineAll mean(mean(pic_tmp.xlim(xXlineAll(end)+[-0.1 0.1]).zlim(zXlineAll(end)+[-0.1 0.1]).Ey))];
    scatter(pic_tmp.twci,xXlineAll(end),abs(EyXlineAll(end))*100+1,iX)
    drawnow
  end
  disp(sprintf('%g',it))
end

dA = diff(Aval(1:pic.nt));
dt = diff(times(1:pic.nt));
dAdt_ = dA./dt;
dAdt = interp1(times(1:end-1)+0.5*dt,dAdt_,times);

% for it = 1:pic.nt
%   pic_tmp = pic.twcilim(times(it));
%   EyXline05(it) = mean(mean(pic_tmp.xlim(xXline(it)+0.5*[-1 1]).zlim(zXline(it)+0.5*[-1 1]).Ey));
% end

%% UB
times = pic.twci;
clear UB
for it = 1:pic.nt
  pic_tmp = pic.twcilim(times(it));

  Bx = pic_tmp.Bx;
  Bz = pic_tmp.Bz;
  By = pic_tmp.By;
  Babs = sqrt(Bx.^2+By.^2+Bz.^2);
  UB(it) = sum(Babs(:).^2)/2;  
  disp(sprintf('%g',it))
end


%% Write attributes
indsave = 1:(pic.nt-3);
indsave = (pic.nt-2):pic.nt;
indsave = 1:pic.nt;
indsave = 2:(pic.nt-1);
%indsave = 1:pic.nt;
h5write_attr(pic.subset('t',indsave),times(indsave),'RE',EyXline(indsave))
%h5write_attr(pic.subset('t',indsave),times(indsave),'RA',dAdt(indsave))
h5write_attr(pic.subset('t',indsave),times(indsave),'UB',UB(indsave))
h5write_attr(pic.subset('t',indsave),times(indsave),'xline_position',[xXline(indsave)' zXline(indsave)'])

%% Write ancillary data (not attributes), for example A
pic = no02m;
timesteps = pic.twpelim([15100:100:15900 16100:100:16900],'exact').twpe;

for time = timesteps
  pic_tmp = pic.twpelim(time);
  Bx = pic_tmp.Bx;
  Bz = pic_tmp.Bz;
  A = vector_potential(pic_tmp.xi,pic_tmp.zi,Bx,Bz);
  imagesc(A'); colorbar; pause(0.1)
  h5write_fields_ancillary(pic_tmp,pic_tmp.twpe,'A',A)  
end