function h5write_fields(data_dir,filePath,timesteps,nSpecies)
% H5WRITE_FIELDS Write Michael's simulation output data to h5 file.
% H5WRITE_FIELDS(dirData,h5FilePath,timesteps,nSpecies)
% 
% dirData - directory of data
% h5FilePath - directory and file name
% timesteps - timesteps to resave, in units of wpe^-1
%             if empty [], check which files are in the data directory and 
%             use those
% nSpecies - number of species: required in order to read the data right
%
% Examples:
%   % Mounted tesla from tiny laptop to Fountain
%   h5write_fields('/Users/cecilia/tesla/cno062/df_cold_protons_n04/data/','/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5',7000:50:8000,6)
%

h5exist = 0;
if exist(filePath,'file')
  h5exist = 1;
  disp(sprintf('File %s exists. Loading file to obtain existing times.',filePath))
  pic = PIC(filePath);
else
  fid = H5F.create(filePath);
  h5writeatt(filePath,'/','software','micPIC')
end

fileList = dir([data_dir '*fields*']);
nFiles = numel(fileList);
if isempty(timesteps)
  timesteps = zeros(nFiles,1);
  for iFile = 1:nFiles
    timesteps(iFile) = str2num(fileList(iFile).name(8:12)); % wpe
  end
end

for itime = 1:numel(timesteps)
  timestep = timesteps(itime);  
  if h5exist && not(isempty(find(timestep==pic.twpe)))
    disp(sprintf('twpe = %g already exists, skipping.',timestep))
    continue
  end
  txtfile = sprintf('%s/fields-%05.0f.dat',data_dir,timestep); 
  disp(sprintf('Reading twpe = %g.',timestep))
  
  
  % read unnormalized data
  if not(exist(txtfile,'file'))
    warning(sprintf('File %s does not exist.',txtfile))
  end
    
  tic; [varstrs,vars] = read_data_no_normalization(txtfile,nSpecies); toc  
  
  nss = numel(vars{find(contains(varstrs,'mass'))}); % number of species
  nnx = vars{find(contains(varstrs,'nnx'))}; % number of grid points in x
  nnz = vars{find(contains(varstrs,'nnz'))}; % number of grid points in z
  
  % Load some vars and remove from vars and varstrs
  iter = vars{find(contains(varstrs,'it'))};   vars(find(contains(varstrs,'it'))) = [];   varstrs(find(contains(varstrs,'it'))) = [];
  time = vars{find(contains(varstrs,'time'))}; vars(find(contains(varstrs,'time'))) = []; varstrs(find(contains(varstrs,'time'))) = [];
  dt   = vars{find(contains(varstrs,'dt'))};   vars(find(contains(varstrs,'dt'))) = [];   varstrs(find(contains(varstrs,'dt'))) = [];
  mass = vars{find(contains(varstrs,'mass'))}; vars(find(contains(varstrs,'mass'))) = []; varstrs(find(contains(varstrs,'mass'))) = [];
  q    = vars{find(contains(varstrs,'q'))};    vars(find(contains(varstrs,'q'))) = [];    varstrs(find(contains(varstrs,'q'))) = [];
  dfac = vars{find(contains(varstrs,'dfac'))}; vars(find(contains(varstrs,'dfac'))) = []; varstrs(find(contains(varstrs,'dfac'))) = [];
  str_iteration = sprintf('%010.0f',iter); % same group format as SMILEI
  
  disp(['time = ' num2str(time)])
  
  % Remaining non-datasize matrices are the same for each time, so only
  % save one time
    
  % loop through variables, and save to h5 file in 'filePath'
  %%
  tic
  nvars = numel(vars);  
  for ivar = 1:nvars
    data = vars{ivar};
  %  varstrs{ivar}
    
    % check data size, and see if it needs splitting up, and where to save
    if not(h5exist) && itime == 1 && not(sum(ismember(size(data),[nnx nnz]))>=2)
      %continue % implement later
      disp(['/simulation_information/' varstrs{ivar}])
      h5create(filePath,['/simulation_information/' varstrs{ivar}], size(data),'ChunkSize',[1 1]);
      h5write(filePath,['/simulation_information/' varstrs{ivar}], data);
      continue % jump to next variable
    end
    
    % From here, we only read data field variables
    % Check if it is a field (E,B) or species data (n,j,vv)
    if ndims(data) == 2 && all(size(data) == [nnx nnz]) % is (E,B)
      dataset_name = ['/data/' str_iteration '/' varstrs{ivar}];
      disp(dataset_name)
      h5create(filePath, dataset_name, size(data),'ChunkSize',[50 50]);
      h5write( filePath, dataset_name, data);
    elseif size(data,3) == nss % is (n,vs,vv)
      for iSpecies = 1:nss
        data_tmp = data(:,:,iSpecies);
        dataset_name = ['/data/' str_iteration '/' varstrs{ivar} '/' num2str(iSpecies)];
        disp(dataset_name)
        h5create(filePath, dataset_name, size(data_tmp),'ChunkSize',[50 50]);
        h5write( filePath, dataset_name, data_tmp);
        % Also write species data as attributes
        h5writeatt(filePath,dataset_name, 'mass',mass(iSpecies)) 
        h5writeatt(filePath,dataset_name, 'charge',q(iSpecies)) 
        h5writeatt(filePath,dataset_name, 'dfac',dfac(iSpecies)) 
        % info.Groups(1).Groups(ig).Datasets(id).Attributes.Name
      end
    else
      warning(sprintf('Variable: %s skipped',varstrs{ivar}))
      continue      
    end   
    % Write attributes for group (iteration)
    h5writeatt(filePath,['/data/' str_iteration '/'], 'time',time)
    h5writeatt(filePath,['/data/' str_iteration '/'], 'dt',dt)
    % info.Groups(1).Groups(ig).Attributes.Name
  end
  toc

  % Write attributes
  % UB
  bx = vars{find(strcmp(varstrs,'bx'))};
  by = vars{find(strcmp(varstrs,'by'))};
  bz = vars{find(strcmp(varstrs,'bz'))};
  UB = (bx.^2+by.^2+bz.^2)/2;
  UB = sum(UB(:));
  h5writeatt(filePath,['/data/' str_iteration '/'],'UB',UB)

  % UT, UK
  if 0
  n = vars{find(strcmp(varstrs,'dns'))};
  jx = vars{find(strcmp(varstrs,'jx'))};
  jy = vars{find(strcmp(varstrs,'jy'))};
  jz = vars{find(strcmp(varstrs,'jz'))};
  pxx = vars{find(strcmp(varstrs,'pxx'))};
  pyy = vars{find(strcmp(varstrs,'pyy'))};
  pzz = vars{find(strcmp(varstrs,'pzz'))};

  for iSpecies = 1:nss
    pdyn = mass(iSpecies)/mass(1)*0.5*(jx(:,:,iSpecies).^2 + jy(:,:,iSpecies).^2 + jz(:,:,iSpecies).^2)./n(:,:,iSpecies);
    p = (pxx(:,:,iSpecies)+pyy(:,:,iSpecies)+pzz(:,:,iSpecies))/3; % scalar pressure
    UT(1,iSpecies) = 3/2*nansum(p(:));
    UK(1,iSpecies) = nansum(pdyn(:));
  end
  h5writeatt(filePath,['/data/' str_iteration '/'],'UK',UK)
  h5writeatt(filePath,['/data/' str_iteration '/'],'UT',UT)
  end
end

disp('Done.')