function fun_pic_h5write(data_dir,filePath,timesteps,nSpecies)
% fun_pic_h5write(data_dir,filePath,timesteps)
% data_dir - directory of data
% h5 file path - directory and file name
% timesteps - timesteps to resave, in units of wpe^-1
%             if empty [], check which files are in the data directory and 
%             use those
h5exist = 0;
if exist(filePath,'file')
  h5exist = 1;
  pic = PIC(filePath);
else
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
  
  
  % read unnormalized data
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
    
  
  % Remaining non-datasize matrices are the same for each time, so only
  % save one time
    
  % loop through variables, and save to h5 file in 'filePath'
  %%
  nvars = numel(vars);  
  for ivar = 1:nvars
    data = vars{ivar};
    
    % check data size, and see if it needs splitting up, and where to save
    if not(h5exist) && itime == 1 && not(sum(ismember(size(data),[nnx nnz]))>=2)
      %continue % implement later
      disp(['/simulation_information/' varstrs{ivar}])
      h5create(filePath,['/simulation_information/' varstrs{ivar}], size(data));
      h5write(filePath,['/simulation_information/' varstrs{ivar}], data);
      continue % jump to next variable
    end
    
    % From here, we only read data field variables
    % Check if it is a field (E,B) or species data (n,j,vv)
    if ndims(data) == 2 && all(size(data) == [nnx nnz]) % is (E,B)
      dataset_name = ['/data/' str_iteration '/' varstrs{ivar}];
      disp(dataset_name)
      h5create(filePath, dataset_name, size(data));
      h5write( filePath, dataset_name, data);
    elseif size(data,3) == nss % is (n,vs,vv)
      for iSpecies = 1:nss
        data_tmp = data(:,:,iSpecies);
        dataset_name = ['/data/' str_iteration '/' varstrs{ivar} '/' num2str(iSpecies)];
        disp(dataset_name)
        h5create(filePath, dataset_name, size(data_tmp));
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
end

disp('Done.')