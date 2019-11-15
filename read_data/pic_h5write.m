% write pic fortran data ino h5 format
data_dir = '/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data/';
filePath = [data_dir 'fields.h5'];

datasize = size(E.x);
timesteps = 200:200:12000;
fields = {'Ex','Ey','Ez','Bx','By','Bz'};
for itime = 1:numel(timesteps)  
  for ifield = [1:numel(fields)]
    timestep = timesteps(itime);
    str_timestep = sprintf('%010.0f',timestep);
    h5create([data_dir 'fields.h5'],['/data/' str_timestep '/' fields{ifield}],datasize);
  end
end

%h5write([data_dir 'fields.h5'], ['/data/' str_timestep '/Ey'], E.y);
%% For comparison
filePath = '/Users/cno062/Data/SMILEI/GEMchallange/Fields0.h5';
info = h5info(filePath);

%% '/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data_separated/';
data_dir_h5 = '/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data_h5/';
data_dir_separated = '/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data_separated/';
filePath = [data_dir_h5 'fields.h5'];
dirs = dir(data_dir_separated);
datasets = dirs(find(not(contains({dirs.name},{'.','same_for_all_times'})))); % remove some dirs
datasets = {datasets.name};

%dirs = dir([data_dir_separated datasets{1}]);
%dataset_str = dirs(find(not(ismember({dirs.name},{'.','..'})))); % remove some dirs
%dataset_str = {dataset_str.name};
timesteps = 200:200:10800;

%datasize = [6400,1600];


for itime = 9%4:numel(timesteps)  
  for ifield = 23%1:numel(datasets)
    timestep = timesteps(itime);
    str_timestep = sprintf('%010.0f',timestep);
    
    % load data
    vardir = [data_dir_separated datasets{ifield}];
    varstr_reload = sprintf('%s/%s-%05.0f.mat',vardir,datasets{ifield},timestep);
    data_tmp  = load(varstr_reload,datasets{ifield}); 
    data_tmp = eval(sprintf('data_tmp.%s',datasets{ifield}));
    % check for components, i saved them as structures, with field components
    if isa(data_tmp,'numeric') % no field components
      disp(['/data/' str_timestep '/' datasets{ifield}])
      h5create([data_dir_h5 'fields.h5'],['/data/' str_timestep '/' datasets{ifield}], size(data_tmp));
      h5write([data_dir_h5 'fields.h5'], ['/data/' str_timestep '/' datasets{ifield}], data_tmp);
    elseif isa(data_tmp,'struct')
      components = fieldnames(data_tmp);
      for icomp = 1:numel(components)
        disp(['/data/' str_timestep '/' datasets{ifield} components{icomp}])
        datasize = size(subsref(data_tmp,substruct('.',components{icomp})));
        h5create([data_dir_h5 'fields.h5'],['/data/' str_timestep '/' datasets{ifield} components{icomp}], datasize);
        h5write([data_dir_h5 'fields.h5'], ['/data/' str_timestep '/' datasets{ifield} components{icomp}], subsref(data_tmp,substruct('.',components{icomp})));        
      end
    else      
      error('Unknown data type.')
    %  h5write([data_dir_h5 'fields.h5'], ['/data/' str_timestep '/' datasets{ifield}], data_tmp);
    end
  end
end

%h5write([data_dir 'fields.h5'], ['/data/' str_timestep '/Ey'], E.y);

%% '/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data/';
% Better to save original data, and only necessary quantities
data_dir    = '/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data/';
data_dir_h5 = '/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/';
filePath = [data_dir_h5 'fields.h5'];

timesteps = 000:200:12000;
% I forgot vxs,vys,vzs, time attribute, dt attribute, up until 48
for itime = 51:numel(timesteps)
  timestep = timesteps(itime);  
  txtfile = sprintf('%s/fields-%05.0f.dat',data_dir,timestep); 
  
  
  % read unnormalized data
  tic; [varstrs,vars] = read_data_no_normalization(txtfile,6); toc
  nss = numel(vars{find(contains(varstrs,'mass'))}); % number of species
  nnx = vars{find(contains(varstrs,'nnx'))}; % number of grid points in x
  nnz = vars{find(contains(varstrs,'nnz'))}; % number of grid points in z
  
  % Load some vars and remove from vars and varstrs
  iter = vars{find(contains(varstrs,'it'))};   vars(find(contains(varstrs,'it'))) = [];   varstrs(find(contains(varstrs,'it'))) = [];
  time = vars{find(contains(varstrs,'time'))}; vars(find(contains(varstrs,'time'))) = []; varstrs(find(contains(varstrs,'time'))) = [];
  dt = vars{find(contains(varstrs,'dt'))}; vars(find(contains(varstrs,'dt'))) = []; varstrs(find(contains(varstrs,'dt'))) = [];
  mass = vars{find(contains(varstrs,'mass'))}; vars(find(contains(varstrs,'mass'))) = []; varstrs(find(contains(varstrs,'mass'))) = [];
  q = vars{find(contains(varstrs,'q'))}; vars(find(contains(varstrs,'q'))) = []; varstrs(find(contains(varstrs,'q'))) = [];
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
    if itime == 1 && not(sum(ismember(size(data),[nnx nnz]))>=2)
      %continue % implement later
      disp(['/simulation_information/' varstrs{ivar}])
      h5create([data_dir_h5 'fields.h5'],['/simulation_information/' varstrs{ivar}], size(data));
      h5write([data_dir_h5 'fields.h5'],['/simulation_information/' varstrs{ivar}], data);
      continue % jump to next variable
    end
    
    % From here, we only read data field variables
    % Check if it is a field (E,B) or species data (n,j,vv)
    if ndims(data) == 2 && all(size(data) == [nnx nnz]) % is (E,B)
      dataset_name = ['/data/' str_iteration '/' varstrs{ivar}];
      disp(dataset_name)
      h5create([data_dir_h5 'fields.h5'], dataset_name, size(data));
      h5write( [data_dir_h5 'fields.h5'], dataset_name, data);
    elseif size(data,3) == nss % is (n,vs,vv)
      for iSpecies = 1:nss
        data_tmp = data(:,:,nss);
        dataset_name = ['/data/' str_iteration '/' varstrs{ivar} '_' num2str(iSpecies)];
        disp(dataset_name)
        h5create([data_dir_h5 'fields.h5'], dataset_name, size(data_tmp));
        h5write( [data_dir_h5 'fields.h5'], dataset_name, data_tmp);
        % Also write species data as attributes
        h5writeatt([data_dir_h5 'fields.h5'],dataset_name, 'mass',mass(iSpecies)) 
        h5writeatt([data_dir_h5 'fields.h5'],dataset_name, 'charge',q(iSpecies)) 
        h5writeatt([data_dir_h5 'fields.h5'],dataset_name, 'dfac',dfac(iSpecies)) 
        % info.Groups(1).Groups(ig).Datasets(id).Attributes.Name
      end
    else
      warning(sprintf('Variable: %s skipped',varstrs{ivar}))
      continue      
    end   
    % Write attributes for group (iteration)
    h5writeatt([data_dir_h5 'fields.h5'],['/data/' str_iteration '/'], 'time',time)
    h5writeatt([data_dir_h5 'fields.h5'],['/data/' str_iteration '/'], 'dt',dt)
    % info.Groups(1).Groups(ig).Attributes.Name
  end
end

if 0 % Usage examples
  info = h5info(filePath);
  {info.Groups(2).Datasets.Name} % simulation information
  x = h5read(filePath,'/simulation_information/xe');
end


%% '/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data/';
% Better to save original data, and only necessary quantities
data_dir    = '/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data/';
data_dir_h5 = '/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data_h5/';
filePath = [data_dir_h5 'fields.h5'];

timesteps = 200:200:11000;
% I forgot vxs,vys,vzs, time attribute, dt attribute, up until 48
for itime = 1:numel(timesteps)
  timestep = timesteps(itime);  
  txtfile = sprintf('%s/fields-%05.0f.dat',data_dir,timestep); 
  
  
  % read unnormalized data
  tic; [varstrs,vars] = read_data_no_normalization(txtfile,4); toc
  nss = numel(vars{find(contains(varstrs,'mass'))}); % number of species
  nnx = vars{find(contains(varstrs,'nnx'))}; % number of grid points in x
  nnz = vars{find(contains(varstrs,'nnz'))}; % number of grid points in z
  
  % Load some vars and remove from vars and varstrs
  iter = vars{find(contains(varstrs,'it'))};   vars(find(contains(varstrs,'it'))) = [];   varstrs(find(contains(varstrs,'it'))) = [];
  time = vars{find(contains(varstrs,'time'))}; vars(find(contains(varstrs,'time'))) = []; varstrs(find(contains(varstrs,'time'))) = [];
  dt = vars{find(contains(varstrs,'dt'))}; vars(find(contains(varstrs,'dt'))) = []; varstrs(find(contains(varstrs,'dt'))) = [];
  mass = vars{find(contains(varstrs,'mass'))}; vars(find(contains(varstrs,'mass'))) = []; varstrs(find(contains(varstrs,'mass'))) = [];
  q = vars{find(contains(varstrs,'q'))}; vars(find(contains(varstrs,'q'))) = []; varstrs(find(contains(varstrs,'q'))) = [];
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
    if itime == 1 && not(sum(ismember(size(data),[nnx nnz]))>=2)
      %continue % implement later
      disp(['/simulation_information/' varstrs{ivar}])
      h5create([data_dir_h5 'fields.h5'],['/simulation_information/' varstrs{ivar}], size(data));
      h5write([data_dir_h5 'fields.h5'],['/simulation_information/' varstrs{ivar}], data);
      continue % jump to next variable
    end
    
    % From here, we only read data field variables
    % Check if it is a field (E,B) or species data (n,j,vv)
    if ndims(data) == 2 && all(size(data) == [nnx nnz]) % is (E,B)
      dataset_name = ['/data/' str_iteration '/' varstrs{ivar}];
      disp(dataset_name)
      h5create([data_dir_h5 'fields.h5'], dataset_name, size(data));
      h5write( [data_dir_h5 'fields.h5'], dataset_name, data);
    elseif size(data,3) == nss % is (n,vs,vv)
      for iSpecies = 1:nss
        data_tmp = data(:,:,nss);
        dataset_name = ['/data/' str_iteration '/' varstrs{ivar} '_' num2str(iSpecies)];
        disp(dataset_name)
        h5create([data_dir_h5 'fields.h5'], dataset_name, size(data_tmp));
        h5write( [data_dir_h5 'fields.h5'], dataset_name, data_tmp);
        % Also write species data as attributes
        h5writeatt([data_dir_h5 'fields.h5'],dataset_name, 'mass',mass(iSpecies)) 
        h5writeatt([data_dir_h5 'fields.h5'],dataset_name, 'charge',q(iSpecies)) 
        h5writeatt([data_dir_h5 'fields.h5'],dataset_name, 'dfac',dfac(iSpecies)) 
        % info.Groups(1).Groups(ig).Datasets(id).Attributes.Name
      end
    else
      warning(sprintf('Variable: %s skipped',varstrs{ivar}))
      continue      
    end   
    % Write attributes for group (iteration)
    h5writeatt([data_dir_h5 'fields.h5'],['/data/' str_iteration '/'], 'time',time)
    h5writeatt([data_dir_h5 'fields.h5'],['/data/' str_iteration '/'], 'dt',dt)
    % info.Groups(1).Groups(ig).Attributes.Name
  end
end

if 0 % Usage examples
  info = h5info(filePath);
  {info.Groups(2).Datasets.Name} % simulation information
  x = h5read(filePath,'/simulation_information/xe');
end
