% Loads data and saves the different variables individually, this way one
% can for example easily load any variable for all times.

%% Define times
timesteps = 00200:200:05000;%00200:200:08800;
times = timesteps/50;

timesteps = 05978:1:06000;%00200:200:08800;
times = timesteps/200;

ntimes = numel(timesteps);
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/';
data_dir = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data/';
data_dir_resave = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data_separated/';


savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_04/';
data_dir = '/Volumes/pic/in_progress/df_cold_protons_04/data/';
data_dir_resave = '/Volumes/pic/in_progress/df_cold_protons_04/data_separated/';

savedir_root = '/Users/cno062/Research/PIC/michael_run/';
data_dir = '/Volumes/Fountain/Data/PIC/michael_run/data/';
data_dir_resave = '/Volumes/Fountain/Data/PIC/michael_run/data_separated/';

savedir_root = '/Users/cno062/Research/PIC/michael_run/';
data_dir = '/Volumes/pic/finished_runs/turbulencerun/data/';
data_dir_resave = '/Volumes/pic/finished_runs/turbulencerun/data_separated/';

%% Variables to load
% only saved once
varstrs_same = {'x','z','dfac','teti','nnx','nnz','wpewce','mass','it','time','dt','xmax','zmax','q'};
vardir_same = [data_dir_resave 'same_for_all_times'];
%sim_info = load([vardir_same '/sim_info.mat']);

% saved for each time step
varstrs = {'A','E','B',...
        'ni1','ne1','ni2','ne2',...
        'vi1','ve1','vi2','ve2',...
        'ji1','je1','ji2','je2',...
        'pi1','pe1','pi2','pe2',...
        'ti1','te1','ti2','te2'...
        };
% saved for each time step
varstrs = {'A','E','B',...
        'ni1','ne1','ni2','ne2','ni3','ne3',...
        'vi1','ve1','vi2','ve2','vi3','ve3',...
        'ji1','je1','ji2','je2','ji3','je3'...
        };      
varstrs = {'E','B'...
        'ni1','ne1','ni2','ne2',...        
        'vi1','ve1','vi2','ve2'};
varstrs = {'vi1','ve1','vi2','ve2'};
%varstrs = {'A'...
%        };

varstrs = {'ni12'};

if strcmp(varstrs,'all');
  
end
nvars = numel(varstrs);

%% Loop over times, load data then do whatever
for ivar = 1:nvars
  vardir = [data_dir_resave varstrs{ivar}];
  disp(sprintf('Loading %s*.dat...',vardir))
  tic;
  isInitialized = 0;
  for itime = 1:ntimes
    %% Load data
    timestep = timesteps(itime);    
    %timestep = 09000;
    varstr_reload = sprintf('%s/%s-%05.0f.mat',vardir,varstrs{ivar},timestep);
    data_tmp  = load(varstr_reload,varstrs{ivar}); 
    data_tmp = eval(sprintf('data_tmp.%s',varstrs{ivar}));    
    if isnumeric(data_tmp)
      if not(isInitialized)
        data = zeros([ntimes size(data_tmp)]); 
        isInitialized = 1;
      end
      data(itime,:,:,:,:,:) = data_tmp;
    elseif isstruct(data_tmp)
      var_fields = fields(data_tmp);
      nfields = numel(var_fields);
      vec_fields = {'x';'y';'z'};
      tens_fields = {'xx';'xy';'xz';'yy';'yz';'zz'};
      if numel(var_fields)==3 && all(cellfun(@isequal,var_fields,vec_fields))
        datasize = size(data_tmp.x);
        if not(isInitialized)
          data = zeros([ntimes datasize numel(vec_fields)]); 
          isInitialized = 1;
        end
        for ifield = 1:nfields
          %data_field = eval([]);
          data(itime,:,:,ifield) = eval(['data_tmp.' var_fields{ifield}]);
        end
        %data(itime,:,:,) = data_tmp;
      elseif numel(var_fields)==6 && all(cellfun(@isequal,var_fields,tens_fields))        
      end
    end
  end
  eval(sprintf('%s = data;',varstrs{ivar}))
  toc
end
