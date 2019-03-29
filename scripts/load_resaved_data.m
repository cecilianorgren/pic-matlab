% Loads data and saves the different variables individually, this way one
% can for example easily load any variable for all times.

%% Define times
timesteps = 00200:200:08800;
times = timesteps/50;
ntimes = numel(timesteps);
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/';
data_dir = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data/';
data_dir_resave = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data_separated/';


%% Variables to load
% only save once
varstrs_same = {'x','z','dfac','teti','nnx','nnz','wpewce','mass','it','time','dt','xmax','zmax','q'};
vardir_same = [data_dir_resave 'same_for_all_times'];
if not(exist(vardir_same,'dir'))
  mkdir(vardir_same)
end
  
% save for each time step
varstrs = {'A','E','B',...
        'ni1','ne1','ni2','ne2',...
        'vi1','ve1','vi2','ve2',...
        'ji1','je1','ji2','je2',...
        'pi1','pe1','pi2','pe2',...
        'ti1','te1','ti2','te2'...
        };
nvars = numel(varstrs);
for ivar = 1:nvars
  vardir = [data_dir_resave varstrs{ivar}];
  if not(exist(vardir,'dir'))
    mkdir(vardir)
  end
end

%% Loop over times, load data then resave
for ivar = 1%:nvars
  vardir = [data_dir_resave varstrs{ivar}];
  disp(sprintf('Loading  %s...',vardir))
  for itime = 1:ntimes
    %% Load data
    timestep = timesteps(itime);    
    timestep = 09000;
    varstr_reload = sprintf('%s/%s-%05.0f.mat',vardir,varstrs{ivar},timestep);
    data_tmp  = load(varstr_reload,varstrs{ivar}); 
    data_tmp = eval('data_tmp.%s',varstrs{ivar});
    cat_dim = 4;
    data = cat(cat_dim,data,data_tmp);
    1;
  end
end
