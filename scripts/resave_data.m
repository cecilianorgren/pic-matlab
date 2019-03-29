% Loads data and saves the different variables individually, this way one
% can for example easily load any variable for all times.

%% Define times
timesteps = 00200:200:10800;
times = timesteps/50;
ntimes = numel(timesteps);
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/';
data_dir = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data/';
data_dir_resave = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data_separated/';


%% Variables to save
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
for itime = 2:ntimes
  %% Load data
  timestep = timesteps(itime);
  disp(sprintf('timestep = %05.0f/%05.0f',timestep,timesteps(end)))
  txtfile = sprintf('%s/fields-%05.0f.dat',data_dir,timestep); % michael's perturbation
  disp('Loading data...')
  tic; [x,z,E,B,...
        ni1,ne1,ni2,ne2,...
        vi1,ve1,vi2,ve2,...
        ji1,je1,ji2,je2,...
        pi1,pe1,pi2,pe2,...
        ti1,te1,ti2,te2,...
        dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q]... 
        = read_fields_(txtfile); toc
  A = vector_potential(x,z,B.x,B.z); % vector potential  
  
  disp('Resaving data into separate files... \n')
  tic
  for ivar = 1:nvars    
    fprintf('%s ',varstrs{ivar})
    vardir = [data_dir_resave varstrs{ivar}];
    varstr_resave = sprintf('%s/%s-%05.0f',vardir,varstrs{ivar},timestep);
    save(varstr_resave,varstrs{ivar})    
  end
  toc
end
%%
if 1 % same basic simulation data that is the same for all runs
  vardir = [data_dir_resave vardir_same];
  varstr_resave = sprintf('%s/%s',vardir_same,'sim_info');
  save(varstr_resave,varstrs_same{:})
end
