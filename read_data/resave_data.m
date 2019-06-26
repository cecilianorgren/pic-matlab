% Loads data and saves the different variables individually, this way one
% can for example easily load any variable for all times.

%% Define times
timesteps = 00200:200:12000;
times = timesteps/50;
ntimes = numel(timesteps);
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/';
data_dir = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data/';
data_dir_resave = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data_separated/';
%data_dir = '/Volumes/pic/in_progress/df_cold_protons_04/data/';
%data_dir_resave = '/Volumes/pic/in_progress/df_cold_protons_04/data_separated/';


%% Variables to save
% only save once
varstrs_same = {'x','z','dfac','teti','nnx','nnz','wpewce','mass','it','time','dt','xmax','zmax','q'};
vardir_same = [data_dir_resave 'same_for_all_times'];
if not(exist(vardir_same,'dir'))
  mkdir(vardir_same)
end
  
% save for each time step
varstrs = {'A','E','B',...
        'ni1','ne1','ni2','ne2','ni3','ne3',...
        'vi1','ve1','vi2','ve2','vi3','ve3',...
        'ji1','je1','ji2','je2','ji3','je3',...
        'pi1','pe1','pi2','pe2','pi3','pe3',...
        'ti1','te1','ti2','te2','ti3','te3'...
        };
      varstrs = {};
nvars = numel(varstrs);
for ivar = 1:nvars
  vardir = [data_dir_resave varstrs{ivar}];
  if not(exist(vardir,'dir'))
    mkdir(vardir)
  end
end

%% Loop over times, load data then resave
nss = 4;
for itime = 1
  %% Load data
  timestep = timesteps(itime);
  disp(sprintf('timestep = %05.0f/%05.0f',timestep,timesteps(end)))
  txtfile = sprintf('%s/fields-%05.0f.dat',data_dir,timestep); % michael's perturbation
  disp('Loading data...')
  tic;     
    [x,z,E,B,...
    dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q,...
    ni1,ne1,ni2,ne2,...
    vi1,ve1,vi2,ve2,...
    ji1,je1,ji2,je2,...
    pi1,pe1,pi2,pe2,...
    ti1,te1,ti2,te2...
    ] = read_fields_(txtfile,'nss',nss); 
  toc
        %= read_fields_(txtfile); toc
  A = vector_potential(x,z,B.x,B.z); % vector potential  
  
  disp('Resaving data into separate files...')
  % nss = 6, 200s/timestep
  tic
  for ivar = 1:nvars    
    fprintf('%s ',varstrs{ivar})
    vardir = [data_dir_resave varstrs{ivar}];
    varstr_resave = sprintf('%s/%s-%05.0f',vardir,varstrs{ivar},timestep);
    save(varstr_resave,varstrs{ivar})    
  end
  fprintf('\n')
  toc
end
%%
if 1 % same basic simulation data that is the same for all runs
  vardir = [data_dir_resave vardir_same];
  varstr_resave = sprintf('%s/%s',vardir_same,'sim_info');
  save(varstr_resave,varstrs_same{:})
end
