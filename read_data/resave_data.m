% Loads data and saves the different variables individually, this way one
% can for example easily load any variable for all times.

savedir_root = '/Users/cno062/Research/PIC/michael_run/';
data_dir = '/Volumes/Fountain/Data/PIC/michael_run/data/';
data_dir_resave = '/Volumes/Fountain/Data/PIC/michael_run/data_separated/';
data_dir = '/Volumes/pic/finished_runs/turbulencerun/data/';
data_dir_resave = '/Volumes/pic/finished_runs/turbulencerun/data_separated/';

%% Define times
if 0 % manually
  timesteps = 05978:1:06000;
  times = timesteps/200;
else % all in folder
  fileList = dir([data_dir '*fields*']);
  nFiles = numel(fileList);
  timesteps = zeros(nFiles,1);
  for iFile = 1:nFiles
    timesteps(iFile) = str2num(fileList(iFile).name(8:12));
  end
end
timesteps = 05200:200:06000;
ntimes = numel(timesteps);
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
varstrs = {'A','E','B',...
        'ni1','ne1','ni2','ne2',...
        'vi1','ve1','vi2','ve2',...
        };
varstrs = {'E','B',...
        'ni12','ne12',...
        'vi12','ve12',...
        };
nvars = numel(varstrs);
for ivar = 1:nvars
  vardir = [data_dir_resave varstrs{ivar}];
  if not(exist(vardir,'dir'))
    mkdir(vardir)
  end
end

%% Loop over times, load data then resave
nss = 4;
for itime = 1:ntimes % 141
  %% Load data
  timestep = timesteps(itime);
  disp(sprintf('timestep = %05.0f/%05.0f (%.0f/%.0f)',timestep,timesteps(end),itime,ntimes))
  txtfile = sprintf('%s/fields-%05.0f.dat',data_dir,timestep); % michael's perturbation
  disp('Loading data...')
  tic;     
  if 0
    [x,z,E,B,...
    dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q,...
    ni1,ne1,ni2,ne2,...
    vi1,ve1,vi2,ve2,...
    ji1,je1,ji2,je2,...
    pi1,pe1,pi2,pe2,...
    ti1,te1,ti2,te2...
    ] = read_fields_(txtfile,'nss',nss); 
  else
    txtfile = sprintf('%s/fields-%05.0f.dat',data_dir,timestep); % michael's perturbation
    all_data = read_fields_adaptive(txtfile,'nss',nss,'group',{[1 3],[2 4]});
    ndata = size(all_data,1);
    disp('Loaded: ')
    for idata = 1:ndata
      eval([all_data{idata,1} '= all_data{idata,2};' ]);
      fprintf('%s ',all_data{idata,1})
    end
    fprintf('\n') 
  end
  toc
  
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
