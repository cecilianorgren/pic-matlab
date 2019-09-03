function out = fun_load_resaved_data(data_dir_resave,varstrs,timesteps)
% FUN_LOAD_RESAVED_DATA Loads single variable.
%   out = FUN_LOAD_RESAVED_DATA(data_dir_resave,varstr,timesteps)
%
%   Examples:
%   timesteps = 00200:200:05000;
%   data_dir_resave = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data_separated/';
%   ve1_ts = FUN_LOAD_RESAVED_DATA(data_dir_resave,'ve1',timesteps)
% 
%   timesteps = 00200:200:05000;
%   data_dir_resave = '/Volumes/pic/in_progress/df_cold_protons_04/data_separated/';
%   ve1_ts = FUN_LOAD_RESAVED_DATA(data_dir_resave,'ve1',timesteps)
% 
%   data_dir_resave = '/Volumes/Fountain/Data/PIC/michael_run/data_separated/';
%   timesteps = 05978:1:06000;
%   ve1_ts = FUN_LOAD_RESAVED_DATA(data_dir_resave,'ve1',timesteps)

%% Variables to load
% only saved once
varstrs_same = {'x','z','dfac','teti','nnx','nnz','wpewce','mass','it','time','dt','xmax','zmax','q'};
vardir_same = [data_dir_resave 'same_for_all_times'];
sim_info = load([vardir_same '/sim_info.mat']);

% saved for each time step
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
  out = data;  
  toc
end
