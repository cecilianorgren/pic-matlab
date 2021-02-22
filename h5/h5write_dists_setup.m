timestep = 10000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_test/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_test/data_h5/dists.h5';
distIndRead = 463:667;
distIndRead = 668:918;
nSpecies = 6;
iteration = nobg.twpelim(timestep).iteration;
mass = [25 1 25 1 25 1];
charge = [1 -1 1 -1 1 -1];
tag = 'line4';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tag)
%%
timestep = 5000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/data_h5/dists.h5';
distIndRead = 131:243;
nSpecies = 6;
iteration = df04n.twpelim(timestep).iteration;
mass = [25  1 25 1 25 1];
charge = [1 -1 1 -1 1 -1];
tag = 'idr vertical';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tag)

%% 
timestep = 24000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 990:910;
tags = arrayfun(@(s)sprintf('A=%g',s),all_tags,'UniformOutput',false);
nSpecies = 6;
iteration = no02m.twpelim(timestep).iteration;
mass = [100 1 100 1 100 1];
charge = [1 -1 1 -1 1 -1];
tag = 'line vertical';
for ic = 1:2000, tags{ic} = tag; end
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tags)

%% 
timestep = 20000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 301:418;
tags = arrayfun(@(s)sprintf('A=%g',s),all_tags,'UniformOutput',false);
nSpecies = 6;
iteration = no02m.twpelim(timestep).iteration;
mass = [100 1 100 1 100 1];
charge = [1 -1 1 -1 1 -1];
%tag = 'idr vertical';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tags)
%% 
timestep = 24000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 1601:1738;1445:1600;1738;%1738
tags = arrayfun(@(s)sprintf('A=7.5',s),1:2000,'UniformOutput',false);
%tags = arrayfun(@(s)sprintf('',s),1:2000,'UniformOutput',false);
% for itag = 1:numel(all_tags)
%   tags{distIndRead(itag)} = sprintf('A=%g',all_tags(itag));
% end

%tags = arrayfun(@(s)sprintf('A=%g',s),all_tags,'UniformOutput',false);
nSpecies = 6;
iteration = no02m.twpelim(timestep).iteration;
mass = [100 1 100 1 100 1];
charge = [1 -1 1 -1 1 -1];
%tag = 'idr vertical';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tags)

%% 
timestep = 23000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 1:280;%1738
tags = arrayfun(@(s)sprintf('A=7.5',s),1:280,'UniformOutput',false);
%tags = arrayfun(@(s)sprintf('',s),1:2000,'UniformOutput',false);
% for itag = 1:numel(all_tags)
%   tags{distIndRead(itag)} = sprintf('A=%g',all_tags(itag));
% end

%tags = arrayfun(@(s)sprintf('A=%g',s),all_tags,'UniformOutput',false);
nSpecies = 6;
iteration = no02m.twpelim(timestep).iteration;
mass = [100 1 100 1 100 1];
charge = [1 -1 1 -1 1 -1];
%tag = 'idr vertical';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tags)

%% Add attributes
