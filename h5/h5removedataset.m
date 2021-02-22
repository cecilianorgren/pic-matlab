filePath = sim.file;
fid = H5F.open(filePath,'H5F_ACC_RDWR','H5P_DEFAULT');
H5L.delete(fid,'/scalar_timeseries/U/T/1','H5P_DEFAULT');
H5F.close(fid);

%% Remove trajectories

filePath = '/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5';
fid = H5F.open(filePath,'H5F_ACC_RDWR','H5P_DEFAULT');
iTr0 = 410;
for iTr = 2:23  
  iTr_str = sprintf('%06.0f',iTr+iTr0);
  group_name = ['/traj/' iTr_str '/'];
  H5L.delete(fid,[group_name],'H5P_DEFAULT'); 
%   H5L.delete(fid,[group_name 'y'],'H5P_DEFAULT'); 
%   H5L.delete(fid,[group_name 'z'],'H5P_DEFAULT'); 
%   H5L.delete(fid,[group_name 'vx'],'H5P_DEFAULT'); 
%   H5L.delete(fid,[group_name 'vy'],'H5P_DEFAULT'); 
%   H5L.delete(fid,[group_name 'vz'],'H5P_DEFAULT'); 
%   H5L.delete(fid,[group_name 'Ex'],'H5P_DEFAULT'); 
%   H5L.delete(fid,[group_name 'Ey'],'H5P_DEFAULT'); 
%   H5L.delete(fid,[group_name 'Ez'],'H5P_DEFAULT'); 
%   H5L.delete(fid,[group_name 'Bx'],'H5P_DEFAULT'); 
%   H5L.delete(fid,[group_name 'By'],'H5P_DEFAULT'); 
%   H5L.delete(fid,[group_name 'Bz'],'H5P_DEFAULT'); 
end
H5F.close(fid);


%% Remove distributions
filePath = ds01.file;
fid = H5F.open(filePath,'H5F_ACC_RDWR','H5P_DEFAULT');
iDist = 276:283;
for iTr = iDist 
  iTr_str = sprintf('%05.0f',iTr);
  group_name = ['data/0000023000/' iTr_str '/'];
  H5L.delete(fid,[group_name],'H5P_DEFAULT');
end
H5F.close(fid);

%% Remove distributions
filePath = ds100.file;
fid = H5F.open(filePath,'H5F_ACC_RDWR','H5P_DEFAULT');
iDist = 1;
for iTr = iDist 
  iTr_str = sprintf('%05.0f',iTr);
  %group_name = ['data/0000046000/' iTr_str '/']; % single distribution
  group_name = ['data/0000046000/']; % entire group
  H5L.delete(fid,[group_name],'H5P_DEFAULT');
end
H5F.close(fid);