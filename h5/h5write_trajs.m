function h5write_trajs(filePath,tr_arr,fpeaks)
% H5WRITE_TRAJS Write trajectories to h5 file.
% H5WRITE_TRAJS(h5FilePath,structure_trajectories)
% 
% dirData - directory of data
% h5FilePath - directory and file name
% structure_trajectories - structure array with needed trajectory
%     information

h5exist = 0;
if exist(filePath,'file')
  h5exist = 1;
  disp(sprintf('File %s exists. Loading file to obtain existing times.',filePath))
  traj = PICTraj(filePath);
else
  % do i need to initialize something? dont remember
  %h5writeatt(filePath,'/','software','micPIC')
end

doAppend = 1; % for future

m = 1; % hardcoded for now, need to change
q = 1;

iTr0 = traj.ntr; % add number to append correctly.

for iTr = 1:numel(tr_arr)
  % TO IMPLEMENT: check if the trajectory is already there
  
  [iPeak,iDist] = ind2sub(size(tr_arr),iTr);
  % Get trajectory structure from array
  tr = tr_arr(iTr);
  
  % Get fields at location
  Ex = tr.Ex;
  Ey = tr.Ey;
  Ez = tr.Ez;
  Bx = tr.Bx;
  By = tr.By;
  Bz = tr.Bz;
  %vxB = cross_product(tr.vx,tr.vy,tr.vz,Bx,By,Bz,'components',1); % calculate force
  
  % Information about how initial r0, v0 were chosen  
  fpeaks_info.iteration = fpeaks(iTr).timeIteration; % need to fix time for distributions class
  fpeaks_info.x = (fpeaks(iTr).dist_x1+fpeaks(iTr).dist_x2)/2;
  fpeaks_info.z = (fpeaks(iTr).dist_z1+fpeaks(iTr).dist_z2)/2;
  fpeaks_info.x1 = fpeaks(iTr).dist_x1;
  fpeaks_info.x2 = fpeaks(iTr).dist_x2;
  fpeaks_info.z1 = fpeaks(iTr).dist_z1;
  fpeaks_info.z2 = fpeaks(iTr).dist_z2;
  fpeaks_info.spacingPeaks = fpeaks(iTr).spacingPeaks;
  fpeaks_info.nPeaks = fpeaks(iTr).nPeaks;  
  fpeaks_info.iPeak = fpeaks(iTr).iPeak;
  fpeaks_info.iSpecies = fpeaks(iTr).iSpecies;
  
  %tags = {'outflow','agu'};
  fpeaks_ind = -1;
  
  iTr_str = sprintf('%06.0f',iTr+iTr0);
  group_name = ['/traj/' iTr_str '/'];
  
  %continue
  %
  
  h5create(filePath, [group_name 't'], size(tr.t));
  %catch
  %warning('h5 structure %s already exists',[group_name 't'])      
  %end
  
  
  h5create(filePath, [group_name 'x'], size(tr.t));  
  h5create(filePath, [group_name 'y'], size(tr.t));
  h5create(filePath, [group_name 'z'], size(tr.t));
  h5create(filePath, [group_name 'vx'], size(tr.t));
  h5create(filePath, [group_name 'vy'], size(tr.t));
  h5create(filePath, [group_name 'vz'], size(tr.t));
  h5create(filePath, [group_name 'Ex'], size(tr.t));  
  h5create(filePath, [group_name 'Ey'], size(tr.t));
  h5create(filePath, [group_name 'Ez'], size(tr.t));
  h5create(filePath, [group_name 'Bx'], size(tr.t));
  h5create(filePath, [group_name 'By'], size(tr.t));
  h5create(filePath, [group_name 'Bz'], size(tr.t));
  
  h5write(filePath, [group_name 't'], tr.t);
  h5write(filePath, [group_name 'x'], tr.x);
  h5write(filePath, [group_name 'y'], tr.y);
  h5write(filePath, [group_name 'z'], tr.z);
  h5write(filePath, [group_name 'vx'], tr.vx);
  h5write(filePath, [group_name 'vy'], tr.vy);
  h5write(filePath, [group_name 'vz'], tr.vz);
  h5write(filePath, [group_name 'Ex'], Ex);
  h5write(filePath, [group_name 'Ey'], Ey);
  h5write(filePath, [group_name 'Ez'], Ez);
  h5write(filePath, [group_name 'Bx'], Bx);
  h5write(filePath, [group_name 'By'], By);
  h5write(filePath, [group_name 'Bz'], Bz);
    
  h5writeatt(filePath, group_name,'t0', 160); % tr.t0
  h5writeatt(filePath, group_name,'x0', tr.x0);
  h5writeatt(filePath, group_name,'y0', tr.y0);
  h5writeatt(filePath, group_name,'z0', tr.z0);
  h5writeatt(filePath, group_name,'vx0', tr.vx0);
  h5writeatt(filePath, group_name,'vy0', tr.vy0);
  h5writeatt(filePath, group_name,'vz0', tr.vz0);  
  h5writeatt(filePath, group_name,'m', m);
  h5writeatt(filePath, group_name,'q', q);
  %h5writeatt(filePath, group_name,'tags', tags);
  h5writeatt(filePath, group_name,'fpeaks_x', fpeaks_info.x);
  h5writeatt(filePath, group_name,'fpeaks_z', fpeaks_info.z);
  h5writeatt(filePath, group_name,'fpeaks_x1', fpeaks_info.x1);
  h5writeatt(filePath, group_name,'fpeaks_z1', fpeaks_info.z1);  
  h5writeatt(filePath, group_name,'fpeaks_x2', fpeaks_info.x2);
  h5writeatt(filePath, group_name,'fpeaks_z2', fpeaks_info.z2);
  h5writeatt(filePath, group_name,'fpeaks_spacingPeaks', fpeaks_info.spacingPeaks);
  h5writeatt(filePath, group_name,'fpeaks_nPeaks', fpeaks_info.nPeaks);
  h5writeatt(filePath, group_name,'fpeaks_iPeak', fpeaks_info.iPeak);
  h5writeatt(filePath, group_name,'fpeaks_simIteration', fpeaks_info.iteration);
  h5writeatt(filePath, group_name,'RelTol', tr.options.RelTol);
  h5writeatt(filePath, group_name,'AbsTol', tr.options.AbsTol);

  h5disp(filePath,group_name)
    
end
disp('Done.')