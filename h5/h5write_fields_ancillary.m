function h5write_fields_ancillary(pic_orig,timesteps,varstr,data)
% H5WRITE_FIELDS Write ancillaty field to h5 file. For example A.
% H5WRITE_FIELDS(dirData,h5FilePath,timesteps,nSpecies)
% 
% dirData - directory of data
% h5FilePath - directory and file name
% timesteps - timesteps to resave, in units of wpe^-1
%             if empty [], check which files are in the data directory and 
%             use those
% nSpecies - number of species: required in order to read the data right
%
% Examples:
%   % Mounted tesla from tiny laptop to Fountain
%   ...


pic = pic_orig.twcilim(timesteps,'exact');
if not(isequal(timesteps,pic.twpe))
  error('Check timesteps.')  
end

for itime = 1:numel(timesteps)
  timestep = timesteps(itime);  
  pic_tmp = pic.twpelim(timestep,'exact');    
  
  % read unnormalized data    
  nnx = pic_tmp.nx; % number of grid points in x
  nnz = pic_tmp.nz; % number of grid points in z
  
  % Load iteration info
  iter = pic_tmp.iteration;
  str_iteration = sprintf('%010.0f',iter); % same group format as SMILEI                 
        
  dataset_name = ['/data/' str_iteration '/' varstr];
  disp(dataset_name)
  try
    h5create(pic_tmp.file, dataset_name, size(data));
  catch
    disp('Dataset already exists, overwriting.')
  end
  h5write( pic_tmp.file, dataset_name, data);
       
end
  

disp('Done.')