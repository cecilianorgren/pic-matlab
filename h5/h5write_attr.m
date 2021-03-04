function h5write_attr(pic_orig,timesteps,attr_str,attr_data)
%function h5write_attr(picobj,timesteps,attr_str,attr_data)
% H5WRITE_ATTR Write some auxcilliary data to h5 file.
% H5WRITE_ATTR(picobj,timesteps,attr_str,attr_data)
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
%   h5write_fields('/Users/cecilia/tesla/cno062/df_cold_protons_n04/data/','/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5',7000:50:8000,6)
%

pic = pic_orig.twcilim(timesteps,'exact');
if not(isequal(timesteps,pic.twci))
  error('Check timesteps.')  
end
if not(size(attr_data,1)==pic.nt) 
  attr_data = attr_data';
end

for itime = 1:pic.nt
  pic_tmp = pic.twcilim(pic.twci(itime));    
  str_iteration = sprintf('%010.0f',pic_tmp.iteration); % same group format as SMILEI      
  dataset_name = ['/data/' str_iteration '/'];
  disp(dataset_name)    
  % Write attributes
  
  h5writeatt(pic.file,dataset_name, attr_str ,attr_data(itime,:))
end

disp('Done.')