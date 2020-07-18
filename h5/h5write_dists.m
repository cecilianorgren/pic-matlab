function h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tag)
% H5WRITE_DISTS Write Michael's simulation output distribution data to h5 file.
% H5WRITE_DISTS(dirData,h5filePath,distIndices,nSpecies,iteration)
%
% dirData - directory of data
% h5filepath - directory and file name
% distIndices - distribution numbers to resave, if empty [], go through 
%   all in directory and check them against existing ones (if any)
% nSpecies - number of species: required in order to read the data right
% iteration - time iteration (iteration is used in h5 data structure)
%
% Example
%   timestep = 10000;
%   dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_test/distributions/%05.0f/',timestep);
%   h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_test/data_h5/dists.h5';
%   distIndRead = 1:150;
%   nSpecies = 6;
%   iteration = nobg.twpelim(timestep).iteration;
%   mass = [25 1 25 1 25 1];
%   charge = [1 -1 1 -1 1 -1];
%   tag = 'line';
%
%   % h5write_dists(dirData,h5filepath,1:300,nSpecies,timestep,iteration)
%   h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tag)

h5exist = 0;
newh5file = 0;

if exist('tag','var') && not(isempty(tag))
  tag = tag;
else
  tag = '';
end

str_iteration = sprintf('%010.0f',iteration);
  
if exist(h5FilePath,'file')
  h5exist = 1;
  disp(sprintf('File %s exists. Loading file to obtain existing distributions.',h5FilePath))
  dist = PICDist(h5FilePath);
  newh5file = 0;   % if I want to add some ovearching information later
else
  H5F.create(h5FilePath); % does not work in older matlab versions
  info = h5info(h5FilePath);
  newh5file = 1;
  %h5writeatt(h5FilePath,'/','test',1)
  h5writeatt(h5FilePath,'/','software','micPIC')    
end

% data_dir can both be overarching directory with subdirectories for each
% time (or something else), or the directory for a given time
% first try for data files
fileList = dir([dirData '*.dat']);
if isempty(fileList)
  warning('No .dat files in ''%s'', please check path.',pathFile)
  return
end

% if isempty(fileList)
%   dirList = dir([data_dir '*00*']);
% end
% nDirs = numel(dirList);
% for 

nFiles = numel(fileList);

if newh5file
  id = 0;
else
  twpe_exist = dist.twpe;
  if find(twpe_exist==timestep)
    ds = dist.twpelim(timestep);
    if not(isempty(ds))
      id = ds.nd{1};
    else
      id = 0;
    end
  else
    id = 0;
  end  
end
id = distIndRead(1)-1;

for iFile = distIndRead
  id = id + 1;
  distFilePath = sprintf('%s%04.0f.dat',dirData,iFile);
  h5write_dists_single(distFilePath,h5FilePath,id,nSpecies,iteration)
  
  %if h5exist && not(isempty(find(timestep==pic.twpe)))
  %  disp(sprintf('twpe = %g already exists, skipping.',timestep))
  %  continue
  %end
  %txtfile = sprintf('%s/fields-%05.0f.dat',data_dir,timestep); 
  
  %h5write_dists_single(datFilePath,h5FilePath,nSpecies,iteration)
  

end
% Embedded function
function h5write_dists_single(datFilePath,h5FilePath,distnumber,nSpecies,iteration)
  % datFilePath - path to -dat file
  % h5FilePath - path to h5 file
  %   
  
  

  % Read distributions  
  nv = 101;
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] = read_distributions(datFilePath,nSpecies,nv);

  % Write data to file
  %for isp = 1:nss
    %dataset_name = ['/data/' str_iteration '/' num2str(distnumber),'/'];
    %h5create(filePath, dataset_name, size(fxyz));
    %h5write(filePath, dataset_name, fxyz);
    
    dataset_name = ['/data/' str_iteration '/' num2str(distnumber,'%05.0f'),'/fxyz'];
    try
      h5create(h5FilePath, dataset_name, size(fxyz));
    catch
      warning('h5 structure %s already exists, overwriting.',dataset_name)      
    end
    h5write(h5FilePath, dataset_name, fxyz);
    h5writeatt(h5FilePath, dataset_name,'x', [xlo xhi]);
    h5writeatt(h5FilePath, dataset_name,'z', [zlo zhi]);
    h5writeatt(h5FilePath, dataset_name,'ic', ic);
    h5writeatt(h5FilePath, dataset_name,'vxa', vxa);
    h5writeatt(h5FilePath, dataset_name,'vya', vya);
    h5writeatt(h5FilePath, dataset_name,'vza', vza);
    h5writeatt(h5FilePath, dataset_name,'axes', axes);
    h5writeatt(h5FilePath, dataset_name,'tag', tag);

    h5writeatt(h5FilePath, ['/data/' str_iteration],'twpe', timestep);    
    h5writeatt(h5FilePath, ['/data/' str_iteration],'twci', timestep/50);
  
    h5disp(h5FilePath,dataset_name(1:(end-4)))
    %pause
end
end