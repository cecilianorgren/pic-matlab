function h5write_dists_add_attr(h5FilePath,attr_str,attr_data)

disp(sprintf('Loading file to obtain existing distributions.',h5FilePath))
dist = PICDist(h5FilePath);

% if isempty(distIndRead)
%   distIndRead = dist.indices{1};
% end

if not(dist.nd{1}==numel(attr_str))
  error(sprintf('Number of attributes does not match the number of distributions.'))
  return;
end
iteration = dist.iteration;

distInd = dist.indices{1};
id = 0;
for iFile = distInd
  id = id + 1;  
  str_iteration = sprintf('%010.0f',iteration);  
  dataset_name = ['/data/',str_iteration,'/',ds100.dists{1}{id},'/fxyz'];    
  h5writeatt(h5FilePath, dataset_name, attr_str{id}, attr_data{id});
  %if h5exist && not(isempty(find(timestep==pic.twpe)))
  %  disp(sprintf('twpe = %g already exists, skipping.',timestep))
  %  continue
  %end
  %txtfile = sprintf('%s/fields-%05.0f.dat',data_dir,timestep); 
  
  %h5write_dists_single(datFilePath,h5FilePath,nSpecies,iteration)
  

end
dataset_name = ['/data/' str_iteration '/' num2str(distnumber,'%05.0f'),'/fxyz'];