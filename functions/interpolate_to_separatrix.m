function [new_x,new_z,new_data] = interpolate_to_separatrix(data,x,z,xq,zq,A,varargin)
% [] = interpolate_to_separatrix(D,x,z,xq,zq,A);
%   
%   First interpolates to separatrix obtained from contourcs, then again
%   interpolates to values xq (or zq if xq is empty). If both are empty,
%   this second step is skipped, but then the new x's can in theory be
%   different for each time step, and a cell matrix is returned for data.
%
%   Empty xq not implemented yet.
%   
%   plot(new_x',new_z') % plots separatrix contour
%   plot(new_x',new_data')
%   plot(new_x',new_data(:,:,1)') % plots first component

  new_nz = [];
  return_mat = 1;
  [ntimes,nx,nz,ncomp] = get_size(x,z,data);  
  doInterpX = 1;
  if isempty(xq)
    xq = x;      
  end
  new_nx = numel(xq);
  new_x = nan(ntimes,new_nx);
  new_z = nan(ntimes,new_nx);
  new_data = nan(ntimes,new_nx,ncomp);
  [X,Z] = meshgrid(x,z);
    
  for it = 1:ntimes
    it
    A_tmp = squeeze(A(it,:,:));
    [saddle_locations,saddle_values] = saddle(A_tmp,'sort');
    AX = saddle_values(1);
    S = contourcs(x,z,A_tmp',AX*[1 1]*0.999);
    new_x(it,:) = xq;
    new_z(it,:) = interp1(S(1).X,S(1).Y,xq); % only need to do this one ? no
    

    for icomp = 1:ncomp      
      variable = squeeze(data(it,:,:,icomp));
      Vq = interp2(X,Z,variable',S(1).X,S(1).Y);
      if doInterpX
        Vq = interp1(S(1).X,Vq,xq);
        new_data(it,:,icomp) = Vq;
      end
      
    end
  end

end
function [nt,nx,nz,ncomp] = get_size(x,z,D)
  nx = numel(x);
  nz = numel(z);
  size_data = size(D);
  size_tmp = size_data;
  ind_x = find(size_tmp == nx,1,'first'); size_tmp = size_tmp(setdiff(1:numel(size_tmp),ind_x));
  ind_z = find(size_tmp == nz,1,'first'); size_tmp = size_tmp(setdiff(1:numel(size_tmp),ind_z));
  
  if ind_x == 1
    nt = 1;
  elseif ind_x == 2
    nt = size_data(1);
  end
  ind_t = find(size_tmp == nt,1,'first'); size_tmp = size_tmp(setdiff(1:numel(size_tmp),ind_t));
  
  if not(isempty(size_tmp))
    ncomp = size_tmp(1);
  else
    ncomp = 1; 
  end    
end