function varargout = fluxtube_volume(A,levels,varargin)
% FLUXTUBE_VOLUME Find relative volume within given levels of A.
%   [A_volume,A_map] = FLUXTUBE_VOLUME(A,-30:1:1,'interpolate',2);
%   [A_volume,A_map,A_levels] = FLUXTUBE_VOLUME(A,-30:1:1,'interpolate',2);
%       A_volume - relative volume 
%       A_map - map of relative volume
%       A_levels - map of levels - same as contourf
%   
%   Options:
%     'interpolate' - interpolates the matrix using interp2, a value of k=1
%                     results in 2^k-1 interpolated points between sample 
%                     values, used for better accuracy
%     'repolate' - downsamples data to original size (buggy: only works for 
%                  integer values of (1+k))
%
%   Examples:
%     [saddle_locations,saddle_values] = saddle(A,'sort');
%     dA = 0.5;
%     A_edg = saddle_values(1):dA:max(A(:));
%     A_lev = A_edg(1:end-1) + dA;
%     [A_vol,A_map] = FLUXTUBE_VOLUME(A,A_edg);
%     plot(A_lev,A_vol)
%
% See also: vector_potential, interp2

% Defualt options
doInterpolate = 0;
doRepolate = 0;
 
have_options = nargin > 2;
args = varargin;
while have_options
  switch(lower(args{1}))
    case {'interpolate'}
      l = 2;
      interpolate_value = args{2};      
      doInterpolate = 1;
    case {'repolate'}
      l = 1;      
      doRepolate = 1;  
    otherwise 
      warning(sprintf('Input argument unknown. Skipping.'))
  end
  args = args((l+1):end);
  if isempty(args), break, end
end

if any(diff(levels) < 0)
  error(sprintf('Levels must be monotonically increasing. Hint: min(A) = %g, max(A) = %g.',min(A(:)),max(A(:))))
  return;
end

if doInterpolate
  if 0
  [old_nx,old_nz] = size(A);
  old_x = 1:old_nx;
  old_z = 1:old_nz;
  new_x = (1:new_nx)/max(old_nx);
  new_z = (1:new_nz)/max(old_nz);
  [old_X,old_Z] = meshgrid(old_x,old_z);
  [new_X,new_Z] = meshgrid(new_x,new_z);
  new_A = interp2(old_X',old_Z',A,new_X',new_Z');
  A = new_A;
  else
    old_A = A;
    A = interp2(A,interpolate_value);
  end
end

nlevels = numel(levels)-1;
n_points_within_level = nan(nlevels,1);
n_points_total = numel(A);
A_map = nan(size(A));
A_levels = nan(size(A));

if 0
  A_map = A*0;
  [C,H] = contourf(A,levels);
  A_levels = H.ZData;
  for ilevel = 1:nlevels
    C = find(A_levels == levels(ilevel));
    n_points_within_level(ilevel) = numel(C);
    A_map(C) = numel(C)/n_points_total;
  end
else
  for ilevel = 1:nlevels
    I1 = find(A>levels(ilevel));
    I2 = find(A<levels(ilevel+1));
    [C,IA,IB] = intersect(I1,I2);
    n_points_within_level(ilevel) = numel(C);
    A_map(C) = numel(C)/n_points_total;
    A_levels(C) = mean(levels(ilevel:ilevel+1));
  end
end

if doRepolate % sample matrices back down to original size
  %A_map = interp2(A_map,1/(1+interpolate_value));
  A_map = downsample(A_map,1+interpolate_value);
  A_map = downsample(A_map',1+interpolate_value)';  
  A_levels = downsample(A_levels,1+interpolate_value);
  A_levels = downsample(A_levels',1+interpolate_value)';  
end

if nargout == 1
	varargout{1} = n_points_within_level/n_points_total;
elseif nargout == 2
  varargout{1} = n_points_within_level/n_points_total;
  varargout{2} = A_map;
elseif nargout == 3
  varargout{1} = n_points_within_level/n_points_total;
  varargout{2} = A_map;  
  varargout{3} = A_levels;  
end
  