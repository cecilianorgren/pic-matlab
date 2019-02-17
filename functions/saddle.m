function varargout = saddle(A,MinPeakProminence)
% Find saddle point in 2D matrix.
% [inds,vals] = saddle(A,option);
% indxy = saddle(A,option);
% inds - array with: [indx,indy]
%   option: 
%     MinPeakProminence - decides the sensitivity to peaks (findpeaks.m), 
%                         default is 1e-2
%                         inds = saddle(A,1e-2);

if not(exist('MinPeakProminence','var'))
  MinPeakProminence = 1e-2;
  fprintf('Using ''MinPeakProminence'' %g to find saddle points. \n',MinPeakProminence)
end
[nx,ny] = size(A);

indices = nan(nx,1);
vals = nan(nx,1);

for ix = 1:nx
  [val,ind] = min(A(ix,:));
  y_indices_of_min_A(ix) = ind;
  vals_of_min_A(ix) = val;
end
[pks,locs]= findpeaks(vals_of_min_A,'MinPeakProminence',MinPeakProminence);
if numel(locs) == 1
  ix_saddle = locs;
  iy_saddle = y_indices_of_min_A(locs);
else
  ix_saddle = tocolumn(locs);
  iy_saddle = tocolumn(y_indices_of_min_A(locs));
end

if nargout == 1
  varargout(1) = {[ix_saddle,iy_saddle]};
elseif nargout == 2
  varargout(1) = {[ix_saddle,iy_saddle]};
  varargout(2) = {vals_of_min_A(locs)};
end
if 0 % plot for diagnostics
  imagesc(A')
  hold on
  plot(ix_saddle,iy_saddle,'k*')
  contour(A',vals_of_min_A(locs))
  hold off
end