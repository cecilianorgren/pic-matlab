function varargout = saddle(A,varargin)
% Find saddle point in 2D matrix.
%   inds = saddle(A,options);
%   [inds,vals] = saddle(A,options);
%   [inds,vals,widths] = saddle(A,options);
%   [inds,vals,widths,peakprominence] = saddle(A,options);
%   [inds,vals] = saddle(A,'sort','plot');
%   inds - indices of saddle points, array with [indx,indy]
%   vals - saddle values
%   widths - peak widths, see findpeaks.m
%   peakprominence - peak prominence, see findpeaks.m
%
%   options: 
%     MinPeakProminence - decides the sensitivity to peaks (findpeaks.m), 
%                         default is 1e-2
%                         inds = saddle(A,'MinPeakProminence',1e-2);
%     sort              - sorts indices and values in descending order
%                         based on values
%     plot              - plot results, for diagnostic purposes
%
% See also, FINDPEAKS

% Defaults
doPlot = 0;
doSort = 0;
doRefine = 0;
MinPeakProminence = 1e-2;
findpeaks_varargin = {};

args = varargin;
nargs = numel(args);
if nargs > 0; have_options = 1; else have_options = 0; end
while have_options
  l = 1;
  switch(lower(args{1}))   
    case 'minpeakprominence'
      l = 2;
      MinPeakProminence = args{2};   
      % add to findpeaks_varargin below
    case 'plot'
      l = 1;      
      doPlot = 1;
    case {'interp','refine'}
      l = 1;
      doRefine = 1;
    case 'sort'
      l = 1;
      doSort = 0;
      findpeaks_varargin{end+1} = 'SortStr';
      findpeaks_varargin{end+1} = 'descend';
    otherwise
      disp(sprintf('Unknown argument: %s',args{1}))
  end
  args = args(l+1:end);  
  if isempty(args), break, end    
end
%disp(sprintf('MinPeakProminence = %g',MinPeakProminence))
findpeaks_varargin{end+1} = 'MinPeakProminence';
findpeaks_varargin{end+1} = MinPeakProminence;
% if not(exist('MinPeakProminence','var'))
%   MinPeakProminence = 1e-2;
%   fprintf('Using ''MinPeakProminence'' %g to find saddle points. \n',MinPeakProminence)
% end

if doRefine
  
end
[nx,ny] = size(A);

indices = nan(nx,1);
vals = nan(nx,1);
if A(1,1) > 0
end

for ix = 1:nx
  [val,ind] = max(A(ix,:));
  y_indices_of_max_A(ix) = ind;
  vals_of_max_A(ix) = val;
end
%[pks,locs,W,P] = findpeaks(vals_of_min_A,'MinPeakProminence',MinPeakProminence);
[pks,locs,W,P] = findpeaks(-vals_of_max_A,findpeaks_varargin{:});
A_saddle = -pks;
if isempty(A_saddle)
  warning('No saddle points found, consider changing ''MinPeakProminence''.')  
end



if numel(locs) == 1
  ix_saddle = locs;
  iy_saddle = y_indices_of_max_A(locs);
else
  ix_saddle = tocolumn(locs);
  iy_saddle = tocolumn(y_indices_of_max_A(locs));
end

if doSort % not used, insert arguments directly into findpeaks
  [~,ind_sort] = sort(A_saddle,'descend');
  ix_saddle = ix_saddle(ind_sort);
  iy_saddle = iy_saddle(ind_sort);
  A_saddle = A_saddle(ind_sort);
end
if doPlot % plot for diagnostics
  hca = subplot(2,1,1);
  ax = plotyy(hca,1:nx,vals_of_max_A,1:nx,y_indices_of_max_A);  
  hold(hca,'on')
  h_saddle = plot(hca,ix_saddle,A_saddle,'cx');
  h_prominence = errorbar(ix_saddle,A_saddle,P,'m.');
  h_min_prominence = errorbar(ix_saddle,A_saddle,ones(size(A_saddle))*MinPeakProminence,'k.');
  hold(hca,'off')
  ax(1).YLim = [floor(min(vals_of_max_A(:))-MinPeakProminence) ceil(max(vals_of_max_A(:))+MinPeakProminence)];
  ax(1).XLabel.String = 'x index';
  ax(1).YLabel.String = 'minimum A';
  ax(2).YLabel.String = 'y index where A is minimum';
  ax(1).XLim = [1 nx];
  ax(2).XLim = [1 nx];  
  legend([h_saddle,h_min_prominence,h_prominence],{'chosen saddle points','MinPeakProminence','PeakProminence'},'location','northwest')
  ax(1).Title.String = sprintf('MeanPeakProminence = %g',MinPeakProminence);
  
  hca = subplot(2,1,2);
  imagesc(hca,A')
  hca.YDir = 'normal';
  hold(hca,'on')
  h_all = plot(hca,1:nx,y_indices_of_max_A,'k-');
  h_saddle = plot(hca,ix_saddle,iy_saddle,'cx');
  if 0 % plot corresponding contours
    if numel(A_saddle) == 1
      hc = contour(hca,A',A_saddle*[1 1]);
    else
      hc = contour(hca,A',A_saddle);
    end	
  end
  hold(hca,'off')
  hca.XLabel.String = 'x index';
  hca.YLabel.String = 'y index';
  legend([h_all,h_saddle],{'y index where A is minimum','chosen saddle points'},'location','northwest')
  hca.Title.String = sprintf('MeanPeakProminence = %g',MinPeakProminence);
  hcf=gcf; all_axes = findobj(hcf.Children,'type','axes'); linkprop(all_axes,'XLim');
end

if nargout == 1
  varargout(1) = {[ix_saddle,iy_saddle]};
elseif nargout == 2
  varargout(1) = {[ix_saddle,iy_saddle]};
  varargout(2) = {vals_of_max_A(locs)};
elseif nargout == 3
  varargout(1) = {[ix_saddle,iy_saddle]};
  varargout(2) = {vals_of_max_A(locs)};
  varargout(3) = {W};
elseif nargout == 4
  varargout(1) = {[ix_saddle,iy_saddle]};
  varargout(2) = {vals_of_max_A(locs)};
  varargout(3) = {W};
  varargout(4) = {P};
end

% Help functions
function col = tocolumn(vector)
% TOCOLUMN - outputs a column vector whatever the format of
%   the input vector. Returns zero if input is rank higher
%   than 1 (matrices etc), or input if input is empty or
%   scalar.      
%
% See also TOROW

% Anders.Eriksson@irfu.se 021208

if isempty(vector)                     % input is empty
  col = vector;
else
  s = size(vector);
  if (max(size(s)) < 3) && find(s == 1) % input is a vector (or scalar)
    if s(1) > 1                        % input is column vector
      col = vector;
    else                               % input is row vector (or scalar)
      col = vector';
    end
  else                                 % input is matrix or higher rank
    col = 0;
  end
end  
end
end