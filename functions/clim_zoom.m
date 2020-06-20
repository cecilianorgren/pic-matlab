function clim_zoom(h)
% CLIM_ZOOM
% Takes the current view of a plot and zooms in the caxis to fit the data

if not(exist('h','var'))
  fig = gcf; % get current figure
  h = findobj(fig.Children,'type','axes');
end

objs = findobj(h,'type','Image');
no = numel(objs);
for io = 1:no
  xdata = objs(io).XData;
  ydata = objs(io).YData;
  xlim = objs(io).Parent.XLim;
  ylim = objs(io).Parent.YLim;
  
  ix = ind_from_lim(xdata,xlim);
  iy = ind_from_lim(ydata,ylim);
  
  cdata_orig = objs(io).CData';
  cdata = cdata_orig(ix,iy);
  clim = [min(cdata(:)) max(cdata(:))];
  
  objs(io).Parent.CLim = clim;
end

function out = ind_from_lim(var,value,varargin)
  % method is the same for xlim, zlim ilim, i, twpelim, twcilim

  % Defaults
  doExact = 0; % find all exact matches, can be any number
  doBounding = 0; % find a given number of values bounding the given value
  nBounding = 0;
  doClosest = 0;
  nClosest = 1; % only the closest index, can be one or many, for example to reduce cadence

  if numel(value) == 1
    doClosest = 1;        
  end

  % Check input
  have_options = 0;
  nargs = numel(varargin);      
  if nargs > 0, have_options = 1; args = varargin(:); end      
  while have_options
    l = 1;
    switch(lower(args{1}))
      case 'closest'
        l = 2;
        doClosest = 1;
        nClosest = args{2};            
        args = args(l+1:end);    
      case 'bounding'
        l = 2;
        doBounding = 1;
        nBounding = args{2};
        args = args(l+1:end);
      case 'exact'
        l = 1;
        doExact = 1;
        args = args(l+1:end);
      otherwise
        warning(sprintf('Input ''%s'' not recognized.',args{1}))
        args = args(l+1:end);
    end        
    if isempty(args), break, end    
  end

  % Find indices
  if doBounding
    i0 = find(abs(var-value(1)) == min(abs(var-value(1))));
    i1 = i0 - nBounding;
    i2 = i0 + nBounding;
    % Check so that indices are not outside range
    if i1 < 1
      i1 = 1; 
    end
    if i2 > numel(var) 
      i2 = numel(var); 
    end
    inds = i1:i2;
  elseif doClosest        
    ii = abs(var-value(1));
    [is, index] = sort(abs(var-value(1)));
    inds = sort(index(1:nClosest)); % why did i just not take the first closest index?     
  elseif doExact
    [~,inds,~] = intersect(var,value);
  else        
    i1 = find(var >= value(1),1,'first'); 
    i2 = find(var <= value(2),1,'last'); 
    inds = i1:i2;
  end
  out = inds;
end
end