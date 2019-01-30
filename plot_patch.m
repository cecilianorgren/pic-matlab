function varargout = plot_patch(varargin)

[ax,args,nargs] = axescheck(varargin{:});
if isempty(ax), ax = gca; end

doSort = 0;
doColors = 0;
doStacked = 1;
doBase = 0; % add zero level

x = args{1}; args = args(2:end);
y = args{1}; args = args(2:end);

while ~isempty(args)
  switch args{1}
    case 'sort'
      doSort = 1;
      sortOption = args{2};
      l = 2;
    case {'color','colors'}
      doColors = 1;
      colors = args{2};
      if size(colors,1) < (size(y,1)-1)
        warning(sprintf('Insufficient number of colors provided. Recycling %g times.',ceil((size(y,1)-1)/size(colors,1))))
        colors = repmat(colors,ceil((size(y,1)-1)/size(colors,1)),1);
      end
      l = 2;
    case 'base'
      doBase = 1;
      base = args{2};
      l = 2;     
  end
  args = args((l+1):end);
end
  
if doSort
  [~,ind_sort] = sort(y(:,1),1,sortOption);
  y = y(ind_sort,:);
end
y = cumsum(y,1);
%y = y-repmat(y(:,1),1,size(y,2));

if doBase
  yBase = repmat(base,1,size(y,2));
  yBase(isnan(y(1,:))) = NaN;  
  y = [yBase; y];
end

nComps = size(y,1);
hl = plot(ax,x,y(2:end,:));

if doColors
  for iLine = 1:numel(hl)
    hl(iLine).Color = colors(iLine,:);
  end
end

hold(ax,'on')
for iComp = 1:(nComps-1)
  xPatch = [x, NaN, x(end:-1:1)];
  yPatch = [y(iComp,:), NaN, y(iComp+1,end:-1:1)];
  isNan = isnan(yPatch);
  xPatch(isNan) = [];
  yPatch(isNan) = [];

  hp = patch(ax,xPatch,yPatch,'k','Parent', ax);
  hp.FaceAlpha = 0.5;
  hp.EdgeColor = hl(iComp).Color*0;
  hp.FaceColor = hl(iComp).Color;
  hp_out(iComp,1) = hp;
end
hold(ax,'off')

switch nargout
  case 1
    varargout(1) = {hp_out};
  case 2
    varargout(1) = {hp_out};
    varargout(2) = {hl};
end