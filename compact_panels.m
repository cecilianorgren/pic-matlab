function compact_panels(varargin)

if nargin == 1
  space = varargin{1};
else
  space = 0.01;
end

fig = gcf;
ax = findobj(fig.Children,'Type','Axes');
nax = numel(ax);
for iax = 1:nax
  positions(iax,:) = ax(iax).Position;
end
[~,sortind] = sort(positions(:,2),'ascend'); % highest axes last
ax = ax(sortind);
positions = positions(sortind,:);

yminus = positions(:,2);
yplus =  positions(:,2) + positions(:,4);
old_space = yminus(2:end)-yplus(1:(end-1));

ylowbound = min(yminus);
yhighbound = max(yplus);
miny = min(positions(:,2));

new_yheight = (yhighbound-ylowbound-(nax-1)*space)/nax;
for iax = 1:nax
  ax(iax).Position(4) = new_yheight;
  if iax == 1 % keep as is
    
  else
    ax(iax).Position(2) = ax(iax-1).Position(2) + ax(iax-1).Position(4) + space;
  end
end
for iax = 2:nax
  ax(iax).XTickLabels = [];
end
1;