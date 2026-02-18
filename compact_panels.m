function compact_panels(varargin)

doRemoveLabels = 1;

[ax,args,nargs] = irf.axescheck(varargin{:});
nargin = nargs;
if isempty(ax)
  fig = gcf;
  ax = findobj(fig.Children,'Type','Axes');  
end
nax = numel(ax);
for iax = 1:nax
  positions(iax,:) = ax(iax).Position;
end
mean_width = mean(positions(:,3),1);
mean_height = mean(positions(:,4),1);

doY = 1;
doX = 0;
keep_label = 0;
if nargin == 0
  space = 0.01;
elseif nargin == 1
  space = args{1};
elseif nargin == 2
  space = args{1};
  space_x = args{2};
  doX = 1;
elseif nargin == 3
  space = args{1};
  space_x = args{2};
  doX = 1;
  doRemoveLabels = ~args{3};
end


% find unique y-positions

unique_xpos = sort(uniquetol(positions(:,1),0.2*mean_width),'descend');
unique_ypos = sort(uniquetol(positions(:,2),0.2*mean_height),'descend');
n_rows = numel(unique_ypos);
n_cols = numel(unique_xpos);

for iax = 1:nax % find row and col
  [~,irow] = min(abs(positions(iax,2)-unique_ypos));
  row(iax) = irow;
  
  [~,icol] = min(abs(positions(iax,1)-unique_xpos));
  col(iax) = icol;
  %row(iax) = find(positions(iax,2)==unique_ypos);
  %col(iax) = find(positions(iax,1)==unique_xpos);
  iax_rowcol(iax,:) = [row(iax) col(iax)];
end

%[~,sortind] = sort(positions(:,2),'ascend'); % highest axes last
%ax = ax(sortind);
%positions = positions(sortind,:);

yminus = positions(:,2);
yplus =  positions(:,2) + positions(:,4);

ylowbound = min(yminus);
yhighbound = max(yplus);

new_yheight = (yhighbound-ylowbound-(n_rows-1)*space)/n_rows;
ybottom = min(unique_ypos);


for iax = 1:nax  
  ax(iax).Position(4) = new_yheight;  
  ax(iax).Position(2) = ybottom + (new_yheight + space)*(n_rows-row(iax));
end

if doX  
  xminus = positions(:,1);
  xplus =  positions(:,1) + positions(:,3);
  %xlowbound = min([space min(xminus)]);
  %xhighbound = max([1+space max(xplus)]);
  xlowbound = 0.5*space_x;
  xhighbound = 1-space_x;
  xlowbound = min(xminus);
  xhighbound = max(xplus);
  new_xheight = (xhighbound-xlowbound-(n_cols-1)*space_x)/n_cols;
  xbottom = min(unique_xpos);
  xbottom = xlowbound;

  for iax = 1:nax  
    ax(iax).Position(3) = new_xheight;  
    ax(iax).Position(1) = xbottom + (new_xheight + space_x)*(n_cols-col(iax));
  end
  xshift = 0.0;
  for iax = 1:nax  
    ax(iax).Position(1) = ax(iax).Position(1)-xshift;    
  end
end

if doRemoveLabels
for iax = 1:nax
  if not(ybottom == ax(iax).Position(2)) % not bottom row
    ax(iax).XTickLabels = [];
    ax(iax).XLabel.String = '';
  end
end
if doX
  for iax = 1:nax
    if not(xbottom == ax(iax).Position(1)) % not left row
      ax(iax).YTickLabels = [];
      ax(iax).YLabel.String = '';
    end
  end
end
end
1;