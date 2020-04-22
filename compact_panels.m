function compact_panels(varargin)

doY = 1;
doX = 0;
if nargin == 0
  space = 0.01;
elseif nargin == 1
  space = varargin{1};
elseif nargin == 2
  space = varargin{1};
  space_x = varargin{2};
  doX = 1;
end


fig = gcf;
ax = findobj(fig.Children,'Type','Axes');
nax = numel(ax);
for iax = 1:nax
  positions(iax,:) = ax(iax).Position;
end
% find unique y-positions

unique_xpos = sort(unique(positions(:,1)),'descend');
unique_ypos = sort(unique(positions(:,2)),'descend');
n_rows = numel(unique_ypos);
n_cols = numel(unique_xpos);

for iax = 1:nax % find row and col
  row(iax) = find(positions(iax,2)==unique_ypos);
  col(iax) = find(positions(iax,1)==unique_xpos);
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
  if iax == 1 % keep as is
    
  else
    %ax(iax).Position(2) = ax(iax-1).Position(2) + ax(iax-1).Position(4) + space;
   % ax(iax).Position(4) = new_yheight;
   
  end
end

if doX  
  xminus = positions(:,1);
  xplus =  positions(:,1) + positions(:,3);
  %xlowbound = min([space min(xminus)]);
  %xhighbound = max([1+space max(xplus)]);
  xlowbound = 0.5*space_x;
  xhighbound = 1-space_x;
  new_xheight = (xhighbound-xlowbound-(n_cols-1)*space_x)/n_cols;
  xbottom = min(unique_xpos);
  xbottom = xlowbound;

  for iax = 1:nax  
    ax(iax).Position(3) = new_xheight;  
    ax(iax).Position(1) = xbottom + (new_xheight + space_x)*(n_cols-col(iax));
    if iax == 1 % keep as is

    else
      %ax(iax).Position(2) = ax(iax-1).Position(2) + ax(iax-1).Position(4) + space;
     % ax(iax).Position(4) = new_yheight;

    end
  end
  xshift = 0.0;
  for iax = 1:nax  
    ax(iax).Position(1) = ax(iax).Position(1)-xshift;    
  end
end

for iax = 1:nax
  if not(ybottom == ax(iax).Position(2)) % not bottom row
    ax(iax).XTickLabels = [];
    ax(iax).XLabel.String = '';
  end
end
1;