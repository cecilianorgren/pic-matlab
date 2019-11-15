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
%old_space = yminus(2:end)-yplus(1:(end-1));

ylowbound = min(yminus);
yhighbound = max(yplus);
%miny = min(positions(:,2));

%new_yheight = (yhighbound-ylowbound-(nax-1)*space)/nax;
new_yheight = (yhighbound-ylowbound-(n_rows-1)*space)/n_rows;

ybottom = min(unique_ypos);
%new_y_positions = ybottom + (new_yheight + space)*(0:(n_rows-1));


for iax = 1:nax
  
  ax(iax).Position(4) = new_yheight;
  ax(iax).Position(2) = ybottom + (new_yheight + space)*(n_rows-row(iax));
  if iax == 1 % keep as is
    
  else
    %ax(iax).Position(2) = ax(iax-1).Position(2) + ax(iax-1).Position(4) + space;
   % ax(iax).Position(4) = new_yheight;
   
  end
end
for iax = 2:nax
  ax(iax).XTickLabels = [];
end
1;