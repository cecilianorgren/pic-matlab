function [hax,hlab] = label_panels(h)
% Sort labels from top left -> top right, etc

fontsize = 12;
color = [0 0 0];

nh = numel(h);
hid = 1:nh;
for ih = 1:nh
  x(ih) = h(ih).Position(1);
  y(ih) = h(ih).Position(2);
end

[y_1,I] = sort(y,'descend');
hid_1 = hid(I);

unique_y = unique(y);
nrows = numel(unique_y);


hid_end = [];
for irow = 1:nrows
  irow_tmp = find(y==unique_y(irow));
  [col_1,I] = sort(x(irow_tmp));
  hid_end = [hid_end hid_1(I)];
  hid_1(I) = [];
end

for ih = 1:nh  
  hca = h(hid_end(ih));
  hleg(ih) = irf_legend(hca,sprintf('%g',ih),[0.98 0.98],'color',color,'fontsize',fontsize);
end

hax = h(hid_end);
hlab = hleg;
