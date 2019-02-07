function out = find_closest_ind(x,xval)

nind = numel(xval);
all_ind = nan(nind,1);
for iind = 1:nind
  all_ind(iind,1) = find(abs(x-xval(iind))==min(abs(x-xval(iind))));
end
out = all_ind;