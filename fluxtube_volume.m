function out = fluxtube_volume(A,levels)

if any(diff(levels)) < 0
  error(sprintf('Levels must be monotonically increasing.'))
  return;
end

nlevels = numel(levels)-1;
n_points_within_level = nan(nlevels,1);

for ilevel = 1:nlevels
  I1 = find(A>levels(ilevel));
  I2 = find(A<levels(ilevel+1));
  [C,IA,IB] = intersect(I1,I2);
  n_points_within_level(ilevel) = numel(C);
end

out = n_points_within_level;