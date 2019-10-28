function varargout = calc_fluxtube_content(x,z,A,edgesA,vars)

if any(diff(edgesA) < 0)
  error(sprintf('Levels must be monotonically increasing. Hint: min(A) = %g, max(A) = %g.',min(A(:)),max(A(:))))
  return;
end


nlevels = numel(edgesA)-1;
n_points_within_level = nan(nlevels,1);
n_points_total = numel(A);
A_map = nan(size(A));
A_levels = nan(size(A));

% Calculate which grid point belong to which A level intervals
nVars = numel(vars);
for ilevel = 1:nlevels
  I1 = find(A>edgesA(ilevel));
  I2 = find(A<edgesA(ilevel+1));
  [C,IA,IB] = intersect(I1,I2);
  n_points_within_level(ilevel) = numel(C);
  A_map(C) = numel(C)/n_points_total;
  A_levels(C) = mean(edgesA(ilevel:ilevel+1));  
  for iVar = 1:nVars
    variable = vars{iVar};
    var_sum{iVar}(ilevel) = nanmean(variable(I1(IA)));
  end
end
A_vol_rel = n_points_within_level/n_points_total;

varargout{1} = A_vol_rel;
varargout{2} = var_sum;





