levels_edges = -30:0.25:1;
levels_centers = levels_edges(1:end-1) + 0.5*diff(levels_edges(1:end));
tic; [A_volume,A_map,A_levels] = fluxtube_volume(A,levels_edges); toc
[saddle_locations,saddle_values] = saddle(A);

%% Find x line
xline_ind = find(saddle_values == max(saddle_values));
xline_A = saddle_values(xline_ind);

ind_inner = find(A<xline_A);
ind_outer = find(A>xline_A);

A_outer = A; A_outer(ind_inner) = NaN;
A_inner = A; A_inner(ind_outer) = NaN;

A_map_outer = A_map; A_map_outer(ind_inner) = NaN;
A_map_inner = A_map; A_map_inner(ind_outer) = NaN;

A_levels_outer = A_levels; A_levels_outer(ind_inner) = NaN;
A_levels_inner = A_levels; A_levels_inner(ind_outer) = NaN;

%% Initialize figure
npanels = nvars;
nrows = 4;
ncols = 2;
npanels = nrows*ncols;
isub = 1; 
for ipanel = 1:npanels
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
end

isub = 1;
if 1 % A_volume vs levels
  hca = h(isub); isub = isub + 1;
  hline = plot(hca,levels_centers,A_volume);
  hline.Marker = '.';
  hca.XLabel.String = 'A_level';
  hca.XLabel.Interpreter = 'none';  
  hca.YLabel.String = 'A_volume';
  hca.YLabel.Interpreter = 'none';  
end
if 1 % A
  hca = h(isub); isub = isub + 1;
  imagesc(hca,x,z,A')
  hcb = colorbar('peer',hca);
  hca.Title.String = 'A';
  hca.Title.Interpreter = 'none';
end
if 1 % A_levels
  hca = h(isub); isub = isub + 1;
  imagesc(hca,x,z,A_levels')
  hcb = colorbar('peer',hca);
  hca.Title.String = 'A_levels';
  hca.Title.Interpreter = 'none';
end
if 1 % A_volume, automatic caxis
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,A_map');
  hcb = colorbar('peer',hca);  
  hca.Title.String = 'A_map';
  hca.Title.Interpreter = 'none';
end
if 1 % A_volume, outside outermost saddle point (main x line)
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,A_map_outer');
  hcb = colorbar('peer',hca);
end
if 1 % A_volume, inside outermost saddle point (main x line)
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,A_map_inner');
  hcb = colorbar('peer',hca);
end
if 1 % A_volume, outside outermost saddle point (main x line)
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,A_levels_outer');
  hcb = colorbar('peer',hca);
end
if 1 % A_volume, inside outermost saddle point (main x line)
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,x,z,A_levels_inner');
  hcb = colorbar('peer',hca);
end


for ipanel = 1:npanels
  h(ipanel).YDir = 'normal';
end