function out = cmap2color(cmap,range,val)

range = [min(range) max(range)];
ncmap = size(cmap,1);
val = val - range(1);
range = range - range(1);
% method 'linear' can't handle end points
color = interp1((1:ncmap)',cmap,val/diff(range)*ncmap,'spline');
out = color;