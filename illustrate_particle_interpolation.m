% check interpolation
x1 = 0; x2 = 1;
y1 = 0; y2 = 1;

nx = 100;
ny = 100;
xedges = linspace(x1,x2,nx);
yedges = linspace(y1,y2,ny);
dx = xedges(2) - xedges(1);
dy = yedges(2) - yedges(1);
xcenters = xedges(1:end-1)+dx;
ycenters = yedges(1:end-1)+dy;

np = 10000;
xp = x1 + (x2-x1)*rand(np,1);
yp = y1 + (y2-y1)*rand(np,1);

[N,XEDGES,YEDGES] = histcounts2(xp,yp,xedges,yedges);

nrows = 1;
ncols = 3;
h = setup_subplots(nrows,ncols);
isub = 1;

hca = h(isub); isub = isub + 1;
imagesc(hca,XEDGES,YEDGES,N/np);
hb = colorbar('peer',hca);

hca = h(isub); isub = isub + 1;
imagesc(hca,XEDGES,YEDGES,smooth2(N,1));
hb = colorbar('peer',hca);

hca = h(isub); isub = isub + 1;
imagesc(hca,XEDGES,YEDGES,smooth2(N,2));
hb = colorbar('peer',hca);



function [x,y] = interpolate_particles(xp_orig,yp_orig,dx,dy)
  
  square = dx*dy; % square
  triangle = 1;
  
  x = 1;
end
