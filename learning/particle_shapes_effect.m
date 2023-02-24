% Specify grid area
xmin = 0; xmax = 1;
ymin = 0; ymax = 1;

% Simulation grid
dx = 0.1;
dy = 1;
xg = xmin:dx:xmax;
yg = ymin:dy:ymax;
[XG,YG] = meshgrid(xg,yg);

% Particle shape grid
% dx_ps = dx/10;
% dy_ps = dy/10;
% xg_ps = xmin:dx_ps:xmax;
% yg_ps = ymin:dy_ps:ymax;
% [XG_ps,YG_ps] = meshgrid(xg_ps,yg_ps);
b_xmin = -4*dx; b_xmax = -4*dx;
b_ymin = -4*dy; b_ymax = -4*dy;

dx_ps = dx/10;
dy_ps = dy/10;
xg_ps = b_xmin:dx_ps:b_xmax;
yg_ps = b_ymin:dy_ps:b_ymax;
[XG_ps,YG_ps] = meshgrid(xg_ps,yg_ps);


% Initialize particles in an area
N = 10000;
xp = xmin-0.5*dx + (xmax-xmin+dx)*rand(N,1);
yp = ymin-0.5*dy + (ymax-ymin+dy)*rand(N,1);

%xp = 0.5*(xmax-xmin) + (xmax-xmin+dx)*randn(N,1)*0.3;
%yp = 0.5*(ymax-ymin) + (ymax-ymin+dy)*randn(N,1)*0.3;

% Define particle function, done below at end of file
%f1 = @(xg,xp,dx,yg,yp,dy) (heaviside((xg-xp)+0.5*dx) - heaviside((xg-xp)-0.5*dx)).*(heaviside((yg-yp)+0.5*dy) - heaviside((yg-yp)-0.5*dy));

% Collect density into grid
n1 = zeros(numel(yg),numel(xg)); % initialize density matrix
n2 = zeros(numel(yg),numel(xg)); % initialize density matrix
tic
for ip = 1:N
  % Particle shape function
  b1x = b1(XG_ps,xp(ip),dx);
  b1y = b1(YG_ps,yp(ip),dy);
  b1xy = b1x.*b1y;

  b2x = b2(XG_ps,xp(ip),dx);
  b2y = b2(YG_ps,yp(ip),dy);
  b2xy = b2x.*b2y;
  % imagesc(b2xy)

  % Integrate partcle shape function for each grid area
  intb1 = intb(XG_ps,YG_ps,b1xy,XG,YG,dx,dy);
  intb2 = intb(XG_ps,YG_ps,b2xy,XG,YG,dx,dy);
  n1 = n1 + intb1/N;
  n2 = n2 + intb2/N;
end

toc


%
nrows = 3;
ncols = 1;
h = gobjects(nrows,ncols);
ipanel = 0;
for irow = 1:nrows
  for icol = 1:ncols
    ipanel = ipanel + 1;
    h(irow,icol) = subplot(nrows,ncols,ipanel);
  end
end
isub = 1;

hca = h(isub) ; isub = isub +1;
scatter(hca,xp,yp,50,'filled')
%shading(hca,'flat')
hca.XLim = [xmin-0.5*dx xmax+0.5*dx];
hca.YLim = [ymin-0.5*dy ymax+0.5*dy];
hca.XTick = xg;
hca.YTick = yg;
hca.XGrid = 'on';
hca.YGrid = 'on';

hca = h(isub) ; isub = isub +1;
imagesc(hca,xg,yg,n1)
hca.YDir = 'normal';
%shading(hca,'flat')
hcb = colorbar('peer',hca);

hca = h(isub) ; isub = isub +1;
imagesc(hca,xg,yg,n2)
hca.YDir = 'normal';
%shading(hca,'flat')
hcb = colorbar('peer',hca);

for ip = 1:nrows*ncols
  h(ip).Position(3) = h(2).Position(3);
end

hlinks = linkprop(h([2 3]),{'CLim'});


% Particle shape functions
function b = b1_(xg,xp,dx,pm)
% Defines function on a grid finer than the 'simulation grid', only include
% points close to particle, defined by parameter pm
  b = zeros(pm*2+1,pm*2+1);
  
  % Define function on subgrid grid
  for ig = 1:numel(b)
    xx = xg(ig)-xp;
    if xx > -0.5*dx && xx < 0.5*dx
      b(ig) = 1;    
    else 
      b(ig) = 0;
    end     
  end
end
function b = b2(xg,xp,dx)
% Defines function on a grid finer than the 'simulation grid'
  b = zeros(size(xg));
  
  % Define function on subgrid grid
  for ig = 1:numel(b)
    xx = xg(ig)-xp;
    if xx > -dx && xx <= 0    
      b(ig) = 1 + xx/(dx);
    elseif xx < dx && xx >= 0    
      b(ig) = 1 + -xx/(dx);
    else 
      b(ig) = 0;
    end     
  end
end

function fout = intb(xb,yb,b,xg,yg,dx,dy)
% Integrates particle shape function on 'simulation grid'
  fout = zeros(size(xg));
  
  bsum = 1;sum(b,'all');
  % Define function on subgrid grid
  for ig = 1:numel(fout)
    % Extract part of b belonging to this grid
    %indx = find(all([(xb>xg(ig)-0.5*dx);(xb<xg(ig)+0.5*dx)],1));
    %indy = find(all([(yb>yg(ig)-0.5*dy);(yb<yg(ig)+0.5*dy)],1));
    ind = find(all(cat(3,xb>xg(ig)-0.5*dx,xb<xg(ig)+0.5*dx,yb>yg(ig)-0.5*dy,yb<yg(ig)+0.5*dy),3));
    fout(ig) = sum(b(ind),'all')/bsum;    

  end
end

function b = b1(xg,xp,dx)
% Defines function on a grid finer than the 'simulation grid'
  b = zeros(size(xg));
  
  % Define function on subgrid grid
  for ig = 1:numel(b)
    xx = xg(ig)-xp;
    if xx > -0.5*dx && xx < 0.5*dx
      b(ig) = 1;    
    else 
      b(ig) = 0;
    end     
  end
end
function b = b2(xg,xp,dx)
% Defines function on a grid finer than the 'simulation grid'
  b = zeros(size(xg));
  
  % Define function on subgrid grid
  for ig = 1:numel(b)
    xx = xg(ig)-xp;
    if xx > -dx && xx <= 0    
      b(ig) = 1 + xx/(dx);
    elseif xx < dx && xx >= 0    
      b(ig) = 1 + -xx/(dx);
    else 
      b(ig) = 0;
    end     
  end
end

function fout = intb(xb,yb,b,xg,yg,dx,dy)
% Integrates particle shape function on 'simulation grid'
  fout = zeros(size(xg));
  
  bsum = 1;sum(b,'all');
  % Define function on subgrid grid
  for ig = 1:numel(fout)
    % Extract part of b belonging to this grid
    %indx = find(all([(xb>xg(ig)-0.5*dx);(xb<xg(ig)+0.5*dx)],1));
    %indy = find(all([(yb>yg(ig)-0.5*dy);(yb<yg(ig)+0.5*dy)],1));
    ind = find(all(cat(3,xb>xg(ig)-0.5*dx,xb<xg(ig)+0.5*dx,yb>yg(ig)-0.5*dy,yb<yg(ig)+0.5*dy),3));
    fout(ig) = sum(b(ind),'all')/bsum;    

  end
end