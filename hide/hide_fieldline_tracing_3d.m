%% 3D with Hide's data
% openw,lun,'epar.dat', /get_lun
% for iz=0,50 do begin
% for iy=0,50 do begin
% for ix=0,50 do begin
% printf,lun,epar(ix,iy,iz)
% end
% end
% end
% close,lun
% free_lun,lun
% 
% openw,lun,'b.dat', /get_lun
% for iz=0,50 do begin
% for iy=0,50 do begin
% for ix=0,50 do begin
% printf,lun,bx(ix,iy,iz), by(ix,iy,iz), bz(ix,iy,iz)
% end
% end
% end
% close,lun
% free_lun,lun
directory = '/Users/cecilia/Data/Hide/ellip0.10eta79/';
% See content
%s = whos('-file',directory);
%varnames = {s.name};

b = load([directory 'b.dat']);
epar = load([directory 'epar.dat']);
nx = 51;
ny = 51;
nz = 51;
bx = nan(nx,ny,nz);
by = nan(nx,ny,nz);
bz = nan(nx,ny,nz);
ic = 0;
x = linspace(0,2,nx);
y = linspace(0,2,ny);
z = linspace(0,2,nz);

Bx = reshape(b(:,1),nx,ny,nz);
By = reshape(b(:,2),nx,ny,nz);
Bz = reshape(b(:,3),nx,ny,nz);
Epar = reshape(epar,nx,ny,nz);

%% Plot Hide's data
ix = 2;
%Bx = Bx*0;
x0 = x(ix);
y0 = 0.5;
z0 = 0.9;
ds = 0.005;
nsteps = 600; % if nsteps is not an even number, integration will stop 
                % when the fieldline arclength is above nsteps
tic;
%ix = 25;
%Bx_ = squeeze(Bx(ix,:,:));
%By_ = squeeze(By(ix,:,:));
%Bz_ = squeeze(Bz(ix,:,:));
[linearc,linex,liney,linez,linebx,lineby,linebz] = fieldline3(x0,y0,z0,x,y,z,Bx,By,Bz,ds,nsteps);
Epar_along_line = interpfield3(x,y,z,Epar,linex,liney,linez);
intEpar_along_line = -cumtrapz(linearc,Epar_along_line);
toc;

nrows = 2;
ncols = 3;
h = setup_subplots(nrows,ncols);
isub = 1;
if 1  % Bx
  hca = h(isub); isub = isub + 1;
  pcolor(hca,y,z,squeeze(Bx(ix,:,:))')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('B_x(x = %g)',x(ix));
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,liney,linez,'k','linewidth',1.5)  
  plot(hca,liney(1),linez(1),'g.',liney(end),linez(end),'r.')
  hold(hca,'off')
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 1  % By
  hca = h(isub); isub = isub + 1;
  pcolor(hca,y,z,squeeze(By(ix,:,:))')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('B_y(x = %g)',x(ix));
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,liney,linez,'k','linewidth',1.5)
  plot(hca,liney(1),linez(1),'g.',liney(end),linez(end),'r.')
  hold(hca,'off')
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 1  % Bz
  hca = h(isub); isub = isub + 1;
  pcolor(hca,y,z,squeeze(Bz(ix,:,:))')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('B_z(x = %g)',x(ix));
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,liney,linez,'k','linewidth',1.5)
  plot(hca,liney(1),linez(1),'g.',liney(end),linez(end),'r.')
  hold(hca,'off')
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 1  
  hca = h(isub); isub = isub + 1;  
  plot3(hca,linex,liney,linez,linex(1),liney(1),linez(1),'go',linex(end),liney(end),linez(end),'rx')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';
  hca.ZLabel.String = 'z';
end
if 1  
  hca = h(isub); isub = isub + 1;  
  plot(hca,linearc,Epar_along_line)
  hca.XLabel.String = 'l_{||}';
  hca.YLabel.String = 'E_{||}'; 
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1  
  hca = h(isub); isub = isub + 1;  
  plot(hca,linearc,intEpar_along_line)  
  hca.XLabel.String = 'l_{||}';
  hca.YLabel.String = '-\int E_{||} dl_{||}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end

%% Plot Hide's data, 2D map of int(Epar)
intEpar = zeros(ny,nz);
ystart = zeros(ny,nz);
ystop = zeros(ny,nz);
zstart = zeros(ny,nz);
zstop = zeros(ny,nz);

ix = 2;
tic
for iy = 2:ny-1
  disp(sprintf('iy = %g',iy))
  for iz = 2:nz-1
    %disp(sprintf('[iy,iz] = [%g,%g]',iy,iz))
    x0 = x(ix);
    y0 = y(iy);
    z0 = z(iz);
    ds = 0.005; 
    
    nsteps = inf; % if nsteps is not an even number, integration will stop 
                    % when the fieldline arclength is above nsteps    
    [linearc,linex,liney,linez,linebx,lineby,linebz] = fieldline3(x0,y0,z0,x,y,z,Bx,By,Bz,ds,nsteps);
    Epar_along_line = interpfield3(x,y,z,Epar,linex,liney,linez);
    intEpar_along_line = -trapz(linearc,Epar_along_line);
    intEpar(iy,iz) = intEpar_along_line;
    ystart(iy,iz) = liney(1);
    ystop(iy,iz) = liney(end);
    zstart(iy,iz) = linez(1);
    zstop(iy,iz) = linez(end);
  end
end
toc
%%
hca = subplot(1,1,1);
imagesc(hca,y(2:end-1),z(2:end-1),intEpar')
hca.YDir = 'normal';
hb = colorbar('peer',hca);
hb.YLabel.String = '\int E_{||} dl';
hca.XLabel.String = 'y';
hca.ZLabel.String = 'z';

hold(hca,'on')
for iy = 2:ny-1  
  for iz = 2:nz-1
    plot(hca,[ystart(iy,iz) ystop(iy,iz)],[zstart(iy,iz) zstop(iy,iz)],'k')
    plot(hca,[ystart(iy,iz)],[zstart(iy,iz)],'k.')
  end
end
%%
nrows = 2;
ncols = 3;
h = setup_subplots(nrows,ncols);
isub = 1;
if 1  % Bx
  hca = h(isub); isub = isub + 1;
  pcolor(hca,y,z,squeeze(Bx(ix,:,:))')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('B_x(x = %g)',x(ix));
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,liney,linez,'k','linewidth',1.5)  
  plot(hca,liney(1),linez(1),'g.',liney(end),linez(end),'r.')
  hold(hca,'off')
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 1  % By
  hca = h(isub); isub = isub + 1;
  pcolor(hca,y,z,squeeze(By(ix,:,:))')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('B_y(x = %g)',x(ix));
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,liney,linez,'k','linewidth',1.5)
  plot(hca,liney(1),linez(1),'g.',liney(end),linez(end),'r.')
  hold(hca,'off')
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 1  % Bz
  hca = h(isub); isub = isub + 1;
  pcolor(hca,y,z,squeeze(Bz(ix,:,:))')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('B_z(x = %g)',x(ix));
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,liney,linez,'k','linewidth',1.5)
  plot(hca,liney(1),linez(1),'g.',liney(end),linez(end),'r.')
  hold(hca,'off')
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 1  
  hca = h(isub); isub = isub + 1;  
  plot3(hca,linex,liney,linez,linex(1),liney(1),linez(1),'go',linex(end),liney(end),linez(end),'rx')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';
  hca.ZLabel.String = 'z';
end
if 1  
  hca = h(isub); isub = isub + 1;  
  plot(hca,linearc,Epar_along_line)
  hca.XLabel.String = 'l_{||}';
  hca.YLabel.String = 'E_{||}'; 
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1  
  hca = h(isub); isub = isub + 1;  
  plot(hca,linearc,intEpar_along_line)  
  hca.XLabel.String = 'l_{||}';
  hca.YLabel.String = '-\int E_{||} dl_{||}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end


%% Get selection of fieldlines, to make sketch figure
intEpar = zeros(ny,nz);
ystart = zeros(ny,nz);
ystop = zeros(ny,nz);
zstart = zeros(ny,nz);
zstop = zeros(ny,nz);

%clear lines;
ix = 2;
tic
for iy = 2:1:ny-1
  disp(sprintf('iy = %g',iy))
  for iz = 2:1:nz-1
    %disp(sprintf('[iy,iz] = [%g,%g]',iy,iz))
    x0 = x(ix);
    y0 = y(iy);
    z0 = z(iz);
    ds = 0.005; 
    
    nsteps = inf; % if nsteps is not an even number, integration will stop 
                    % when the fieldline arclength is above nsteps    
    [linearc,linex,liney,linez,linebx,lineby,linebz] = fieldline3(x0,y0,z0,x,y,z,Bx,By,Bz,ds,nsteps);    
    Epar_along_line = interpfield3(x,y,z,Epar,linex,liney,linez);
    intEpar_along_line = -trapz(linearc,Epar_along_line);
    intEpar_along_line_cum = -cumtrapz(linearc,Epar_along_line);
    intEpar(iy,iz) = intEpar_along_line;
    ystart(iy,iz) = liney(1);
    ystop(iy,iz) = liney(end);
    zstart(iy,iz) = linez(1);
    zstop(iy,iz) = linez(end);
    
    lines(iy,iz).s = linearc;
    lines(iy,iz).x = linex;
    lines(iy,iz).y = liney;
    lines(iy,iz).z = linez;
    lines(iy,iz).epar  = Epar_along_line;
    lines(iy,iz).epar_int  = intEpar_along_line_cum;
  end
end
toc

%% Find equicontour lines of Epar
ix = 2;
tic
phi_lev = linspace(0,13e-4,10);
for iy = 2:1:ny-1
  disp(sprintf('iy = %g',iy))
  for iz = 2:1:nz-1
    phi_par = lines(iy,iz).epar_int;    
    % find point closest 9to each level
    for ip = 1:numel(phi_lev)
      [lev,ind] = min(abs(phi_par-phi_lev(ip)));     
      
      phi(iy,iz,ip) = phi_par(ind);
      x_phi(iy,iz,ip) = lines(iy,iz).x(ind);
      y_phi(iy,iz,ip) = lines(iy,iz).y(ind);
      z_phi(iy,iz,ip) = lines(iy,iz).z(ind);
      if 0%max(phi_par)>5e-4        
        %%
        plot(lines(iy,iz).s,phi_par,[0 2],phi_lev(ip),'-o',lines(iy,iz).s(ind), phi_par(ind),'*')
        title(sprintf('ind = %g, philev=%g',ind,phi_lev(ip)))
        1;
      end
    end
  end
end

%% Plot
nlines = numel(lines);
holdon = 0;
hca = subplot(1,1,1);
hca.XLim = [0 2];
hca.YLim = [-1 1];
hca.ZLim = [-1 1];

for iline = 1:nlines
  l = lines(iline);
  plot3(hca,l.x,l.y,l.z)
  if not(holdon), hold(hca,'on'); end
end

%[C,H] = contour3(...)
  
hold(hca,'off')

%% Plot
nplanes = numel(phi_lev);
holdon = 0;
hca = subplot(1,1,1);
hca.XLim = [0 2];

for ip = 2:8
  l = lines(ip);
  xx = x_phi(:,:,ip); 
  yy = y_phi(:,:,ip); 
  zz = z_phi(:,:,ip); 
  scatter3(hca,xx(:),yy(:),zz(:),'.')
  if not(holdon), hold(hca,'on'); end
end

%[C,H] = contour3(...)
  
hold(hca,'off')
hca.ZLim = [0.5 1.5];
hca.YLim = [0 2];
hca.XLim = [0.5 2];

%% Plot
nlines = numel(lines);
holdon = 0;
hca = subplot(1,1,1);
%hca.XLim = [0 2];
%hca.YLim = [-1 1];
%hca.ZLim = [-1 1];

for iline = 1:nlines
  l = lines(iline);
  plot(hca,l.s,-l.epar_int)
  if not(holdon), hold(hca,'on'); end
  drawnow
  pause(0.01)
end

%[C,H] = contour3(...)
  
hold(hca,'off')