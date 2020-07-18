%% 
pic = nobg.twpelim(10000).xlim([130 170]).zlim([-15 15]);
pic = no02.twpelim(3000).xlim([130 200]).zlim([-15 15]);
% Load data
x = pic.xi;
z = pic.zi;
Bx = pic.Bx;
By = pic.By;
Bz = pic.Bz;
A = pic.A; 
Ex = pic.Ex;
Ez = pic.Ez;
Epar = pic.Epar;
ni = pic.ni;
ne = pic.ne;
%%
% This is the line integration
x0 = 150;
z0 = 5;
dx = 0.05;
dy = 0.05;
dz = 0.05;
nsteps = 40.1; % if nsteps is not an even number, integration will stop 
                % when the fieldline arclength is above nsteps
tic;
[linearc,linex,liney,linez,linebx,lineby,linebz] = fieldline(x0,z0,x,z,Bx,By,Bz,dx,dy,dz,nsteps);
toc;
if 1 % Interpolate line to tighter values
  
end
% Get field values at the points of the field line
By_along_fieldline = interpfield(x,z,By,linex,linez); 
Ex_along_fieldline = interpfield(x,z,Ex,linex,linez); 
Epar_along_fieldline = interpfield(x,z,Epar,linex,linez); 
ni_along_fieldline = interpfield(x,z,ni,linex,linez); 
ne_along_fieldline = interpfield(x,z,ne,linex,linez); 
% These fields can be integrated using some integration function like trapz
% ...

% Plotting
doA = 1;
nrows = 4; 
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;
isArc = [];
isMap = [];

if 0 % By
  isMap(end+1) = isub;
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,By')
  hca.YDir = 'normal';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'B_y';
  
  if doA    
    hold(hca,'on')
    clim = hca.CLim;
    levA = floor(min(A(:))):1:ceil(max(A(:)));
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;
    contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(hca,'off')
    hca.CLim = clim;
  end
  if 1 % plot magnetic field line
    hold(hca,'on')
    plot(hca,linex,linez,'.k',linex(1),linez(1),'go',linex(end),linez(end),'rx')
    hold(hca,'off')    
  end
end
if 1 % Epar
  isMap(end+1) = isub;
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,Epar')
  hca.YDir = 'normal';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_{||}';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  if doA    
    hold(hca,'on')
    clim = hca.CLim;
    levA = floor(min(A(:))):1:ceil(max(A(:)));
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;
    contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(hca,'off')
    hca.CLim = clim;
  end
  if 1 % plot magnetic field line
    hold(hca,'on')
    plot(hca,linex,linez,'.k',linex(1),linez(1),'go',linex(end),linez(end),'rx')
    hold(hca,'off')    
  end
end
if 1 % ni
  isMap(end+1) = isub;
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,ni')
  hca.YDir = 'normal';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_{i}';
  colormap(hca,pic_colors('waterfall'))
  
  if doA
    hold(hca,'on')
    clim = hca.CLim;
    levA = floor(min(A(:))):1:ceil(max(A(:)));
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;
    contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(hca,'off')
    hca.CLim = clim;
  end
  if 1 % plot magnetic field line
    hold(hca,'on')
    plot(hca,linex,linez,'.k',linex(1),linez(1),'go',linex(end),linez(end),'rx')
    hold(hca,'off')    
  end
end
if 0 % Plot field line in 2D
  hca = h(isub); isub = isub + 1;  
  plot(hca,linex,liney,linex(1),liney(1),'go',linex(end),liney(end),'rx')
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
end
if 0 % Plot field line in 3D
  hca = h(isub); isub = isub + 1;  
  plot3(hca,linex,liney,linez,linex(1),liney(1),linez(1),'go',linex(end),liney(end),linez(end),'rx')
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';
  hca.ZLabel.String = 'z';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
end
if 0 % By along that line
  isArc(end+1) = isub;
  hca = h(isub); isub = isub + 1;  
  plot(hca,linearc,By_along_fieldline)
  hca.XLabel.String = 'Arclength';
  hca.YLabel.String = 'B_y';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % Epar along that line + int Epar
  isArc(end+1) = isub;
  hca = h(isub); isub = isub + 1;  
  plot(hca,linearc,Epar_along_fieldline,linearc,-cumtrapz(linearc,Epar_along_fieldline))
  hca.XLabel.String = 'Arclength';
  hca.YLabel.String = 'E_{||}, \int E_{||}dl_{||}'; 
  legend(hca,{'E_{||}','-\int E_{||}dl_{||}'},'location','best') 
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLim = linearc([1 end]);
end
if 0 % Epar along that line
  isArc(end+1) = isub;
  hca = h(isub); isub = isub + 1;  
  plotyy(hca,linearc,Epar_along_fieldline)
  hca.XLabel.String = 'Arclength';
  hca.YLabel.String = 'E_{||}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLim = linearc([1 end]);
end
if 1 % n along that line
  isArc(end+1) = isub;
  hca = h(isub); isub = isub + 1;  
  plot(hca,linearc,ni_along_fieldline,linearc,ne_along_fieldline)
  hca.XLabel.String = 'Arclength';
  hca.YLabel.String = 'n_{i}';  
  legend(hca,{'n_i','n_e'},'location','best')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLim = linearc([1 end]);
end

hlinksArc = linkprop(h(isArc),{'XLim'});
hlinksMap = linkprop(h(isMap),{'XLim','YLim'});

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
dx = 0.005; 
dy = 0.005;
dz = 0.005;
nsteps = 600; % if nsteps is not an even number, integration will stop 
                % when the fieldline arclength is above nsteps
tic;
%ix = 25;
%Bx_ = squeeze(Bx(ix,:,:));
%By_ = squeeze(By(ix,:,:));
%Bz_ = squeeze(Bz(ix,:,:));
[linearc,linex,liney,linez,linebx,lineby,linebz] = fieldline3(x0,y0,z0,x,y,z,Bx,By,Bz,dx,dy,dz,nsteps);
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


