
pic = df04.twpelim(8000);
% Load data
x = pic.xi;
z = pic.zi;
Bx = pic.Bx;
By = pic.By;
Bz = pic.Bz;
A = pic.A; 
Ex = pic.Ex;
Ez = pic.Ez;

% This is the line integration
x0 = 180;
z0 = 3;
dx = 0.1;
dy = 0.1;
dz = 0.1;
nsteps = 30.1; % if nsteps is not an even number, integration will stop 
                % when the fieldline arclength is above nsteps
[linearc,linex,liney,linez,linebx,lineby,linebz] = fieldline(x0,z0,x,z,Bx,By,Bz,dx,dy,dz,nsteps);
% Get field values at the points of the field line
By_along_fieldline = interpfield(x,z,By,linex,linez); 
Ex_along_fieldline = interpfield(x,z,Ex,linex,linez); 
% These fields can be integrated using some integration function like trapz
% ...

% Plotting
doA = 1;
nrows = 3; 
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % By
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
    levA = -30:1:0;
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
if 1 % Plot field line in 2D
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
if 1 % By along that line
  hca = h(isub); isub = isub + 1;  
  plot(linearc,By_along_fieldline)
  hca.XLabel.String = 'Arclength';
  hca.YLabel.String = 'B_y';  
end



