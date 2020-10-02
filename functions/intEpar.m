tag = 'nobg';
twpe = 1000;
xlim = [99 211];
zlim = [-1 20];
pic = nobg.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');
field = 'B';
% Load field and make unit vectors
Fx = pic.Bx;
Fy = pic.By;
Fz = pic.Bz;
fabs = sqrt(Fx.^2 + Fy.^2 + Fz.^2);
fx = Fx./fabs;
fy = Fy./fabs;
fz = Fz./fabs;
A = pic.xlim(xlim).zlim(zlim).A;
Epar = pic.xlim(xlim).zlim(zlim).Epar;
Ex = pic.xlim(xlim).zlim(zlim).Ex;
Ey = pic.xlim(xlim).zlim(zlim).Ey;
Ez = pic.xlim(xlim).zlim(zlim).Ez;

%% First make afew lines and plot them as scatter lines, to see how the 
% overall structure is
x = pic.xi;
z = pic.zi;
ds = 0.02;
dspl = 10;
z0_all = [2:0.2:9 9.1:0.2:14];
z0_all = [2:0.2:19];
%z0_all = [1:0.2:1.8];
x0_all = 100+zeros(numel(z0_all),1);
nlines = numel(z0_all);
h = setup_subplots(4,1);
contour(h(1),x(2:3:end),z(2:3:end),A(2:3:end,2:3:end)',-50:1:50,'color',[0 0 0]+0.4)
contour(h(2),x(2:3:end),z(2:3:end),A(2:3:end,2:3:end)',-50:1:50,'color',[0 0 0]+0.4)
contour(h(3),x(2:3:end),z(2:3:end),A(2:3:end,2:3:end)',-50:1:50,'color',[0 0 0]+0.4)
contour(h(4),x(2:3:end),z(2:3:end),A(2:3:end,2:3:end)',-50:1:50,'color',[0 0 0]+0.4)
hb(1) = colorbar('peer',h(1));
hb(2) = colorbar('peer',h(2));
hb(3) = colorbar('peer',h(3));
hb(4) = colorbar('peer',h(4));
clim{1} = [0 1e-3];
clim{2} = [0 1e-3];
clim{3} = [0 1e-3];
clim{4} = [0 1e-3];
hold(h(1),'on')
hold(h(2),'on')
hold(h(3),'on')
hold(h(4),'on')
colormap(pic_colors('blue_red'))
for iline = 1:nlines
  tic;
  [linearc,linex,liney,linez,linebx,lineby,linebz] = fieldline(x0_all(iline),z0_all(iline),x,z,Fx,Fy,Fz,ds,inf);
  toc;
  Epar_along_fieldline = interpfield(x,z,Epar,linex,linez); 
  intEpar_along_fieldline = -cumtrapz(linearc,Epar_along_fieldline);
  Ex_along_fieldline = interpfield(x,z,Ex,linex,linez); 
  Ey_along_fieldline = interpfield(x,z,Ey,linex,linez); 
  Ez_along_fieldline = interpfield(x,z,Ez,linex,linez); 
  intEpar_along_fieldline_x = -cumtrapz(linex,Ex_along_fieldline);
  intEpar_along_fieldline_y = -cumtrapz(liney,Ey_along_fieldline);
  intEpar_along_fieldline_z = -cumtrapz(linez,Ez_along_fieldline);
  
  isub = 0;
  if 1
    isub = isub + 1;
    hca = h(isub);
    clim{isub} = [min([clim{isub}(1) min(intEpar_along_fieldline)]) max([clim{isub}(2) intEpar_along_fieldline])];
    scatter(hca,linex(1:dspl:end),linez(1:dspl:end),5,intEpar_along_fieldline(1:dspl:end))
    hca.CLim = clim{isub};
  end
  if 1
    isub = isub + 1;
    hca = h(isub);
    clim{isub} = [min([clim{isub}(1) min(intEpar_along_fieldline_x)]) max([clim{isub}(2) intEpar_along_fieldline_x])];
    scatter(hca,linex(1:dspl:end),linez(1:dspl:end),5,intEpar_along_fieldline_x(1:dspl:end))
    hca.CLim = clim{isub};
  end
  if 1
    isub = isub + 1;
    hca = h(isub);
    clim{isub} = [min([clim{isub}(1) min(intEpar_along_fieldline_y)]) max([clim{isub}(2) intEpar_along_fieldline_y])];
    scatter(hca,linex(1:dspl:end),linez(1:dspl:end),5,intEpar_along_fieldline_y(1:dspl:end))
    hca.CLim = clim{isub};
  end
  if 1
    isub = isub + 1;
    hca = h(isub);
    clim{isub} = [min([clim{isub}(1) min(intEpar_along_fieldline_z)]) max([clim{isub}(2) intEpar_along_fieldline_z])];
    scatter(hca,linex(1:dspl:end),linez(1:dspl:end),5,intEpar_along_fieldline_z(1:dspl:end))
    hca.CLim = clim{isub};
  end
  drawnow
  pause(0.1)  
end
contour(h(1),x(2:3:end),z(2:3:end),A(2:3:end,2:3:end)',-50:1:50,'color',[0 0 0]+0.4)
contour(h(2),x(2:3:end),z(2:3:end),A(2:3:end,2:3:end)',-50:1:50,'color',[0 0 0]+0.4)
contour(h(3),x(2:3:end),z(2:3:end),A(2:3:end,2:3:end)',-50:1:50,'color',[0 0 0]+0.4)
contour(h(4),x(2:3:end),z(2:3:end),A(2:3:end,2:3:end)',-50:1:50,'color',[0 0 0]+0.4)
hb(1).YLabel.String = '-\int E_{||}dl_{||}';
hb(2).YLabel.String = '-\int E_{x}dl_{x}';
hb(3).YLabel.String = '-\int E_{y}dl_{y}';
hb(4).YLabel.String = '-\int E_{z}dl_{z}';
%hold(hca,'off') 9
c_eval('h(?).XLabel.String = ''x'';',1:4);
c_eval('h(?).YLabel.String = ''z'';',1:4);
hca.YLabel.String = 'z';
colormap(pic_colors('blue_red'))
h(1).Title.String = sprintf('%s, tw_{pe}= %g,',tag,pic.twpe);
hlink = linkprop(h,{'XLim','YLim','CLim'});
hca.CLim = [-1.6 1.6];
compact_panels(0.01)
h(1).XLim(1) = 110;

%%
% We need to fill the grid with cumulative interpolated values of -intEds
dx = (pic.xi(2)-pic.xi(1));
dz = (pic.zi(2)-pic.zi(1));
ds = dz/2;
      
      for ix = 2:pic.nx  
        1;
        disp(sprintf('ix = %g/%g',ix,pic.nx))
        for iz = 1:pic.nz          
          try            
            out = interp(obj,x,z,t,field,varargin);
          % interpolate from current to previous x          
          bx_ = bx(ix,iz);%*sign(bx(ix,iz));
          by_ = by(ix,iz);%*sign(bx(ix,iz));
          bz_ = bz(ix,iz);%*sign(bx(ix,iz));          
          % make interpolation in z-plane
          znew = pic.zi(iz)+dz;
          if znew > pic.zi(end) || znew < pic.zi(1)
            continue
          end
          iz1 = find(pic.zi<znew,1,'last');
          iz2 = find(pic.zi>znew,1,'first');
          
          % distance to points
          dz1 = (pic.zi(iz1)-znew)/dz;
          dz2 = (pic.zi(iz2)-znew)/dz;
          % Values at the bounding pounts
%           Bx1 = Bx(ix,iz1); Bx2 = Bx(ix,iz2);
%           By1 = By(ix,iz1); By2 = By(ix,iz2);
%           Bz1 = Bz(ix,iz1); Bz2 = Bz(ix,iz2);
          Epar1 = Epar(ix,iz1); Epar2 = Epar(ix,iz2);
          intEpar1 = intEpar(ix,iz1); intEpar2 = intEpar(ix,iz2);
          
%           Bxprev = dz2 * Bx1 + dz1 * Bx2;
%           Bxprev = dz2 * By1 + dz1 * By2;
%           Bxprev = dz2 * Bz1 + dz1 * Bz2;
          Eparprev = dz2 * Epar1 + dz1 * Epar2;
          intEparprev = dz2 * intEpar1 + dz1 * intEpar2;
          intEpar(ix,iz) = intEparprev-Eparprev*ds;
          catch
            1;
          end
        end
      end
      out = intEpar;