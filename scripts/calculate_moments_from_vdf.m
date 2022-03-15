% Load PICDist object
%dist = PICDist('/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');
%dist = PICDist('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');
no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');


% Load 3D phase space distribution
iSpecies = 1;
it = 1;
id = 1;
fxyz = dist.f(it,id,iSpecies); 

% Load moments for comparison
%pxy_all = no02m.twpelim(fxyz.twpe).pxy(iSpecies);
%pxy_at_f = no02m.twpelim(fxyz.twpe).xlim(fxyz.x).zlim(fxyz.z).pxy(iSpecies);

dx = diff(fxyz.x); % not used
dz = diff(fxyz.z);

vx = fxyz.v;
vy = fxyz.v;
vz = fxyz.v;

dvx = vx(2)-vx(1);
dvy = vy(2)-vy(1);
dvz = vz(2)-vz(1);

% 2D reduced distributions, for plotting
fxy = squeeze(sum(fxyz.f,3))*dvz;
fxz = squeeze(sum(fxyz.f,2))*dvy;
fyz = squeeze(sum(fxyz.f,1))*dvx;

% Calculate density
% n = int(f)*dvx*dvy*dvz

n = sum(fxyz.f(:))*dvx*dvy*dvz;

% Calculate bulk velocity v = j/n, j = flux
% n*vx = int(f*vx)dvx*dvy*dvz
% n*vy = int(f*vy)dvx*dvy*dvz
% n*vz = int(f*vz)dvx*dvy*dvz

[VX,VY,VZ] = ndgrid(vx,vy,vz);

fvx = fxyz.f.*VX;
fvy = fxyz.f.*VY;
fvz = fxyz.f.*VZ;

vx_bulk = sum(fvx(:))*dvx*dvy*dvz/n;
vy_bulk = sum(fvy(:))*dvx*dvy*dvz/n;
vz_bulk = sum(fvz(:))*dvx*dvy*dvz/n;

Pxx_arg = fxyz.f.*(VX-vx_bulk).*(VX-vx_bulk);
Pxy_arg = fxyz.f.*(VX-vx_bulk).*(VY-vy_bulk);
Pxz_arg = fxyz.f.*(VZ-vz_bulk).*(VZ-vz_bulk);


% Plot
nrows = 1;
ncols = 3;
isub = 1;

colormap(pic_colors('candy4'))

if 1 % pxy_all
  hca = subplot(nrows,ncols,isub); isub = isub + 1;
  pcolor(hca,no02m.xi,no02m.xi,pxy_all')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'p_{xy}';
  hold(hca,'on')
  plot(hca,[f.x f.x([2 1])],[f.z(1) f.z(2) f.z(2) f.z(1)],'k')
  hold(hca,'off')
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end

if 1 % fxy
  hca = subplot(nrows,ncols,isub); isub = isub + 1;
  pcolor(hca,vx,vy,fxy')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f(v_x,v_y)';
  hold(hca,'on')
  plot(hca,vx_bulk,vy_bulk,'ko')
  hold(hca,'off')
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 0
  hca = subplot(nrows,ncols,isub); isub = isub + 1;
  pcolor(hca,vx,vz,fxz')
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,vx_bulk,vz_bulk,'ko')
  hold(hca,'off')
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end
if 0
  hca = subplot(nrows,ncols,isub); isub = isub + 1;
  pcolor(hca,vy,vz,fyz')
  shading(hca,'flat')
  hold(hca,'on')
  plot(hca,vy_bulk,vz_bulk,'ko')
  hold(hca,'off')
  hca.XLabel.String = 'v_y';
  hca.YLabel.String = 'v_z';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
end

if 1 % fxy*(vx-ux)*(vy-uy)
  hca = subplot(nrows,ncols,isub); isub = isub + 1;
  pcolor(hca,vx,vy,squeeze(sum(Pxy_arg,3))')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hcb.YLabel.String = 'f(v_x,v_y)(v_x-v_{x,bulk})(v_y-v_{y,bulk})';
  hold(hca,'on')
  plot(hca,vx_bulk,vy_bulk,'ko')
  hold(hca,'off')
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  colormap(hca,pic_colors('blue_red'))
end