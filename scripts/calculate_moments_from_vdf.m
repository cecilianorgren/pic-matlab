% Load PICDist object
dist = PICDist('/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');

% Load 3D phase space distribution
fxyz = dist.f(1,1,1); 

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

n = sum(fxyz.f,'all')*dvx*dvy*dvz;

% Calculate bulk velocity v = j/n, j = flux
% n*vx = int(f*vx)dvx*dvy*dvz
% n*vy = int(f*vy)dvx*dvy*dvz
% n*vz = int(f*vz)dvx*dvy*dvz

[VX,VY,VZ] = ndgrid(vx,vy,vz);

fvx = fxyz.f.*VX;
fvy = fxyz.f.*VY;
fvz = fxyz.f.*VZ;

vx_bulk = sum(fvx,'all')*dvx*dvy*dvz/n;
vy_bulk = sum(fvy,'all')*dvx*dvy*dvz/n;
vz_bulk = sum(fvz,'all')*dvx*dvy*dvz/n;


% Plot
nrows = 3;
ncols = 1;
isub = 1;

colormap(pic_colors('candy4'))

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,vx,vy,fxy')
shading(hca,'flat')
hold(hca,'on')
plot(hca,vx_bulk,vy_bulk,'ko')
hold(hca,'off')
hca.XLabel.String = 'v_x';
hca.YLabel.String = 'v_y';
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.Layer = 'top';

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
