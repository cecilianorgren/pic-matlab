% Grid
x = xe;
z = ze;
dx = x(2)-x(1);
dz = z(2)-z(1);
[X,Z] = meshgrid(x,z);

% Densities
ne = squeeze(dns(:,:,2)+dns(:,:,4));
ni = squeeze(dns(:,:,1)+dns(:,:,3));
nq = (ni - ne); % charge density

% Velocities
vex = squeeze(vxs(:,:,2) + vxs(:,:,4));
vey = squeeze(vys(:,:,2) + vys(:,:,4));
vez = squeeze(vzs(:,:,2) + vzs(:,:,4));
vix = squeeze(vxs(:,:,1) + vxs(:,:,3));
viy = squeeze(vys(:,:,1) + vys(:,:,3));
viz = squeeze(vzs(:,:,1) + vzs(:,:,3));

% Vector potential, for plotting xz-plane magnetic field lines
Ay1 = cumsum(bz,1)*dx; - cumsum(bx,2)*dz;

% Pauls, looks better
ixm = 10;
Ay2(ixm,1) = 0;
for j=2:nnz 
  Ay2(ixm,j) = Ay2(ixm,j-1) + dz*bx(ixm,j-1);
end
for ix=ixm+1:nnx
  Ay2(ix,:) = Ay2(ix-1,:) - bz(ix-1,:)*dx; % ;     advance to the right
end
for ix=ixm-1:-1:1
  Ay2(ix,:) = Ay2(ix+1,:)+ bz(ix,:)*dx; % ;     advance to the left
end

Ay = Ay2;

% Field aligned quantitieas, ex: EdotB/absB
Babs = sqrt(bx.^2 + by.^2 + bz.^2);
Epar = (ex.*bx + ey.*by + ez.*bz)./Babs;
vepar = (vex.*bx + vey.*by + vez.*bz)./Babs;
vipar = (vix.*bx + viy.*by + viz.*bz)./Babs;

%%
hca = subplot(1,1,1);
imagesc(hca,x,z,Babs')
hcb = colorbar('peer',hca);
clim = hca.CLim;
hold(hca,'on')
contour(hca,x,z,A',20,'color',[0.3 0.3 0.3])
hold(hca,'off')
hca.CLim = clim;
