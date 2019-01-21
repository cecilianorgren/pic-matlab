% Grid
x = xe;
z = ze;
dx = x(2)-x(1);
dz = z(2)-z(1);
[X,Z] = meshgrid(x,z);


ie = [2 4];
ii = [2 4];

% Densities
ne = squeeze(dns(:,:,2)+dns(:,:,4));
ni = squeeze(dns(:,:,1)+dns(:,:,3));
nq = (ni - ne); % charge density

% Velocities
vex = squeeze(sum(vxs(:,:,ie),3));
vey = squeeze(sum(vys(:,:,ie),3));
vez = squeeze(sum(vxs(:,:,ie),3));
vix = squeeze(sum(vxs(:,:,ii),3));
viy = squeeze(sum(vys(:,:,ii),3));
viz = squeeze(sum(vxs(:,:,ii),3));

% Field aligned quantitieas, ex: EdotB/absB
Babs = sqrt(bx.^2 + by.^2 + bz.^2);
Epar = (ex.*bx + ey.*by + ez.*bz)./Babs;
vepar = (vex.*bx + vey.*by + vez.*bz)./Babs;
vipar = (vix.*bx + viy.*by + viz.*bz)./Babs;