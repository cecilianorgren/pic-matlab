nz = 100; 
nvz = 101;
z = linspace(-5,5,nz); 
vz = linspace(-1.5,1.5,nvz);
[Z,VZ] = meshgrid(z,vz);

phi = -0.3*exp(-z.^2);
PHI = repmat(phi,[nvz,1]);
phi0 = 0;
v0 = VZ;
m = 1;
VZ_ =  (v0.^2 + 2*PHI/m).^0.5;

hca = subplot(2,1,1);
pcolor(hca,z,vz,PHI); shading(hca,'flat');
hca.YDir = 'normal';
hca.XLabel.String = 'z';
hca.YLabel.String = 'v_z';
hb = colorbar('peer',hca);
hb.YLabel.String = '\phi(z)';

hca = subplot(2,1,2);
contourf(hca,z,vz,real(VZ_),-10:0.1:10); shading(hca,'flat');
hca.YDir = 'normal';
hca.XLabel.String = 'z';
hca.YLabel.String = 'v_z';
hb = colorbar('peer',hca);
hb.YLabel.String = 'v_z(z,v_{z,0},\phi(z))';

%V0(iabove) = vph + ((V(iabove)-vph).^2 - 2*units.e*(PHI(iabove))/units.me).^0.5; % same as Schamel free streaming
%V0(ibelow) = vph - ((V(ibelow)-vph).^2 - 2*units.e*(PHI(ibelow))/units.me).^0.5;
