% PIC_CALC_SCRIPT
% Calculates various auxillary from quantities provided by read_fields:
%
%  [x,z,E,B,...
%   ni1,ne1,ni2,ne2,...
%   vi1,ve1,vi2,ve2,...
%   ji1,je1,ji2,je2,...
%   pi1,pe1,pi2,pe2,...
%   ti1,te1,ti2,te2,...
%   dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q] = read_fields(txtfile);

% Stream functions, not correct
c_eval('Se?.xz = vector_potential(x,z,ve?.x,ve?.z);',1:2) % stream function
c_eval('Si?.xz = vector_potential(x,z,vi?.x,vi?.z);',1:2) % stream function

pB = B.abs.^2/2; % magnetic pressure
KB = magnetic_field_curvature(x,z,B.x,B.y,B.z); % magnetic curvature
BB = tensor_product(B.x,B.y,B.z,B.x,B.y,B.z);

c_eval('vve? = tensor_product(ve?.x,ve?.y,ve?.z,ve?.x,ve?.y,ve?.z);',1:2)
c_eval('vvi? = tensor_product(vi?.x,vi?.y,vi?.z,vi?.x,vi?.y,vi?.z);',1:2)

c_eval('nmvve?.xx = mass(2)/mass(1)*ne?.*vve?.xx;',1:2)
c_eval('nmvve?.xy = mass(2)/mass(1)*ne?.*vve?.xy;',1:2)
c_eval('nmvve?.xz = mass(2)/mass(1)*ne?.*vve?.xz;',1:2)
c_eval('nmvve?.yy = mass(2)/mass(1)*ne?.*vve?.yy;',1:2)
c_eval('nmvve?.yz = mass(2)/mass(1)*ne?.*vve?.yz;',1:2)
c_eval('nmvve?.zz = mass(2)/mass(1)*ne?.*vve?.zz;',1:2)
c_eval('nmvvi?.xx = mass(1)/mass(1)*ni?.*vvi?.xx;',1:2)
c_eval('nmvvi?.xy = mass(1)/mass(1)*ni?.*vvi?.xy;',1:2)
c_eval('nmvvi?.xz = mass(1)/mass(1)*ni?.*vvi?.xz;',1:2)
c_eval('nmvvi?.yy = mass(1)/mass(1)*ni?.*vvi?.yy;',1:2)
c_eval('nmvvi?.yz = mass(1)/mass(1)*ni?.*vvi?.yz;',1:2)
c_eval('nmvvi?.zz = mass(1)/mass(1)*ni?.*vvi?.zz;',1:2)

% Magnetic field unit vector
b.x = B.x./B.abs;
b.y = B.y./B.abs;
b.z = B.z./B.abs;

% Motional electric field
c_eval('ve?xB = cross_product(ve?.x,ve?.y,ve?.z,B.x,B.y,B.z);',1:2) % electron motional electric field
c_eval('vi?xB = cross_product(vi?.x,vi?.y,vi?.z,B.x,B.y,B.z);',1:2) % ion motional electric field

% Electric field in differnt frames
c_eval('E_ve?xB.x = E.x + ve?xB.x; E_ve?xB.y = E.y + ve?xB.y; E_ve?xB.z = E.z + ve?xB.z;',1:2) % electron motional electric field
c_eval('E_vi?xB.x = E.x + vi?xB.x; E_vi?xB.y = E.y + vi?xB.y; E_vi?xB.z = E.z + vi?xB.z;',1:2) % electron motional electric field

% Energy dissipation in lab frame
c_eval('je?E = je?.x.*E.x + je?.y.*E.y + je?.z.*E.z;',1:2)
c_eval('ji?E = ji?.x.*E.x + ji?.y.*E.y + ji?.z.*E.z;',1:2)

% Poynting flux
ExB = cross_product(E.x,E.y,E.z,B.x,B.y,B.z); 

% Pressure gradients
c_eval('gradpe? = grad_scalar(x,z,pe?.scalar);',1:2) %
c_eval('gradpi? = grad_scalar(x,z,pi?.scalar);',1:2) %
c_eval('gradpe?_smooth = grad_scalar(x,z,smooth2(pe?.scalar,1));',1:2)
c_eval('gradpi?_smooth = grad_scalar(x,z,smooth2(pi?.scalar,1));',1:2)

% Diamagnetic drifts
c_eval('vDe? = cross_product(gradpe?_smooth.x,gradpe?_smooth.y,gradpe?_smooth.z,B.x./B.abs./B.abs./ne?,B.y./B.abs./B.abs./ne?,B.z./B.abs./B.abs./ne?);',1:2)
c_eval('vDi? = cross_product(-gradpi?_smooth.x,-gradpi?_smooth.y,-gradpi?_smooth.z,B.x./B.abs./B.abs./ni?,B.y./B.abs./B.abs./ni?,B.z./B.abs./B.abs./ni?);',1:2)

% ExB drift
vExB = cross_product(E.x,E.y,E.z,B.x./B.abs./B.abs,B.y./B.abs./B.abs,B.z./B.abs./B.abs); % Poynting flux
vExB.x = ExB.x./B.abs./B.abs; vExB.x(B.abs<0.1) = 0;
vExB.y = ExB.y./B.abs./B.abs; vExB.y(B.abs<0.1) = 0;
vExB.z = ExB.z./B.abs./B.abs; vExB.z(B.abs<0.1) = 0;

% Thermal velocities
c_eval('vte? = 2*sqrt(te?.scalar/(mass(2)/mass(1))); vte? = real(vte?);',1:2) % obs, fix temperature instead
c_eval('vti? = 2*sqrt(ti?.scalar/(mass(1)/mass(1))); vti? = real(vti?);',1:2) % obs, fix temperature instead

% Frequencies
c_eval('wce? = B.abs/(mass(2)/mass(1));',1:2)
c_eval('wci? = B.abs/(mass(1)/mass(1));',1:2)

% Length scales
c_eval('re? = vte?./wce?;',1:2)
c_eval('ri? = vti?./wci?;',1:2)

% Energy densities
UB.tot = 0.5*B.abs.^2;
UB.x = 0.5*B.x.^2;
UB.y = 0.5*B.y.^2;
UB.z = 0.5*B.z.^2;
c_eval('Uke? = mass(2)/mass(1)*0.5*ne?.*(ve?.x.^2 + ve?.y.^2 + ve?.z.^2);',1:2)
c_eval('Uki? = mass(1)/mass(1)*0.5*ni?.*(vi?.x.^2 + vi?.y.^2 + vi?.z.^2);',1:2)
c_eval('Ute? = pe?.scalar;',1:2)
c_eval('Uti? = pi?.scalar;',1:2)
Uktot = Uki1 + Uki2 + Uke1 + Uke2;
Uttot = Uti1 + Uti2 + Ute1 + Ute2;

% Fluxes
je.x = je1.x + je2.x;
je.y = je1.y + je2.y;
je.z = je1.z + je2.z;
ji.x = ji1.x + ji2.x;
ji.y = ji1.y + ji2.y;
ji.z = ji1.z + ji2.z;

% Current, previously called jtot
J.x = ji1.x + ji2.x - je1.x - je2.x;
J.y = ji1.y + ji2.y - je1.y - je2.y;
J.z = ji1.z + ji2.z - je1.z - je2.z;
J.abs = sqrt(J.x.^2+J.y.^2+J.z.^2);

% Perpendicular velocities
c_eval('ve?.perp.x = ve?.x-ve?.par.*B.x./B.abs;',1:2)
c_eval('ve?.perp.y = ve?.y-ve?.par.*B.y./B.abs;',1:2)
c_eval('ve?.perp.z = ve?.z-ve?.par.*B.z./B.abs;',1:2)
c_eval('vi?.perp.x = vi?.x-vi?.par.*B.x./B.abs;',1:2)
c_eval('vi?.perp.y = vi?.y-vi?.par.*B.y./B.abs;',1:2)
c_eval('vi?.perp.z = vi?.z-vi?.par.*B.z./B.abs;',1:2)

% Angles, pitch angles
c_eval('pae? = acosd(ve?.par./sqrt(ve?.x.^2+ve?.y.^2+ve?.z.^2));',1:2)
c_eval('pai? = acosd(vi?.par./sqrt(vi?.x.^2+vi?.y.^2+vi?.z.^2));',1:2)

c_eval('angle_Bve? = angle_vec(B,ve?);',1:2)
c_eval('angle_Bvi? = angle_vec(B,vi?);',1:2)

angles_to_do = {'vi1','vi2','ve1','ve2'};
for ia = 1:3
  for ib = (ia+1):4
    aa = angles_to_do{ia};
    bb = angles_to_do{ib};
    eval(sprintf('angle_%s%s = angle_vec(%s,%s);',aa,bb,aa,bb))
  end
end
angle_vi1vi2 = angle_vec(vi1,vi2);
angle_ve1ve2 = angle_vec(ve1,ve2);

angle_vi1vi2 = angle_vec(vi1,vi2);
angle_ve1ve2 = angle_vec(ve1,ve2);
angle_vi1vi2 = angle_vec(vi1,vi2);
angle_ve1ve2 = angle_vec(ve1,ve2);
angle_vi1vi2 = angle_vec(vi1,vi2);
angle_ve1ve2 = angle_vec(ve1,ve2);

angle_jije = angle_vec(ji,je);

% Combine densities
ne = ne1 + ne2;
ni = ni1 + ni2;

% Combined velocities 
ve.x = je.x./(ne1+ne2);
ve.y = je.y./(ne1+ne2);
ve.z = je.z./(ne1+ne2);
vi.x = ji.x./(ni1+ni2);
vi.y = ji.y./(ni1+ni2);
vi.z = ji.z./(ni1+ni2);

% Remove some moments with very low densities
nlim = 0.05;
indrem = find(ne1<nlim);
ve2.x(ne2<nlim) = 0;  
ve2.y(ne2<nlim) = 0;  
ve2.z(ne2<nlim) = 0;  
vi2.x(ni2<nlim) = 0;  
vi2.y(ni2<nlim) = 0;  
vi2.z(ni2<nlim) = 0; 
nlim = 0.05;
indrem = find(ne1<nlim);
te2.scalar(ne2<nlim) = NaN;  
te2.scalar(ne2<nlim) = NaN;  
te2.scalar(ne2<nlim) = NaN;  
ti2.scalar(ni2<nlim) = NaN;  
ti2.scalar(ni2<nlim) = NaN;  
ti2.scalar(ni2<nlim) = NaN;  