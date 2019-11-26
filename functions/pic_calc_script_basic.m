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
iss = 1:2;

% Magnetic field unit vector
B.abs = sqrt(B.x.^2 + B.y.^2 + B.z.^2);
b.x = B.x./B.abs;
b.y = B.y./B.abs;
b.z = B.z./B.abs;

% Parallel and perpendicular electric field
E.abs = sqrt(E.x.^2 + E.y.^2 + E.z.^2);
E.par = (E.x.*b.x + E.y.*b.y + E.z.*b.z);
E.perp.x = E.x - E.par.*b.x;
E.perp.y = E.y - E.par.*b.y;
E.perp.z = E.z - E.par.*b.z;

E_smooth.x = smooth2(E.x,1);
E_smooth.y = smooth2(E.y,1);
E_smooth.z = smooth2(E.z,1);

E_smooth2.x = smooth2(E_smooth.x,1);
E_smooth2.y = smooth2(E_smooth.y,1);
E_smooth2.z = smooth2(E_smooth.z,1);

if not(isfield(B,'abs')) B.abs = sqrt(B.x.^2 + B.y.^2 + B.z.^2); end
pB = B.abs.^2/2; % magnetic pressure

KB = magnetic_field_curvature(x,z,B.x,B.y,B.z); % magnetic curvature
BB = tensor_product(B.x,B.y,B.z,B.x,B.y,B.z);

% Scalar pressure and temperatures
c_eval('pe?.scalar = (pe?.xx + pe?.yy + pe?.zz)/3;',iss)
c_eval('pi?.scalar = (pi?.xx + pi?.yy + pi?.zz)/3;',iss)
c_eval('te?.scalar = (te?.xx + te?.yy + te?.zz)/3;',iss)
c_eval('ti?.scalar = (ti?.xx + ti?.yy + ti?.zz)/3;',iss)


% Force terms
min_ne = 0.02;
min_ni = 0.02;
c_eval('gradpene?.x = gradpe?_smooth.x./ne?; gradpene?.x(ne?<min_ne) = 0;',iss)
c_eval('gradpene?.y = gradpe?_smooth.y./ne?; gradpene?.y(ne?<min_ne) = 0;',iss)
c_eval('gradpene?.z = gradpe?_smooth.z./ne?; gradpene?.z(ne?<min_ne) = 0;',iss)
c_eval('gradpini?.x = gradpi?_smooth.x./ni?; gradpini?.x(ni?<min_ni) = 0;',iss)
c_eval('gradpini?.y = gradpi?_smooth.y./ni?; gradpini?.y(ni?<min_ni) = 0;',iss)
c_eval('gradpini?.z = gradpi?_smooth.z./ni?; gradpini?.z(ni?<min_ni) = 0;',iss)

% Motional electric field
c_eval('ve?xB = cross_product(ve?.x,ve?.y,ve?.z,B.x,B.y,B.z);',iss) % electron motional electric field
c_eval('vi?xB = cross_product(vi?.x,vi?.y,vi?.z,B.x,B.y,B.z);',iss) % ion motional electric field



c_eval('vve? = tensor_product(ve?.x,ve?.y,ve?.z,ve?.x,ve?.y,ve?.z);',iss)
c_eval('vvi? = tensor_product(vi?.x,vi?.y,vi?.z,vi?.x,vi?.y,vi?.z);',iss)

c_eval('nmvve?.xx = mass(2)/mass(1)*ne?.*vve?.xx;',iss)
c_eval('nmvve?.xy = mass(2)/mass(1)*ne?.*vve?.xy;',iss)
c_eval('nmvve?.xz = mass(2)/mass(1)*ne?.*vve?.xz;',iss)
c_eval('nmvve?.yy = mass(2)/mass(1)*ne?.*vve?.yy;',iss)
c_eval('nmvve?.yz = mass(2)/mass(1)*ne?.*vve?.yz;',iss)
c_eval('nmvve?.zz = mass(2)/mass(1)*ne?.*vve?.zz;',iss)
c_eval('nmvvi?.xx = mass(1)/mass(1)*ni?.*vvi?.xx;',iss)
c_eval('nmvvi?.xy = mass(1)/mass(1)*ni?.*vvi?.xy;',iss)
c_eval('nmvvi?.xz = mass(1)/mass(1)*ni?.*vvi?.xz;',iss)
c_eval('nmvvi?.yy = mass(1)/mass(1)*ni?.*vvi?.yy;',iss)
c_eval('nmvvi?.yz = mass(1)/mass(1)*ni?.*vvi?.yz;',iss)
c_eval('nmvvi?.zz = mass(1)/mass(1)*ni?.*vvi?.zz;',iss)


% Electric field in different frames
c_eval('E_ve?xB.x = E.x + ve?xB.x; E_ve?xB.y = E.y + ve?xB.y; E_ve?xB.z = E.z + ve?xB.z;',iss) % electron motional electric field
c_eval('E_vi?xB.x = E.x + vi?xB.x; E_vi?xB.y = E.y + vi?xB.y; E_vi?xB.z = E.z + vi?xB.z;',iss) % electron motional electric field

% Energy dissipation in lab frame
%c_eval('je?E = je?.x.*E.x + je?.y.*E.y + je?.z.*E.z;',iss)
%c_eval('ji?E = ji?.x.*E.x + ji?.y.*E.y + ji?.z.*E.z;',iss)

% Poynting flux
ExB = cross_product(E.x,E.y,E.z,B.x,B.y,B.z); 

% Pressure gradients
c_eval('gradpe? = grad_scalar(x,z,pe?.scalar);',iss) %
c_eval('gradpi? = grad_scalar(x,z,pi?.scalar);',iss) %
c_eval('gradpe?_smooth = grad_scalar(x,z,smooth2(pe?.scalar,1));',iss)
c_eval('gradpi?_smooth = grad_scalar(x,z,smooth2(pi?.scalar,1));',iss)

% Inertia gradients
c_eval('gradx_nmvve?xx = grad_scalar(x,z,nmvve?.xx);',iss) %
c_eval('gradx_nmvvi?xx = grad_scalar(x,z,nmvvi?.xx);',iss) %

c_eval('div_nmvvi? = div_tensor(x,z,nmvvi?);',iss) %
c_eval('div_nmvve? = div_tensor(x,z,nmvve?);',iss) %


% Diamagnetic drifts
c_eval('vDe? = cross_product(gradpe?_smooth.x,gradpe?_smooth.y,gradpe?_smooth.z,B.x./B.abs./B.abs./ne?,B.y./B.abs./B.abs./ne?,B.z./B.abs./B.abs./ne?);',iss)
c_eval('vDi? = cross_product(-gradpi?_smooth.x,-gradpi?_smooth.y,-gradpi?_smooth.z,B.x./B.abs./B.abs./ni?,B.y./B.abs./B.abs./ni?,B.z./B.abs./B.abs./ni?);',iss)

% ExB drift
vExB = cross_product(E.x,E.y,E.z,B.x./B.abs./B.abs,B.y./B.abs./B.abs,B.z./B.abs./B.abs); % Poynting flux
vExB.x = ExB.x./B.abs./B.abs; vExB.x(B.abs<0.1) = 0;
vExB.y = ExB.y./B.abs./B.abs; vExB.y(B.abs<0.1) = 0;
vExB.z = ExB.z./B.abs./B.abs; vExB.z(B.abs<0.1) = 0;

% Thermal velocities
c_eval('vte? = 2*sqrt(te?.scalar/(mass(2)/mass(1))); vte? = real(vte?);',iss) % obs, fix temperature instead
c_eval('vti? = 2*sqrt(ti?.scalar/(mass(1)/mass(1))); vti? = real(vti?);',iss) % obs, fix temperature instead

% Frequencies
c_eval('wce? = B.abs/(mass(2)/mass(1));',iss)
c_eval('wci? = B.abs/(mass(1)/mass(1));',iss)

% Length scales
c_eval('re? = vte?./wce?;',iss)
c_eval('ri? = vti?./wci?;',iss)

% Energy densities
UB.tot = 0.5*B.abs.^2;
UB.x = 0.5*B.x.^2;
UB.y = 0.5*B.y.^2;
UB.z = 0.5*B.z.^2;
c_eval('Uke? = mass(2)/mass(1)*0.5*ne?.*(ve?.x.^2 + ve?.y.^2 + ve?.z.^2);',iss)
c_eval('Uki? = mass(1)/mass(1)*0.5*ni?.*(vi?.x.^2 + vi?.y.^2 + vi?.z.^2);',iss)
c_eval('Ute? = 3/2*pe?.scalar;',iss)
c_eval('Uti? = 3/2*pi?.scalar;',iss)
Uktot = Uki1 + Uki2 + Uki3 + Uke1 + Uke2 + Uke3;
Uttot = Uti1 + Uti2 + Uti3 + Ute1 + Ute2 + Ute3;

% Mean energy densitites, multiply with the box size to get integrated
% quantities
U.B_mean = mean(UB.tot(:));
c_eval('U.Uke?_mean = mean(Uke?(:));',iss)
c_eval('U.Uki?_mean = mean(Uki?(:));',iss)
c_eval('U.Ute?_mean = mean(Ute?(:));',iss)
c_eval('U.Uti?_mean = mean(Uti?(:));',iss)
U.B_sum = sum(UB.tot(:));
c_eval('U.Uke?_sum = sum(Uke?(:));',iss)
c_eval('U.Uki?_sum = sum(Uki?(:));',iss)
c_eval('U.Ute?_sum = sum(Ute?(:));',iss)
c_eval('U.Uti?_sum = sum(Uti?(:));',iss)

% Fluxes
je.x = je1.x + je2.x + je3.x;
je.y = je1.y + je2.y + je3.y;
je.z = je1.z + je2.z + je3.z;
ji.x = ji1.x + ji2.x + ji3.x;
ji.y = ji1.y + ji2.y + ji3.y;
ji.z = ji1.z + ji2.z + ji3.z;

ratio.jije.x = ji.x./je.x;
ratio.jije.y = ji.y./je.y;
ratio.jije.z = ji.z./je.z;

% Perpendicular fluxes
c_eval('je?.perp.x = je?.x-je?.par.*B.x./B.abs;',iss)
c_eval('je?.perp.y = je?.y-je?.par.*B.y./B.abs;',iss)
c_eval('je?.perp.z = je?.z-je?.par.*B.z./B.abs;',iss)
c_eval('ji?.perp.x = ji?.x-ji?.par.*B.x./B.abs;',iss)
c_eval('ji?.perp.y = ji?.y-ji?.par.*B.y./B.abs;',iss)
c_eval('ji?.perp.z = ji?.z-ji?.par.*B.z./B.abs;',iss)
je.perp.x = je1.perp.x + je2.perp.x + je3.perp.x;
je.perp.y = je1.perp.y + je2.perp.y + je3.perp.y;
je.perp.z = je1.perp.z + je2.perp.z + je3.perp.z;
ji.perp.x = ji1.perp.x + ji2.perp.x + ji3.perp.x;
ji.perp.y = ji1.perp.y + ji2.perp.y + ji3.perp.y;
ji.perp.z = ji1.perp.z + ji2.perp.z + ji3.perp.z;

c_eval('je?.par = ne?.*ve?.par;',iss)
c_eval('ji?.par = ni?.*vi?.par;',iss)
ji.par = ji1.par + ji2.par + ji3.par;
je.par = je1.par + je2.par + je3.par;


je.abs = sqrt(je.x.^2 + je.y.^2 + je.z.^2);
ji.abs = sqrt(ji.x.^2 + ji.y.^2 + ji.z.^2);

% Current, previously called jtot
J.x = ji1.x + ji2.x + ji3.x - je1.x - je2.x - je3.x;
J.y = ji1.y + ji2.y + ji3.y - je1.y - je2.y - je3.y;
J.z = ji1.z + ji2.z + ji3.z - je1.z - je2.z - je3.z;
J.abs = sqrt(J.x.^2+J.y.^2+J.z.^2);
J.perp.x = ji1.perp.x + ji2.perp.x - je1.perp.x - je2.perp.x;
J.perp.y = ji1.perp.y + ji2.perp.y - je1.perp.y - je2.perp.y;
J.perp.z = ji1.perp.z + ji2.perp.z - je1.perp.z - je2.perp.z;
J.par = ji1.par + ji2.par - je1.par - je2.par;


% Perpendicular velocities
c_eval('ve?.perp.x = ve?.x-ve?.par.*B.x./B.abs;',iss)
c_eval('ve?.perp.y = ve?.y-ve?.par.*B.y./B.abs;',iss)
c_eval('ve?.perp.z = ve?.z-ve?.par.*B.z./B.abs;',iss)
c_eval('vi?.perp.x = vi?.x-vi?.par.*B.x./B.abs;',iss)
c_eval('vi?.perp.y = vi?.y-vi?.par.*B.y./B.abs;',iss)
c_eval('vi?.perp.z = vi?.z-vi?.par.*B.z./B.abs;',iss)

% Angles, pitch angles
if 0
c_eval('pae? = acosd(ve?.par./sqrt(ve?.x.^2+ve?.y.^2+ve?.z.^2));',iss)
c_eval('pai? = acosd(vi?.par./sqrt(vi?.x.^2+vi?.y.^2+vi?.z.^2));',iss)

c_eval('angle_Bve? = angle_vec(B,ve?);',iss)
c_eval('angle_Bvi? = angle_vec(B,vi?);',iss)

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
end

% Combine densities
ne = ne1 + ne2 + ne3;
ni = ni1 + ni2 + ni3;

% Combined velocities 
ve.x = je.x./(ne);
ve.y = je.y./(ne);
ve.z = je.z./(ne);
vi.x = ji.x./(ni);
vi.y = ji.y./(ni);
vi.z = ji.z./(ni);

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