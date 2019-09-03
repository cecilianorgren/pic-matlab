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
iss = 1:3;

% Motional electric field
c_eval('ve?xB = cross_product(ve?.x,ve?.y,ve?.z,B.x,B.y,B.z);',iss) % electron motional electric field
c_eval('vi?xB = cross_product(vi?.x,vi?.y,vi?.z,B.x,B.y,B.z);',iss) % ion motional electric field

% Scalar pressure and temperatures
c_eval('div_pe? = div_tensor(x,z,pe?);',iss) %
c_eval('div_pi? = div_tensor(x,z,pi?);',iss) %


c_eval('pe?.scalar = (pe?.xx + pe?.yy + pe?.zz)/3;',iss)
c_eval('pi?.scalar = (pi?.xx + pi?.yy + pi?.zz)/3;',iss)
c_eval('te?.scalar = (te?.xx + te?.yy + te?.zz)/3;',iss)
c_eval('ti?.scalar = (ti?.xx + ti?.yy + ti?.zz)/3;',iss)

% Inertia, convective derivative
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

c_eval('div_nmvvi? = div_tensor(x,z,nmvvi?);',iss) %
c_eval('div_nmvve? = div_tensor(x,z,nmvve?);',iss) %



% Force terms
min_ne = 0.02;
min_ni = 0.02;
c_eval('gradpene?.x = gradpe?_smooth.x./ne?; gradpene?.x(ne?<min_ne) = 0;',iss)
c_eval('gradpene?.y = gradpe?_smooth.y./ne?; gradpene?.y(ne?<min_ne) = 0;',iss)
c_eval('gradpene?.z = gradpe?_smooth.z./ne?; gradpene?.z(ne?<min_ne) = 0;',iss)
c_eval('gradpini?.x = gradpi?_smooth.x./ni?; gradpini?.x(ni?<min_ni) = 0;',iss)
c_eval('gradpini?.y = gradpi?_smooth.y./ni?; gradpini?.y(ni?<min_ni) = 0;',iss)
c_eval('gradpini?.z = gradpi?_smooth.z./ni?; gradpini?.z(ni?<min_ni) = 0;',iss)


% Electric field in different frames
c_eval('E_ve?xB.x = E.x + ve?xB.x; E_ve?xB.y = E.y + ve?xB.y; E_ve?xB.z = E.z + ve?xB.z;',iss) % electron motional electric field
c_eval('E_vi?xB.x = E.x + vi?xB.x; E_vi?xB.y = E.y + vi?xB.y; E_vi?xB.z = E.z + vi?xB.z;',iss) % electron motional electric field


% Pressure gradients
c_eval('gradpe? = grad_scalar(x,z,pe?.scalar);',iss) %
c_eval('gradpi? = grad_scalar(x,z,pi?.scalar);',iss) %
c_eval('gradpe?_smooth = grad_scalar(x,z,smooth2(pe?.scalar,1));',iss)
c_eval('gradpi?_smooth = grad_scalar(x,z,smooth2(pi?.scalar,1));',iss)

% Inertia gradients
c_eval('gradx_nmvve?xx = grad_scalar(x,z,nmvve?.xx);',iss) %
c_eval('gradx_nmvvi?xx = grad_scalar(x,z,nmvvi?.xx);',iss) %

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