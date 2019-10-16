function out = vector_par_perp(V,B,varargin)
% add aprallel and eprpendicular components

doAdd = 0;
if nargin == 3 && strcmp(varargin{1},'add')
  doAdd = 1;
end

Babs = sqrt(B.x.^2 + B.y.^2 + B.z.^2);
b.x = B.x./Babs;
b.y = B.y./Babs;
b.z = B.z./Babs;


Vpar = (V.x.*b.x + V.y.*b.y + V.z.*b.z);
Vperp.x = V.x - Vpar.*b.x;
Vperp.y = V.y - Vpar.*b.y;
Vperp.z = V.z - Vpar.*b.z;

if doAdd
  out = V;
  out.par = Vpar;
  out.perp = Vperp;
end
