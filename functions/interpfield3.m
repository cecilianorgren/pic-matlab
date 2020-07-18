function fout = interpfield3(x,y,z,f,xq,yq,zq)
% INTERPFIELD Linear interpolation of 3D field.
%   bout = interpfield(x,y,z,f,xq,yq,zq)
%   x, y, z - grid
%   f - field to interploate (defined on grid x,y,z)  
%   xq, yq, zq - point to interpolate to, can be scalar values or vectors
%
nq = numel(xq);
for iq = 1:nq  
  % find two closest grid points in x and z
  %   11 ------- 21
  %    |  o      |      o is point (xq,zq)
  %    |         |
  %    |         |
  %   12---------22  
  
  % First do two 2D interpolations, then one 1D interpolation
  ix1 = find(x<xq(iq),1,'last');
  ix2 = find(x>xq(iq),1,'first');
  iy1 = find(y<yq(iq),1,'last');
  iy2 = find(y>yq(iq),1,'first');
  iz1 = find(z<zq(iq),1,'last');
  iz2 = find(z>zq(iq),1,'first');
  % pick out x and z values for these points
  x1 = x(ix1);
  x2 = x(ix2);
  y1 = y(iy1);
  y2 = y(iy2); 
  z1 = z(iz1);
  z2 = z(iz2);  
  % pick out field values at the eight corners ()
  f111 = f(ix1,iy1,iz1);
  f112 = f(ix1,iy1,iz2);
  f121 = f(ix1,iy2,iz1);
  f122 = f(ix1,iy2,iz2);
  f211 = f(ix2,iy1,iz1);
  f212 = f(ix2,iy1,iz2);
  f221 = f(ix2,iy2,iz1);
  f222 = f(ix2,iy2,iz2);  
  % distances between point 'o' and edges of box (x1,x2,z1,z2)
  xx = (xq(iq) - x1)/(x2-x1);
  yy = (yq(iq) - y1)/(y2-y1);  
  zz = (zq(iq) - z1)/(z2-z1);
  % bilinear interpolation
  % you can solve for coefficients by matrix operations, but for now I just
  % used the expression directly: https://en.wikipedia.org/wiki/Bilinear_interpolation
  % Here are instructions for 3D: https://en.wikipedia.org/wiki/Trilinear_interpolation    
  % First do two 2D interpolations in yz-plane, then one 1D interpolation
  % along x (which order you do the planes in doesn't matter)
  % Bilinear
  f1 = f111*(1-yy)*(1-zz) + f121*yy*(1-zz) + f112*(1-yy)*zz + f122*yy*zz;
  f2 = f211*(1-yy)*(1-zz) + f221*yy*(1-zz) + f212*(1-yy)*zz + f222*yy*zz;
  % Linear
  f12(iq) = f1*(1-xx) + f2*xx;
      
  % debug
  %disp(sprintf('istep = %g, [x1,x1,z1,z1] = [%g,%g,%g,%g],  [b11,b12,b21,b22] = [%g,%g,%g,%g], binterp = %g',istep,x1,x2,z1,z2,b11,b12,b21,b22,bout))
end
fout = f12;