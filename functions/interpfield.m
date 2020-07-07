function [bout] = interpfield(x,z,b,xq,zq)
% INTERPFIELD Linear interpolation of 2D field.
%   bout = interpfield(x,z,b,xq,zq)
%   x, z - grid
%   b - field to interploate (defined on grid x,z)  
%   xq, zq - point to interpolate to, can be scalar values or two vectors
%
nq = numel(xq);
for iq = 1:nq  
  % find two closest grid points in x and z
  %   11 ------- 21
  %    |  o      |      o is point (xq,zq)
  %    |         |
  %    |         |
  %   12---------22  
  ix1 = find(x<xq(iq),1,'last');
  ix2 = find(x>xq(iq),1,'first');
  iz1 = find(z<zq(iq),1,'last');
  iz2 = find(z>zq(iq),1,'first');
  % pick out x and z values for these points
  x1 = x(ix1);
  x2 = x(ix2);
  z1 = z(iz1);
  z2 = z(iz2);  
  % pick out field values at the four corners
  b11 = b(ix1,iz1);
  b12 = b(ix1,iz2);
  b21 = b(ix2,iz1);
  b22 = b(ix2,iz2);
  % distances between point 'o' and edges of box (x1,x2,z1,z2)
  xx = (xq(iq) - x1)/(x2-x1);
  zz = (zq(iq) - z1)/(z2-z1);  
  % bilinear interpolation
  % you can solve for coefficients by matrix operations, but for now I just
  % used the expression directly: https://en.wikipedia.org/wiki/Bilinear_interpolation
  % Here are instructions for 3D: https://en.wikipedia.org/wiki/Trilinear_interpolation
  bout(iq) = b11*(1-xx)*(1-zz) + b21*xx*(1-zz) + b12*(1-xx)*zz + b22*xx*zz;
  
  % debug
  %disp(sprintf('istep = %g, [x1,x1,z1,z1] = [%g,%g,%g,%g],  [b11,b12,b21,b22] = [%g,%g,%g,%g], binterp = %g',istep,x1,x2,z1,z2,b11,b12,b21,b22,bout))
end