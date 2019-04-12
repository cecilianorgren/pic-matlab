function [r1,r2,r3] = new_csys(r1_orig,r2_orig)

if isstruct(r1_orig)
  r1 = r1_orig;
elseif isnumeric(r1_orig) && numel(r1_orig) == 3
  r1.x = r1_orig(1);
  r1.y = r1_orig(2);
  r1.z = r1_orig(3);
end
if isstruct(r2_orig)
  r2 = r2_orig;
elseif isnumeric(r2_orig) && numel(r2_orig) == 3
  r2.x = r2_orig(1);
  r2.y = r2_orig(2);
  r2.z = r2_orig(3);
end

r1.abs = sqrt(r1.x.^2 + r1.y.^2 + r1.z.^2);
r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);

r1.x = r1.x./r1.abs;
r1.y = r1.y./r1.abs;
r1.z = r1.z./r1.abs;


r3 = cross_product(r1.x,r1.y,r1.z,r2.x,r2.y,r2.z);
r2 = cross_product(r3.x,r3.y,r3.z,r1.x,r1.y,r1.z);
r3 = cross_product(r1.x,r1.y,r1.z,r2.x,r2.y,r2.z);

r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);
r2.x = r2.x./r2.abs;
r2.y = r2.y./r2.abs;
r2.z = r2.z./r2.abs;

r3.abs = sqrt(r3.x.^2 + r3.y.^2 + r3.z.^2);
r3.x = r3.x./r3.abs;
r3.y = r3.y./r3.abs;
r3.z = r3.z./r3.abs;
