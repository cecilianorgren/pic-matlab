function [r1_out,r2_out,r3_out] = construct_coordsystem(a,b,c)


r1 = a;
% Make sure its a unit vector
r1_abs = sqrt(r1.x.^2 + r1.y.^2 + r1.z.^2);
r1.x = r1.x./r1_abs;
r1.y = r1.y./r1_abs;
r1.z = r1.z./r1_abs;

% ±only supports scalar vector for now
r2 = cross_product(r1.x,r1.y,r1.z,b(1),b(2),b(3));
r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);
r2.x = r2.x./r2.abs;
r2.y = r2.y./r2.abs;
r2.z = r2.z./r2.abs;
r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);
r2 = cross_product(r2.x,r2.y,r2.z,r1.x,r1.y,r1.z);


r3 = cross_product(r1.x,r1.y,r1.z,r2.x,r2.y,r2.z);
r3.abs = sqrt(r3.x.^2 + r3.y.^2 + r3.z.^2);

r1_out.x = r1.x;
r1_out.y = r1.y;
r1_out.z = r1.z;
r2_out.x = r2.x;
r2_out.y = r2.y;
r2_out.z = r2.z;
r3_out.x = r3.x;
r3_out.y = r3.y;
r3_out.z = r3.z;