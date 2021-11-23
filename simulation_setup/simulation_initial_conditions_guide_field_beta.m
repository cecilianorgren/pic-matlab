syms z z1 z2 z3 z4 a ...
  Bx0(z,z1,z2,z3,z4,a) ...
  By0(z,z1,z2,z3,z4,a) ...
  Bmag0(z,z1,z2,z3,z4,a) ...
  f(z,z1,z2,z3,z4,a) ...
  Bx(z,z1,z2,z3,z4,a) ...
  By(z,z1,z2,z3,z4,a) ...
  Bmag(z,z1,z2,z3,z4,a) 

B0 = 1;
L = 1;
lin = 1;

Bx0 = B0*tanh(z/L);
By0 = B0;
theta = atan(By0,Bx0);
Bmag0 = sqrt(Bx0^2 + By0^2);

f1 = +a*(1 + tanh((z-z1)/lin));
f2 = -a*(1 + tanh((z-z2)/lin));
f3 = +a*(1 + tanh((z-z3)/lin));
f4 = -a*(1 + tanh((z-z4)/lin));
f = a + f1+f2+f3+f4;

%ff = symfun(f,[z z1 z2 z3 z4 a]);
mff = matlabFunction(f);

Bx = Bmag0*ff*cos(theta);
By = Bmag0*ff*sin(theta);
Bmag = sqrt(Bx^2 + By^2);

zvec = linspace(-25,25,1000);
z1 = 2;
z2 = 5;
z3 = 7;
z4 = 10;
a = 0.1;

mfBx0 = matlabFunction(Bx0);
mfBy0 = matlabFunction(By0);
mfBmag0 = matlabFunction(Bmag0);

mfBx = matlabFunction(Bx);
mfBy = matlabFunction(By);
mfBmag = matlabFunction(Bmag);


h = setup_subplots(3,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,zvec,mff(zvec,z1,z2,z3,z4,a))

hca = h(isub); isub = isub + 1;
plot(hca,zvec,mfBx0(zvec,z1,z2,z3,z4,a),zvec,mfBy0(zvec,z1,z2,z3,z4,a),zvec,mfBmag0(zvec,z1,z2,z3,z4,a))

hca = h(isub); isub = isub + 1;
plot(hca,zvec,mfBx(zvec,z1,z2,z3,z4,a),zvec,mfBy(zvec,z1,z2,z3,z4,a),zvec,mfBmag(zvec,z1,z2,z3,z4,a))
