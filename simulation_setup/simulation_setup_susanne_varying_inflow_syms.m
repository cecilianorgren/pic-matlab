syms x y z % l lin z1 z2 a
% if you want to see the expressions with the variables, just uncomment
% above, however, then the plotting breaks down, because now I only plot it
% as a function of z, you would have to feed all the other variables into
% to the mfBtot(z,l,lin,z1,z1,a) etc (not sure if the order is right
% there). Type mfBtot in command window to verify order.
% If I include >> l lin z1 z2 a << as syms above I get for example:
% >> Jy
%  
% Jy =
%  
% ((a*(tanh((z - z1)/lin) + 1) - a*(tanh((z - z2)/lin) + 1))*((a*(tanh((z - z1)/lin)^2 - 1))/lin - (a*(tanh((z - z2)/lin)^2 - 1))/lin))/(1 - (a*(tanh((z - z1)/lin) + 1) - a*(tanh((z - z2)/lin) + 1))^2)^(1/2)
%
% instead of:
% >> Jy
%  
% Jy =
%  
% ((tanh(2*z - 6)^2/2 - tanh(2*z - 12)^2/2)*(tanh(2*z - 6)/4 - tanh(2*z - 12)/4))/(1 - (tanh(2*z - 6)/4 - tanh(2*z - 12)/4)^2)^(1/2)

% Specify some variables
zvec = linspace(-10,10,500);
a = 0.25;
z1 = 3;
z2 = 6;
lin = 0.1;
l = 1;

B0 = 1;
R = [x; y; z];
% You can also try it with a constant B0, and later add the tangent field.
%Btot = B0*tanh(z/l);
Btot = B0;
f1 = a*(1 + tanh((z-z1)/lin));
f2 = -a*(1 + tanh((z-z2)/lin));
f = f1;
By = B0*f;
% If I add the sign(Btot) here, the Jy looks better but gets even more 
% complicted (with some dirac function)
Bx = sqrt(Btot^2 - By^2);%*sign(Btot);
Bz = 0;
B = [Bx;By;Bz];
J = curl(B,R);
Jx = J(1);
Jy = J(2);
Jz = J(3);


ff = symfun(f,[z]);
fBtot = symfun(Btot,[z]);
fBx = symfun(Bx,[z]);
fBy = symfun(By,[z]);
fJx = symfun(Jx,[z]);
fJy = symfun(Jy,[z]);

% Make into matlab handle functions (for plotting)
mff = matlabFunction(ff);
mfBtot = matlabFunction(fBtot);
mfBy = matlabFunction(fBy);
mfBx = matlabFunction(fBx);
mfJy = matlabFunction(fJy);
mfJx = matlabFunction(fJx);

% Plot
nrows = 3;
for irow = 1:nrows
  h(irow) = subplot(nrows,1,irow);
end
isub = 1;

if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,zvec,mff(zvec))
  legend(hca,{'f'})
end
if 1
  hca = h(isub); isub = isub + 1;
  hlines = plot(hca,zvec,mfBtot(zvec),zvec,mfBx(zvec),zvec,mfBy(zvec),zvec,sqrt(mfBx(zvec).^2+mfBy(zvec).^2));
  hlines(4).LineStyle = '--';
  legend(hca,{'|B|','B_x','B_y','(B_x^2+B_y^2)^{1/2}'},'location','best')
end
if 1
  hca = h(isub); isub = isub + 1;
  hlines = plot(hca,zvec,mfJx(zvec),zvec,mfJy(zvec));  
  legend(hca,{'J_x','J_y'},'location','best')
  %legend(hca,{'|B|','B_x','B_y','(B_x^2+B_y^2)^{1/2}'},'location','best')
end
