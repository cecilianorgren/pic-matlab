n0 = 1;
nc = 0.2;
L = 1;
z0 = 5;
n1 = @(z,L) n0*cosh((z-z0)/L).^(-2);
n2 = @(z,L) nc*0.5*(1+tanh((abs(z-z0)-2*L)/(0.5*L)));
n3 = @(z,L) nc*0.5*(1+tanh(((z-z0)-2*L)/(0.5*L)));
n4 = @(z,L) nc*0.5*(1+tanh((-(z-z0)-2*L)/(0.5*L)));
%n3 = @(z,L) nc*tanh(z/(3*L)).^2;

z = linspace(0,10,100);

plot(z,n1(z,L),z,n3(z,L),z,n4(z,L))

%%
n0 = 1;
nc = 0.1;
L = 1;
Lf = 1;
z0 = 40;
n1_ = @(z,L) n0*cosh((z-z0)/L).^(-2);
n2_ = @(z,L) nc*0.5*(1+tanh((abs(z-z0)-2*L)/(0.5*L)));
n3_ = @(z,L) nc*0.5*(1+tanh(((z-z0)-2*L)/(0.5*L))).*(3+cos((2*pi/(Lf))*(z-z0)))/4;
%            nb*0.5*(1+tanh(((y-y0)-2*L)/(0.5*L)))*(3+cos((2*pi/(Lf))*(y-y0)))/4
n4_ = @(z,L) 2*nc*0.5*(1+tanh((-(z-z0)-2*L)/(0.5*L))).*(3+cos((2*pi/(Lf))*(z-z0)))/4;

n5_ = @(z,L) cos((2*pi/2/L)*(z-z0));
%n3 = @(z,L) nc*tanh(z/(3*L)).^2;

z = linspace(0,80,2000);

plot(z,n1_(z,L),z,n3_(z,L),z,n4_(z,L))
