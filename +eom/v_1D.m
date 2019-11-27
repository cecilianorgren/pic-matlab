function  x_res = vlasov(x,Ex_)
% x_res = eom_general(t,x_vect,Bx_,By_,Bz_,Ex_,Ey_,Ez_)
% q = -e


% physical constants
e = 1.6022e-19;
me = 9.10939999e-31;

v = fv(2);

E = Ex_(x)

x_res = zeros(1,1);
x_res(2) = (-e/me)*(E); % dvx/dt = ax;               