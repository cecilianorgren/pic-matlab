syms t Ey ts E0 B0 z L zmin
Bx = B0*tanh(z/L);
% indefinite: int(bx) = B0*L*log(cosh(z/L))
%FH = int(abs(Bx),z,zmin,zmax);
FH = 2*B0*L*log(cosh(zmax/L));
Ey = E0*t*exp(-t/ts);

F = int(Ey,t);



mfEy = matlabFunction(Ey);
mfF = matlabFunction(F);
mfF0 = matlabFunction(FH);

mfF_norm = @(E0,t,ts)-E0.*ts.*exp(-t./ts).*(t+0*ts);

tvec = linspace(0,200,1000);


nrows = 2;
ncols = 2;
ip = 1;
h(ip) = subplot(nrows,ncols,ip); ip = ip + 1;
h(ip) = subplot(nrows,ncols,ip); ip = ip + 1;
h(ip) = subplot(nrows,ncols,ip); ip = ip + 1;
h(ip) = subplot(nrows,ncols,ip); ip = ip + 1;


E0 = 1;
ts = 30;

isub = 1;
hca = h(isub); isub = isub + 1;
plot(hca,tvec,mfEy(E0,tvec,ts))

hca = h(isub); isub = isub + 1;
plot(hca,tvec,mfF(E0,tvec,ts),tvec,mfF_norm(E0,tvec,ts),tvec,cumtrapz(tvec,mfEy(E0,tvec,ts)),'--')


hca = h(isub); isub = isub + 1;
plot(hca,tvec,cumtrapz(tvec,mfEy(E0,tvec,ts)))

%%
Ey = @(E0,t,ts) E0.*(t./ts).*exp(-t./ts);


t = 0:400;
ts = 10:50;
E0 = 1;
[T,TS] = meshgrid(t,ts);
EY = Ey(E0*TS.^(-1),T,TS);

nrows = 2;
ncols = 1;
ip = 1;
for ip = 1:(nrows*ncols)
h(ip) = subplot(nrows,ncols,ip); ip = ip + 1;
end

isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,T,TS,EY)
shading(hca,'flat')
hca.XLabel.String = 't';
hca.YLabel.String = 't_s';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'E_y for E0 = E_0/t_s';

hca = h(isub); isub = isub + 1;
pcolor(hca,T,TS,cumtrapz(t,EY,2))
shading(hca,'flat')
hca.XLabel.String = 't';
hca.YLabel.String = 't_s';
hcb = colorbar('peer',hca);
hcb.YLabel.String = '\int E_y dt';

