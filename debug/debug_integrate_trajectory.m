t0 = 120;
tspan = [60,t0,200];
%tspan = [t0-1,t0,t0+1];
tic;
tr1 = df04.integrate_trajectory_old(r0,v0,tspan,m,q);    
[Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB(tr1.x,tr1.z,tr1.t);  % interpolate
tr1.t0 = t0;
tr1.Ex = Ex; 
tr1.Ey = Ey;
tr1.Ez = Ez;
tr1.Bx = Bx;
tr1.By = By;
tr1.Bz = Bz;
toc

tic
tr2 = df04.integrate_trajectory(r0,v0,tspan,m,q); 
toc

%% the interpolation is different
h = setup_subplots(9,1);
isub = 1;
hca = h(isub); isub = isub + 1; plot(hca,tr1.x,tr1.z,'o',tr2.x,tr2.z,'x')
hca = h(isub); isub = isub + 1; plot(hca,tr1.t,tr1.x,'o',tr2.t,tr2.x,'x')
hca = h(isub); isub = isub + 1; plot(hca,tr1.t,tr1.z,'o',tr2.t,tr2.z,'x')
hca = h(isub); isub = isub + 1; plot(hca,tr1.t,tr1.Ex,'-',tr2.t,tr2.Ex)
hca = h(isub); isub = isub + 1; plot(hca,tr1.t,tr1.Ey,'-',tr2.t,tr2.Ey)
hca = h(isub); isub = isub + 1; plot(hca,tr1.t,tr1.Ez,'-',tr2.t,tr2.Ez)
hca = h(isub); isub = isub + 1; plot(hca,tr1.t,tr1.Bx,'-',tr2.t,tr2.Bx)
hca = h(isub); isub = isub + 1; plot(hca,tr1.t,tr1.By,'-',tr2.t,tr2.By)
hca = h(isub); isub = isub + 1; plot(hca,tr1.t,tr1.Bz,'-',tr2.t,tr2.Bz)