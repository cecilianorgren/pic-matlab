%% Non-maxwellianity
pic = df04;
ds = ds04;
% Pick distribution
iSpecies = 1;
xval = 197;
zval = 2;
twpe = 7000;
ds = ds.twpelim(twpe).dxlim([0 0.21]).zfind(zval).xfind(xval);

n = mean(mean(pic.twpelim(twpe).xlim([ds.xi1{1} ds.xi2{1}]).zlim([ds.zi1{1} ds.zi2{1}]).n(iSpecies)));
vx = mean(mean(pic.twpelim(twpe).xlim([ds.xi1{1} ds.xi2{1}]).zlim([ds.zi1{1} ds.zi2{1}]).vx(iSpecies)));
vy = mean(mean(pic.twpelim(twpe).xlim([ds.xi1{1} ds.xi2{1}]).zlim([ds.zi1{1} ds.zi2{1}]).vy(iSpecies)));
vz = mean(mean(pic.twpelim(twpe).xlim([ds.xi1{1} ds.xi2{1}]).zlim([ds.zi1{1} ds.zi2{1}]).vz(iSpecies)));
P = mean(mean(pic.twpelim(twpe).xlim([ds.xi1{1} ds.xi2{1}]).zlim([ds.zi1{1} ds.zi2{1}]).p(iSpecies)));
T = P/n;

f = ds.f(1,1,3); % 3d and 2d
fxyz = ds.fxyz(1,1,iSpecies); % 3d: fout = fxyz(obj, it, id, iss, sumdim)
% units should be #/(vA^3*di^3), and we want to integrate over velocities,
% such that the units will be #/(di^3)

% Need for normalization of f
% vpost = vcode*(wpe/wce)*(mi/me)^0.5
v_norm = pic.wpewce*(pic.mime)^0.5;
x_norm = 1/sqrt(pic.mime);

v3x3_norm = v_norm^3*x_norm^2;
fxyz_norm = fxyz.f/v3x3_norm;


dv = fxyz.v(2)-fxyz.v(1);
dx = fxyz.x(2)-fxyz.x(1);
dz = fxyz.z(2)-fxyz.z(1);

%f_n = sum(fxyz_norm(:))*dv*dv*dv;

disp(sprintf('moments:      n = %8.5f, vx = %8.5f, vy = %8.5f, vz = %8.5f, T = %8.5f, P = %8.5f',n,vx,vy,vz,T,P))


[fn,fv,fT] = calculate_moments(fxyz); % Moments that go into a maxwellian distribution
disp(sprintf('distribution: n = %8.5f, vx = %8.5f, vy = %8.5f, vz = %8.5f, T = %8.5f',fn,fv.x,fv.y,fv.z,fT))
f = construct_maxwellian_from_f(fxyz); % make a maxwellian based on moments

% Calculate moments and compare to averaged moments over the same domain

