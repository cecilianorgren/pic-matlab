% sort out distributions etc
% first time, first distribution, 3rd species
iTime = 2;
iDist = 2;
iSpecies = 1;

fxyz = ds100.fxyz(iTime,iDist,iSpecies);
f = ds100.f(iTime,iDist,iSpecies);

% density from moments
n_fxyz_mom = mean(mean(no02m.twpelim(fxyz.twpe).xlim(fxyz.x).zlim(fxyz.z).n(iSpecies)));
n_f_mom = mean(mean(no02m.twpelim(f.twpe).xlim(f.x).zlim(f.z).n(iSpecies)));

% density from distributions
n_fxyz_dist = sum(fxyz.f(:))*fxyz.dv^3; % n ~ int f*dv3
n_f_dist = sum(f.f(:))*f.dv^3; % n ~ int f*dv3
n_fyz_dist = sum(f.fyz(:))*f.dv^2; % n ~ int fyz*dv2
n_fxz_dist = sum(f.fxz(:))*f.dv^2; % n ~ int fxz*dv2
n_fxy_dist = sum(f.fxy(:))*f.dv^2; % n ~ int fyz*dv2

[n_fxyz_mom n_fxyz_dist n_fxyz_mom/n_fxyz_dist] % factor of 1.2 off (why?)
% Can it have something to do with the size of the sample box and range of 
% velocities? i.e. v's and x'z, z's ... PROBABLY NOT
% x2-x1 = 0.06
% (f.v(end)-f.v(1))*diff(f.x)
%
tmpind = cell(ds100.nt,1);
tmpind{iTime} = iDist;
vaxes = f.v;
vaxes = linspace(f.v(1),f.v(end),101);
fred = ds100.update_inds(tmpind).reduce_1d_new('x',iSpecies,vaxes);

n_fred_dist = sum(f.f(:))*f.dv^3; % n ~ int f*dv3
n_fredx_dist = sum(fred.fvx(:))*(fred.v(2)-fred.v(1)); % n ~ int fyz*dv2
n_fredy_dist = sum(fred.fvy(:))*(fred.v(2)-fred.v(1)); % n ~ int fxz*dv2
n_fredz_dist = sum(fred.fvz(:))*(fred.v(2)-fred.v(1)); % n ~ int fyz*dv2

ds = ds100.update_inds(tmpind);
xdist = (ds.xi1{1}+ds.xi2{1})/2;
zdist = (ds.zi1{1}+ds.zi2{1})/2;
dxdist = ds.xi1{1}-ds.xi2{1};
dzdist = ds.zi1{1}-ds.zi2{1};
tdist = repmat(twpe,size(xdist));
Bx = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'Bx');
By = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'By');
Bz = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'Bz');
% 
% Bx = 1;
% By = 0.5;
% Bz = 0;

fred_par = ds.reduce_1d_new('x',[iSpecies],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});

n_fredpar_dist = sum(fred_par.fvpar(:))*(fred_par.vpar_edges(2)-fred_par.vpar_edges(1)); % n ~ int fyz*dv2

[n_fxyz_mom n_fxyz_dist n_fredz_dist n_fredpar_dist]
1;

h = setup_subplots(3,3);
isub = 1;
if 1 % fxy
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,f.v,f.v,f.fxy')
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
end
if 1 % fxz
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,f.v,f.v,f.fxy')
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
end
if 1 % fyz
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,f.v,f.v,f.fyz')
  hca.XLabel.String = 'v_y';
  hca.YLabel.String = 'v_z';
end
if 1 % fx,fy,fz
  hca = h(isub); isub = isub + 1;  
  plot(hca,fred.v,fred.fvx,fred.v,fred.fvy,fred.v,fred.fvz)
  hca.XLabel.String = 'v';
  hca.YLabel.String = 'f';
end
if 1 % f(E)
  hca = h(isub); isub = isub + 1;  
  loglog(hca,fred_par.Epitch_center,squeeze(sum(fred_par.fpitchE,3)))
  hca.XLabel.String = 'E';
  hca.YLabel.String = 'f';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
















