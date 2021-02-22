 function fout = rotate_dist(f,rx,ry,rz)
% f is a structure with fields v and 3D f(vx,vy,vz)

doInterp = 0;
if doInterp
  f_old = f.f;
  v_old = f.v;
  v_new = v_old; % this cuts of the corners, but there's never many particles out there anyway.
  
  nv_old = numel(v_old);
  nv_interp = 2*nv_old;
  v_old_interp = linspace(v_old(1),v_old(end),nv_interp);
  % interpolate to finer grid
  [X,Y,Z] = meshgrid(v_old,v_old,v_old);
  [Xq,Yq,Zq] = meshgrid(v_old_interp,v_old_interp,v_old_interp);
  f_old_interp = interp3(X,Y,Z,f_old,Xq,Yq,Zq);
                
  v_old =  v_old_interp;
  f_old = f_old_interp;
  
  dv_old = v_old(2) - v_old(1);
  v_new = v_old;
  
  dv = v_new(2) - v_new(1);
  v_new_edges = [v_new-0.5*dv v_new(end)+0.5*dv];
  % make sure new coordinate system is orthogonal and made up of unit vectors
  rx = rx./sqrt(sum(rx.^2));
  ry = cross(rx,cross(ry,rx)); ry = ry./sqrt(sum(ry.^2));
  rz = cross(rx,ry); rz = rz./sqrt(sum(rz.^2));

  [VX_old,VY_old,VZ_old] = ndgrid(v_old,v_old,v_old);

  [VX_new,VY_new,VZ_new] = rotate_xyz(VX_old,VY_old,VZ_old,rx,ry,rz);

  % bin

  V_new = [VX_new(:),VY_new(:),VZ_new(:)];
  
  [count edges mid loc] = histcn(V_new,v_new_edges,v_new_edges,v_new_edges,'AccumData', f_old(:)*dv_old^3);
  f_new = count/(dv^3);

  fout.f = f_new;
  fout.v = v_new;
  fout.dv = dv;
  fout.r1 = rx;
  fout.r2 = ry;
  fout.r3 = rz;
else
  f_old = f.f;
  v_old = f.v;
  dv_old = v_old(2) - v_old(1);
  v_new = v_old; % this cuts of the corners, but there's never many particles out there anyway.

  dv = v_new(2) - v_new(1);
  v_new_edges = [v_new-0.5*dv v_new(end)+0.5*dv];
  % make sure new coordinate system is orthogonal and made up of unit vectors
  rx = rx./sqrt(sum(rx.^2));
  ry = cross(rx,cross(ry,rx)); ry = ry./sqrt(sum(ry.^2));
  rz = cross(rx,ry); rz = rz./sqrt(sum(rz.^2));

  [VX_old,VY_old,VZ_old] = ndgrid(v_old,v_old,v_old);

  [VX_new,VY_new,VZ_new] = rotate_xyz(VX_old,VY_old,VZ_old,rx,ry,rz);

  % bin

  V_new = [VX_new(:),VY_new(:),VZ_new(:)];

  [count edges mid loc] = histcn(V_new,v_new_edges,v_new_edges,v_new_edges,'AccumData', f_old(:)*dv_old^3);
  f_new = count/(dv^3);

  fout.f = f_new;
  fout.v = v_new;
  fout.dv = dv;
  fout.r1 = rx;
  fout.r2 = ry;
  fout.r3 = rz;

end
doPlot = 0;
if doPlot
  %%
  nrows = 3;
  ncols = 2;
  npanels = nrows*ncols;
  h = setup_subplots(3,2,'vertical');
  isub = 1;

  % old
  hca = h(isub); isub = isub + 1;
  imagesc(hca,v_old,v_old,squeeze(sum(f_old,3))')
  hca.XLabel.String = 'v_{x}^{old}';
  hca.YLabel.String = 'v_{y}^{old}';

  hca = h(isub); isub = isub + 1;
  imagesc(hca,v_old,v_old,squeeze(sum(f_old,2))')
  hca.XLabel.String = 'v_{x}^{old}';
  hca.YLabel.String = 'v_{z}^{old}';

  hca = h(isub); isub = isub + 1;
  imagesc(hca,v_old,v_old,squeeze(sum(f_old,1))')
  hca.XLabel.String = 'v_{y}^{old}';
  hca.YLabel.String = 'v_{z}^{old}';

  % new
  hca = h(isub); isub = isub + 1;
  imagesc(hca,v_old,v_old,squeeze(sum(f_new,3))')
  hca.XLabel.String = 'v_{x}^{new}';
  hca.YLabel.String = 'v_{y}^{new}';

  hca = h(isub); isub = isub + 1;
  imagesc(hca,v_old,v_old,squeeze(sum(f_new,2))')
  hca.XLabel.String = 'v_{x}^{new}';
  hca.YLabel.String = 'v_{z}^{new}';

  hca = h(isub); isub = isub + 1;
  imagesc(hca,v_old,v_old,squeeze(sum(f_new,1))')
  hca.XLabel.String = 'v_{y}^{new}';
  hca.YLabel.String = 'v_{z}^{new}';

  for ip = 1:npanels
    h(ip).YDir = 'normal';
    axis(h(ip),'square')
  end
end