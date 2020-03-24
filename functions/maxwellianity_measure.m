function out = maxwellianity_measure(v_inp,f_inp)
fin.v = v_inp;
fin.f = f_inp;
[n,v,p] = calculate_moments(fin);
t = p/n;
vt = sqrt(t);
vdx = v.x;
vdy = v.y;
vdz = v.z;

%fun_max = get_f_maxwellian()
%fun_fmax = f_maxwellian(n,[v.x v.y v.z],[vt vt vt],3);

nPop = numel(n);
f0_str = ['f0 = @(vx,vy,vz) ' sprintf('n(%g)*(1/pi./vt(%g).^2)^(3/2)*exp(-(vx-vdx(%g)).^2./vt(%g).^2-(vy-vdy(%g)).^2./vt(%g).^2-(vz-vdz(%g)).^2./vt(%g).^2)+',repmat((1:nPop),8,1))];
f0_str = [f0_str(1:end-1) ';'];
eval(f0_str)


[VX,VY,VZ] = ndgrid(v_inp,v_inp,v_inp);

f_max = f0(VX,VY,VZ);

eps = sqrt(sum(abs(f_max(:)-f_inp(:))));

out = eps;
if 1
  %%
  h = setup_subplots(3,2);
  isub = 1;
  
  hca = h(isub); isub = isub + 1;
  pcolor(hca,v_inp,v_inp,squeeze(sum(f_inp,3))')
  hb = colorbar('peer',hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  hca.Title.String = 'Simulation distribution, f_{sim}';
  irf_legend(hca,{sprintf('n = %.2f',n);sprintf('v = [%.2f, %.2f, %.2f]',vdx,vdy,vdz);sprintf('t = %.2f',t)},[0.02 0.98],'color',[0 0 0])
  
  hca = h(isub); isub = isub + 1;
  pcolor(hca,v_inp,v_inp,squeeze(sum(f_inp,2))') 
  hb = colorbar('peer',hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
  
  hca = h(isub); isub = isub + 1;
  pcolor(hca,v_inp,v_inp,squeeze(sum(f_max,3))') 
  hb = colorbar('peer',hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  hca.Title.String = {'Maxwellian distribution based on','moments from simulation distribution, f_{max}'};
  
  hca = h(isub); isub = isub + 1;
  pcolor(hca,v_inp,v_inp,squeeze(sum(f_max,2))') 
  hb = colorbar('peer',hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
  
  hca = h(isub); isub = isub + 1;
  pcolor(hca,v_inp,v_inp,squeeze(sum(f_max-f_inp,3))') 
  hb = colorbar('peer',hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  hca.Title.String = 'Difference: f_{max}-f_{sim}';
  irf_legend(hca,{sprintf('sqrt(sum(|f_{max}-f_{sim}|)) = %.0f',eps)},[0.02 0.98],'color',[0 0 0])
  
  hca = h(isub); isub = isub + 1;
  pcolor(hca,v_inp,v_inp,squeeze(sum(f_max-f_inp,2))') 
  hb = colorbar('peer',hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
  
  hlinks = linkprop(h,{'XLim','YLim','CLim'});
  for ip = 1:numel(h)
    shading(h(ip),'flat')
    h(ip).XGrid = 'on';
    h(ip).YGrid = 'on';
    h(ip).Layer = 'top';
    colormap(h(ip),pic_colors('blue_red'))
  end
  
  h(1).CLim = 3*max(max(squeeze(sum(f_max,2))))*[-1 1];
end