twpe = 3000:1000:7000;

for it = 1:numel(twpe)
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] = read_distributions(sprintf('/Users/cno062/Data/PIC/varying_guide_field/distributions/entire_box/twpe%05.f.dat',twpe(it)),4,101);

ispecies = 2;

[VX,VY,VZ] = ndgrid(axes(:,ispecies),axes(:,ispecies),axes(:,ispecies));
V2 = VX.^2 + VY.^2 + VZ.^2;
%hist(V2(:))
v2_edges = logspace(-3,3,30);

Fxyz = fxyz(:,:,:,ispecies);


[fv2 edges mid loc] = histcn(V2(:), v2_edges, 'AccumData', Fxyz(:), 'Fun', @mean);

loglog(mid{1},fv2)
hold on
end

%%
paths = {sprintf('/Users/cno062/Data/PIC/varying_guide_field/distributions/entire_box/twpe%05.f.dat',7000),...
         sprintf('/Users/cno062/Data/PIC/varying_guide_field/distributions/entire_box/twpe%05.f_vmax_i7_e20_nv101.dat',7000)};

ispecies = 1;
for ip = 1:2
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] = read_distributions(paths{ip},4,101);
  dv = axes(2,ispecies) - axes(1,ispecies);
  [VX,VY,VZ] = ndgrid(axes(:,ispecies),axes(:,ispecies),axes(:,ispecies));
  V2 = VX.^2 + VY.^2 + VZ.^2;
  v2_edges = logspace(-3,0.5,30);
  %v2_edges = linspace(0,100,100);
  Fxyz = fxyz(:,:,:,ispecies)/(dv*dv*dv);
  [fv2 edges mid loc] = histcn(V2(:), v2_edges, 'AccumData', Fxyz(:), 'Fun', @mean);

  res(ip).v = axes(:,ispecies);
  res(ip).fv2 = fv2;
  res(ip).v_center = mid{1};
  res(ip).v_edge = edges{1};
end

loglog(res(1).v_center,res(1).fv2,res(2).v_center,res(2).fv2);

%% compare runs
v = linspace(-20,20,100);
dv = v(2)-v(1);
ffun = @(vx,vy,vz) exp((-vx.^2-vy.^2-vz.^2)/10);

[VX,VY,VZ] = ndgrid(v,v,v);

F = ffun(VX,VY,VZ);
V2 = VX.^2 + VY.^2 + VZ.^2;

V2 = V2(:);
F = F(:);

v2edges_lin = linspace(0,0.1*max(V2),40);
v2edges_log = logspace(-3,log10(0.1*max(V2)),30);

vedges_lin = sqrt(v2edges_lin); % v = sqrt(v^2)
vedges_log = sqrt(v2edges_log);


vmid_lin = vedges_lin(2:end)-0.5*(vedges_lin(2:end)-vedges_lin(1:end-1));
vmid_log = vedges_log(2:end)-0.5*(vedges_log(2:end)-vedges_log(1:end-1));
dv_lin = diff(vedges_lin);
dv_log = diff(vedges_log);
vol_lin = vmid_lin.^2.*dv_lin;
vol_log = vmid_log.^2.*dv_log;

%[NV edges mid loc] = histcn(V2,v2edges);
[NF_lin edges mid_lin loc] = histcn(V2,v2edges_lin,'AccumData',F*dv.^3);
[NF_log edges mid_log loc] = histcn(V2,v2edges_log,'AccumData',F*dv.^3);
%[N edges mid loc] = histcn(V,vedges);

loglog(mid_lin{1},NF_lin./vol_lin',mid_log{1},NF_log./vol_log')
legend('lin','log')

%loglog(mid_lin{1},vol_lin,mid_log{1},NF_log./1')


%% compare times
twpe = 3000:1000:7000;
for ip = 1:numel(paths)  
  spl = strsplit(paths{ip},'/'); 
  files{ip} = spl{end};
end

clear res;
ispecies = [1 3];
%ispecies = [2 4];
for it = 1:numel(twpe)
  paths = {sprintf('/Users/cno062/Data/PIC/varying_guide_field/distributions_cn/entire_box/baseline/twpe%05.f_baseline.dat',twpe(it)),...
          sprintf('/Users/cno062/Data/PIC/varying_guide_field/distributions_cn/entire_box/twpe%05.f_vmax_i7_e20_nv101.dat',twpe(it))};
  for ip = 1:2
    [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] = read_distributions(paths{ip},4,101);
   
    v = axes(:,ispecies(1));
    vmax = max(v);
  
    v2edges = linspace(0,3*vmax^2,40);
    v2edges = logspace(-1,log10(3*vmax^2),100);
  
    f3d = sum(fxyz(:,:,:,ispecies),4)/numel(ispecies)/(v(2)-v(1))^3;
    [f,v2mid,vol] = v2_binning(v,f3d,v2edges);
  
    res(ip,it).v = v;
    res(ip,it).fv2 = f;
    res(ip,it).v_center = v2mid;  
    res(ip,it).vol = vol;  
    %res(ip,it)
  end
end

res(ip,it)
ip = 2;
hca = subplot(1,1,1);
loglog(hca,res(ip,1).v_center,res(ip,1).fv2,'-o',...
           res(ip,2).v_center,res(ip,2).fv2,'-o',...
           res(ip,3).v_center,res(ip,3).fv2,'-o',...
           res(ip,4).v_center,res(ip,4).fv2,'-o',...
           res(ip,5).v_center,res(ip,5).fv2,'-o');
hold(hca,'on')
plot(hca,[1 1]*vmax^2,hca.YLim,'--k')
plot(hca,[2 2]*vmax^2,hca.YLim,'--k')
plot(hca,[3 3]*vmax^2,hca.YLim,'--k')
hold(hca,'off')

ylim = hca.YLim;

hold(hca,'on')

switch ispecies(1)
  case 1
    plot(hca,v2edges,3*1e-1*exp(-v2edges/0.9),'k-')
    %plot(hca,v2edges,2*0.1e-3*exp(-v2edges/2.1),'k-')
  case 2
    plot(hca,v2edges,2.7*1e-2*exp(-v2edges/4.5),'k-')
    %plot(hca,v2edges,2*0.1e-4*exp(-v2edges/20),'k-')
end

hold(hca,'off')
hca.YLim = ylim;


hca.Title.String = [files{ip} ', species = [' num2str(ispecies) ']'];
hca.Title.Interpreter = 'none';

%legend(hca,files,'interpreter','none')

legend(hca,arrayfun(@(s) sprintf('twpe = %g ',s),twpe,'UniformOutput',false)','interpreter','none')


%% compare runs
twpe = 6000;
paths = {sprintf('/Users/cno062/Data/PIC/varying_guide_field/distributions_cn/entire_box/baseline/twpe%05.f_baseline.dat',twpe),...
         sprintf('/Users/cno062/Data/PIC/varying_guide_field/distributions_cn/entire_box/twpe%05.f_vmax_i7_e20_nv101.dat',twpe)};
for ip = 1:numel(paths)  
  spl = strsplit(paths{ip},'/'); 
  files{ip} = spl{end};
end




ispecies = [1 3];
%ispecies = [2 4];
for ip = 1:2
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] = read_distributions(paths{ip},4,101);
 
  v = axes(:,ispecies(1));
  vmax = max(v);

  v2edges = linspace(0,3*vmax^2,40);
  v2edges = logspace(-1,log10(3*vmax^2),100);

  f3d = sum(fxyz(:,:,:,ispecies),4)/numel(ispecies)/(v(2)-v(1))^3;
  [f,v2mid,vol] = v2_binning(v,f3d,v2edges);

  res(ip).v = v;
  res(ip).fv2 = f;
  res(ip).v_center = v2mid;  
  res(ip).vol = vol;  
end

hca = subplot(1,1,1);
loglog(hca,res(1).v_center,res(1).fv2,'-o',res(2).v_center,res(2).fv2,'-s');
hold(hca,'on')
plot(hca,[1 1]*vmax^2,hca.YLim,'--k')
plot(hca,[2 2]*vmax^2,hca.YLim,'--k')
plot(hca,[3 3]*vmax^2,hca.YLim,'--k')
hold(hca,'off')

ylim = hca.YLim;
%
hold(hca,'on')

switch ispecies(1)
  case 1
    plot(hca,v2edges,3*1e-1*exp(-v2edges/0.9),'k-')
    %plot(hca,v2edges,2*0.1e-3*exp(-v2edges/2.1),'k-')
  case 2
    plot(hca,v2edges,2.7*1e-2*exp(-v2edges/4.5),'k-')
    %plot(hca,v2edges,2*0.1e-4*exp(-v2edges/20),'k-')
end

hold(hca,'off')
hca.YLim = ylim;


hca.Title.String = ['t\omega_{pe} = ' num2str(twpe) ', species = [' num2str(ispecies) ']'];

legend(hca,files,'interpreter','none')
