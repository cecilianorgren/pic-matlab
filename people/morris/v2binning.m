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