% compare loading of data
% paul, oxygen
[xe,ze,ex,ey,ez,bx,by,bz,dni,dne,jix,jiy,jiz,...
jex,jey,jez,vix,viy,viz,ti,te,nnx,nnz,...
dni_h,dne_h,jix_h,jiy_h,jiz_h,jex_h,jey_h,jez_h,vix_h,viy_h,viz_h,ti_h,te_h,a,...
wpewce,mass,pxxi_h,pyyi_h,pzzi_h,pxxe_h,pyye_h,pzze_h,vex_h,vey_h,vez_h,vex,vey,vez]...
= func_read_file_fields_oxygen([pathFields fileName]);

%% paul
[xe,ze,ex,ey,ez,bx,by,bz,dni,dne,jix,jiy,jiz,...
jex,jey,jez,vix,viy,viz,ti,te,nnx,nnz,...
dni_h,dne_h,jix_h,jiy_h,jiz_h,jex_h,jey_h,jez_h,vix_h,viy_h,viz_h,ti_h,te_h,a,...
wpewce,mass,pxxi_h,pyyi_h,pzzi_h,pxxe_h,pyye_h,pzze_h,vex_h,vey_h,vez_h,vex,vey,vez]...
= read_fields([pathFields fileName]);

%% mine
tic
[t,r,my_e,my_b,my_ni,my_ne,my_ve,my_vi,my_je,my_ji,...
  my_pe,my_pi,my_dfac,my_teti,my_nnx,my_nnz,my_wpewce,my_mass,my_it,my_dt,my_xmax,my_zmax,my_q] = read_fields_basic([pathFields fileName]);
toc
my_pb = (my_b.x.^2 + my_b.y.^2 + my_b.z.^2)/2;
my_pp = my_pe.scalar + my_pi.scalar;
my_ptot = my_pb + my_pp;

babs = sqrt(my_b.x.^2 + my_b.y.^2 + my_b.z.^2);
bmax = max(max(babs));

epar = (my_e.x.*bx + my_e.y.*by + my_e.z.*bz)./babs;
eabs = sqrt(my_e.x.^2 + my_e.y.^2 + my_e.z.^2);
emax = max(max(eabs));
%% plot for comparison
ncols = 4;
nrows = 3;
npanels = nrows*ncols;
for ipanel = 1:npanels
  h(ipanel) = subplot(nrows,ncols,ipanel);
end

isub = 1;
if 0 % bx
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,bx'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'bx';
  colorbar('peer',hca);
end
if 0 % bx my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_b.x'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my bx';
  colorbar('peer',hca);
end
if 0 % ez
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,ez'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'ez';
  colorbar('peer',hca);
end
if 0 % ez my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_e.z'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my ez';
  colorbar('peer',hca);
end
if 0 % dni
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,dni' + dni_h'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'dni + dni h';
  colorbar('peer',hca);
end
if 0 % ni my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_ni.tot'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my ni tot';
  colorbar('peer',hca);
end
if 0 % dne
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,dne' + dne_h'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'dne + dne h';
  colorbar('peer',hca);
end
if 0 % vex
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,vex'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'vex';
  colorbar('peer',hca);
  %hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.CLim = 2*[-1 1]; 
end
if 1 % vex my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_ve.x'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my vex';
  colorbar('peer',hca);  
  %hca.CLim = max(abs(hca.CLim))*[-1 1]; 
  hca.CLim = 10*[-1 1]; 
end
if 0 % viy
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,viy'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'viy';
  colorbar('peer',hca);
  %hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.CLim = max(abs(hca.CLim))*[-1 1]; 
end
if 0 % viy my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_vi.y'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my viy';
  colorbar('peer',hca);
  %hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.CLim = 10*[-1 1]; 
end
if 0 % jy
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,viy' - vey' + viy_h' - vey_h'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'jy + jy h';
  colorbar('peer',hca);
  %hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.CLim = max(abs(hca.CLim))*[-1 1]; 
  hca.CLim = 1.5*[-1 1]; 
end
if 1 % jy my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_ji.y' + my_je.y'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my jy';
  hcb = colorbar('peer',hca);
  %hca.CLim = max(abs(hca.CLim))*[-1 1];
  hcb.YLim = [-0.5 1.5];
  hca.CLim = [-1.5 1.5]; 
end
if 1 % e par my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,epar'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my e par';
  hcb = colorbar('peer',hca);
  %hca.CLim = max(abs(hca.CLim))*[-1 1];
  %hcb.YLim = [-0.5 1.5];
  %hca.CLim = [-1.5 1.5]; 
end
if 1 % b abs my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,babs'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my b abs';
  hcb = colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hcb.YLim(1) = 0;
  %hcb.YLim = [-0.5 1.5];
  %hca.CLim = [-1.5 1.5]; 
end
if 1 % pb my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_pb'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my pb';
  hcb = colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hcb.YLim(1) = 0;
end
if 1 % ne my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_ne.tot'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my ne tot';
  hcb = colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hcb.YLim(1) = 0;
end
if 1 % ni my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_ni.tot'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my ni tot';
  hcb = colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hcb.YLim(1) = 0;
end
if 0 % jiy
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,jiy' + jiy_h'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'jiy + jiy h';
  colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
if 0 % jiy my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_ji.y'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my jiy';
  colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
if 0 % jiy
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,jey'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'jey';
  colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
if 0 % jiy my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_je.y'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my jey';
  colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
if 0 % pizz
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,pzzi_h'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'pizz';
  colorbar('peer',hca);
end
if 0 % pizz my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_pi.zz1'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my pizz';
  colorbar('peer',hca);
end
if 0 % pezz
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,pzze_h'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'pezz';
  colorbar('peer',hca);
end
if 0 % pezz my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_pe.zz1'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my pezz';
  colorbar('peer',hca);
end
if 1 % pe my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,my_pe.scalar'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'pe';
  hcb = colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
if 1 % pi my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,my_pi.scalar'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'pi';
  hcb = colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
if 1 % pp my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,my_pp'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'pp';
  hcb = colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
if 1 % ptot my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,r.x,r.z,my_ptot'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'my pp tot';
  hcb = colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
if 1 % te my
  hca = h(isub); isub = isub + 1;
  imagesc(hca,xe,ze,my_pe.scalar'./my_ne.tot'); 
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca.Title.String = 'te';
  hcb = colorbar('peer',hca);
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
colormap(cn.cmap('blue_red'))