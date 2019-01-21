%run_try.load_path
time = 00000;
fileName = sprintf('fields-%05.0f.dat',time);
%fileName = txtfile;

%[xe,ze,ex,ey,ez,bx,by,bz,n,vx,vy,vz,dfac,teti,nnx,nnz,wpewce,mass,it,dt,xmax,zmax,q,pxx,pyy,pzz,pxy,pxz,pyz] ...
%[xe,ze,ex,ey,ez,bx,by,bz,ni,ne,ve,vi,je,ji,pe,pi,dfac,teti,nnx,nnz,wpewce,mass,it,dt,xmax,zmax,q] = read_fields_basic([pathFields fileName]);  
[time,r,e,b,ni,ne,ve,vi,je,ji,pe,pi,dfac,teti,nnx,nnz,wpewce,mass,it,dt,xmax,zmax,q] = read_fields_basic(txtfile);  

if 0
 [xe,ze,ex,ey,ez,bx,by,bz,dni,dne,jix,jiy,jiz,...
    jex,jey,jez,vix,viy,viz,ti,te,nnx,nnz, ...
    dni_h,dne_h,jix_h,jiy_h,jiz_h,jex_h,jey_h,jez_h,vix_h,viy_h,viz_h,ti_h,te_h,a, ...
    wpewce,mass,pxxi_h,pyyi_h,pzzi_h,pxxe_h,pyye_h,pzze_h,vex_h,vey_h,vez_h,vex,vey,vez] ...
    = read_fields([pathFields fileName]);
end
%[axes,xlo,xhi,zlo,zhi,ic,fxyz ,fxy,fxz,fyz,vxa,vya,vza] ...
%  = read_parts([pathFields fileName]);

Ay = vector_potential(xe,ze,bx,bz);
nA = 10;
cA = [0.8 0.8 0.8];

babs = sqrt(bx.^2 + by.^2 + bz.^2);
bmax = max(max(babs));
b.abs = babs;

epar = (ex.*bx + ey.*by + ez.*bz)./babs;
eabs = sqrt(ex.^2 + ey.^2 + ez.^2);
emax = max(max(eabs));

ve.par = (ve.x.*bx + ve.y.*by + ve.z.*bz)./babs;
ve.abs = sqrt(ve.x.^2 + ve.y.^2 + ve.z.^2);
vemax = max(max(ve.abs))*0.5;

vi.par = (vi.x.*bx + vi.y.*by + vi.z.*bz)./babs;
vi.abs = sqrt(vi.x.^2 + vi.y.^2 + vi.z.^2);
vimax = max(max(vi.abs))*0.5;

%pdiag = pxx + pyy + pzz;
%pmax = max(max(max(pdiag)));

if emax == 0 
  emax = 1;
end

%% plot, initial current 
fig = figure(30);
nrows = 3;
ncols = 3;
npanels = nrows*ncols;
isub = 1; 
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel,'Parent',fig); isub = isub + 1;  
end

isub = 1;
if 0 % babs
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,babs');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '|B|';
 
  hca.CLim = bmax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  hcb.YLim = bmax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % bx
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,bx');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'B_x';

  hca.CLim = bmax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % dbx/dx
  hca = h(isub); isub = isub + 1;
  
  himag = imagesc(hca,xe,ze,bx');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'B_x';

  %hca.CLim = bmax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % vex
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ve.x');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'v_{e,x}';

  hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % vey
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ve.y');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'v_{e,y}';

  hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % vez
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ve.z');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'v_{e,z}';

  hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % jey
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,je.y');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'j_{e,y}';

  %hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % jiy
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ji.y');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'j_{i,y}';

  %hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % jey + jiy
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,je.y'+ji.y');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'j_{e,y}+j_{i,y}';

  %hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % babs^2
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,0.5*babs.^2');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'B^2';
 
  hca.CLim = bmax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  hcb.YLim = bmax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % pe scalar
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,pe.scalar');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'P_e';
 
  %hca.CLim = bmax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  %hcb.YLim = bmax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % pi scalar
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,pi.scalar');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'P_i';
 
  %hca.CLim = bmax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  %hcb.YLim = bmax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % pscalar
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,pe.scalar'+pi.scalar');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'P_e+P_i';
 
  %hca.CLim = bmax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  %hcb.YLim = bmax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % pe+pi+b^2
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,pe.scalar'+pi.scalar'+babs.^2');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'P_e+P_i+B^2';
 
  %hca.CLim = bmax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  %hcb.YLim = bmax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end

%% basic plot
nrows = 2;
ncols = 2;
npanels = nrows*ncols;
isub = 1; 
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
end

isub = 1;
if 0 % babs
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,babs');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '|B|';
 
  hca.CLim = bmax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  hcb.YLim = bmax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % bx
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,bx');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'B_x';

  hca.CLim = bmax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % by
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,by');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'B_y';

  hca.CLim = bmax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % bz
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,bz');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'B_z';

  hca.CLim = bmax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % eabs
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,eabs');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '|E|';
 
  hca.CLim = emax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  hcb.YLim = emax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % ex
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ex');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'E_x';

  hca.CLim = emax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % ey
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ey');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'E_y';

  hca.CLim = bmax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % ez
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ez');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'E_z';

  hca.CLim = emax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % vex
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ve.x');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'v_{e,x}';

  hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % vey
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ve.y');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'v_{e,y}';

  hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % vez
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ve.z');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'v_{e,z}';

  hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % jex
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,vj.x');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'j_{e,x}';

  hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % jey
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,je.y');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'j_{e,y}';

  %hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % jez
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,je.z');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'j_{e,z}';

  hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % jiy
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ji.y');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'j_{i,y}';

  %hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % vepar
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ve.par');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'v_{||}';

  hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % veabs
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ve.abs');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '|v_e|';

  hca.CLim = vemax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % babs^2
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,babs.^2');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'B^2';
 
  hca.CLim = bmax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  hcb.YLim = bmax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % pscalar
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,pe.scalar'+pi.scalar');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'P_e+P_i';
 
  %hca.CLim = bmax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  %hcb.YLim = bmax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % pe+pi+b^2
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,pe.scalar'+pi.scalar'+babs.^2');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'P_e+P_i+B^2';
 
  %hca.CLim = bmax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  %hcb.YLim = bmax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end

%% diagnostics plots
nrows = 3;
ncols = 4;
npanels = nrows*ncols;
isub = 1; 
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
end

isub = 1;
if 1 % babs
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,babs');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '|B|';
 
  hca.CLim = bmax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  hcb.YLim = bmax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % bx
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,bx');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'B_x';

  hca.CLim = bmax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % by
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,by');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'B_y';

  hca.CLim = bmax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 1 % bz
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,bz');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'B_z';

  hca.CLim = bmax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % eabs
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,eabs');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = '|E|';
 
  hca.CLim = emax*[-1 1];  
  
  hcb = colorbar('peer',hca);
  hcb.YLim = emax*[0 1];
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % ex
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ex');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'E_x';

  hca.CLim = emax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % ey
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ey');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'E_y';

  hca.CLim = bmax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % ez
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ez');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'E_z';

  hca.CLim = emax*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')  
end
if 0 % n1
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,squeeze(dns(:,:,1))');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'n_1';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')
end
if 0 % n2
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,squeeze(dns(:,:,2))');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'n_2';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')
end
if 0 % n3
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,squeeze(dns(:,:,3))');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'n_3';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')
end
if 0 % n4
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,squeeze(dns(:,:,4))');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'n_4';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')
end
if 0 % pdiag1 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,squeeze(pdiag(:,:,1))');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'pdiag_1';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = pmax*[-1 1];
  hcb.YLim = pmax*[0 1];
  
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')
end
if 0 % pdiag2 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,squeeze(pdiag(:,:,2))');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'pdiag_2';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = pmax*[-1 1];
  hcb.YLim = pmax*[0 1];
  
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')
end
if 0 % pdiag3 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,squeeze(pdiag(:,:,3))');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'pdiag_3';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = pmax*[-1 1];
  hcb.YLim = pmax*[0 1];
  
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')
end
if 0 % pdiag4 
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,squeeze(pdiag(:,:,4))');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'pdiag_4';

  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
  hca.CLim = pmax*[-1 1];
  hcb.YLim = pmax*[0 1];
  
  hold(hca,'on')
  hcont = contour(hca,xe,ze,Ay',10,'color',cA,'linewidth',1.0);  
  hold(hca,'off')
end
