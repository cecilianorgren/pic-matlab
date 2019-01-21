txtfile = '/Users/cno062/Data/PIC/df_cold_protons/data/fields-02200.dat';
%%
[xe,ze,ex,ey,ez,bx,by,bz,dni,dne,jix,jiy,jiz,...
    jex,jey,jez,vix,viy,viz,ti,te,nnx,nnz, ...
    dni_h,dne_h,jix_h,jiy_h,jiz_h,jex_h,jey_h,jez_h,vix_h,viy_h,viz_h,ti_h,te_h,a, ...
    wpewce,mass,pxxi_h,pyyi_h,pzzi_h,pxxe_h,pyye_h,pzze_h,vex_h,vey_h,vez_h,vex,vey,vez] ...
    = func_read_file_fields_oxygen(txtfile);  
%%
[time,r,e,b,ni,ne,ve,vi,je,ji,pe,pi,dfac,teti,nnx,nnz,wpewce,mass,it,dt,xmax,zmax,q] = read_fields_basic(txtfile);  
  
%% Plot
nrows = 2;
ncols = 2;
npanels = nrows*ncols;
for ipanel = 1:npanels  
  h(ipanel) = subplot(nrows,ncols,ipanel);
end
isub = 1;

cA = [0 0 0];
doA = 0;
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
if 0 % bx
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
if 1 % vex
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ve.x2');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'v_{e,x}'; 

  hca.CLim = max(abs(ve.x(:)))*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,a',10,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vex
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,vex');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'v_{e,x}';

  hca.CLim = max(abs(vex_h(:)))*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,a',10,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vex
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,ve.x1');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'v_{e,x}';

  hca.CLim = max(abs(ve.x(:)))*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,a',10,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
end
if 1 % vex
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,xe,ze,vex_h');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.Title.String = 'v_{e,x}';

  hca.CLim = max(abs(vex_h(:)))*[-1 1];
  
  hcb = colorbar('peer',hca);
  colormap(hca,cn.cmap('blue_red'));
    
  if doA
    hold(hca,'on')
    hcont = contour(hca,xe,ze,a',10,'color',cA,'linewidth',1.0);  
    hold(hca,'off')  
  end
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