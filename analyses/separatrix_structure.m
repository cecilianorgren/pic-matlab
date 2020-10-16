% separatrix_structure
%df04n = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/data_h5/fields.h5');
tlim = [4 200];
zlim = [-2,10];
xlim = [130,280];
pic = df04n.twcilim(tlim).zlim(zlim).xlim(xlim);
%pic = df04n;
%sep_xz = pic.separatrix_location; 
Ey = squeeze(mean(pic.zlim([-0.2 0.2]).Ey,2));
ni = squeeze(mean(pic.zlim([-0.2 0.2]).ni,2));
nihot = squeeze(mean(pic.zlim([-0.2 0.2]).n(1),2));
nicold = ni-nihot;

%%
nrows = 5;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

zlevels = 0:0.5:15;
clim_re = [0 0.3];
cmap_re = pic_colors('waterfall');

if 1 % z(t,x)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep_xz.x(:,1)',sep_xz.twci,sep_xz.z'); 
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'z_{sep} (d_{i0})';
  hold(hca,'on')  
  contour(hca,sep_xz.x(:,1)',sep_xz.twci,sep_xz.z',zlevels,'k');
  hold(hca,'off')
  hold(hca,'on')  
  plot(hca,sep_xz.xline_x,sep_xz.twci,'k','linewidth',1.5);
  hold(hca,'off')  
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 't (\omega_{ci})';
  colormap(hca,pic_colors('waterfall'))
end
if 1 % z(t,x)/x(t,x)
  hca = h(isub); isub = isub + 1;
  xline_mat = repmat(tocolumn(sep_xz.xline_x),1,size(sep_xz.x,2));
  xline_mat = repmat(torow(sep_xz.xline_x),size(sep_xz.x,1),1);
  z_over_x = abs(sep_xz.z./(sep_xz.x-xline_mat));
  pcolor(hca,sep_xz.x(:,1),sep_xz.twci,z_over_x'); 
  drawnow;
  clim = hca.CLim;
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'z_{sep}/(x_{sep}-x_{xline})';
  hold(hca,'on')  
  contour(hca,sep_xz.x(:,1)',sep_xz.twci,sep_xz.z',zlevels,'k');
  hold(hca,'off')
  hold(hca,'on')  
  plot(hca,sep_xz.xline_x,sep_xz.twci,'k','linewidth',1.5);
  hold(hca,'off')  
  hca.CLim = clim_re;%prctile(z_over_x(:),99)*[-1 1];
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 't (\omega_{ci})';
  colormap(hca,cmap_re)
end
if 1 % atand(z(t,x)/x(t,x))
  hca = h(isub); isub = isub + 1;
  xline_mat = repmat(torow(sep_xz.xline_x),size(sep_xz.x,1),1);
  opening_angle = 2*atand(abs(sep_xz.z./(sep_xz.x-xline_mat)));
  pcolor(hca,sep_xz.x(:,1),sep_xz.twci,opening_angle'); 
  drawnow;
  clim = hca.CLim;
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = '\phi (deg)';
  hold(hca,'on')  
  contour(hca,sep_xz.x(:,1)',sep_xz.twci,sep_xz.z',zlevels,'k');
  hold(hca,'off')
  hold(hca,'on')  
  plot(hca,sep_xz.xline_x,sep_xz.twci,'k','linewidth',1.5);
  hold(hca,'off')  
  hca.CLim = [0 60];
  %hca.CLim = clim_re;%prctile(z_over_x(:),99)*[-1 1];
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 't (\omega_{ci})';
  colormap(hca,cmap_re)
end
if 1 % Ey(z=0)
  hca = h(isub); isub = isub + 1;
  %pcolor(hca,sep_xz.x(:,1),sep_xz.twci,nicold'); 
  pcolor(hca,pic.xi,pic.twci,Ey'); 
  drawnow;
  clim = hca.CLim;
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_y (v_{A0}B_0)';
  hold(hca,'on')  
  contour(hca,sep_xz.x(:,1)',sep_xz.twci,sep_xz.z',zlevels,'k');
  hold(hca,'off')
  hold(hca,'on')  
  plot(hca,sep_xz.xline_x,sep_xz.twci,'k','linewidth',1.5);
  hold(hca,'off')  
  hca.CLim = clim_re;
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 't (\omega_{ci})';
  colormap(hca,cmap_re)
end
if 1 % ni(z=0)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep_xz.x(:,1),sep_xz.twci,nicold'); 
  drawnow;
  clim = hca.CLim;
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_{i,cold} (n_0)';
  hold(hca,'on')  
  contour(hca,sep_xz.x(:,1)',sep_xz.twci,sep_xz.z',zlevels,'k');
  hold(hca,'off')
  hold(hca,'on')  
  plot(hca,sep_xz.xline_x,sep_xz.twci,'k','linewidth',1.5);
  hold(hca,'off')  
  hca.CLim = [0 1];
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 't (\omega_{ci})';
  colormap(hca,cmap_re)
end

hlinks = linkprop(h,{'XLim','YLim'});
compact_panels(0.01)
h(1).Title.String = pic.file;
h(1).Title.String = pic.file; h(1).Title.Interpreter = 'none';
%hlinks = linkprop(h,{'YLim'});

for ip = 1:npanels
  h(ip).XTick = 0:10:max(pic.xi);
  h(ip).YTick = 0:20:500;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).Layer = 'top';
  %h(ip).CLim = clim;
  h(ip).FontSize = 14;
end
