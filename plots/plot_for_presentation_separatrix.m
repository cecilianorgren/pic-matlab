% plot for presentation
[saddle_locations,saddle_values] = saddle(A,'sort');
AX = saddle_values(1);
S = contourcs(x,z,A',AX*[1 1]*0.999);
[X,Z] = meshgrid(x,z);

xlim = [x(1) x(end)]+[3 -3]; 
zlim = [z(1) z(end)]; zlim = [-10 10];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:1:ipx2;
ipz = ipz1:1:ipz2;

% plot
nrows = 6;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1
  hca = h(isub); isub = isub + 1;  
  variable = ve12.x;
  imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)')
  %hca.Title.String = varstrs{ivar};
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_{e,x}/v_{A0}';
  hca.XLabel.String = 'x/d_{di0}';
  hca.YLabel.String = 'z/d_{di0}';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = 6*[-1 1];
  hold(hca,'on')
  contour(hca,x(1:5:end),z(1:5:end),A(1:5:end,1:5:end)',-25:2:0,'k')
  plot(hca,S(1).X,S(1).Y,'linewidth',3,'color',0*pic_colors('3')+[0 0 0])
  hold(hca,'off')
  hca.YDir = 'normal';
end

if 1
  hca = h(isub); isub = isub + 1;
  variable = ve12.par/vte12.par(100,100);
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);      
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.Title.String = varstrs{ivar};
  hca.YLabel.String = 'v_{e,||}/v_{te,0}';
  hca.XLabel.String = 'x/d_{di0}';
end
if 1
  hca = h(isub); isub = isub + 1;
  variable = ne12/ne12(100,100);
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);      
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.Title.String = varstrs{ivar};
  hca.YLabel.String = 'n_{e}/n_{0}';
  hca.XLabel.String = 'x/d_{di0}';
end
if 1 % flux/flux0
  hca = h(isub); isub = isub + 1;
  variable = ne12.*ve12.par/ne12(100,100)/vte12.par(100,100);
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);      
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.Title.String = varstrs{ivar};
  hca.YLabel.String = 'n_ev_{e,||}/n_0v_{te,0}';
  hca.XLabel.String = 'x/d_{di0}';
end
if 1 % vtepar/vt0
  hca = h(isub); isub = isub + 1;
  variable = vte12.par/vte12.par(100,100);
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);      
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.Title.String = varstrs{ivar};
  hca.YLabel.String = 'v_{te,||}/v_{te,0}';
  hca.XLabel.String = 'x/d_{di0}';
end
if 1 % Epar
  hca = h(isub); isub = isub + 1;
  variable = E.par;
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);      
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.Title.String = varstrs{ivar};
  hca.YLabel.String = 'E_{||}/B_0v_{A0}';
  hca.XLabel.String = 'x/d_{di0}';
end

irf_plot_axis_align
compact_panels(0.01)
hlink = linkprop(h,{'XLim'});
%%
varstrs = {'smooth2(vte12.par,np_smooth)/vte12.par(100,100)','ve12.par/vte12.par(100,100)','ne12/ne12(100,100)','ne12.*ve12.par/ne12(100,100)/vte12.par(100,100)'};




nvars = numel(varstrs);
vars = cell(nvars,1);

h = setup_subplots(2,nvars,'vertical');
isub = 1;
for ivar = 1:nvars
  variable = eval(varstrs{ivar});
  Vq = interp2(X,Z,variable',S(1).X,S(1).Y);
  vars{ivar} = Vq;
  hca = h(isub); isub = isub + 1;
  plot(hca,S(1).X,Vq)
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Title.String = varstrs{ivar};
  
  hca = h(isub); isub = isub + 1;
  imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)')
  hold(hca,'on')
  plot(hca,S(1).X,S(1).Y,'linewidth',2,'color',[0 0 0])
  hca.YDir = 'normal';
  hold(hca,'off')
  colormap(hca,pic_colors('blue_red'))
  hca.Title.String = varstrs{ivar};
  
end

hlink = linkprop(h,{'XLim'});