np_smooth = 2;
varstrs = {'smooth2(vte12.par,np_smooth)/vte12.par(100,100)','ve12.par/vte12.par(100,100)','ne12/ne12(100,100)','ne12.*ve12.par/ne12(100,100)/vte12.par(100,100)'};

%varstrs = {'ve12.z','vi12.z','ne12.*ve12.z','ni12.*vi12.z','ni12.*vi12.z-ne12.*ve12.z'};

[saddle_locations,saddle_values] = saddle(A,'sort');
AX = saddle_values(1);

S = contourcs(x,z,A',AX*[1 1]*0.999);

[X,Z] = meshgrid(x,z);

xlim = [x(1) x(end)]; 
zlim = [z(1) z(end)]; zlim = [-10 10];
ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:1:ipx2;
ipz = ipz1:1:ipz2;

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