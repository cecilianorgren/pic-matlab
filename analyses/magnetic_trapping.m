%ds100 = PICDist('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');
%no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
iIon = [1 3 5];
iEle = [2 4 6];
iIonCold = [3 5];
iEleCold = [4 6];
iIonHot = [1];
iEleHot = [2];

%%
twpe = 23000;
% Field region
xlim = [20 180]; 
zlim = [-12 12];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
nSpecies = numel(pic.mass);

%clear allA
strA = {'A=-4','A=-5','A=-6','A=-7','A=-8','A=-9'};
nA = numel(strA);
for iA = 1:nA % save data in structure table
  ds = ds100.twpelim(twpe).findtag(strA(iA));  
  allA(iA).ds = ds;
  allA(iA).x = (ds.xi1{1} + ds.xi2{1})/2;
  allA(iA).z = (ds.zi1{1} + ds.zi2{1})/2;
  allA(iA).dx = ds.xi2{1} - ds.xi1{1};
  allA(iA).dz = ds.zi2{1} - ds.zi1{1};
  allA(iA).arclength = [0 cumsum(sqrt(diff(allA(iA).x).^2 + diff(allA(iA).z).^2))];
  allA(iA).arccenter = (arclength(end)-arclength(1))/2;
  allA(iA).arc_z0 = arclength(find(allA(iA).z==min(allA(iA).z)));
  allA(iA).darc = allA(iA).arclength(2)-allA(iA).arclength(1);
  allA(iA).arcedges = [allA(iA).arclength(1)-0.5*allA(iA).darc, allA(iA).arclength+0.5*allA(iA).darc];
  range = allA(iA).dx(1)*0.5*[-1 1];  
  if 0
  allA(iA).Bx = pic.get_points(allA(iA).x,allA(iA).z,ds.twci,range,'Bx');
  allA(iA).By = pic.get_points(allA(iA).x,allA(iA).z,ds.twci,range,'By');
  allA(iA).Bz = pic.get_points(allA(iA).x,allA(iA).z,ds.twci,range,'Bz');
  allA(iA).Ex = pic.get_points(allA(iA).x,allA(iA).z,ds.twci,range,'Ex');
  allA(iA).Ey = pic.get_points(allA(iA).x,allA(iA).z,ds.twci,range,'Ey');
  allA(iA).Ez = pic.get_points(allA(iA).x,allA(iA).z,ds.twci,range,'Ez');
  cellB = {allA(iA).Bx,allA(iA).By,allA(iA).Bz};
  allA(iA).fred = ds.reduce_1d_new('x',[3 5],[],'vpar',cellB,'pitch',cellB);
  end
end
 
%% Plot
twpe = 24000;
xlim = [50 120]; 
zlim = [-12 12];
pic_lim = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);

nrows = nA+1;
ncols = 2;
h = setup_subplots(nrows,ncols,'horizontal');
isMap = [];
isArc = [];
isPitch = [];
isVpar = [];

doE = 0; colorE = [0.0 0.0 0.0];
doV = 0; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;

isub = 1;
if 1 % line position on map, ntop
  isMap(end+1) = isub; 
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic_lim.xi,pic_lim.zi,pic_lim.n(3)');
  colormap(hca,pic_colors('thermal'));
  hcb = colorbar('peer',hca);
  hca.CLim = max(max(get(findobj(hca.Children,'type','Image'),'CData')))*[-1 1];
  hcb.YLabel.String = 'n';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.CLim = [0 0.5];
  if 1 % plot_boxes
    hold(hca,'on')
    for iA = 1:nA
      allA(iA).ds.plot_boxes(hca);
    end
    hold(hca,'off')
  end
end
if 1 % line position on map, nbot
  isMap(end+1) = isub; 
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic_lim.xi,pic_lim.zi,pic_lim.n(5)');
  colormap(hca,pic_colors('thermal'));
  hcb = colorbar('peer',hca);
  hca.CLim = max(max(get(findobj(hca.Children,'type','Image'),'CData')))*[-1 1];
  hcb.YLabel.String = 'n';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.CLim = [0 0.5];
  if 1 % plot_boxes
    hold(hca,'on')
    for iA = 1:nA
      allA(iA).ds.plot_boxes(hca);
    end
    hold(hca,'off')
  end
end
for iA = 1:nA
  if 1 % log 10 fi35(vpar)
    hca = h(isub); isub = isub + 1;    
    isArc = [isArc isub-1];
    isVpar = [isVpar isub-1];
    %isDistEle = [isDistEle isub-1];
    
    fred = allA(iA).fred;
    arc = allA(iA).arclength;
    arcedges = allA(iA).arcedges;
    
    %surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',log10(fred.fvpar)')
    surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fred.fvpar')
    view(hca,[0 0 1]); 
    shading(hca,'flat')
    hca.XLabel.String = 'arclength (d_i)';
    hca.YLabel.String = 'v_{||}';
    %colormap(hca,pic_colors('candy4'))  
    colormap(hca,pic_colors('thermal'))  
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = ['log_{10}f_{i,cold}(l_{||},v_{||})'];
    %hca.CLim(2) = prctile(fred.fvpar(:),99);
    %hca.CLim = fi_clim;
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    %hca.YLim = [-2.5 2.5];
    %hca.CLim = [-6 -1];
    if 0 %doE
      hold(hca,'on')
      %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
      plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
      hold(hca,'off')
    end
    if 0% doV
      hold(hca,'on')
      plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
      hold(hca,'off')
    end
  end
  if 1 % log 10 fi35(pitchangle)
    hca = h(isub); isub = isub + 1;    
    isArc = [isArc isub-1];
    isPitch = [isPitch isub-1];
    
    fred = allA(iA).fred;
    arc = allA(iA).arclength;
    arcedges = allA(iA).arcedges;
        
    iE = 5;
    surf(hca,arcedges,fred.pitch,zeros(numel(arcedges),numel(fred.pitch))',log10(squeeze(mean(fred.fpitch(:,iE,:),2)))')
    
    view(hca,[0 0 1]); 
    shading(hca,'flat')
    hca.XLabel.String = 'arclength (d_i)';
    hca.YLabel.String = '\theta (\circ)';
    %colormap(hca,pic_colors('candy4'))  
    colormap(hca,pic_colors('thermal'))  
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = ['log_{10}f_{i,cold}(l_{||},v_{||})'];
    %hca.CLim(2) = prctile(fred.fvpar(:),99);
    %hca.CLim = fi_clim;
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    %hca.YLim = [-2.5 2.5];
    %hca.CLim = [-6 -1];
    if 0 %doE
      hold(hca,'on')
      %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
      plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
      hold(hca,'off')
    end
    if 0% doV
      hold(hca,'on')
      plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
      hold(hca,'off')
    end
  end
end

hlinks = linkprop(h(isVpar),{'YLim'});