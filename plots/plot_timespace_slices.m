%% Bz
zlim = 0.1*[-1 1];

Bz_z0 = mean(df04.zlim(zlim).Bz,3); % mean over zlim

h = setup_subplots(2,2);
isub = 1;


if 1
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).Bz,3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'B_z';
  hca.YLabel.String = 't\omega_{ci}';
  hca.Title.String = 'n_c = 0.4 n_0';
  hca.XLabel.String = 'x/d_i';
end
if 1
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).Bz,3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'B_z';
  hca.YLabel.String = 't\omega_{ci}';
  hca.Title.String = 'n_c = 0.8 n_0';
  hca.XLabel.String = 'x/d_i';
end
if 1
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.UB/pic(1).UB,mean(pic.zlim(zlim).Bz,3))
  hb = colorbar('peer',hca);
  hca.YDir = 'reverse';
  hb.YLabel.String = 'B_z';
  hca.YLabel.String = 'U_B/U_{B,t=0}';
  hca.Title.String = 'n_c = 0.4 n_0';
  hca.XLabel.String = 'x/d_i';
end
if 1
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.UB/pic(1).UB,mean(pic.zlim(zlim).Bz,3))
  hb = colorbar('peer',hca);
  hca.YDir = 'reverse';
  hb.YLabel.String = 'B_z';
  hca.YLabel.String = 'U_B/U_{B,t=0}';
  hca.Title.String = 'n_c = 0.8 n_0';
  hca.XLabel.String = 'x/d_i';
end
colormap(pic_colors('blue_red'))
hlinks = linkprop(h,{'CLim'}); hlinks.Targets(1).CLim = [-1 1];
hlinks = linkprop(h(1:2),{'YLim'}); hlinks.Targets(1).YLim = [0 max([df04.twci, df08.twci])];
hlinks = linkprop(h(3:4),{'YLim'}); hlinks.Targets(1).YLim = [min([df04.UB/df04(1).UB, df08.UB/df08(1).UB]) 1];

%% n
zlim = 0.0+0.1*[-1 1];

h = setup_subplots(2,2);
isub = 1;
levA = -25:1:0;
doA = 0;
if 1
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).n(1),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_p';
  hca.YLabel.String = 't\omega_{ci}';
  hca.Title.String = 'n_c = 0.4 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
if 1
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).n(1),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_p';
  hca.YLabel.String = 't\omega_{ci}';
  hca.Title.String = 'n_c = 0.8 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
if 1
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.UB/pic(1).UB,mean(pic.zlim(zlim).n(1),3))
  hb = colorbar('peer',hca);
  hca.YDir = 'reverse';
  hb.YLabel.String = 'B_z';
  hca.YLabel.String = 'U_B/U_{B,t=0}';
  hca.Title.String = 'n_c = 0.4 n_0';
  hca.XLabel.String = 'x/d_i';
  %hold(hca,'on')
  %contour(hca,pic.xi,pic.UB/pic(1).UB,mean(pic.zlim(zlim).A,3),levA,'k')
  %hold(hca,'off')
end
if 1
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.UB/pic(1).UB,mean(pic.zlim(zlim).n(1),3))
  hb = colorbar('peer',hca);
  hca.YDir = 'reverse';
  hb.YLabel.String = 'B_z';
  hca.YLabel.String = 'U_B/U_{B,t=0}';
  hca.Title.String = 'n_c = 0.8 n_0';
  hca.XLabel.String = 'x/d_i';
  %hold(hca,'on')
  %contour(hca,pic.xi,pic.UB/pic(1).UB,mean(pic.zlim(zlim).A,3),levA,'k')
  %hold(hca,'off')
end
colormap(pic_colors('candy'))

hlinks = linkprop(h,{'CLim'}); hlinks.Targets(1).CLim = 1 + 1*[-1 1];
hlinks = linkprop(h(1:2),{'YLim'}); hlinks.Targets(1).YLim = [0 max([df04.twci, df08.twci])];
hlinks = linkprop(h(3:4),{'YLim'}); hlinks.Targets(1).YLim = [min([df04.UB/df04(1).UB, df08.UB/df08(1).UB]) 1];

%% p
zlim = 0.0+0.1*[-1 1];

h = setup_subplots(2,2);
isub = 1;
levA = -25:1:0;
doA = 0;
if 1
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).p(1),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p_p';
  hca.YLabel.String = 't\omega_{ci}';
  hca.Title.String = 'n_c = 0.4 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
if 1
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).p(1),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p_p';
  hca.YLabel.String = 't\omega_{ci}';
  hca.Title.String = 'n_c = 0.8 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
if 1
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.UB/pic(1).UB,mean(pic.zlim(zlim).p(1),3))
  hb = colorbar('peer',hca);
  hca.YDir = 'reverse';
  hb.YLabel.String = 'p_p';
  hca.YLabel.String = 'U_B/U_{B,t=0}';
  hca.Title.String = 'n_c = 0.4 n_0';
  hca.XLabel.String = 'x/d_i';
  %hold(hca,'on')
  %contour(hca,pic.xi,pic.UB/pic(1).UB,mean(pic.zlim(zlim).A,3),levA,'k')
  %hold(hca,'off')
end
if 1
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.UB/pic(1).UB,mean(pic.zlim(zlim).p(1),3))
  hb = colorbar('peer',hca);
  hca.YDir = 'reverse';
  hb.YLabel.String = 'p_p';
  hca.YLabel.String = 'U_B/U_{B,t=0}';
  hca.Title.String = 'n_c = 0.8 n_0';
  hca.XLabel.String = 'x/d_i';
  %hold(hca,'on')
  %contour(hca,pic.xi,pic.UB/pic(1).UB,mean(pic.zlim(zlim).A,3),levA,'k')
  %hold(hca,'off')
end
colormap(pic_colors('candy'))

hlinks = linkprop(h,{'CLim'}); hlinks.Targets(1).CLim = 1 + 1*[-1 1];
hlinks = linkprop(h(1:2),{'YLim'}); hlinks.Targets(1).YLim = [0 max([df04.twci, df08.twci])];
hlinks = linkprop(h(3:4),{'YLim'}); hlinks.Targets(1).YLim = [min([df04.UB/df04(1).UB, df08.UB/df08(1).UB]) 1];

%% pB, pP, pD
zlim = 0.0+0.1*[-1 1];

h = setup_subplots(3,2);
sp04 = [1 3 5];
sp08 = [1 3];
isub = 1;
levA = -25:1:0;
doA = 1;
if 1 % pB
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).PB,3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p_B';
  hca.YLabel.String = 't\omega_{ci}';
  hca.Title.String = sprintf('n_c = 0.4 n_0, z = [%g,%g]',zlim(1),zlim(2));
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
if 1 % pB
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).PB,3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p_B';
  hca.YLabel.String = 't\omega_{ci}';
  hca.Title.String = 'n_c = 0.8 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
if 1 % pP
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).p(sp04),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p_{th}';
  hca.YLabel.String = 't\omega_{ci}';
  %hca.Title.String = 'n_c = 0.4 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
if 1 % pP
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).p(sp08),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p_{th}';
  hca.YLabel.String = 't\omega_{ci}';
  %hca.Title.String = 'n_c = 0.8 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
if 1 % pDyn
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).pD(sp04),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p_{dyn}';
  hca.YLabel.String = 't\omega_{ci}';
  %hca.Title.String = 'n_c = 0.4 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
if 1 % pDyn
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).pD(sp08),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p_{dyn}';
  hca.YLabel.String = 't\omega_{ci}';
  %hca.Title.String = 'n_c = 0.8 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA    
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
colormap(pic_colors('candy'))

c_eval('h(?).CLim = [0 0.3];',1:2)
c_eval('h(?).CLim = [0 1];',3:4)
c_eval('h(?).CLim = [0 0.4];',5:6)
hlinks = linkprop(h,{'XLim','YLim'});
%hlinks = linkprop(h(1:2),{'YLim'}); hlinks.Targets(1).YLim = [0 max([df04.twci, df08.twci])];
%hlinks = linkprop(h(3:4),{'YLim'}); hlinks.Targets(1).YLim = [min([df04.UB/df04(1).UB, df08.UB/df08(1).UB]) 1];
compact_panels(0.02)

%% spreading of plasma
zlim = 0.0+0.1*[-1 1];

h = setup_subplots(3,2);
sp04 = [1 3 5];
sp08 = [1 3];
isub = 1;
levA = -25:1:0;
doA = 1;
doVA = 1;
if doVA
  x0 = mean(df04.xi);
  xx = 0:10:x0;
  fun_t = @(x,t0) t0-1*x; 
end
if 1 % pB
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).PB,3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p_B';
  hca.YLabel.String = 't\omega_{ci}';
  hca.Title.String = sprintf('n_c = 0.4 n_0, z = [%g,%g] d_i',zlim(1),zlim(2));
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
if 1 % pB
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).PB,3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p_B';
  hca.YLabel.String = 't\omega_{ci}';
  hca.Title.String = sprintf('n_c = 0.8 n_0, z = [%g,%g] d_i',zlim(1),zlim(2));
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
end
if 1 % vix df04
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).vx(sp04),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'vx_{}';
  hca.YLabel.String = 't\omega_{ci}';
  %hca.Title.String = 'n_c = 0.4 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
  if doVA
    hold(hca,'on')
    for t0 = 0:50:550
      plot(hca,xx,fun_t(xx,t0),'--','color',0.6*[1 1 1])
      plot(hca,xx+x0,-fun_t(xx,t0-550+250),'--','color',0.6*[1 1 1])      
    end
    hold(hca,'off')    
  end
  irf_legend(hca,'ions',[0.02 0.98])
end
if 1 % vix df08
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).vx(sp08),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'vx_{}';
  hca.YLabel.String = 't\omega_{ci}';
  %hca.Title.String = 'n_c = 0.8 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
  if doVA
    hold(hca,'on')
    for t0 = 0:50:550
      plot(hca,xx,fun_t(xx,t0),'--','color',0.6*[1 1 1])
      plot(hca,xx+x0,-fun_t(xx,t0-550+250),'--','color',0.6*[1 1 1])      
    end
    hold(hca,'off')    
  end
  irf_legend(hca,'ions',[0.02 0.98])
end
if 1 % vex df04
  hca = h(isub); isub = isub + 1;
  pic = df04;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).vx(sp04+1),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'vx_{}';
  hca.YLabel.String = 't\omega_{ci}';
  %hca.Title.String = 'n_c = 0.4 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
  if doVA
    hold(hca,'on')
    for t0 = 0:50:550
      plot(hca,xx,fun_t(xx,t0),'--','color',0.6*[1 1 1])
      plot(hca,xx+x0,-fun_t(xx,t0-550+250),'--','color',0.6*[1 1 1])      
    end
    hold(hca,'off')    
  end
  irf_legend(hca,'electrons',[0.02 0.98])
end
if 1 % vex df08
  hca = h(isub); isub = isub + 1;
  pic = df08;
  imagesc(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).vx(sp08+1),3))
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'vx_{}';
  hca.YLabel.String = 't\omega_{ci}';
  %hca.Title.String = 'n_c = 0.8 n_0';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    contour(hca,pic.xi,pic.twci,mean(pic.zlim(zlim).A,3),levA,'k')
    hold(hca,'off')
  end
  if doVA
    hold(hca,'on')
    for t0 = 0:50:550
      plot(hca,xx,fun_t(xx,t0),'--','color',0.6*[1 1 1])
      plot(hca,xx+x0,-fun_t(xx,t0-550+250),'--','color',0.6*[1 1 1])      
    end
    hold(hca,'off')    
  end
  irf_legend(hca,'electrons',[0.02 0.98])
end
colormap(pic_colors('candy'))

c_eval('h(?).CLim = [0 0.3];',1:2)
c_eval('h(?).CLim = [-1 1]; colormap(h(?),pic_colors(''blue_red''))',3:4)
c_eval('h(?).CLim = 1.5*[-1 1]; colormap(h(?),pic_colors(''blue_red''))',5:6)
hlinks = linkprop(h,{'XLim','YLim'});
%hlinks = linkprop(h(1:2),{'YLim'}); hlinks.Targets(1).YLim = [0 max([df04.twci, df08.twci])];
%hlinks = linkprop(h(3:4),{'YLim'}); hlinks.Targets(1).YLim = [min([df04.UB/df04(1).UB, df08.UB/df08(1).UB]) 1];
compact_panels(0.02)

