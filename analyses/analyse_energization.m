
pic = df04.zlim(2+[-0.2 0.2]);

vv1 = mean(pic.vv_diag(1),3);
vv3 = mean(pic.vv_diag(3),3);
vv4 = mean(pic.vv_diag(4),3);
vv6 = mean(pic.vv_diag(6),3);

%%
vv1_04 = nanmean(df04.zlim(0+[-0.2 0.2]).njp(1),3);
vv1_08 = nanmean(df08.zlim(0+[-0.2 0.2]).njp(1),3);
%%
h = setup_subplots(4,1);
isub = 1;

if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df04.xi,df04.twci,vv1_04)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'm_i\Sigma v_{ii}^2/3 (hot ions)';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 't\omega_{ci0}';
  hca.YDir = 'normal';
  hca.Title.String = 'n_c = 0.4 n_0';
end
if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df08.xi,df08.twci,vv1_08)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'm_i\Sigma v_{ii}^2/3 (hot ions)';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 't\omega_{ci0}';
  hca.YDir = 'normal';
  hca.Title.String = 'n_c = 0.8 n_0';
end
if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df04.xi,df04.UB/df04.i(1).UB,vv1_04)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'm_i\Sigma v_{ii}^2/3 (hot ions)';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'U_B/U_{B,t=0}';
  hca.YDir = 'reverse';
  hca.Title.String = 'n_c = 0.4 n_0';
end
if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df08.xi,df08.UB/df08.i(1).UB,vv1_08)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'm_i\Sigma v_{ii}^2/3 (hot ions)';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'U_B/U_{B,t=0}';
  hca.YDir = 'reverse';
  hca.Title.String = 'n_c = 0.8 n_0';
end
hlink = linkprop(h(1:4),{'CLim'});

%%
tic; [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = df04.zlim(0+[-0.1 0.1]).njp(1); toc;
p1_04 = mean((pxx.^2 + pyy.^2 + pzz.^2)/3,3);
t1_04 = mean((pxx.^2 + pyy.^2 + pzz.^2)/3./n,3);

tic; [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = df08.zlim(0+[-0.1 0.1]).njp(1); toc;
p1_08 = mean((pxx.^2 + pyy.^2 + pzz.^2)/3,3);
t1_08 = mean((pxx.^2 + pyy.^2 + pzz.^2)/3./n,3);

pB_04 = mean(df04.zlim(0+[-0.1 0.1]).PB,3);
pB_08 = mean(df08.zlim(0+[-0.1 0.1]).PB,3);
%%
h = setup_subplots(6,1);
isub = 1;

if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df04.xi,df04.twci,p1_04)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p (hot ions)';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 't\omega_{ci0}';
  hca.YDir = 'normal';
  hca.Title.String = 'n_c = 0.4 n_0';
end
if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df08.xi,df08.twci,p1_08)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p (hot ions)';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 't\omega_{ci0}';
  hca.YDir = 'normal';
  hca.Title.String = 'n_c = 0.8 n_0';
end
if 0
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df04.xi,df04.UB/df04.i(1).UB,p1_04)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p (hot ions)';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'U_B/U_{B,t=0}';
  hca.YDir = 'reverse';
  hca.Title.String = 'n_c = 0.4 n_0';
end
if 0
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df08.xi,df08.UB/df08.i(1).UB,p1_08)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p (hot ions)';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'U_B/U_{B,t=0}';
  hca.YDir = 'reverse';
  hca.Title.String = 'n_c = 0.8 n_0';
end
if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df04.xi,df04.twci,t1_04)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'T (hot ions)';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 't\omega_{ci0}';
  hca.YDir = 'normal';
  hca.Title.String = 'n_c = 0.4 n_0';
end
if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df08.xi,df08.twci,t1_08)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'T (hot ions)';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 't\omega_{ci0}';
  hca.YDir = 'normal';
  hca.Title.String = 'n_c = 0.8 n_0';
end
if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df04.xi,df04.twci,pB_04)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'P_B';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 't\omega_{ci0}';
  hca.YDir = 'normal';
  hca.Title.String = 'n_c = 0.4 n_0';
end
if 1
  hca = h(isub); isub = isub + 1;
  imagesc(hca,df08.xi,df08.twci,pB_08)
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'P_B';
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 't\omega_{ci0}';
  hca.YDir = 'normal';
  hca.Title.String = 'n_c = 0.8 n_0';
end
hlink_p = linkprop(h(1:2),{'CLim'});
%hlink_p.Targets{1}.CLim = [0 1]; 
hlink_t = linkprop(h(3:4),{'CLim'});
hlink_pb = linkprop(h(5:6),{'CLim'});
colormap(pic_colors('candy2'))


%%
it = 42;
tic; [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = df04.i(it).zlim(0+[-0.1 0.1]).njp(1); toc;
p1_04 = mean((pxx + pyy + pzz)/3,3);
t1_04 = mean((pxx + pyy + pzz)/3./n,3);
n1_04 = mean(n,3);

tic; [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = df04.i(it).zlim(0+[-0.1 0.1]).njp([3 5]); toc;
p35_04 = mean((pxx + pyy + pzz)/3,3);
t35_04 = mean((pxx + pyy + pzz)/3./n,3);
n35_04 = mean(n,3);

tic; [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = df08.i(it).zlim(0+[-0.1 0.1]).njp(1); toc;
p1_08 = mean((pxx + pyy + pzz)/3,3);
t1_08 = mean((pxx + pyy + pzz)/3./n,3);
n1_08 = mean(n,3);

tic; [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = df08.i(it).zlim(0+[-0.1 0.1]).njp(3); toc;
p3_08 = mean((pxx + pyy + pzz)/3,3);
t3_08 = mean((pxx + pyy + pzz)/3./n,3);
n3_08 = mean(n,3);


pB_04 = mean(df04.i(it).zlim(0+[-0.1 0.1]).PB,3);
pB_08 = mean(df08.i(it).zlim(0+[-0.1 0.1]).PB,3);

%% one time, cut in x
nrows = 6; 
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(npanels,1);
isub = 1;

if 1 % p 04
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.xi,pB_04,df04.xi,p1_04-0.6,df04.xi,p35_04)
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'p';
  hca.Title.String = 'n_c = 0.4 n_0';
end
if 1 % p 08
  hca = h(isub); isub = isub + 1;
  plot(hca,df08.xi,pB_08,df08.xi,p1_08-0.6,df08.xi,p3_08)
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'p';
  hca.Title.String = 'n_c = 0.8 n_0';
end
if 1 % t 04
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.xi,t1_04,df04.xi,t35_04)
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'T';
  hca.Title.String = 'n_c = 0.4 n_0';
end
if 1 % t 08
  hca = h(isub); isub = isub + 1;
  plot(hca,df08.xi,t1_08,df08.xi,t3_08)
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'T';
  hca.Title.String = 'n_c = 0.8 n_0';
end
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.xi,n1_04+n35_04,df04.xi,n1_04,df04.xi,n35_04)
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'n';
  hca.Title.String = 'n_c = 0.4 n_0';
end
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,df08.xi,n1_08+n3_08,df08.xi,n1_08,df08.xi,n3_08)
  hca.XLabel.String = 'x/d_{i0}';
  hca.YLabel.String = 'n';
  hca.Title.String = 'n_c = 0.8 n_0';
end

hlink = linkprop(h,{'XLim'});
h(1).XLim = [100 300];
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:npanels)


%% thermalization var(x,z), ions
for it = 5:5:60
%it = 40;
zlim = 5.99*[-1 1];
xlim = mean(df04.xi) + 100*[-1 1];
pic = df04.i(it).xlim(xlim).zlim(zlim);

tic; [n1,jx1,jy1,jz1,pxx1,pxy1,pxz1,pyy1,pyz1,pzz1] = pic.njp(1); toc;
p1 = (pxx1 + pyy1 + pzz1)/3;
t1 = p1./n1;
%uk1 = ;
ut1 = (3/2)*p1;

tic; [n3,jx3,jy3,jz3,pxx3,pxy3,pxz3,pyy3,pyz3,pzz3] = pic.njp([3]); toc;
p3 = (pxx3 + pyy3 + pzz3)/3;
t3 = p3./n3;
ut3 = (3/2)*p3;
%n3 = mean(n3,3);

tic; [n35,jx35,jy35,jz35,pxx35,pxy35,pxz35,pyy35,pyz35,pzz35] = pic.njp([3 5]); toc;
p35 = (pxx35 + pyy35 + pzz35)/3;
t35 = p35./n35;
ut35 = (3/2)*p35;

% Plot
nrows = 3;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 %  n1
  hca = h(isub); isub = isub + 1;
  var = squeeze(n1);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n1';
  hca.CLim(1) = 0;
end
if 1 %  n3
  hca = h(isub); isub = isub + 1;
  var = squeeze(n3);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n3';
  hca.CLim(1) = 0;
end
if 1 %  n35
  hca = h(isub); isub = isub + 1;
  var = squeeze(n35);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n35';
  hca.CLim(1) = 0;
end
if 1 %  p1
  hca = h(isub); isub = isub + 1;
  var = squeeze(p1);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p1';
  hca.CLim(1) = 0;
end
if 1 %  p3
  hca = h(isub); isub = isub + 1;
  var = squeeze(p3);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p3';
  hca.CLim(1) = 0;
end
if 1 %  p35
  hca = h(isub); isub = isub + 1;
  var = squeeze(p35);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p35';
  hca.CLim(1) = 0;
end
if 1 %  t1
  hca = h(isub); isub = isub + 1;
  var = squeeze(t1);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 't1';
  hca.CLim(1) = 0;
end
if 1 %  t3
  hca = h(isub); isub = isub + 1;
  var = squeeze(t3);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 't3';
  hca.CLim(1) = 0;
end
if 1 %  t35
  hca = h(isub); isub = isub + 1;
  var = squeeze(t35);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 't35';
  hca.CLim(1) = 0;
end

h(1).Title.String = sprintf('twci = %g',pic.twci);

doA = 1;
for ipanel = 1:npanels
  if doA
    levA = -25:1:0;
    hold(h(ipanel),'on')
    iAx = 1:5:pic.nx;
    iAz = 1:5:pic.nz;
    A = squeeze(pic.A);
    contour(h(ipanel),pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(h(ipanel),'off')
  end
end
hlink_n135 = linkprop(h([1 2 3]),{'CLim'}); hlink_n135.Targets(1).CLim = [0 1.8];
hlink_p35 = linkprop(h([5 6]),{'CLim'}); hlink_p35.Targets(1).CLim = [0 0.1];
hlink_t35 = linkprop(h([8 9]),{'CLim'}); hlink_t35.Targets(1).CLim = [0 0.2];
h(4).CLim = [0 1];
h(7).CLim = [0 1];
compact_panels(0.02,0.05)
cn.print(sprintf('npt_1_3_35_it=%g_ion',it))
end
%% thermalization var(x,z), electrons
for it = 5:5:60

%it = 40;
zlim = 7.99*[-1 1];
xlim = mean(df04.xi) + 100*[-1 1];
pic = df04.i(it).xlim(xlim).zlim(zlim);

tic; [n2,jx2,jy2,jz2,pxx2,pxy2,pxz2,pyy2,pyz2,pzz2] = pic.njp(2); toc;
p2 = (pxx2 + pyy2 + pzz2)/3;
t2 = p2./n2;
%uk2 = ;
ut2 = (3/2)*p2;

tic; [n4,jx4,jy4,jz4,pxx4,pxy4,pxz4,pyy4,pyz4,pzz4] = pic.njp([4]); toc;
p4 = (pxx4 + pyy4 + pzz4)/3;
t4 = p4./n4;
ut4 = (3/2)*p4;
%n4 = mean(n3,3);

tic; [n46,jx46,jy46,jz46,pxx46,pxy46,pxz46,pyy46,pyz46,pzz46] = pic.njp([4 6]); toc;
p46 = (pxx46 + pyy46 + pzz46)/3;
t46 = p46./n46;
ut46 = (3/2)*p46;

% % Plot
nrows = 3;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 %  n1
  hca = h(isub); isub = isub + 1;
  var = squeeze(n2);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n1';
  hca.CLim(1) = 0;
end
if 1 %  n3
  hca = h(isub); isub = isub + 1;
  var = squeeze(n4);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n3';
  hca.CLim(1) = 0;
end
if 1 %  n35
  hca = h(isub); isub = isub + 1;
  var = squeeze(n46);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n35';
  hca.CLim(1) = 0;
end
if 1 %  p1
  hca = h(isub); isub = isub + 1;
  var = squeeze(p2);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p1';
  hca.CLim(1) = 0;
end
if 1 %  p3
  hca = h(isub); isub = isub + 1;
  var = squeeze(p4);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p3';
  hca.CLim(1) = 0;
end
if 1 %  p35
  hca = h(isub); isub = isub + 1;
  var = squeeze(p46);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'p35';
  hca.CLim(1) = 0;
end
if 1 %  t1
  hca = h(isub); isub = isub + 1;
  var = squeeze(t2);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 't1';
  hca.CLim(1) = 0;
end
if 1 %  t3
  hca = h(isub); isub = isub + 1;
  var = squeeze(t4);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 't3';
  hca.CLim(1) = 0;
end
if 1 %  t35
  hca = h(isub); isub = isub + 1;
  var = squeeze(t46);
  imagesc(hca,pic.xi,pic.zi,var')
  hca.XLabel.String = 'x/di';
  hca.YLabel.String = 'z/di';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 't35';
  hca.CLim(1) = 0;
end

h(1).Title.String = sprintf('twci = %g',pic.twci);

doA = 1;
for ipanel = 1:npanels
  if doA
    levA = -25:1:0;
    hold(h(ipanel),'on')
    iAx = 1:5:pic.nx;
    iAz = 1:5:pic.nz;
    A = squeeze(pic.A);
    contour(h(ipanel),pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(h(ipanel),'off')
  end
end
hlink_n135 = linkprop(h([1 2 3]),{'CLim'}); hlink_n135.Targets(1).CLim = [0 1.8];
hlink_p35 = linkprop(h([5 6]),{'CLim'}); hlink_p35.Targets(1).CLim = [0 0.15];
hlink_t35 = linkprop(h([8 9]),{'CLim'}); hlink_t35.Targets(1).CLim = [0 0.15];
h(4).CLim = [0 1]/5;
h(7).CLim = [0 1]/5;
compact_panels(0.02,0.05)
cn.print(sprintf('npt_1_3_35_it=%g_ele',it))
end