%% Reconnection rate
pic = gf05.xlim([-90 110]).zlim([-2 2]);
for it = 1:pic.nt
  A = pic(it).A;  
  [saddle_locations,saddle_values] = saddle(A,'sort');
  xX = pic.xi(saddle_locations(1,1));
  zX = pic.zi(saddle_locations(1,2));
  Ey_1di = pic(it).xlim(xX+0.5*[-1 1]).zlim(zX+0.5*[-1 1]).Ey;
  Ey_2di = pic(it).xlim(xX+1*[-1 1]).zlim(zX+1*[-1 1]).Ey;  
  R.twpe(it) = pic(it).twpe;
  R.xXline(it) = xX;
  R.zXline(it) = zX;
  R.AXline(it) = saddle_values(1);
  R.E_1di(it) = mean(Ey_1di(:));
  R.E_2di(it) = mean(Ey_2di(:));  
end
%R.A = interp1(R.twpe(2:end)-diff(R.twpe(1:end)),diff(R.AXline),R.twpe)/20;

[AX,H1,H2] = plotyy(R.twpe,[R.E_1di;R.E_2di],R.twpe,R.xXline);
AX(1).YLabel.String = 'R_E (v_{A0}B_0)';
hl = legend(AX(1),{'1x1 d_{i0}','2x2 d_{i0}'});
hl.Title.String = 'Box size';
AX(2).YLabel.String = 'x_{x line} (d_{i0})';
AX(2).YLabel.String = 'x_{xline} (d_{i0})';
AX(1).XLabel.String = 'tw_{pe}';

%% Some plot
twpe = [1400 2200 3000 3800 4600 5400];
twpe = [3000 3800 4600];
%twpe = gf05.twpe;
pic = gf05.twpelim(twpe,'exact').xlim([55 145]).zlim([-7 7]);

nrows = pic.nt;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;

doA = 1;
doAx = 0; % plot separatrix, NOT impl.
if doA % Flux function, set parameters
  cA = 0*[0.8 0.8 0.8];
  nA = 20;
  nA = [-50:1:50];
  ipxA = 1:10:pic.nx;
  ipzA = 1:10:pic.nz;     
end
  
links = cell(10,1); 
hbs = gobjects(0);

for it = 1:pic.nt
  pic_tmp = pic(it);
  varcount = 0;
  if doA
    A_tmp = pic_tmp.A;
  end
  if 0 % By
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;
    links{varcount} = [links{varcount} isub - 1]; 
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.By');   
    hb = colorbar('peer',hca);
    hca.CLim = 0.5 + 0.5*[-1 1];
    colormap(hca,pic_colors('blue_red'));
    hb.YLabel.String = 'B_{y}';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end  
    irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
  end
  if 0 % Ey
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;
    links{varcount} = [links{varcount} isub - 1]; 
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.Ey');   
    hb = colorbar('peer',hca);
    hca.CLim = [0 0.2];
    colormap(hca,pic_colors('candy'));
    hb.YLabel.String = 'B_{y}';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end  
    irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
  end
  if 0 % vabs hot e
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;    
    links{varcount} = [links{varcount} isub - 1];     
    vx = pic_tmp.vx(2);
    vy = pic_tmp.vy(2);
    vz = pic_tmp.vz(2);
    vabs = sqrt(vx.^2 + vy.^2 +vz.^2);
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,vabs');
    hb = colorbar('peer',hca);
    hca.CLim = [0 5];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = '|v_{eh}|';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
    irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
  end
  if 0 % vabs cold e
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;    
    links{varcount} = [links{varcount} isub - 1];     
    vx = pic_tmp.vx(4);
    vy = pic_tmp.vy(4);
    vz = pic_tmp.vz(4);
    vabs = sqrt(vx.^2 + vy.^2 +vz.^2);
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,vabs');
    hb = colorbar('peer',hca); hbs(end+1) = hb;
    hca.CLim = [0 5];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = '|v_{ec}|';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
    irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
  end
  if 0 % vabs/wce
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;    
    links{varcount} = [links{varcount} isub - 1]; 
    Bx = pic_tmp.Bx;
    By = pic_tmp.By;
    Bz = pic_tmp.Bz;
    Babs = sqrt(Bx.^2 + By.^2 +Bz.^2);
    vx = pic_tmp.vx(4);
    vy = pic_tmp.vy(4);
    vz = pic_tmp.vz(4);
    vabs = sqrt(vx.^2 + vy.^2 +vz.^2);
    me = 1;
    wce = Babs/me;
    rs = vabs./wce;                    
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,rs');
    hb = colorbar('peer',hca); hbs(end+1) = hb;     
    hca.CLim = [0 10];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = '|v_{ec}|/\omega_{ce} (d_{i0})';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
    irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
  end
  if 0 % 1/KB
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;    
    links{varcount} = [links{varcount} isub - 1]; 
    Bx = pic_tmp.Bx;
    By = pic_tmp.By;
    Bz = pic_tmp.Bz;
    KB = magnetic_field_curvature(pic_tmp.xi,pic_tmp.zi,Bx,By,Bz); % magnetic curvature
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,1./KB.abs');
    hb = colorbar('peer',hca); hbs(end+1) = hb;     
    hca.CLim = [0 10];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = sprintf('1/|K| (d_{i0})',KB.units);
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
    irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
  end
  if 0 % KB
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;    
    links{varcount} = [links{varcount} isub - 1]; 
    Bx = pic_tmp.Bx;
    By = pic_tmp.By;
    Bz = pic_tmp.Bz;
    KB = magnetic_field_curvature(pic_tmp.xi,pic_tmp.zi,Bx,By,Bz); % magnetic curvature
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,KB.abs');
    hb = colorbar('peer',hca);
    hca.CLim = [0 0.4];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = sprintf('|K| (d_{i0}^{-1})',KB.units);
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
    irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
  end
  if 0 % KBvabs/wce
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;    
    links{varcount} = [links{varcount} isub - 1]; 
    Bx = pic_tmp.Bx;
    By = pic_tmp.By;
    Bz = pic_tmp.Bz;
    Babs = sqrt(Bx.^2 + By.^2 +Bz.^2);
    vx = pic_tmp.vx(4);
    vy = pic_tmp.vy(4);
    vz = pic_tmp.vz(4);
    vabs = sqrt(vx.^2 + vy.^2 +vz.^2);
    me = 1;
    wce = Babs/me;
    rs = vabs./wce;            
    KB = magnetic_field_curvature(pic_tmp.xi,pic_tmp.zi,Bx,By,Bz); % magnetic curvature
    rsKB = rs.*KB.abs;
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,rsKB');
    hb = colorbar('peer',hca); hbs(end+1) = hb;
    hca.CLim = [0 1];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = '|K||v_e|/wce ';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
    irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
  end
  if 1 % nec
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;
    links{varcount} = [links{varcount} isub - 1]; 
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.n(4)');
    hb = colorbar('peer',hca); hbs(end+1) = hb;
    hca.CLim = [0 0.4];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = 'n_{e,c}';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
    irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
  end
  if 0 % nic
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;
    links{varcount} = [links{varcount} isub - 1]; 
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.n(3)');
    hb = colorbar('peer',hca);
    hca.CLim = [0 0.5];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = 'n_{i,c}';
  end
  if 0 % neh
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;
    links{varcount} = [links{varcount} isub - 1]; 
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.n(2)');
    hb = colorbar('peer',hca);
    hca.CLim = [0 1.5];
    colormap(hca,pic_colors('waterfall'));
  end
  if 0 % pexx
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;
    links{varcount} = [links{varcount} isub - 1];     
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.pxx(4)');
    hb = colorbar('peer',hca);
    hca.CLim = [0 0.01];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = 'p_{ec,xx}';
  end
  if 0 % peyy
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;
    links{varcount} = [links{varcount} isub - 1];     
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.pyy(4)');
    hb = colorbar('peer',hca);
    hca.CLim = [0 0.01];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = 'p_{ec,yy}';
  end
  if 0 % pezz
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;
    links{varcount} = [links{varcount} isub - 1];     
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.pzz(4)');
    hb = colorbar('peer',hca);
    hca.CLim = [0 0.01];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = 'p_{ec,zz}';
  end
  if 1 % pec
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;
    links{varcount} = [links{varcount} isub - 1];     
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.p(4)');
    hb = colorbar('peer',hca); hbs(end+1) = hb;
    hca.CLim = [0 0.05];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = 'p_{ec}';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
  end
  if 0 % pic
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;
    links{varcount} = [links{varcount} isub - 1];     
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.p(3)');
    hb = colorbar('peer',hca); hbs(end+1) = hb;
    hca.CLim = [0 0.5];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = 'p_{ic}';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
  end
  if 1 % tec
    hca = h(isub); isub = isub + 1;    
    varcount = varcount + 1;
    links{varcount} = [links{varcount} isub - 1];     
    p = pic_tmp.p(3);
    n = pic_tmp.n(3);
    t = p./n;
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,t');
    hb = colorbar('peer',hca); hbs(end+1) = hb;
    hca.CLim = [0 0.2];
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = 't_{ic}';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
  end
  drawnow
end


for ilink = 1:varcount
  linkprop(h(links{ilink}),{'CLim'});
end

linkprop(h,{'XLim','YLim'});
drawnow

if 1
  %%
for ilink = 1:varcount
  htop = hbs(links{ilink}(1));
  htop.Location = 'northoutside';
  for it = 2:pic.nt
    delete(hbs(links{ilink}(it)));    
  end
end
end
%h(1).CLim = [0 0.4];
%h(2).CLim = 0.5 + 0.5*[-1 1];

for ip = 1:npanels
  h(ip).YDir = 'normal';
  h(ip).XLabel.String = 'x (di)';
  h(ip).YLabel.String = 'z (di)';
end
drawnow
compact_panels(0.01,0.01)

for ivar = 2:ncols
  for it = 1:pic.nt
    h(links{ivar}(it)).YLabel.String = '';
    h(links{ivar}(it)).YTickLabel = []';
  end
end
%% Get X line location

for it = 1:pic.nt
  pic_tmp = pic(it);
  [saddle_locations,saddle_values] = saddle(pic_tmp.A,'sort');
  xXline(it) = pic_tmp.xi(saddle_locations(1));
  zXline(it) = pic_tmp.xi(saddle_locations(2));
end

%% Check potential numerical heating 

%% Plot for VR application, single time
twpe = [5000];
pic = gf05.twpelim(twpe,'exact').xlim([65 135]).zlim([-6 6]);
pic_tmp = pic;

nrows = 3;
ncols = 2;
npanels = nrows*ncols;
%h = setup_subplots(nrows,ncols,'horizontal');
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

doA = 1;
doAx = 0; % plot separatrix, NOT impl.
if doA % Flux function, set parameters
  cA = 0*[0.8 0.8 0.8];
  nA = 20;
  nA = [-50:1:50];
  ipxA = 1:10:pic.nx;
  ipzA = 1:10:pic.nz;     
  A_tmp = pic.A;
end
  
links = cell(10,1); 
hbs = gobjects(0);

%
KB = pic.magnetic_curvature;
Bx = pic_tmp.Bx;
By = pic_tmp.By;
Bz = pic_tmp.Bz;
Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
b.x =  Bx./Babs;
b.y =  By./Babs;
b.z =  Bz./Babs;

wce = Babs/(pic.mass(4)/pic.mime);

Ex = pic_tmp.Ex;
Ey = pic_tmp.Ey;
Ez = pic_tmp.Ez;
Epar = (Ex.*b.x + Ey.*b.y + Ez.*b.z);


if 0 % make par/perp temperatures
  %%
  tic
  clear p t
  p.xx = pic_tmp.pxx(4);
  p.yy = pic_tmp.pyy(4);
  p.zz = pic_tmp.pzz(4);
  n = pic_tmp.n(4);
  t.xx = p.xx./n;
  t.yy = p.yy./n;
  t.zz = p.zz./n;
  
  t.xy = t.xx*0;
  t.xz = t.xx*0;
  t.yz = t.xx*0;
  
  r1 = b; % magnetic field unit vector
  r2 = cross_product(r1.x,r1.y,r1.z,0,1,0);
  r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);
  r2.x = r2.x./r2.abs;
  r2.y = r2.y./r2.abs;
  r2.z = r2.z./r2.abs;
  r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);
  r2 = cross_product(r2.x,r2.y,r2.z,r1.x,r1.y,r1.z);
  r3 = cross_product(r1.x,r1.y,r1.z,r2.x,r2.y,r2.z);
  r3.abs = sqrt(r3.x.^2 + r3.y.^2 + r3.z.^2);
  tic;te1_fac = rotate_tens(t,r1,r2,r3); toc
  
  t_fac = rotate_tens(t,r1,r2,r3);
  t_perp = 0.5*(t_fac.yy + t_fac.zz);
  t_par = t_fac.xx;  
  t_scal = (t_fac.xx + t_fac.yy + t_fac.zz)/3;
  vt_scal = real(sqrt(2*t_scal/(1/25)));
  anis = (t_par./t_perp);
  anis(anis<0) = NaN;
  toc
end



if doA
  A_tmp = pic_tmp.A;
end
if 0 % Babs
  hca = h(isub); isub = isub + 1;            
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,Babs');   
  hb = colorbar('peer',hca);
  hca.CLim = 0.5 + 0.5*[-1 1];
  colormap(hca,pic_colors('blue_red'));
  hb.YLabel.String = '|B|';
  if doA
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end  
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end  
if 0 % By
  hca = h(isub); isub = isub + 1;            
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.By');   
  hb = colorbar('peer',hca);
  hca.CLim = 0.5 + 0.5*[-1 1]*0.99;
  colormap(hca,pic_colors('blue_red'));
  hb.YLabel.String = 'B_{y} (B_0)';
  if doA
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end  
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end 
if 0 % Ez
  hca = h(isub); isub = isub + 1;            
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.Ez');   
  hb = colorbar('peer',hca);
  hca.CLim = 1*[-1 1]*0.99;
  colormap(hca,pic_colors('blue_red'));
  hb.YLabel.String = 'E_{z} (v_{A0}B_0)';
  if doA
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end  
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end 
if 1 % KB
  hca = h(isub); isub = isub + 1;    
  Bx = pic_tmp.Bx;
  By = pic_tmp.By;
  Bz = pic_tmp.Bz;
  KB = magnetic_field_curvature(pic_tmp.xi,pic_tmp.zi,Bx,By,Bz); % magnetic curvature
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,KB.abs');
  hb = colorbar('peer',hca);
  hca.CLim = [0 0.8];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = sprintf('|K| (d_{i0}^{-1})',KB.units);
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
 % irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end
if 1 % 1/KB
  hca = h(isub); isub = isub + 1;    
  Bx = pic_tmp.Bx;
  By = pic_tmp.By;
  Bz = pic_tmp.Bz;
  KB = magnetic_field_curvature(pic_tmp.xi,pic_tmp.zi,Bx,By,Bz); % magnetic curvature
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,1./KB.abs');
  hb = colorbar('peer',hca);
  hca.CLim = [0 20];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = sprintf('1/|K| (d_{i0})',KB.units);
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
 % irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end  
if 0 % vex,cold
  hca = h(isub); isub = isub + 1;            
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.vx([2 4])');   
  hb = colorbar('peer',hca);
  hca.CLim = 5*[-1 1]*0.99;
  colormap(hca,pic_colors('blue_red'));
  hb.YLabel.String = 'v_{ec,x} (v_{A0})';
  if doA
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end  
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end 
if 0 % Epar
  hca = h(isub); isub = isub + 1;    
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,Epar');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = 0.2*[-1 1];
  colormap(hca,pic_colors('blue_red'));
  hb.YLabel.String = 'E_{||} (v_{A0}B_0)';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end
if 1 % vt cold e
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,vt_scal');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 5];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = 'v_{t,ec} = (2T/m_e)^{1/2}';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end 
if 1 % wce
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,wce');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 30];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = '\omega_{ce} (\omega_{ci})';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end 
if 1 % vtec/wce
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,(vt_scal./wce)');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 0.2];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = 'v_{t,ec}/\omega_{ce}';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end 
if 1 % (vtec/wce)/(1/KB)
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,KB.abs'.*(vt_scal./wce)');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 0.1];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = 'K_Bv_{t,ec}/\omega_{ce}';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end 
if 0 % vpar cold e
  hca = h(isub); isub = isub + 1;    
  vx = pic_tmp.vx(4);
  vy = pic_tmp.vy(4);
  vz = pic_tmp.vz(4);
  vpar = (vx.*b.x + vy.*b.y + vz.*b.z);
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,vpar');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = 10*[-1 1];
  colormap(hca,pic_colors('blue_red'));
  hb.YLabel.String = 'v_{||} (v_{A0})';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end
if 0 % vabs cold e
  hca = h(isub); isub = isub + 1;    
  vx = pic_tmp.vx(4);
  vy = pic_tmp.vy(4);
  vz = pic_tmp.vz(4);
  vabs = sqrt(vx.^2 + vy.^2 +vz.^2);
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,vabs');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 5];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = '|v_{ec}| (v_{A0})';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end
if 0 % 1/KB
  hca = h(isub); isub = isub + 1;    
  Bx = pic_tmp.Bx;
  By = pic_tmp.By;
  Bz = pic_tmp.Bz;
  KB = magnetic_field_curvature(pic_tmp.xi,pic_tmp.zi,Bx,By,Bz); % magnetic curvature
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,1./KB.abs');
  hb = colorbar('peer',hca); hbs(end+1) = hb;     
  hca.CLim = [0 10];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = sprintf('1/|K| (d_{i0})',KB.units);
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
end
if 0 % neh
  hca = h(isub); isub = isub + 1;    
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.n(2)');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 1];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = 'n_{e,h}';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end
if 0 % nic
  hca = h(isub); isub = isub + 1;    
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.n(3)');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 0.4];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = 'n_{i,c}';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end
if 0 % nec
  hca = h(isub); isub = isub + 1;    
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.n(4)');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 0.4];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = 'n_{e,c} (n_0)';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
  %irf_legend(hca,sprintf('twpe = %05.0f',pic_tmp.twpe),[0.02 0.98],'color','k')
end
if 0 % pec
  hca = h(isub); isub = isub + 1;    
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.p(4)');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 0.05];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = 'p_{ec}';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
end
if 0 % pic
  hca = h(isub); isub = isub + 1;    
  varcount = varcount + 1;
  links{varcount} = [links{varcount} isub - 1];     
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,pic_tmp.p(3)');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 0.5];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = 'p_{ic}';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
end
if 0 % tec
  hca = h(isub); isub = isub + 1;    
  p = pic_tmp.p(3);
  n = pic_tmp.n(3);
  t = p./n;
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,t');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 0.2];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = 't_{ec}';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
end
if 0 % tec par
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,t_par');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 0.2];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = 't_{ec,||}';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
end
if 0 % tec perp
  hca = h(isub); isub = isub + 1; 
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,t_perp');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [0 0.2];
  colormap(hca,pic_colors('waterfall'));
  hb.YLabel.String = 't_{ec,\perp} (m_i v_{A0}^2)';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
end
if 0 % tec par/perp
  hca = h(isub); isub = isub + 1;      
  imagesc(hca,pic_tmp.xi,pic_tmp.zi,log10(real(anis))');
  hb = colorbar('peer',hca); hbs(end+1) = hb;
  hca.CLim = [-1 1];
  colormap(hca,pic_colors('blue_red'));
  hb.YLabel.String = 'log_{10}(t_{ec,||}/t_{ec,\perp})';
  if doA      
    hold(hca,'on')
    hcont = contour(hca,pic_tmp.xi(ipxA),pic_tmp.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off') 
  end
end
drawnow


linkprop(h,{'XLim','YLim'});
drawnow

for ip = 1:npanels
h(ip).YDir = 'normal';
h(ip).XLabel.String = 'x (di)';
h(ip).YLabel.String = 'z (di)';
end
drawnow
compact_panels(0.01)


