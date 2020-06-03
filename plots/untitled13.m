% movie for interview
gf05 = PIC('/Volumes/Fountain/Data/PIC/guide_field_05_cold_ions/data_h5/fields.h5');
pic = gf05;
times = pic.twci;
ntimes = pic.nt;
xlim = [50 160];
zlim = [-11 11];
clim = [0 0.4];
hbs = gobjects(0);
doVideo = 0;
localuser = datastore('local','user');
savedir  = ['/Users/' localuser '/GoogleDrive/Research/PIC/gf05/'];
movieName = 'n_cold_hot_diff_df04_frameoffront_clim';
if doVideo
  vidObj = VideoWriter([savedir movieName '.mp4'],'MPEG-4');
  open(vidObj);
end

for itime = 14%:ntimes
  pic_tmp = pic.twcilim(times(itime)).xlim(xlim).zlim(zlim);

  Bx = pic_tmp.Bx;
  By = pic_tmp.By;
  Bz = pic_tmp.Bz;
  Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
  b.x =  Bx./Babs;
  b.y =  By./Babs;
  b.z =  Bz./Babs;
  if 1 % make par/perp temperatures
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
  
  %% Set up figure
  nrows = 2;
  ncols = 1; 
  npanels = nrows*ncols;
  %h = setup_subplots(nrows,ncols,'horizontal');
  h = setup_subplots(nrows,ncols,'vertical');
  isub = 1;

  doA = 1;
  doAx = 0; % plot separatrix, NOT impl.
  if doA % Flux function, set parameters
    pic_A = pic.twcilim(times(itime)).xlim(xlim+5*[-1 1]).zlim(zlim+1*[-1 1]);
    cA = 0*[0.8 0.8 0.8];
    nA = 20;
    nA = [-50:1:50];
    ipxA = 1:5:pic_A.nx;
    ipzA = 1:5:pic_A.nz;     
    A_tmp = pic_A.A;
  end
  
  % Panels
  if 1 % tec par
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,t_par');
    hb = colorbar('peer',hca); hbs(isub-1) = hb;
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.CLim = clim;
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = 't_{ec,||} (m_i v_{A0}^2)';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_A.xi(ipxA),pic_A.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
  end
  if 1 % tec perp
    hca = h(isub); isub = isub + 1; 
    imagesc(hca,pic_tmp.xi,pic_tmp.zi,t_perp');
    hb = colorbar('peer',hca); hbs(isub -1) = hb;
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.CLim = clim;    
    colormap(hca,pic_colors('waterfall'));
    hb.YLabel.String = 't_{ec,\perp} (m_i v_{A0}^2)';
    if doA      
      hold(hca,'on')
      hcont = contour(hca,pic_A.xi(ipxA),pic_A.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
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
      hcont = contour(hca,pic_A.xi(ipxA),pic_A.zi(ipzA),squeeze(A_tmp(ipxA,ipzA))',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
      hold(hca,'off') 
    end
  end
  compact_panels(0.01)
  hbs(2).Position(1) = hbs(1).Position(1);
  h(1).Title.String = ['t\omega_{ci} = ' sprintf('%g',times(itime))];
  for ip = 1:npanels
    h(ip).YDir = 'normal';
    h(ip).XLim = xlim;
    h(ip).YLim = zlim;
%    h(ip).XLabel.String = 'x (di)';
%    h(ip).YLabel.String = 'z (di)';
  end
  
  pause(1)
  if doVideo
    set(gcf,'color','white');
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
end

if doVideo, close(vidObj); end