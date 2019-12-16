%% Load objects
%tubr = PIC(); % dont have electron data for turbulence run
df04 

%% Get some electron trajectories
it = 2;
iSpecies = [4]; % cold electron from top
idf = 10;
switch idf
  case 1
    ds = ds04(it).xlim([169]+[-0.1 0.1]).zlim([-0.3 0.3]); % top row.
    spacingPeaks = 5; % for ions its 0.2 vA
    nPeaks = 15;
  case 2
    ds = ds04(it).xlim([180]+[-0.1 0.1]).zlim([-0.3 0.3]); % top row.
    spacingPeaks = 7; % for ions its 0.2 vA
    nPeaks = 15;
  case 3
    ds = ds04(it).xlim([180]+[-0.1 0.1]).zlim(2+[-0.25 0.25]); % top row.
    spacingPeaks = 7; % for ions its 0.2 vA
    nPeaks = 15;
  case 4 % done and added to h5 file
    ds = ds04(it).xlim([188]+[-0.1 0.1]).zlim(3+[-0.25 0.25]); % top row.
    spacingPeaks = 3; % for ions its 0.2 vA
    nPeaks = 5;
  case 5
    ds = ds04(it).xlim([174]+[-0.1 0.1]).zlim(3+[-0.25 0.25]); % top row.
    spacingPeaks = 6; % for ions its 0.2 vA
    nPeaks = 10;
  case 6 % done and added to h5 file
    ds = ds04(it).xlim([205]+[-0.1 0.1]).zlim(0.5+[-0.25 0.25]); % top row.
    spacingPeaks = 2; % for ions its 0.2 vA
    nPeaks = 10;
  case 7
    ds = ds04(it).xlim([184]+[-0.1 0.1]).zlim(1.5+[-0.25 0.25]); % top row.
    spacingPeaks = 5; % for ions its 0.2 vA
    nPeaks = 10;
  case 8
    ds = ds04(it).xlim([188]+[-0.1 0.1]).zlim(0+[-0.25 0.25]); % top row.
    spacingPeaks = 5; % for ions its 0.2 vA
    nPeaks = 15;
  case 9
    ds = ds04(it).xlim([168]+[-0.1 0.1]).zlim(1+[-0.25 0.25]); % top row.
    spacingPeaks = 5; % for ions its 0.2 vA
    nPeaks = 15;
  case 10 % done and added to h5 file
    ds = ds04(it).xlim([203]+[-0.1 0.1]).zlim(0+[-0.25 0.25]); % top row.
    spacingPeaks = 2; % for ions its 0.2 vA
    nPeaks = 10;
end
fpeaks = ds.get_peaks(nPeaks,spacingPeaks,iSpecies);
nDists = ds.nd;
doPlot = 1;
if doPlot
  % plot results
  for id = 1:ds.nd{1}
  f = ds.f(1,id,iSpecies);
  figure(27)
  h = setup_subplots(3,1);
  hca = h(1);
  imagesc(hca,f.v,f.v,f.fxy')
  hca.YDir = 'normal';
  colormap(pic_colors('candy'))
  hold(hca,'on')
  plot(hca,[fpeaks(:,id).vx],[fpeaks(:,id).vy],'k.')
  hold(hca,'off')

  hca = h(2);
  imagesc(hca,f.v,f.v,f.fxz')
  hca.YDir = 'normal';
  colormap(pic_colors('candy'))
  hold(hca,'on')
  plot(hca,[fpeaks(:,id).vx],[fpeaks(:,id).vz],'k.')
  hold(hca,'off')

  hca = h(3);
  imagesc(hca,f.v,f.v,f.fyz')
  hca.YDir = 'normal';
  colormap(pic_colors('candy'))
  hold(hca,'on')
  plot(hca,[fpeaks(:,id).vy],[fpeaks(:,id).vz],'k.')
  hold(hca,'off')
  pause(0.5)
  end
end

%% Loop through points, integrate trajectories
pic = df04;
tspan = [140,160,220];
tspan = [150,160];
m = 1/25; 
q = -1;
tt = tic;
for id = 1:ds.nd{1}
  for iPeak = 1:nPeaks
    fprintf('(id/nd,ipeak/npeaks) = (%g/%g,%g/%g)\n',id,ds.nd{1},iPeak,nPeaks)
    r0 = [fpeaks(iPeak,id).x, fpeaks(iPeak,id).y, fpeaks(iPeak,id).z];
    v0 = [fpeaks(iPeak,id).vx, fpeaks(iPeak,id).vy, fpeaks(iPeak,id).vz];
    tic; 
    tr_tmp = df04.integrate_trajectory(r0,v0,tspan,m,q);
    toc
    tic;
    [Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(tr_tmp.x,tr_tmp.z,tr_tmp.t);  % interpolate
    toc
    tr_tmp.Ex = Ex;
    tr_tmp.Ey = Ey;
    tr_tmp.Ez = Ez;
    tr_tmp.Bx = Bx;
    tr_tmp.By = By;
    tr_tmp.Bz = Bz;
    tr(iPeak,id) = tr_tmp;
    toc(tt)
  end
end
traj = tr;

%% Plot what we have, all trajectories
hca = subplot(1,1,1);
hold(hca,'on')
for id = 1:size(tr,2)
  for iPeak = 1:nPeaks
    plot3(hca,tr(iPeak,id).x,tr(iPeak,id).y,tr(iPeak,id).z)
    %plot3(tr_pass(iPeak,id).x,tr_pass(iPeak,id).y,tr_pass(iPeak,id).z)
    %plot(tr(iPeak,id).t,tr(iPeak,id).vz)
  end
end
hold(hca,'off')

%% Plot on top of field, option to make movie, frame of the DF.
doVideo = 1;
twci = [140 220];
doA = 1; Alev = -25:1:0;
pic0 = df04.twcilim(twci);

%[xDF,vDF,aDF,BDF] = df04.xva_df;
xDF = xDF(1,:);

if doVideo
  vidObj = VideoWriter([savedir 'trajectories_e.mp4'],'MPEG-4');
  open(vidObj);
end

nrows = 2;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);

for it = 1:pic0.nt
  if 0
    x0 = xDF(pic0.it(it));  
    xlim = x0 + [-70 70];
  else
    xlim = [140 260];
  end
    
  zlim = [-10 10];
  twci = pic0.twci(it);
  pic = df04.twcilim(twci).xlim(xlim).zlim(zlim);

  
  pc = pic;
  
 
    
  isub = 1;
  if 0 % Ez
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.Ez')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'E_z';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % By
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.By')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'B_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 1 % Ey
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.Ey')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'E_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % Bz
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.Bz')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'B_z';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % magnetic field curvature
    hca = h(isub); isub = isub + 1;
    
    Bx = pc.Bx;
    By = pc.By;
    Bz = pc.Bz;   
    bcurv = magnetic_field_curvature(pc.xi,pc.zi,Bx,By,Bz);
    imagesc(hca,pc.xi,pc.zi,bcurv.abs')
    hca.CLim = [0 2];
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '|B_{curv}|';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
    
  end
  if 0 % v35y
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.vy([3 5])')
    hca.CLim = [-.5 0.5];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'v_{y,i,cold}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % vxB_y
    hca = h(isub); isub = isub + 1;
    vx = pc.vx([3 5]);
    vy = pc.vy([3 5]);
    vz = pc.vz([3 5]);
    Bx = pc.Bx;
    By = pc.By;
    Bz = pc.Bz;
    vxB = cross_product(vx,vy,vz,Bx,By,Bz);
    imagesc(hca,pc.xi,pc.zi,vxB.y')
    hca.CLim = [-.5 0.5];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '(vxB)_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % n vy T
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = njp(obj,iSpecies);
    hca = h(isub); isub = isub + 1;    
    
    imagesc(hca,pc.xi,pc.zi,vxB.y')
    hca.CLim = [-.5 0.5];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '(vxB)_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end 
  if 1 % vex
    
    hca = h(isub); isub = isub + 1;    
    
    imagesc(hca,pc.xi,pc.zi,pc.vx([4 6])')
    hca.CLim = [-5 5];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'v_{x,cold electrons}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end 
  
  h(1).Title.String = sprintf('twci = %g',pc.twci);
  drawnow
  compact_panels(0.02)
  for ip = 1:npanels
    hca = h(ip);
    hold(hca,'on')
    hca.FontSize = 14;
    hca.YDir ='normal';
    if doA
      iAx = 1:4:pic.nx;
      iAz = 1:4:pic.nz;
      A = pc.A;
      contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',Alev,'k');        
    end
    if 1 % all trajectories
    for itr = 1:numel(tr)           
      
      plot(hca,tr(itr).x,tr(itr).z)
      ii = find(abs(tr(itr).t-pc.twci)==min(abs(tr(itr).t-pc.twci)));
      
      plot(hca,tr(itr).x(ii),tr(itr).z(ii),'ko')      
    end
    else %subset
      for itr = 1:size(iPeakDist,1)
        plot(hca,tr(iPeakDist(itr,1),iPeakDist(itr,2)).x,tr(iPeakDist(itr,1),iPeakDist(itr,2)).z)
        ii = find(abs(tr(iPeakDist(itr,1),iPeakDist(itr,2)).t-pc.twci)==min(abs(tr(iPeakDist(itr,1),iPeakDist(itr,2)).t-pc.twci)));
        plot(hca,tr(iPeakDist(itr,1),iPeakDist(itr,2)).x(ii),tr(iPeakDist(itr,1),iPeakDist(itr,2)).z(ii),'ko')
      end
    end
    hold(hca,'off')
  end
  pause(1)
  if doVideo
    set(gcf,'color','white');
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
end
if doVideo, close(vidObj); end

%%
v = struct();


for iTr = 1:tre.ntr
  v(iTr).abs = sqrt(tre(iTr).vx.^2+tre(iTr).vy.^2+tre(iTr).vz.^2);
  babs = sqrt(tre(iTr).Bx.^2+tre(iTr).By.^2+tre(iTr).Bz.^2);
  v(iTr).par = (tre(iTr).vx.*tre(iTr).Bx + tre(iTr).vy.*tre(iTr).By + tre(iTr).vz.*tre(iTr).Bz)./babs;  
  v(iTr).absmax = max(v(iTr).abs);
  v(iTr).parmax = max(v(iTr).par);
end


[val,ind] = sort([v.parmax],'descend');

iplot = 6:15; % ind(1:5), 1:tre.ntr
h = setup_subplots(2,2);
isub = 1;
if 1 % vpar max vs vabs
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  for iTr = iplot
    vabs = sqrt(tre(iTr).vx.^2+tre(iTr).vy.^2+tre(iTr).vz.^2);
    babs = sqrt(tre(iTr).Bx.^2+tre(iTr).By.^2+tre(iTr).Bz.^2);
    vpar = (tre(iTr).vx.*tre(iTr).Bx + tre(iTr).vy.*tre(iTr).By + tre(iTr).vz.*tre(iTr).Bz)./babs;
    scatter(hca,tre(iTr).x,tre(iTr).z,2,vpar);
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'v_{||}';
    colormap(hca,pic_colors('waterfall'))
  end
  hold(hca,'off')
end
if 1 % vpar max vs vabs
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  for iTr = iplot
    vabs = sqrt(tre(iTr).vx.^2+tre(iTr).vy.^2+tre(iTr).vz.^2);
    babs = sqrt(tre(iTr).Bx.^2+tre(iTr).By.^2+tre(iTr).Bz.^2);
    vpar = (tre(iTr).vx.*tre(iTr).Bx + tre(iTr).vy.*tre(iTr).By + tre(iTr).vz.*tre(iTr).Bz)./babs;
    scatter(hca,tre(iTr).x,tre(iTr).z,2,vabs);
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '|v|';
    colormap(hca,pic_colors('waterfall'))
  end
  hold(hca,'off')
end
if 1 % vpar max vs vabs
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  for iTr = iplot
    vabs = sqrt(tre(iTr).vx.^2+tre(iTr).vy.^2+tre(iTr).vz.^2);
    babs = sqrt(tre(iTr).Bx.^2+tre(iTr).By.^2+tre(iTr).Bz.^2);
    vpar = (tre(iTr).vx.*tre(iTr).Bx + tre(iTr).vy.*tre(iTr).By + tre(iTr).vz.*tre(iTr).Bz)./babs;
    plot(hca,tre(iTr).t,abs(vpar));
    hca.YLabel.String = '|v_{||}|';
  end
  hold(hca,'off')
end
if 1 % vpar max vs vabs
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  for iTr = iplot
    vabs = sqrt(tre(iTr).vx.^2+tre(iTr).vy.^2+tre(iTr).vz.^2);
    babs = sqrt(tre(iTr).Bx.^2+tre(iTr).By.^2+tre(iTr).Bz.^2);
    vpar = (tre(iTr).vx.*tre(iTr).Bx + tre(iTr).vy.*tre(iTr).By + tre(iTr).vz.*tre(iTr).Bz)./babs;
    plot(hca,tre(iTr).t,vabs);
    hca.YLabel.String = '|v|';
    
  end
  hold(hca,'off')
end

%% Plot all trajectories
h = setup_subplots(2,2);
isub = 1;
if 1 % vmin s vmax
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  for iTr = 1:tre.ntr
    vabs = sqrt(tre(iTr).vx.^2+tre(iTr).vy.^2+tre(iTr).vz.^2);
    scatter(hca,max(vabs),min(vabs))
  end
  hold(hca,'off')
end
if 1 % vpar max vs vabs
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  for iTr = 1:tre.ntr
    vabs = sqrt(tre(iTr).vx.^2+tre(iTr).vy.^2+tre(iTr).vz.^2);
    babs = sqrt(tre(iTr).Bx.^2+tre(iTr).By.^2+tre(iTr).Bz.^2);
    vpar = (tre(iTr).vx.*tre(iTr).Bx + tre(iTr).vy.*tre(iTr).By + tre(iTr).vz.*tre(iTr).Bz)./babs;
    scatter(hca,max(vabs),max(vpar))
  end
  hold(hca,'off')
end
if 1 % vpar max vs vabs
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  for iTr = 1:tre.ntr
    vabs = sqrt(tre(iTr).vx.^2+tre(iTr).vy.^2+tre(iTr).vz.^2);
    babs = sqrt(tre(iTr).Bx.^2+tre(iTr).By.^2+tre(iTr).Bz.^2);
    vpar = (tre(iTr).vx.*tre(iTr).Bx + tre(iTr).vy.*tre(iTr).By + tre(iTr).vz.*tre(iTr).Bz)./babs;
    scatter(hca,tre(iTr).x,tre(iTr).z,2,vpar);
    hcb = colorbar('peer',hca);
    colormap(hca,pic_colors('waterfall'))
  end
  hold(hca,'off')
end
if 1 % vpar max vs vabs
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  for iTr = 1:tre.ntr
    vabs = sqrt(tre(iTr).vx.^2+tre(iTr).vy.^2+tre(iTr).vz.^2);
    babs = sqrt(tre(iTr).Bx.^2+tre(iTr).By.^2+tre(iTr).Bz.^2);
    vpar = (tre(iTr).vx.*tre(iTr).Bx + tre(iTr).vy.*tre(iTr).By + tre(iTr).vz.*tre(iTr).Bz)./babs;
    scatter(hca,tre(iTr).x,tre(iTr).z,2,vabs);
    hcb = colorbar('peer',hca);
    colormap(hca,pic_colors('waterfall'))
  end
  hold(hca,'off')
end

%% Sketch of what liouville mapping of electrons corresponds to half n_lb*vte_lb
units = irf_units;
m = 1;
fun_f = @(v,vt,n) n/sqrt(vt)*exp(-(v.^2)/(vt^ 2));
fun_v0 = @(v,phi,vph) vph + sign(v-vph).*((v-vph).^2-2*phi./m).^0.5;
fun_v = @(v0,phi,vph) vph + sign(v0-vph).*((v0-vph).^2+2*phi./m).^0.5;


vt = 1;
vph = 2;
phi = 1;
n = 1;
vv = linspace(-10,10,1000);

%plot(v0,v(v0,phi,vph))
%,vv,f(v0(vv,phi,vph),vt,n)
plot(vv,fun_f(vv,vt,n),vv,real(fun_f(fun_v0(vv,phi,vph),vt,n)))

%% Fermi reflection at DF, Liouville mapping
vDF = 1;
vd = 0;
vt = 1;
m = 1;
n = 1;
vbefore = linspace(-4,4,100);
%vafter = 
fun_f = @(v,vd,vt,n) n/sqrt(vt)*exp(-((v-vd).^2)/(vt^ 2));
fun_vafter = @(vbefore,vDF) abs(-vbefore + 2*vDF-vDF)+vDF;
fun_vbefore = @(vbefore,vDF) -vafter + 2*vDF;


fafter_ = fun_f(vbefore,vd,vt,n);
vafter_ = fun_vafter(vbefore,vDF);

vafter = interp1(vafter_,vafter_,vbefore);
%fafter = vbefore*0;
for iind = 1:numel(vbefore)
  ind = find(vafter_(iind) == vbefore);
  fafter(iind) = sum(fafter_(ind));
end


subplot(2,1,1)
plot(vbefore,fun_vafter(vbefore,vDF))

subplot(2,1,2)
plot(vbefore,fun_f(vbefore,vd,vt,n))
hold on
plot(fun_vafter(vbefore,vDF),fun_f(vbefore,vd,vt,n))
plot(vbefore,fafter)
hold off

%% New overview figure, data, 5 panels
mms_20170706_005403.load_data;
mms_20170706_005403.prepare_data_single_sc;

%% Make reduced distribution
tintZoom = irf.tint('2017-07-06T00:55:20.00Z',5);
tintZoom = irf.tint('2017-07-06T00:55:10.00Z',60);
%tintZoom = irf.tint('2017-07-06T08:18:00.00Z',13);
tintZoom = irf.tint('2017-07-06T00:54:05.00Z',40);
strTintZoom = [irf_time(tintZoom(1),'epochtt>utc_yyyymmdd_HHMMSS') '_' irf_time(tintZoom(2),'epochtt>utc_HHMMSS')];

tintZoom = irf.tint('2017-07-06T00:54:12.00Z/2017-07-06T00:54:20.00Z'); % shorter

eint = [000 40000];
vint = [-Inf Inf];
lowerelim = 100;

eDist = ePDist1.tlim(tintZoom).elim(eint);
iDist = iPDist1.tlim(tintZoom).elim(eint);
ve = gseVe1.tlim(eDist.time).resample(eDist);
vi = gseVi1.tlim(iDist.time).resample(iDist);
scpot_margin = 1; % keep in mind that this also affects the velocity at lower energies
scpot_lim = scPot1.resample(eDist)*scpot_margin;
eLine = dmpaB1.resample(eDist).norm;
iLine = dmpaB1.resample(iDist).norm;

tic; ef1D = eDist.reduce('1D',eLine,'vint',vint,'scpot',scpot_lim,'lowerelim',lowerelim); toc % reduced distribution along B
tic; if1D = iDist.reduce('1D',iLine,'vint',vint); toc % reduced distribution along B
lineVe = ve.dot(eLine); % projection of Vi on B
lineVi = vi.dot(iLine); % projection of Vi on B

%% Plot
ic = 1;
npanels = 5;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];

if 0 % B abs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B abs');
  irf_plot(hca,gseB1.abs);
  hca.YLabel.String = 'B (nT)';
end
if 0 % B GSE
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 0 % iDEF omni
  isub = isub + 1;
  hca = irf_panel('iDEF');  
  [hout,hcb] = irf_spectrogram(hca,iDist.convertto('s^3/m^6').omni.specrec,'log');
%   hold(hca,'on')
%   lineScpot = irf_plot(hca,scPot1,'k');
%   lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
%   hold(hca,'off')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 0 % iPDist pa 64
  isub = isub + 1;
  hca = irf_panel('i PA e64 deflux lowe');  
  eint = [100 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];   
end
if 0 % i psd vpar
  isub = isub + 1;
  hca = irf_panel('fi1D');
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));
  hold(hca,'on')
  irf_plot(hca,{lineVi},'comp')
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hca.YLim = if1D.depend{1}(1,[1 end]);
  hca.YLabel.String = 'v_i (km/s)'; 
  hca.CLim = [-5 0];
  %hca.CLim = max(abs(hca.CLim))*[-1 1];  
end
if 0 % i psd vpar
  isub = isub + 1;
  hca = irf_panel('fi1D*v');
  irf_spectrogram(hca,if1D.specrec('v_f1D*v'));
  hold(hca,'on')
  irf_plot(hca,{lineVi},'comp')
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hca.YLim = if1D.depend{1}(1,[1 end]);
  hca.YLabel.String = 'v_i (km/s)'; 
  hca.CLim = max(abs(hca.CLim))*[-1 1]*0.3;
  colormap(hca,cn.cmap('blue_red'))
end
if 1 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.02 0.95],'fontsize',12);
end
if 1 % Ve
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x*1e-3,gseVe?.y*1e-3,gseVe?.z*1e-3},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.02 0.95],'fontsize',12);
end
if 0 % eDEF omni
  isub = isub + 1;
  hca = irf_panel('eDEF');
  [hout,hcb] = irf_spectrogram(hca,eDist.convertto('s^3/m^6').omni.specrec,'log');
  hold(hca,'on')
  lineScpot = irf_plot(hca,scpot_lim,'k');
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_e','(eV)'};   
  irf_legend(hca,[num2str(scpot_margin) 'V_{sc}'],[0.99 0.1],'color',0*[1 1 1])
end
if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('eLine');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  %hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 1 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Te par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % E perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x,gseE?perp.y,gseE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end

if 0 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end

irf_zoom(h,'x',tintZoom+[0 -15])
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%h(5).CLim = [-35 -28]+12
%colormap('jet');

colormap(h(3),[1 1 1;pic_colors('waterfall')])
h(3).CLim = [-5.5 -3];
h(3).YLim = [-80 80];
h(3).YTick = [-100:50:100];
h(4).YLim = [0 0.16];
h(4).YTick = [0 0.05 0.1 0.15];

%colormap(h(4),cn.cmap('blue_red'))
%h(4).CLim = 100*[-1 1];
%h=irf_plot({gseB1,gseVi1,iPDist1.deflux.omni.specrec('energy'),f1D.specrec('velocity_1D')}); h(3).YScale = 'log'; %h(4).YLim = [-1000 1000];

%% Plot showing the gradual accelerationa dn thermaliztion
% check for indices
inds = 280:360;
eff = ef1D(inds);
%eff.specrec('velocity_1D','10^3 km/s')
irf_plot(eff.specrec('velocity_1D','10^3 km/s'))
%%
inds = 280:5:360;
hca = subplot(1,1,1);

vv = ef1D.depend{1}(1,:)*1e-3;
ff = ef1D.data(inds,:);
%ff(find(ff==0)) = nan;
semilogy(hca,vv,ff)
hca.XLim = [-20 20];

%% 2D reduced distribution showing 
tint = irf.tint('2017-07-06T00:54:22.945Z/2017-07-06T00:54:23.480Z');
tint = irf.tint('2017-07-06T00:54:22.945Z/2017-07-06T00:54:23.480Z');
tint = irf.tint('2017-07-06T00:54:22.10Z/2017-07-06T00:54:22.30Z');
tint = tint;% + [0.1 -0.1];

tint_utc = tint.utc;
tint_str = [tint_utc(1,12:23) ' - ' tint_utc(2,18:23)];

dist = ePDist1.tlim(tint+[-11]);
B = gseB1.tlim(tint).resample(dist).norm;
E = gseE1.tlim(tint).resample(dist).norm;
ExB = E.cross(B).norm;
%vpar = mean(B.data,1); vpar = vpar./sqrt(sum(vpar.^2));
%ExB = cross_product(E(1),E(2),E(3),B(1),B(2),B(3));
vg = (-70:2:70)*1e3;
f2D = dist.reduce('2D',B,ExB,'lowerelim',120,'base','cart','vg',vg);
f1D = dist.reduce('1D',B,'lowerelim',120,'vg',vg);


hca = subplot(1,2,1);
f2D.plot_plane(hca)
hca.XLim = 70*[-1 1];
hca.YLim = 70*[-1 1];
colormap(hca,pic_colors('waterfall'))
hca.CLim = [-14 -10.5];
axis(hca,'square')
hca.Position(2) = 0.15;
hca.Position(4) = 0.75;
hca.FontSize = 12;
hca.Title.String = tint_str;
hca.XLabel.String = 'v_{||} (10^3 km/s)';
hca.YLabel.String = 'v_{\perp} (10^3 km/s)';
hca.XTick = -80:20:80;


hca = subplot(1,2,2);
plot(hca,f1D.depend{1}(1,:)*1e-3,mean(f1D.data)*1e3);
hca.XLim = 70*[-1 1];
%axis(hca,'square')
hca.Position(2) = 0.15;
hca.Position(4) = 0.75;
hca.FontSize = 12;
hca.Title.String = tint_str;
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XLabel.String = 'v_{||} (10^3 km/s)';
hca.YLabel.String = 'f (10^{-3} s/m^4)';
hca.XTick = -80:20:80;

%% Distribution and waves
tintZoom = irf.tint('2017-07-06T00:54:12.00Z/2017-07-06T00:54:19.00Z');
doVph = 1;
doVtrap = 1;

% Plot
ic = 1;
npanels = 1;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];

if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('eLine');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  
  if doVph
    hold(hca,'on')
    irf_plot(hca,tsVphpar*1e-3,'*k')
    hold(hca,'off')
  end
  if doVtrap
    hold(hca,'on')    
    set(hca,'colororder',[0 0 0;0 0 0])
    hline = irf_plot(hca,{vmin.tlim(tintZoom)*1e-3,vmax.tlim(tintZoom)*1e-3},'comp');
    %hline(1).LineWidth = 2;
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_e','(10^3 km/s)'};
  hca.YLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.Interpreter = 'tex';
  hca.YLim = [-64 64];
  hca.YTick = [-80:20:80];
  
  %hca.Children(1).LineWidth= 1.5;
  %hca.Children(2).LineWidth= 1.5;
end
if 0 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
  hca.YLabel.String = 'E_{||} (mV/m)';  
end

if 0 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
h(1).CLim = [-6 -2.5];
%colormap('jet');

hca.YLim = [-40 70];

colormap(h(1),[1 1 1;pic_colors('candy')])
c_eval('h(?).Position(2) = h(?).Position(2) + 0.1;',1)
c_eval('h(?).Position(4) = h(?).Position(4) - 0.1;',1)

%% Four sc E field - to illustrate nonlinear character and possiblity of timing
h = irf_plot(1);
hca = h(1);
tintZoom = irf.tint('2017-07-06T00:54:15.985Z/2017-07-06T00:54:16.00Z');

set(hca,'colororder',mms_colors('1234'))
irf_plot(hca,{gseE1par.tlim(tintZoom),gseE2par.tlim(tintZoom),gseE3par.tlim(tintZoom),gseE4par.tlim(tintZoom)},'comp')
hca.YLabel.String = 'E_{||} (mV/m)';
hca.Position(2) = 0.2;
hca.Position(4) = 0.7;
irf_legend(hca,{'MMS 1';'MMS 2';'MMS 3';'MMS 4'},[0.98 0.98])

irf_zoom(h,'x',tintZoom)
