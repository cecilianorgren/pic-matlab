%% Load data to make fit to
ds04 = PICDist('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/dists.h5');
df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
iIon = [1 3 5];
iEle = [2 4 6];
iIonCold = [3 5];
iEleCold = [4 6];
iIonHot = [1];
iEleHot = [2];

%% Make overview plot to show the location of distribution
twpe = 7000;
% Field
xlim = [170 200]; 
zlim = [-1 7];
pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);
nSpecies = numel(pic.mass);

% Distribution locations and fit parameters
% units = irf_units;% Physical constants
% qe = 1.602e-19;
% me = 9.109e-31;
% mi = 1.673e-27;
% mime = 1836;
% eps0 = 8.854e-12;
% kB = 1.381e-23;
% c = 299792458;
doPlotFit = 1;
mime = 25;
dist_param = 3;
f_max = @(v,n,vt,vd) n/sqrt(pi*vt^2)*exp(0.5*(v-vd).^2./(vt.^2));
switch dist_param
  case 1
    xval = 185; zval = 0;
  case 2
    xval = 181; zval = 1;
  case 3
    xval = 178; zval = 3;
    n = [0.2 0.30 0.11 0.25 0.5];
    vt = [0.20 0.05 0.1 0.15 1.7];
    m = [1 1 1 1 1]*1;
    q = [1 1 1 1 1]*qe; 
    vd = [-0.85 -0.04 -0.23 -0.49 -0.5]; % m/s
    iIncl = 1:2;
    wp = sqrt(n./m);
  case 4
    xval = 177; zval = 2;
    n = [0.43 0.15 0.11 0.25 0.5];
    vt = [0.10 0.07 0.1 0.15 1.7];
    m = [1 1 1 1 1]*1;
    q = [1 1 1 1 1]*qe; 
    vd = [-0.65 -0.09 -0.23 -0.49 -0.5]; % m/s
    iIncl = 1:4;
    wp = sqrt(n./m);
  case 5
    xval = 178; zval = 2;
    n = [0.4 0.30 0.2 0.5];
    vt = [0.15 0.07 0.34 1.7];
    m = [1 1 1 1]*1;
    q = [1 1 1 1]*qe; 
    vd = [-0.64 0 -0.35 -0.5]; % m/s
    iIncl = 1:3;
    wp = sqrt(n./m);
    
    k_min = 0.1; k_max = 4;
    x = -0.5;
end
fmax = f_maxwellian(n(iIncl),vd(iIncl),vt(iIncl),1);

ds = ds04.twpelim(twpe).xfind(xval).zfind(zval).dxlim([0 0.21]);
clear f
for iSpecies = 1:nSpecies
  f(iSpecies) = ds.fxyz(1,1,iSpecies);    
end
for iSpecies = 1:nSpecies
  f(iSpecies).n_map = mean(mean(pic.xlim(f(iSpecies).x).zlim(f(iSpecies).z).n(iSpecies)));
  f(iSpecies).n_dist = sum(f(iSpecies).f(:))*f(iSpecies).dv^3;
end
clear f_rot
for iSpecies = 1:nSpecies
  x_tmp = f(iSpecies).x;
  z_tmp = f(iSpecies).z;
  pic_tmp = pic.xlim(x_tmp).zlim(z_tmp);
  Bx = mean(mean(pic_tmp.Bx));
  By = mean(mean(pic_tmp.By));
  Bz = mean(mean(pic_tmp.Bz));
  r1 = [Bx,By,Bz]; r1 = r1/sqrt(sum(r1.^2));
  r2 = cross(r1,cross([0 1 0],r1)); r2 = r2/sqrt(sum(r2.^2));
  r3 = cross(r1,r2); r3 = r3/sqrt(sum(r3.^2));
  f_rot(iSpecies) = rotate_dist(f(iSpecies),r1,r2,r3);
end

% Figure
nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
isMap = [];
isDist = [];
isDistIon = [];
isDistEle = [];

if 0 % Ex
  hca = h(isub); isub = isub + 1;
  isMap = [isMap isub-1];
  pcolor(hca,pic.xi,pic.zi,pic.Ex')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 'z (d_{i0})';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_x (v_{A0}B_0)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = 0.5*[-1 1];
  ds.plot_boxes(hca);
end
if 0 % ni cold
  hca = h(isub); isub = isub + 1;
  isMap = [isMap isub-1];
  pcolor(hca,pic.xi,pic.zi,pic.n(iIonCold)')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 'z (d_{i0})';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_{i,cold} (v_{A0}B_0)';
  colormap(hca,pic_colors('waterfall'))
  ds.plot_boxes(hca);
  clim = hca.CLim;
  if 1 % A    
    hold(hca,'on')
    A = pic.twpelim(twpe).A;
    iplx = 1:5:pic.nx;
    iplz = 1:5:pic.nz;
    contour(hca,pic.xi(iplx),pic.zi(iplz),A(iplx,iplz)',-25:1:0,'k')
    hold(hca,'off')
  end
  hca.CLim = clim;
end

if 0 % f(vx), all ions
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];  
  isDistIon = [isDistIon isub-1];  
  isHoldOn = 0;
  for iSpecies = iIon
    plot(hca,f(iSpecies).v,sum(sum(f(iSpecies).f,3),2)*f(iSpecies).dv^2)
    if not(isHoldOn)
      hold(hca,'on')
      isHoldOn = 1;
    end    
  end
  % Add up cold populations
  vcold = f(iIonCold(1)).v;
  fcold = vcold*0;  
  for iSpecies = iIonCold
    fcold = fcold + sum(sum(f(iSpecies).f,3),2)*f(iSpecies).dv^2;    
  end
  plot(hca,vcold,fcold)
  
  hold(hca,'off')
  hca.XLabel.String = 'v_x (v_{A0})';
  hca.YLabel.String = 'f_i (d_{i0}^3v_{A0})';
end
if 0 % f(vx), all electrons
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];  
  isDistEle = [isDistEle isub-1];  
  isHoldOn = 0;
  for iSpecies = iEle
    plot(hca,f(iSpecies).v,sum(sum(f(iSpecies).f,3),2)*f(iSpecies).dv^2)
    if not(isHoldOn)
      hold(hca,'on')
      isHoldOn = 1;
    end    
  end
  hold(hca,'off')
  hca.XLabel.String = 'v_x (v_{A0})';
  hca.YLabel.String = 'f_e(v_x) (d_{i0}^3v_{A0})';
end
if 0 % f(vy), all ions
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1]; 
  isDistIon = [isDistIon isub-1];   
  isHoldOn = 0;
  for iSpecies = iIon
    plot(hca,f(iSpecies).v,sum(sum(f(iSpecies).f,3),1)*f(iSpecies).dv^2)
    if not(isHoldOn)
      hold(hca,'on')
      isHoldOn = 1;
    end    
  end
  % Add up cold populations
  vcold = f(iIonCold(1)).v;
  fcold = vcold*0;  
  for iSpecies = iIonCold
    fcold = fcold + sum(sum(f(iSpecies).f,3),1)*f(iSpecies).dv^2;    
  end
  plot(hca,vcold,fcold)
  hold(hca,'off')
  hca.XLabel.String = 'v_y (v_{A0})';
  hca.YLabel.String = 'f_i(v_y) (d_{i0}^3v_{A0})';
end
if 0 % f(vy), all electrons
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];  
  isDistEle = [isDistEle isub-1];  
  isHoldOn = 0;
  for iSpecies = iEle
    plot(hca,f(iSpecies).v,sum(sum(f(iSpecies).f,3),1)*f(iSpecies).dv^2)
    if not(isHoldOn)
      hold(hca,'on')
      isHoldOn = 1;
    end    
  end
  hold(hca,'off')
  hca.XLabel.String = 'v_y (v_{A0})';
  hca.YLabel.String = 'f_e(v_y) (d_{i0}^3v_{A0})';
end
if 0 % f(vz), all ions
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];
  isDistIon = [isDistIon isub-1];  
  isHoldOn = 0;
  for iSpecies = iIon
    plot(hca,f(iSpecies).v,squeeze(sum(sum(f(iSpecies).f,2),1))*f(iSpecies).dv^2)
    if not(isHoldOn)
      hold(hca,'on')
      isHoldOn = 1;
    end    
  end
  % Add up cold populations
  vcold = f(iIonCold(1)).v;
  fcold = vcold*0;  
  for iSpecies = iIonCold
    fcold = fcold + squeeze(sum(sum(f(iSpecies).f,2),1))*f(iSpecies).dv^2;    
  end
  plot(hca,vcold,fcold)
  hold(hca,'off')
  hca.XLabel.String = 'v_z (v_{A0})';
  hca.YLabel.String = 'f_i(v_z) (d_{i0}^3v_{A0})';
end
if 0 % f(vz), all electrons
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1]; 
  isDistEle = [isDistEle isub-1];   
  isHoldOn = 0;
  for iSpecies = iEle
    plot(hca,f(iSpecies).v,squeeze(sum(sum(f(iSpecies).f,2),1))*f(iSpecies).dv^2)
    if not(isHoldOn)
      hold(hca,'on')
      isHoldOn = 1;
    end    
  end
  hold(hca,'off')
  hca.XLabel.String = 'v_z (v_{A0})';
  hca.YLabel.String = 'f_e(v_z) (d_{i0}^3v_{A0})';
end
if 0 % f(vx,vz), cold ions
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];  
  isHoldOn = 0;
  % Add up cold populations
  vcold = f(iIonCold(1)).v;
  fcold = zeros(numel(vcold),numel(vcold));  
  for iSpecies = iIonCold
    fcold = fcold + squeeze(sum(f(iSpecies).f,2))*f(iSpecies).dv^1;
  end
  toplot = log10(fcold);
  pcolor(hca,vcold,vcold,toplot')
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'f(v_x,v_z) (d_{i0}^{-3}v_{A0}^{-1})';
  hold(hca,'off')
  hca.XLabel.String = 'v_x (v_{A0})';
  hca.YLabel.String = 'v_z (v_{A0})';
  colormap(hca,pic_colors('waterfall'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  %hca.CLim(1) = prctile(toplot(not(isnan(toplot))),5);
  hca.CLim = [-2 1];
  %axis(hca,'equal')
  axis(hca,'square')
end
if 1 % f(vpar), cold ions
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];
  isDistIon = [isDistIon isub-1];  
  isHoldOn = 0;
  for iSpecies = iIon
    plot(hca,f_rot(iSpecies).v,squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f(iSpecies).dv^2)
    if not(isHoldOn)
      hold(hca,'on')
      isHoldOn = 1;
    end    
  end
  % Add up cold populations
  vcold = f(iIonCold(1)).v;
  fcold = vcold*0;  
  for iSpecies = iIonCold
    fcold = fcold + squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;    
  end
  plot(hca,vcold,fcold)
  hold(hca,'off')  
  if doPlotFit
    hold(hca,'on')  
    plot(hca,vcold,fmax(vcold))
    hold(hca,'off')  
  end
  hca.XLabel.String = 'v_{||} (v_{A0})';
  hca.YLabel.String = 'f_i(v_{||}) (d_{i0}^3v_{A0})';
end
if 1 % f(vpar), fit
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];
  isDistIon = [isDistIon isub-1];  
  isHoldOn = 0;
  % Add up cold populations
  vcold = f(iIonCold(1)).v;  
 
  plot(hca,vcold,fmax(vcold),vcold,fcold,'linewidth',1.5)
  hold(hca,'on')
  for iS = 1:numel(iIncl)
    fmax_tmp = f_maxwellian(n(iIncl(iS)),vd(iIncl(iS)),vt(iIncl(iS)),1);
    plot(hca,vcold,fmax_tmp(vcold))
  end
  hold(hca,'off')
    
  hca.XLabel.String = 'v_{||} (v_{A0})';
  hca.YLabel.String = 'f_i(v_{||}) (d_{i0}^3v_{A0})';
end
if 1 % f(vpar,vperp1), cold ions
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];  
  isHoldOn = 0;
  % Add up cold populations
  vcold = f_rot(iIonCold(1)).v;
  fcold = zeros(numel(vcold),numel(vcold));  
  for iSpecies = iIonCold
    fcold = fcold + squeeze(sum(f_rot(iSpecies).f,2))*f_rot(iSpecies).dv^1;
  end
  toplot = log10(fcold);
  pcolor(hca,vcold,vcold,toplot')
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'f(v_x,v_z) (d_{i0}^{-3}v_{A0}^{-1})';
  hold(hca,'off')
  hca.XLabel.String = sprintf('v_{||}^{[%.2f,%.2f,%.2f]} (v_{A0})',f_rot(iSpecies).r1);
  hca.YLabel.String = sprintf('v_{perp,2}^{[%.2f,%.2f,%.2f]} (v_{A0})',f_rot(iSpecies).r3);
  colormap(hca,pic_colors('waterfall'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  %hca.CLim(1) = prctile(toplot(not(isnan(toplot))),5);
  hca.CLim = [-2 1];
  axis(hca,'square')
end


hlinksMap = linkprop(h(isMap),{'XLim','YLim'});
hlinksMap = linkprop(h(isMap),{'XLim','YLim'});
hlinksDistIon = linkprop(h(isDistIon),{'XLim','YLim'});
hlinksDistEle = linkprop(h(isDistEle),{'XLim','YLim'});

hlinksDistIon.Targets(1).XLim = 1.5*[-1 1];
%hlinksDistEle.Targets(1).XLim = 7*[-1 1];

for ip = isDist
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end

%% Solver
% Multi-species streaming instability
% Solver for 1D unmagnetized plasma
% If solution is not as expected, change initial guess: x = 0; 
% i-e instability typically low ~wpi, 
% e-e bump-on-tail instability typically higher ~vd*kvec(ik);
% e-e two-stream instablity typically vbulk*k(1)
% If solution seems ok but out of k range change knorm or k_min, k_max.

% Dispersion solver, one surface
% Wavenumber vector
nk = 100;
lguess = [1 0.1];
k_min = 2*pi/lguess(1);
k_max = 2*pi/lguess(2);
k_min = 0.1; k_max = 4;
knorm = 1;  % length
knorm_str = sprintf('L_{d%g}',1);
kvec = linspace(k_min,k_max,nk)/knorm;

wr_store = nan(1,nk);
wi_store = nan(1,nk);
fval_store = nan(1,nk);
x = -0.5;-400;x=k_min*mean(vd);
%close all
for ik = 1:nk  
  xguess = x;
  %xguess = vd(1)*kvec(ik);
  
  af = @(temp) D_streaming(temp,kvec(ik),vt(iIncl),wp(iIncl),vd(iIncl));   
  options = optimset('GradObj','on','display','off','TolFun',1e-4);  
  [x,FVAL,EXITFLAG] = fsolve(af,xguess,options);    
  fprintf('x = %g + %gi \n',real(x),imag(x))
  fval_store(ik) = FVAL; 
  wr_store(ik) = real(x);
  wi_store(ik) = imag(x);  
  plot(kvec(ik),real(x),'o',kvec(ik),imag(x),'*')
  drawnow
  hold('on')
end

vph_store = wr_store./kvec;

rem_ind = nk;350:nk;
wr_store(rem_ind) = NaN;
wi_store(rem_ind) = NaN;
fval_store(rem_ind) = NaN;
vph_store(rem_ind) = NaN;

ikmax = find(wi_store==max(wi_store),1,'first');
kmax = kvec(ikmax);
vphmax = wr_store(ikmax)/kvec(ikmax);
wimax = wi_store(ikmax);
wrmax = wr_store(ikmax);

%% plot solution
pic = df04.twpelim(twpe).xlim(xval + [-3 3]).zlim(zval + [-3 3]);
figure(81)
% Figure
nrows = 2;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
isMap = [];
isDist = [];
isDistIon = [];
isDistEle = [];

doPlotFit = 1;

if 1 % Ex
  hca = h(isub); isub = isub + 1;
  isMap = [isMap isub-1];
  pcolor(hca,pic.xi,pic.zi,pic.Ex')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 'z (d_{i0})';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_x (v_{A0}B_0)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = 0.5*[-1 1];
  ds.plot_boxes(hca);
end
if 1 % ni cold
  hca = h(isub); isub = isub + 1;
  isMap = [isMap isub-1];
  pcolor(hca,pic.xi,pic.zi,pic.n(iIonCold)')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_{i0})';
  hca.YLabel.String = 'z (d_{i0})';
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_{i,cold} (v_{A0}B_0)';
  colormap(hca,pic_colors('waterfall'))
  ds.plot_boxes(hca);
  clim = hca.CLim;
  if 1 % A    
    hold(hca,'on')
    A = pic.twpelim(twpe).A;
    iplx = 1:5:pic.nx;
    iplz = 1:5:pic.nz;
    contour(hca,pic.xi(iplx),pic.zi(iplz),A(iplx,iplz)',-25:1:0,'k')
    hold(hca,'off')
  end
  hca.CLim = clim;
end
if 1 % f(vpar,vperp1), cold ions
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];  
  isHoldOn = 0;
  % Add up cold populations
  vcold = f_rot(iIonCold(1)).v;
  fcold = zeros(numel(vcold),numel(vcold));  
  for iSpecies = iIonCold
    fcold = fcold + squeeze(sum(f_rot(iSpecies).f,2))*f_rot(iSpecies).dv^1;
  end
  toplot = log10(fcold);
  pcolor(hca,vcold,vcold,toplot')
  shading(hca,'flat')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'f(v_x,v_z) (d_{i0}^{-3}v_{A0}^{-1})';
  hold(hca,'off')
  hca.XLabel.String = sprintf('v_{||}^{[%.2f,%.2f,%.2f]} (v_{A0})',f_rot(iSpecies).r1);
  hca.YLabel.String = sprintf('v_{perp,2}^{[%.2f,%.2f,%.2f]} (v_{A0})',f_rot(iSpecies).r3);
  colormap(hca,pic_colors('waterfall'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  %hca.CLim(1) = prctile(toplot(not(isnan(toplot))),5);
  hca.CLim = [-2 1];
  axis(hca,'square')
  hca.XTick = hca.YTick;
end

if 1 % f(vpar), cold ions
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];
  isDistIon = [isDistIon isub-1];  
  isHoldOn = 0;
  for iSpecies = iIon
    plot(hca,f_rot(iSpecies).v,squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f(iSpecies).dv^2,'linewidth',1.5)
    if not(isHoldOn)
      hold(hca,'on')
      isHoldOn = 1;
    end    
  end
  % Add up cold populations
  vcold = f(iIonCold(1)).v;
  fcold = vcold*0;  
  for iSpecies = iIonCold
    fcold = fcold + squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;    
  end
  plot(hca,vcold,fcold,'linewidth',1.5)
  hold(hca,'off')  
  if doPlotFit
    hold(hca,'on')  
    plot(hca,vcold,fmax(vcold),'k--')
    hold(hca,'off')  
  end
  hca.XLabel.String = 'v_{||} (v_{A0})';
  hca.YLabel.String = 'f_i(v_{||}) (d_{i0}^3v_{A0})';
  legend(hca,{'hot ions','cold ions from north','cold ions from south','all cold ions','fit'},'Box','off','location','best')
end
if 1 % f(vpar), fit
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];
  isDistIon = [isDistIon isub-1];  
  isHoldOn = 0;
  % Add up cold populations
  vcold = f(iIonCold(1)).v;
 
  plot(hca,vcold,fmax(vcold),vphmax,fmax(vphmax),'o','linewidth',1.5)
  hold(hca,'on')
  hold(hca,'off')
    
  hca.XLabel.String = 'v_{||} (v_{A0})';
  hca.YLabel.String = 'f_i(v_{||}) (d_{i0}^3v_{A0})';
  legend(hca,{'input to dispersion solver','phase velocity at maximum growth rate'},'Box','off')
end
if 1 % dispersion relation
  hca = h(isub); isub = isub + 1; 
    
  plot(hca,kvec,wr_store,kvec,wi_store,kvec,wr_store./kvec,'linewidth',1.5)
  hold(hca,'on')
  hold(hca,'off')
    
  hca.XLabel.String = 'k (d_{i0}^{-1})';
  hca.YLabel.String = {'\omega_r,\omega_i (\omega_{ci})','v_{ph} (v_{A0})'};
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  legend(hca,{'real frequency','imaginary frequency (growth rate)','phase velocity'},'Box','off','location','best')
end


hlinksMap = linkprop(h(isMap),{'XLim','YLim'});
hlinksMap = linkprop(h(isMap),{'XLim','YLim'});
hlinksDistIon = linkprop(h(isDistIon),{'XLim','YLim'});
hlinksDistEle = linkprop(h(isDistEle),{'XLim','YLim'});

hlinksDistIon.Targets(1).XLim = 1.5*[-1 1];
%hlinksDistEle.Targets(1).XLim = 7*[-1 1];

for ip = isDist
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end

for ip = 1:npanels
  h(ip).FontSize = 13;  
end
