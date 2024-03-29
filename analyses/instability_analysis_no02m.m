%% Load data to make fit to
%ds100 = PICDist('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');
%no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
iIon = [1 3 5];
iEle = [2 4 6];
iIonCold = [3 5];
iEleCold = [4 6];
iIonHot = [1];
iEleHot = [2];

%% Make overview plot to show the location of distribution and fits
twpe = 23000;
% Field
xlim = [65 78]; 
zlim = [-6 6];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
nSpecies = numel(pic.mass);

ds = ds100.twpelim(twpe).findtag({'A=7.5'});
xdist = (ds.xi1{1}+ds.xi2{1})/2;
zdist = (ds.zi1{1}+ds.zi2{1})/2;
arclength = [0 cumsum(sqrt(diff(xdist).^2 + diff(zdist).^2))];
%arclength = arclength(end:-1:1);
arcval = 4; % twpe=24000
arcval = 4.3; % twpe=24000
arcval = 5.2; % twpe=24000
arcval = 4.5; % twpe=24000
arcval = 7; % twpe=23000
%arcval = 1; % twpe=23000, at neutral plane
arccenter = (arclength(end)-arclength(1))/2;
idist = find(abs(arclength-arcval-arccenter)==min(abs(arclength-arcval-arccenter)));
xzarcdist = [xdist zdist arclength];
xval = xdist(idist);
zval = zdist(idist);
ds = ds.twpelim(twpe).xfind(xval).zfind(zval);
% Distribution locations and fit parameters
% units = irf_units;% Physical constants
% qe = 1.602e-19;
% me = 9.109e-31;
% mi = 1.673e-27;
% mime = 1836;
% eps0 = 8.854e-12;
% kB = 1.381e-23;
% c = 299792458;
qe = -1;
doPlotFit = 1;
mime = 25;
dist_param = arcval;
f_max = @(v,n,vt,vd) n/sqrt(pi*vt^2)*exp(0.5*(v-vd).^2./(vt.^2));
switch dist_param
  case 1 % neutral plane
    qe = 1;
    
    n = [0.05 0.11 0.15 0.16 0.16 0.05]; % n0
    vd = [0.2 0.75 0.3 -0.5 1 0.5]; % vA0
    vt = [1.3 0.1 0.7 4.5 5.0 4.0]; % vA0
    
    % move some density across cold ion species
    n = [0.05 0.11 0.22 0.16 0.16 0.05]; % n0
    vd = [0.2 0.75 0.3 -0.5 1 0.5]; % vA0
    vt = [1.3 0.1 0.7 4.5 5.0 4.0]; % vA0
    
    m = [1 1 1 1/100 1/100 1/100]*1; % m0
    q = [1 1 1 -1 -1 -1]*qe; % e0=1?
    iIncl = [1 2 3];
    eIncl = [4 5 6];
    % n = sqrt(ne^2/meps0)
    wpewce = 2;
    mime = 200;
    eps0 = 1/(wpewce^2*mime);
    wp = sqrt(n./m)*sqrt(1/eps0); % normalization? wci? write in units of wci?
    x = 0.1;
    k_min = 0.2;
    
  case 7 % hot ebg
    %%
    qe = 1;
    n = [0.05 0.08 0.15 0.11 0.13 0.04]; 
    vd = [0.5 1.15 0.05 -2.1 3 0.5]; % m/s
    vt = [1.5 0.3 0.07 2.8 3.0 4.0];
    
    n = [0.05 0.08 0.15 0.11 0.13 0.04]; % n0
    vd = [0.5 1.15 0.05 -2.1 2.7 0.5]; % vA0
    vt = [1.5 0.3 0.07 2.8 3.0 4.0]; % vA0
    
    m = [1 1 1 1/100 1/100 1/100]*1; % m0
    m = [1 1 1 1/100 1/100 1/100]*1; % m0
    q = [1 1 1 -1 -1 -1]*qe; % e0=1?
    iIncl = [1 2 3];
    eIncl = [4 5 6];
    % it seems like vph might be correct because the vd,vt normalizations 
    % are correct...?
    % n = sqrt(ne^2/meps0)
    wpewce = 2;
    mime = 200; % WHY IS MIME = 200 ?????
    eps0 = 1/(wpewce^2*mime);
    wp = sqrt(n./m)*sqrt(1/eps0); % normalization? wci? write in units of wci?
    x = 0.1;
    k_min = 0.2;
  case 7000001
    n = [0.05 0.08 0.15 0.14 0.14];
    vd = [0.5 1.15 0.05 -2 3]; % m/s
    vt = [1.5 0.3 0.07 3.5 3.5];
    %vd = [0.5 0.95 0.025 -3 4]; % m/s
    %vt = [1.5 0.20 0.05 6 6];
    m = [1 1 1 1/100 1/100]*1;
    q = [1 1 1 -1 -1]*qe;     
    iIncl = [1 2 3];
    eIncl = [5];
    wp = sqrt(n./m);
    x = 0.1;
    k_min = 0.2;
  case 4.5
    n = [0.05 0.2 0.06 0.155 0.155];
    vd = [0.5 0.900 0.0 -2 3]; % m/s
    vt = [1.5 0.40 0.07 3.8 3.8];
    %vd = [0.5 0.95 0.025 -3 4]; % m/s
    %vt = [1.5 0.20 0.05 6 6];
    m = [1 1 1 1/100 1/100]*1;
    q = [1 1 1 -1 -1]*qe;     
    iIncl = [1 2 3];
    eIncl = [5];
    wp = sqrt(n./m);
    x = 0.1;
    k_min = 0.2;
  case 5.2009
    n = [0.2 0.06 0.3];
    vd = [0.95 0.025 0.5]; % m/s
    vt = [0.20 0.05 5];
    m = [1 1 1/100]*1;
    q = [1 1 -1]*qe;     
    iIncl = 1:2;
    eIncl = [3];
    wp = sqrt(n./m);
    x = 0.1;
    k_min = 0.2;
  case 5.2
    n = [0.05 0.2 0.06 0.155 0.155];
    vd = [0.5 0.95 0.025 -2 3]; % m/s
    vt = [1.5 0.20 0.05 3.8 3.8];
    %vd = [0.5 0.95 0.025 -3 4]; % m/s
    %vt = [1.5 0.20 0.05 6 6];
    m = [1 1 1 1/100 1/100]*1;
    q = [1 1 1 -1 -1]*qe;     
    iIncl = [1 2 3];
    eIncl = [5];
    wp = sqrt(n./m);
    x = 0.1;
    k_min = 0.2;
  case 4
    n = [0.1 0.06 0.07];
    vd = [1.1 -0.08 0.2]; % m/s
    vt = [0.12 0.05 0.1];
    m = [1 1 1]*1;
    q = [1 1 1]*qe;     
    iIncl = 1:2;
    wp = sqrt(n./m);
  case -2
    n = [0.25 0.05 0.07];
    vd = [-0.88 -0.3 0.2]; % m/s
    vt = [0.15 0.2 0.1];
    m = [1 1 1]*1;
    q = [1 1 1]*qe;     
    iIncl = 1:3;
    wp = sqrt(n./m);
  case 2
    n = [0.05 0.17 0.11 0.25 0.5];
    vt = [0.20 0.08 0.1 0.15 1.7];
    m = [1 1 1 1 1]*1;
    q = [1 1 1 1 1]*qe; 
    vd = [1.4 0.06 -0.23 -0.49 -0.5]; % m/s
    iIncl = 1:2;
    wp = sqrt(n./m);
end

% Prepare simulation distributions
clear f
for iSpecies = 1:nSpecies
  f(iSpecies) = ds.fxyz(1,1,iSpecies);    
end
for iSpecies = 1:nSpecies
  f(iSpecies).n_map = mean(mean(pic.xlim(f(iSpecies).x).zlim(f(iSpecies).z).n(iSpecies)));
  f(iSpecies).n_dist = sum(f(iSpecies).f(:))*f(iSpecies).dv^3;
end
if 1 % takes some time, so not necessary to do many times if the dist is the same
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
end
ve = f_rot(iEle(1)).v;    
iSpecies = 2; fehot = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
iSpecies = 4; fecoldtop = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
iSpecies = 6; fecoldbot = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
fetot = fehot + fecoldtop + fecoldbot;
vi = f_rot(iIon(1)).v;   
iSpecies = 1; fihot = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
iSpecies = 3; ficoldtop = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
iSpecies = 5; ficoldbot = squeeze(sum(sum(f_rot(iSpecies).f,3),2))*f_rot(iSpecies).dv^2;
fitot = fehot + ficoldtop + ficoldbot;

fefit = fmax1D(ve,n(eIncl),vd(eIncl),vt(eIncl));
fifit = fmax1D(vi,n(iIncl),vd(iIncl),vt(iIncl));

% Figure
figure(73)
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
if 1 % f(vpar), electrons
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];
  isDistEle = [isDistEle isub-1];
  plot(hca,ve,fehot,ve,fecoldtop,ve,fecoldbot,ve,fetot,ve,fefit,'--')
  
  hca.XLabel.String = 'v_{||} (v_{A0})';
  hca.YLabel.String = 'f_e(v_{||}) (d_{i0}^3v_{A0})';
  legend(hca,{'f_{hot}','f_{cold}^{top}','f_{cold}^{top}','f_{tot}','f_{fit}'},'location','eastoutside')  
end
if 1 % f(vpar), ions
  hca = h(isub); isub = isub + 1;
  isDist = [isDist isub-1];
  isDistIon = [isDistIon isub-1];  
  plot(hca,vi,fihot,vi,ficoldtop,vi,ficoldbot,vi,fitot,vi,fifit,'--')
  
  hca.XLabel.String = 'v_{||} (v_{A0})';
  hca.YLabel.String = 'f_i(v_{||}) (d_{i0}^3v_{A0})';
  legend(hca,{'f_{hot}','f_{cold}^{top}','f_{cold}^{top}','f_{tot}','f_{fit}'},'location','eastoutside')
end
if 1 % f(vpar), electrons and ions fit
  hca = h(isub); isub = isub + 1;
  plotyy(hca,vi,fifit,ve,fefit)
  
  hca.XLabel.String = 'v_{||} (v_{A0})';
  hca.YLabel.String = 'f_i(v_{||}) (d_{i0}^3v_{A0})';
  legend(hca,{'f_{i,fit}','f_{e,fit}'},'location','eastoutside')  
  hca.Position(3) = h(1).Position(3);
end
if 0 % f(vpar,vperp1), cold ions
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

h(1).Title.String = sprintf('x=%g, z=%g',xval,zval);
hlinksMap = linkprop(h(isMap),{'XLim','YLim'});
hlinksMap = linkprop(h(isMap),{'XLim','YLim'});
hlinksDistIon = linkprop(h(isDistIon),{'XLim','YLim'});
%hlinksDistEle = linkprop(h(isDistEle),{'XLim','YLim'});

hlinksDistIon.Targets(1).XLim = 2*[-1 1];
%hlinksDistEle.Targets(1).XLim = 8*[-1 1];

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
figure(77)
nk = 40;
lguess = [10 0.1];
k_min = 2*pi/lguess(1);
k_max = 2*pi/lguess(2);
k_min = 1; k_max = 60;
knorm = 1;  % length
kmin = 1;
knorm_str = sprintf('L_{d%g}',1);
kvec = linspace(k_min,k_max,nk)/knorm;

wr_store = nan(1,nk);
wi_store = nan(1,nk);
fval_store = nan(1,nk);
h = setup_subplots(4,1);
%x = 0.5;-400;x=k_min*mean(vd);
x = 1;
%close all
clear Dsep_all
colors = pic_colors('matlab');
if not(exist('icolor','var')), icolor = 0; end
icolor = mod(icolor,size(colors,1));
icolor = icolor + 1;
holdOn = 0;
for ik = 1:nk
  xguess = x;
  %xguess = 1.5;
  %if xguess>2, xguess = 0.0; end
  %xguess = vd(1)*kvec(ik);
  %allIncl = [1:4];
  allIncl = [1 2 3 4 5 6];
  %allIncl = [1 3 5];
  af = @(temp) D_streaming(temp,kvec(ik),vt(allIncl),wp(allIncl),vd(allIncl));   
  options = optimset('GradObj','on','display','off','TolFun',1e-4);  
  [x,FVAL,EXITFLAG] = fsolve(af,xguess,options);    
  fprintf('ik = %g, x = %g + %gi \n',ik,real(x),imag(x))  
  fval_store(ik) = FVAL; 
  wr_store(ik) = real(x);
  wi_store(ik) = imag(x);  
  %plot(kvec(ik),real(x),'o',kvec(ik),imag(x),'*')
  plot(h(1),kvec(ik),real(x),'*','color',colors(icolor,:))
  plot(h(2),kvec(ik),imag(x),'*','color',colors(icolor,:))
  plot(h(3),kvec(ik),real(x)/kvec(ik),'*','color',colors(icolor,:))
 
  if 0
  [D,Dsep] = D_streaming(x,kvec(ik),vt(allIncl),wp(allIncl),vd(allIncl)); 
  Dsep_all(ik,:) = Dsep;
  hl = plot(h(4),kvec(ik),imag(Dsep)','*');
  c_eval('hl(?).Color = colors(?,:);',1:numel(hl))
  end
  
  c_eval('grid(h(?),''on'');',1:4)
  drawnow
  if not(holdOn)
    c_eval('hold(h(?),''on'')',1:4)
  end
end

vph_store = wr_store./kvec;

rem_ind = [];nk;350:nk;
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

%% Make an all-in-one figure
twpe = 23000;








