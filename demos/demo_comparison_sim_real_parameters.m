units = irf_units;
me = units.me;
mi = units.mp;
tite = 5;

% "Boundary conditions", "initial values"
% Frequency ratio
wpewce = 30;
% Pressure balance: B^2/2mu0 = nkBT

% Mass ratio (not used?)
mime = 25;

% Normalization vs local alfven speed


%wpewce = @(n,B,me) sqrt(n*me/units.eps0./B.^2);
n0 = @(B0,me,wpewce) wpewce.^2*units.eps0*B0.^2/me;
T0 = @(wpewce) wpewce.^-2*me*units.c^2/2/units.kB;
vA0 = @(B0,n0) B0./sqrt(mi.*n0*units.mu0);
cvte0 = @(tite,wpewce) sqrt(tite+1)*wpewce;

nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;

B0_ = 20e-9;
wpewce_ = 1:0.25:10;

hca = h(isub); isub = isub + 1;
plot(hca,wpewce_,vA0(B0_,n0(B0_,me,wpewce_))*1e-3)
hca.XLabel.String = '\omega_{pe0}/\omega_{ce0}';
hca.YLabel.String = 'V_{A0} (km/s)';
irf_legend(hca,sprintf('B_0 = %.0f nT',B0_*1e9),[0.02 0.98])
hca.XGrid = 'on';
hca.YGrid = 'on';

hca = h(isub); isub = isub + 1;
plot(hca,wpewce_,n0(B0_,me,wpewce_)*1e-6)
hca.XLabel.String = '\omega_{pe0}/\omega_{ce0}';
hca.YLabel.String = 'n_0 (cm^{-3})';
irf_legend(hca,sprintf('B_0 = %.0f nT',B0_*1e9),[0.02 0.98])
hca.XGrid = 'on';
hca.YGrid = 'on';

if 0
  hca = h(isub); isub = isub + 1;
  ax = plotyy(hca,wpewce_,T0(wpewce_)*units.kB/units.eV*1e-3,wpewce_,units.kB*T0(wpewce_)/me/units.c^2);
  hca.XLabel.String = '\omega_{pe0}/\omega_{ce0}';
  hca.YLabel.String = 'T_0 (keV)';
  ax(2).YLabel.String = 'k_BT_0 (m_ec^2)';
  if 1
    ax(1).YScale = 'log';
    ax(2).YScale = 'log';
    ax(1).YTick = 10.^[-10:1:10];
    ax(2).YTick = 10.^[-10:1:10];
  end
  hleg = irf_legend(hca,{sprintf('B_0 = %.0f nT',B0_*1e9);sprintf('T_0 = T_{i0} + T_{e0}, T_{i0}/T_{e0} = %g -> T_{e0} = T_0/%g',tite,tite+1)},[0.98 0.98],'k');
  hleg(1).Color = [0 0 0];
  hleg(2).Color = [0 0 0];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1
  hca = h(isub); isub = isub + 1;
  ax = plotyy(hca,wpewce_,T0(wpewce_)*units.kB/units.eV*1e-3,wpewce_,cvte0(tite,wpewce_));
  hca.XLabel.String = '\omega_{pe0}/\omega_{ce0}';
  hca.YLabel.String = 'T_0 (keV)';
  ax(2).YLabel.String = 'c/v_{te0}';
  if 1
    ax(1).YScale = 'log';
    %ax(2).YScale = 'log';
    ax(1).YTick = 10.^[-10:1:10];
    %ax(2).YTick = 10.^[-10:1:10];
    ax(2).YTick = 0:20:100;
    %ax(2).YGrid = 'on';
    ax(2).YLim = [0 70];
  end
  hleg = irf_legend(hca,{sprintf('B_0 = %.0f nT',B0_*1e9);sprintf('T_0 = T_{i0} + T_{e0}, T_{i0}/T_{e0} = %g -> T_{e0} = T_0/%g',tite,tite+1)},[0.1 0.98],'k');
  hleg(1).Color = [0 0 0];
  hleg(2).Color = [0 0 0];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end