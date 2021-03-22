%% Single mass ratio, single temperature ratio
mi = 1;
mime = 100;
me = mi/mime;
vdf = 0.5; % vA

%% Harris sheet
vd = -0.0;
ttot = 0.5;
tite = 5;
ti = ttot/(1+1/tite);
te = ttot-ti;
vti = sqrt(2*ti/mi);
vte = sqrt(2*te/me);
fermi_ih = fun_fermi_acceleration(mi,vti,vd,vdf,1);
fermi_eh = fun_fermi_acceleration(me,vte,vd,vdf,1);
%fun_fermi_acceleration(m,vt,vd,vdf,doPlot)

%% Cold inflow
ttot = 0.02;
tite = 1;
ti = ttot/(1+1/tite);
te = ttot-ti;
vti = sqrt(2*ti/mi);
vte = sqrt(2*te/me);
fermi_ic = fun_fermi_acceleration(mi,vti,0,vdf,1);
fermi_ec = fun_fermi_acceleration(me,vte,-4,vdf,1);

%% One species, one vdf, but many vt
m = 1;
vdf = 1;
vd = -10:0.5:0; nvd = numel(vd);
% nvd = 5; linspace(-2,0,nvd);
nvt = 50; vt = logspace(-2,log10(20),nvt);
[VD,VT] = ndgrid(vd,vt);
vtcomb = zeros(nvd,nvt);
vdcomb = zeros(nvd,nvt);
clear fermi;
doPlot = 0;
nRefl = 1;
for ivt = 1:nvt
  for ivd = 1:nvd
    fermi_tmp = fun_fermi_acceleration(m,sqrt(2)*vt(ivt),vd(ivd),vdf,nRefl,doPlot);      
    vtcomb(ivd,ivt) = fermi_tmp.vtcomb;
    vdcomb(ivd,ivt) = fermi_tmp.vcomb;
    fermi(ivd,ivt) = fermi_tmp;
  end
end


%% Plot
colors = pic_colors('matlab');
nrows = 1;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
  
if 1 
  hca = h(isub); isub = isub + 1;
  hvt = plot(hca,vt/vdf,vtcomb,'-');
  hold(hca,'on')
  hvd = plot(hca,vt/vdf,vdcomb,'-.');
  hold(hca,'off')
  if 0 % cold plasma estimate, NOT WORKING
    % vd/vdf = 1
    % vt/vdf = 1 + |vbef|/vdf;
    % 
    ylim = hca.YLim;    
    hold(hca,'on')
    for ivd = 1:nvd
      vd_tmp = vd(ivd);      
      %ivt = find()
      
      hvd = plot(hca,vt/vdf,1+vd/vdf,'-.');
    end
    hold(hca,'off') 
    hca.YLim = ylim;
  end
  hca.XLabel.String = 'v_{t}^{bef}/v_{DF}';
  hca.YLabel.String = 'v_{t}, v_{d}';
  legend([hvd(1),hvt(1)],{'v_{d}','v_{t}'},'location','northwest')
  hca.Title.String = 'Bulk (v_d) and thermal (v_t) speed of combined populations';
end
if 0 % vt/vd
  hca = h(isub); isub = isub + 1;
  plot(hca,vt/vdf,vtcomb./vdcomb)
  hca.XLabel.String = 'v_{t}^{bef}/v_{DF}';
  hca.YLabel.String = 'v_{t}/v_{d}'; 
  if 1 % add Harris sheet
    
  end
end
if 1 % Ut/Uk lines
  hca = h(isub); isub = isub + 1;
  plot(hca,vt/vdf,(vtcomb./vdcomb).^2)
  hca.YScale = 'lin';
  hca.XLabel.String = 'v_{t}^{bef}/v_{DF}';
  hca.YLabel.String = 'U_t/U_k = (v_{t}/v_{d})^2'; 
  if 0 % Add Harris sheet values
    hold(hca,'on')
    hih = plot(hca,fermi_ih.vtbefore/vdf,(fermi_ih.vtcomb/fermi_ih.vcomb)^2,'*');
    heh = plot(hca,fermi_eh.vtbefore/vdf,(fermi_eh.vtcomb/fermi_eh.vcomb)^2,'*');
    hold(hca,'off')
    legend([hih heh],{'v_{i}^{Harris}','v_{e}^{Harris}'},'location','best')
  end
end
if 1 % Ut/Uk map
  hca = h(isub); isub = isub + 1;
%   pcolor(hca,vt/vdf,vd/vdf,(vtcomb./vdcomb).^2)
%   shading(hca,'flat')
  contourf(hca,vt/vdf,vd/vdf,(vtcomb./vdcomb).^2,[0:10:200])
  hca.XLabel.String = 'v_{t}^{bef}/v_{DF}';
  hca.YLabel.String = 'v_d^{bef}/v_{DF}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'U_t/U_k = (v_{t}/v_{d})^2';   
  colormap(hca,pic_colors('pasteljet'))
  if 0 % Add Harris sheet values
    hold(hca,'on')
    hih = plot(hca,fermi_ih.vtbefore/vdf,(fermi_ih.vtcomb/fermi_ih.vcomb)^2,'*');
    heh = plot(hca,fermi_eh.vtbefore/vdf,(fermi_eh.vtcomb/fermi_eh.vcomb)^2,'*');
    hold(hca,'off')
    legend([hih heh],{'v_{i}^{Harris}','v_{e}^{Harris}'},'location','southeast')
  end
end

for ip = 1:npanels
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).XLim = [0 20];
end

c_eval('h(?).Position(2) = 0.17;',1:3)
c_eval('h(?).Position(4) = 0.75;',1:3)
c_eval('h(?).Position(3) = 0.21;',1:3)
  
  
  