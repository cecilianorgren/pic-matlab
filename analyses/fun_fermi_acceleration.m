function out = fun_fermi_acceleration(m,vt,vd,vdf,nRefl,doPlot)
% out = fun_fermi_acceleration(m,vt,vd,vdf,nRefl,doPlot);
% Estimate the fermi acceleration temperature gain dependence on initial temperature
% v_after = -v_before + 2*v_DF.
% nRefl not really well implemented... Just use 1.

% Have to settle on one single definition of vt, should it be 
% vt = sqrt(T/m) or vt = sqrt(2T/m)?

% Maxwellian distribution
fun_f_max = @(n,v,vd,vt) n.*(1/pi./vt.^2)^(1/2)*exp(-(v-vd).^2./vt.^2);

% TiTe = 1;
% Ttot = 0.01;
% Ti = Ttot/(1 + 1/TiTe);
% Te = Ttot-Ti;
% 
% mime = 100;
% mi = 1;
% me = mi/mime;
% vte = sqrt(Te/me);
% vti = sqrt(Ti/mi);

n = 1;
%nRefl = 2;

% Initialise particles
nv = 20000;
vbefore = linspace(-5*vt,5*vt,nv);
fbefore = fun_f_max(n,vbefore,vd,vt);

% Velocity after reflection
vbefore_orig = vbefore;
for ir = 1:nRefl
  for iv = 1:numel(vbefore)
    if vbefore(iv) >= vdf
      vafter(iv) = vbefore(iv);
      fafter(iv) = 0;
    else
      vafter(iv) = -vbefore(iv) + 2*vdf;
      fafter(iv) = fbefore(iv);
    end
  end
  vbefore = -vafter;
end
vbefore = vbefore_orig;

% Calculate moments
% Initial distribution
dvbefore = [vbefore(2)-vbefore(1) diff(vbefore)];
mom.fbefore = fbefore;
mom.vvecbefore = vbefore;
mom.dvbefore = dvbefore;
mom.nbefore = nansum(fbefore.*dvbefore);
mom.jbefore = nansum(fbefore.*vbefore.*dvbefore);
mom.vbefore = mom.jbefore/mom.nbefore;
mom.pbefore = m*nansum(fbefore.*(vbefore-mom.vbefore).^2.*dvbefore);
mom.tbefore = mom.pbefore/mom.nbefore;
mom.vtbefore = sqrt(1*mom.tbefore/m);

% Transmitted distribution
dvafter = abs([vafter(2)-vafter(1) diff(vafter)]);
mom.fafter = fafter;
mom.vvecafter = vafter;
mom.dvafter = dvafter;
mom.nafter = nansum(fafter.*dvafter);
mom.jafter = nansum(fafter.*vafter.*dvafter);
mom.vafter = mom.jafter/mom.nafter;
mom.pafter = m*nansum(fafter.*(vafter-mom.vafter).^2.*dvafter);
mom.tafter = mom.pafter/mom.nafter;
mom.vtafter = sqrt(1*mom.tafter/m);

% Inital + transmitted distribution
% They are on different grids (vbefore, vafter) and needs to be resampled
nvcomb = 200;
vcomb_edges = linspace(min([vbefore vafter]),max([vbefore vafter]),nvcomb+1);
dvcomb = vcomb_edges(2)-vcomb_edges(1);
vcomb = [vcomb_edges(1:end-1)+0.5*dvcomb];
V = [vbefore vafter];
F = [fbefore fafter];
[dncomb,edges,mid,loc] = histcn([vbefore vafter]',vcomb_edges,'AccumData',[fbefore.*dvbefore fafter.*dvafter]');
fcomb = dncomb'./dvcomb;

mom.fcomb = fcomb;
mom.vveccomb = vcomb;
mom.dvcomb = dvcomb;
mom.ncomb = nansum(fcomb.*dvcomb);
mom.jcomb = nansum(fcomb.*vcomb.*dvcomb);
mom.vcomb = mom.jcomb/mom.ncomb;
mom.pcomb = m*nansum(fcomb.*(vcomb-mom.vcomb).^2.*dvcomb);
mom.tcomb = mom.pcomb/mom.ncomb;
mom.vtcomb = sqrt(1*mom.tcomb/m);

out = mom;

if doPlot
  colors = pic_colors('matlab');
  nrows = 1;
  ncols = 1;
  npanels = nrows*ncols;
  h = setup_subplots(nrows,ncols);
  isub = 1;

  if 1 % Initial and transmitted distribution
    hca = h(isub); isub = isub + 1;
    % vdf
    plot(hca,vdf,0,'*','color',[0 0 0])    
    hold(hca,'on')
    % f    
    plot(hca,vbefore,fbefore,'-','color',colors(1,:),'linewidth',1)    
    plot(hca,vafter,fafter,'-','color',colors(2,:),'linewidth',1)
    plot(hca,vcomb,fcomb,'-','color',colors(3,:),'linewidth',1)
    % vbulk    
    plot(hca,mom.vbefore*[1 1],hca.YLim,'color',colors(1,:),'linestyle','--')
    plot(hca,mom.vafter*[1 1],hca.YLim,'color',colors(2,:),'linestyle','--')
    plot(hca,mom.vcomb*[1 1],hca.YLim,'color',colors(3,:),'linestyle','--')
    hold(hca,'off')
    
    hca.YLabel.String = 'Phase space density';
    hca.XLabel.String = 'v';
    legend(hca,{'v_{DF}','f^{before}','f^{after}','f^{comb}','v_{bulk}^{before}','v_{bulk}^{after}','v_{bulk}^{comb}'})
    irf_legend(hca,{'f^{before}:',sprintf('n = %.3f',mom.nbefore),sprintf('v_d = %.3f',mom.vbefore),sprintf('v_t = %.3f',mom.vtbefore),sprintf('T = %.3f',mom.tbefore),sprintf('P = %.3f',mom.pbefore)}',[0.02 0.98],'color',[0 0 0])
    irf_legend(hca,{'f^{after}:',sprintf('n = %.3f',mom.nafter),sprintf('v_d = %.3f',mom.vafter),sprintf('v_t = %.3f',mom.vtafter),sprintf('T = %.3f',mom.tafter),sprintf('P = %.3f',mom.pafter)}',[0.2 0.98],'color',[0 0 0])
    irf_legend(hca,{'f^{comb}:',sprintf('n = %.3f',mom.ncomb),sprintf('v_d = %.3f',mom.vcomb),sprintf('v_t = %.3f',mom.vtcomb),sprintf('T = %.3f',mom.tcomb),sprintf('P = %.3f',mom.pcomb)}',[0.02 0.02],'color',[0 0 0])
    
    irf_legend(hca,{sprintf('m = %.3f',m)}',[0.98 0.02],'color',[0 0 0])
  end
  if 0 % Flux of initial and transmitted distribution
    hca = h(isub); isub = isub + 1;
    plot(hca,vbefore,fbefore.*vbefore,'o',vafter,fbefore.*vafter,'*')
    hold(hca,'on')
    plot(hca,vdf*[1 1],hca.YLim)
    plot(hca,mom.vbefore*[1 1],hca.YLim)
    plot(hca,mom.vafter*[1 1],hca.YLim)
    hold(hca,'off')
    hca.YLabel.String = 'Flux';
    hca.XLabel.String = 'v';
  end
  if 0
    hca = h(isub); isub = isub + 1;
    hca = subplot(3,1,3);
    plot(hca,vbefore,vafter,'*')
    hold(hca,'on')
    plot(hca,vdf*[1 1],hca.YLim)
    hold(hca,'off')
    hca.XLabel.String = 'v_{before}';
    hca.YLabel.String = 'v_{after}';
  end
end









