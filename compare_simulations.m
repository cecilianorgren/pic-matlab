%% Load PIC object
df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
df08 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data_h5/fields.h5');

%% Scalar timeseries
h = setup_subplots(4,1);
isub = 1;

if 0 % Reconnction rate
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.RA,df08.twci,df08.RA)
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'R_A';
  legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
  hca.YLim(1) = 0;
end
if 0 % Reconnction rate, rescaled 
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.RA*sqrt(0.4/0.2),df08.twci,df08.RA*sqrt(0.8/0.2))
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'R_A(n_c/n_b)^{1/2}';
  legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
  irf_legend(hca,{sprintf('n_b = 0.2 cc',0.2)},[0.02 0.98],'color','k')
  hca.Title.String = 'Alfven scaling';
  hca.YLim(1) = 0;
end
if 0 % Magnetic energy
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.UB,df08.twci,df08.UB)
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_B';
  legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
end
if 0 % Particle energy
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.UT(1),df04.twci,df04.UT(2),df04.twci,df04.UT(3),df04.twci,df04.UT(4),df04.twci,df04.UT(5),df04.twci,df04.UT(6))
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_{T,K}';
  legend(hca,{'hot ions','hot electrons','cold ions','cold electrons','cold ions','cold electrons'},'location','bestoutside','box','off')
end
if 0 % Cold electron thermal energy
  hca = h(isub); isub = isub + 1;
  plot(hca,df08.twci,df08.UT(4),df04.twci,df04.UT(4),df04.twci,df04.UT(6),'--',df04.twci,df04.UT(46))
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_T';
  legend(hca,{'cold electrons, n_c = 0.8 cc','cold electrons top, n_c = 0.4 cc','cold electrons bottom, n_c = 0.4 cc','cold electrons all, n_c = 0.4 cc'},'location','bestoutside','box','off')
end
if 0 % Cold ion thermal energy
  hca = h(isub); isub = isub + 1;
  plot(hca,df08.twci,df08.UT(3),df04.twci,df04.UT(3),df04.twci,df04.UT(5),'--',df04.twci,df04.UT(35))
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_T';
  legend(hca,{'cold ions, n_c = 0.8 cc','cold ions top, n_c = 0.4 cc','cold ions bottom, n_c = 0.4 cc','cold ions all, n_c = 0.4 cc'},'location','bestoutside','box','off')
end
if 0 % Cold electron drift energy
  hca = h(isub); isub = isub + 1;
  plot(hca,df08.twci,df08.UK(4),df04.twci,df04.UK(4),df04.twci,df04.UK(6),'--',df04.twci,df04.UK(46))
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_K';
  legend(hca,{'cold electrons, n_c = 0.8 cc','cold electrons top, n_c = 0.4 cc','cold electrons bottom, n_c = 0.4 cc','cold electrons all, n_c = 0.4 cc'},'location','bestoutside','box','off')
end
if 0 % Cold ion drift energy
  hca = h(isub); isub = isub + 1;
  plot(hca,df08.twci,df08.UK(3),df04.twci,df04.UK(3),df04.twci,df04.UK(5),'--',df04.twci,df04.UK(35))
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_K';
  legend(hca,{'cold ions, n_c = 0.8 cc','cold ions top, n_c = 0.4 cc','cold ions bottom, n_c = 0.4 cc','cold ions all, n_c = 0.4 cc'},'location','bestoutside','box','off')
end
if 0 % Particle energy, ions
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.UT(1)+df04.UK(1),df04.twci,df04.UT(3)+df04.UK(3),df04.twci,df04.UT(5)+df04.UK(5))
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_T + U_K';
  legend(hca,{'hot ions','cold ions 1','cold ions 2'},'location','bestoutside','box','off')
end
if 1 % Particle energy, ions
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.UT(1)+df04.UK(1)+df04.UT(3)+df04.UK(3)+df04.UT(5)+df04.UK(5),df04.twci,df04.UT(135)+df04.UK(135))
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_T + U_K';
  legend(hca,{'separate added 1+3+5','combined 135'},'location','bestoutside','box','off')
end
if 1 % Particle energy, electrons
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.UT(2)+df04.UK(2)+df04.UT(4)+df04.UK(4)+df04.UT(6)+df04.UK(6),df04.twci,df04.UT(246)+df04.UK(246))
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_T + U_K';
  legend(hca,{'separate added 2+4+6','combined 246'},'location','bestoutside','box','off')
end
if 1 % Particle energy, ions cold
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.UT(3)+df04.UK(3)+df04.UT(5)+df04.UK(5),df04.twci,df04.UT(35)+df04.UK(35))
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_T + U_K';
  legend(hca,{'separate added 3+5','combined 35'},'location','bestoutside','box','off')
end
if 1 % Particle energy, cold electrons
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.UT(4)+df04.UK(4)+df04.UT(6)+df04.UK(6),df04.twci,df04.UT(46)+df04.UK(46))
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_T + U_K';
  legend(hca,{'separate added 4+6','combined 46'},'location','bestoutside','box','off')
end
%irf_plot_axis_align


%% Time (fan plots)

h = setup_subplots(4,1);
isub = 1;

zlim = [-0.1 0.1];
if 0 % Reconnction rate
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.RA,df08.twci,df08.RA)
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'R_A';
  legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
  hca.YLim(1) = 0;
end
if 1 % Pressure at midplane
  hca = h(isub); isub = isub + 1;
  sim_tmp = df04.zlim(zlim);
  iSp = 1; % hot ions
  Atmp = mean(sim_tmp.A,3);
  [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp([iSp]); 
  p = (pxx+pyy+pzz)/3;
  var = mean(p,3);
  imagesc(hca,sim_tmp.twci,sim_tmp.xi,var')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hb = colorbar('peer',hca);
  %legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
  %hca.YLim(1) = 0;
end
if 1 % Pressure at midplane
  hca = h(isub); isub = isub + 1;
  sim_tmp = df08.zlim(zlim);
  iSp = 1; % hot ions
  Atmp = mean(sim_tmp.A,3);
  [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = sim_tmp.njp([iSp]); 
  p = (pxx+pyy+pzz)/3;
  var = mean(p,3);
  imagesc(hca,sim_tmp.twci,sim_tmp.xi,var')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hb = colorbar('peer',hca);
  %legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
  %hca.YLim(1) = 0;
end
if 1 % Stress tensor at midplane
  hca = h(isub); isub = isub + 1;
  sim_tmp = df04.zlim(zlim);
  iSp = 1; % hot ions
  Atmp = mean(sim_tmp.A,3);
  [vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv([iSp]); 
  vv = (vxx+vyy+vzz)/3;
  var = mean(vv,3);
  imagesc(hca,sim_tmp.twci,sim_tmp.xi,var')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hb = colorbar('peer',hca);
  %legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
  %hca.YLim(1) = 0;
end
if 1 % Stress tensor at midplane
  hca = h(isub); isub = isub + 1;
  sim_tmp = df08.zlim(zlim);
  iSp = 1; % hot ions
  Atmp = mean(sim_tmp.A,3);
  [vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv([iSp]); 
  vv = (vxx+vyy+vzz)/3;
  var = mean(vv,3);
  imagesc(hca,sim_tmp.twci,sim_tmp.xi,var')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hb = colorbar('peer',hca);
  %legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
  %hca.YLim(1) = 0;
end
if 0 % Stress tensor at midplane
  hca = h(isub); isub = isub + 1;
  sim_tmp = df04.zlim(zlim);
  iSp = 1; % hot ions
  Atmp = mean(sim_tmp.A,3);
  [vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv([iSp]); 
  vv = (vxx+vyy+vzz)/3;
  var = mean(vv,3);imagesc(hca,sim_tmp.twci,sim_tmp.xi,var')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hb = colorbar('peer',hca);
  %legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
  %hca.YLim(1) = 0;
end
if 0 % Stress tensor at midplane
  hca = h(isub); isub = isub + 1;
  sim_tmp = df08.zlim(zlim);
  iSp = 1; % hot ions
  Atmp = mean(sim_tmp.A,3);
  [vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv([iSp]); 
  vv = (vxx+vyy+vzz)/3;
  var = mean(vv,3);
  imagesc(hca,[0 cumsum(diff(sim_tmp.UB))],sim_tmp.xi,var')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  hb = colorbar('peer',hca);
  %legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
  %hca.YLim(1) = 0;
end
if 0 % Stress tensor at midplane as function of spent UB
  hca = h(isub); isub = isub + 1;
  sim_tmp = df04.zlim(zlim);
  iSp = 1; % hot ions
  Atmp = mean(sim_tmp.A,3);
  [vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv([iSp]); 
  vv = (vxx+vyy+vzz)/3;
  var = mean(vv,3);
  imagesc(hca,abs([0 cumsum(diff(sim_tmp.UB))]),sim_tmp.xi,var')
  hca.XLabel.String = 'U_B(t=0)-U_B(t)';
  hca.YLabel.String = 'x/d_i';
  hb = colorbar('peer',hca);
  %legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
  %hca.YLim(1) = 0;
end
if 0 % Stress tensor at midplane as function of spent UB
  hca = h(isub); isub = isub + 1;
  sim_tmp = df08.zlim(zlim);
  iSp = 1; % hot ions
  Atmp = mean(sim_tmp.A,3);
  [vxx,vxy,vxz,vyy,vyz,vzz] = sim_tmp.vv([iSp]); 
  vv = (vxx+vyy+vzz)/3;
  var = mean(vv,3);
  imagesc(hca,abs([0 cumsum(diff(sim_tmp.UB))]),sim_tmp.xi,var')
  hca.XLabel.String = 'U_B(t=0)-U_B(t)';
  hca.YLabel.String = 'x/d_i';
  hb = colorbar('peer',hca);
  %legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'location','best','box','off')
  %hca.YLim(1) = 0;
end

