% compare_energy_partitioning
%% Load simulations 
df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5')
df08 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data_h5/fields.h5')

%% Plots

nrows = 2;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 0 % UB(twci)
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.UB,df08.twci,df08.UB)
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_{B}';
  legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'Box','off','location','best')
end
if 0 % UB(twci)/UB(0)
  hca = h(isub); isub = isub + 1;  
  plot(hca,df04.twci,df04.UB/df04.i(1).UB,df08.twci,df08.UB./df08.i(1).UB)
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_{B}/U_B(t=0)';
  legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'Box','off','location','best')
  hca.YLim(1) = 0;
end
if 0 % dUB(twci)
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.dUB,df08.twci,df08.dUB)
  hca.XLabel.String = 't\omega_{ci}';  
  hca.YLabel.String = '|U_{B}-U_B(t=0)|';
  legend(hca,{'n_c = 0.4 cc','n_c = 0.8 cc'},'Box','off','location','best')
end
if 0 % Partitioning between drift and thermal energy for cold ions in df04 run normalized to UB(1)
  hca = h(isub); isub = isub + 1;
  colors = pic_colors('matlab');
  hold(hca,'on');
  hca.ColorOrder = colors(1:2,:);
  hca.LineStyleOrder = {'-','-.','--'};
  
  plot(hca,df04.twci,df04.UK(3)/df04(1).UB,...
           df04.twci,df04.UT(3)/df04(1).UB,...
           df04.twci,df04.UK(5)/df04(1).UB,...
           df04.twci,df04.UT(5)/df04(1).UB,'linewidth',1.5)
  %hold(hca,'on')
  plot(hca,df04.twci,df04.UK(35)/df04(1).UB,'--',...
           df04.twci,df04.UT(35)/df04(1).UB,'--','linewidth',1.5)  
  hca.ColorOrder = colors(1:2,:);
  hca.LineStyleOrder = {'-','--'};
  hca.XLabel.String = 't\omega_{ci}';  
  hca.YLabel.String = 'U_i^c/U_B(t=0)';
  hleg = legend(hca,{'U_K^{c1}','U_T^{c1}','U_K^{c2}','U_T^{c2}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','on','location','northwest');
  hleg.Title.String = 'n_c = 0.4 n_0';  
  hold(hca,'off')  
end
if 0 % Partitioning between drift and thermal energy for cold ions in df04 run normalized to UB(1)
  hca = h(isub); isub = isub + 1;
  colors = pic_colors('matlab');
  hold(hca,'on');
  hca.ColorOrder = colors(1:2,:);
  hca.LineStyleOrder = {'-','--'};
  
  plot(hca,df04.twci,2*df04.UK(3)/df04(1).UB,...
           df04.twci,2*df04.UT(3)/df04(1).UB,'linewidth',1.5)
  %hold(hca,'on')
  plot(hca,df04.twci,df04.UK(35)/df04(1).UB,'--',...
           df04.twci,df04.UT(35)/df04(1).UB,'--','linewidth',1.5)
  hca.XLabel.String = 't\omega_{ci}';  
  hca.YLabel.String = 'U_i^c/U_B(t=0)';
  hleg = legend(hca,{'2U_K^{c1}','2U_T^{c1}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','on','location','northwest');
  hleg.Title.String = 'n_c = 0.4 n_0';
  hold(hca,'off')
end
if 0 % Partitioning between drift and thermal energy for cold ions in df04 run
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,df04.UK(3),'--',...
           df04.twci,df04.UT(3),'-',...
           df04.twci,df04.UK(5),'--',...
           df04.twci,df04.UT(5),'-')
  hold(hca,'on')
  plot(hca,df04.twci,df04.UK(35),'--',...
           df04.twci,df04.UT(35),'-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 't\omega_{ci}';  
  hca.YLabel.String = 'U_i^c';
  hleg = legend(hca,{'U_K^{c1}','U_T^{c1}','U_K^{c2}','U_T^{c2}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','on','location','northwest');
  hleg.Title.String = 'n_c = 0.4 n_0';
end
if 0 % Partitioning between drift and thermal energy for cold ions in df04 run 
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.twci,2*df04.UK(3),'--',...
           df04.twci,2*df04.UT(3),'-')
  hold(hca,'on')
  plot(hca,df04.twci,df04.UK(35),'--',...
           df04.twci,df04.UT(35),'-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 't\omega_{ci}';  
  hca.YLabel.String = 'U_i^c';
  hleg = legend(hca,{'2*U_K^{c1}','2*U_T^{c1}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','on','location','northwest');
  hleg.Title.String = 'n_c = 0.4 n_0';
end
if 0 % Partitioning between drift and thermal energy for cold ions in df04 run
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.dUB,df04.UK(3),'--',...
           df04.dUB,df04.UT(3),'-',...
           df04.dUB,df04.UK(5),'--',...
           df04.dUB,df04.UT(5),'-')
  hold(hca,'on')
  plot(hca,df04.dUB,df04.UK(35),'--',...
           df04.dUB,df04.UT(35),'-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = '|U_{B}-U_B(t=0)|';  
  hca.YLabel.String = 'U_i^c';
  legend(hca,{'U_K^{c1}','U_T^{c1}','U_K^{c2}','U_T^{c2}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','off','location','best')
end
if 0 % Partitioning between drift and thermal energy for cold ions in df04 run
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.UB/df04.i(1).UB,df04.UK(3),'--',...
           df04.UB/df04.i(1).UB,df04.UT(3),'-',...
           df04.UB/df04.i(1).UB,df04.UK(5),'--',...
           df04.UB/df04.i(1).UB,df04.UT(5),'-')
  hold(hca,'on')
  plot(hca,df04.UB/df04.i(1).UB,df04.UK(35),'--',...
           df04.UB/df04.i(1).UB,df04.UT(35),'-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 'U_{B}/U_B(t=0)';
  hca.YLabel.String = 'U_i^c';
  hca.XDir = 'reverse';
  hca.XLim(2) = 1;
  legend(hca,{'U_K^{c1}','U_T^{c1}','U_K^{c2}','U_T^{c2}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','off','location','best')
end
if 0 % Partitioning between drift and thermal energy for cold ions in df04 run
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.UB/df04.i(1).UB,df04.UK(3)/df04.i(1).UB,'--',...
           df04.UB/df04.i(1).UB,df04.UT(3)/df04.i(1).UB,'-',...
           df04.UB/df04.i(1).UB,df04.UK(5)/df04.i(1).UB,'--',...
           df04.UB/df04.i(1).UB,df04.UT(5)/df04.i(1).UB,'-')
  hold(hca,'on')
  plot(hca,df04.UB/df04.i(1).UB,df04.UK(35)/df04.i(1).UB,'--',...
           df04.UB/df04.i(1).UB,df04.UT(35)/df04.i(1).UB,'-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 'U_{B}/U_B(t=0)';
  hca.YLabel.String = 'U_i^c/U_B(t=0)';
  hca.XDir = 'reverse';
  hca.XLim(2) = 1;
  hleg = legend(hca,{'U_K^{c1}','U_T^{c1}','U_K^{c2}','U_T^{c2}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','off','location','best');
  hleg.Title.String = 'Cold ions';
end
if 0 % Partitioning between drift and thermal energy for cold electrons in df04 run
  hca = h(isub); isub = isub + 1;
  ieref = 20;
  plot(hca,df04.UB/df04.i(1).UB,df04.UK(4)/df04.i(1).UB,'--',...
           df04.UB/df04.i(1).UB,(df04.UT(4)-df04.i(ieref).UT(4))/df04.i(1).UB,'-',...
           df04.UB/df04.i(1).UB,df04.UK(6)/df04.i(1).UB,'--',...
           df04.UB/df04.i(1).UB,(df04.UT(6)-df04.i(ieref).UT(6))/df04.i(1).UB,'-')
  hold(hca,'on')
  plot(hca,df04.UB/df04.i(1).UB,(df04.UK(46)-df04.i(ieref).UK(46))/df04.i(1).UB,'--',...
           df04.UB/df04.i(1).UB,(df04.UT(46)-df04.i(ieref).UT(46))/df04.i(1).UB,'-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 'U_{B}/U_B(t=0)';
  hca.YLabel.String = 'U_e^c/U_B(t=0)';
  hca.XDir = 'reverse';
  hca.XLim(2) = 1;
  hleg = legend(hca,{'U_K^{c1}',...
    sprintf('U_T^{c1} (t_{ref}w_{ci} = %g)',df04.twci(ieref)),...
    'U_K^{c2}',...
    sprintf('U_T^{c2} (t_{ref}w_{ci} = %g)',df04.twci(ieref)),...
    'U_K^{c1+c2}',...
    sprintf('U_T^{c1+c2} (t_{ref}w_{ci} = %g)',df04.twci(ieref))},...
    'Box','off','location','best');
  hleg.Title.String = 'Cold electrons';
end
if 0 % Partitioning between drift and thermal energy for cold and hot electrons in df04 run
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.UB/df04.i(1).UB,df04.UK(4)/df04.i(1).UB,'--',...
           df04.UB/df04.i(1).UB,df04.UT(4)/df04.i(1).UB,'-')
  hca.XLabel.String = 'U_{B}/U_B(t=0)';
  hca.YLabel.String = 'U_e^c/U_B(t=0)';
  hca.XDir = 'reverse';
  hca.XLim(2) = 1;
  %hleg = legend(hca,{'U_K^{c1}','U_T^{c1}','U_K^{c2}','U_T^{c2}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','off','location','best');
  %hleg.Title.String = 'Cold electrons';
end
if 0 % Partitioning between drift and thermal energy for cold and hot electrons in df04 run
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.UB/df04.i(1).UB,df04.UK(1)/df04.i(1).UB,'-',...
           df04.UB/df04.i(1).UB,df04.UT(1)/df04.i(1).UB,'-',...
           df04.UB/df04.i(1).UB,df04.UK(46)/df04.i(1).UB,'--',...
           df04.UB/df04.i(1).UB,df04.UT(46)/df04.i(1).UB,'--')
  hca.XLabel.String = 'U_{B}/U_B(t=0)';
  hca.YLabel.String = 'U_e^c/U_B(t=0)';
  hca.XDir = 'reverse';
  hca.XLim(2) = 1;
  %hleg = legend(hca,{'U_K^{c1}','U_T^{c1}','U_K^{c2}','U_T^{c2}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','off','location','best');
  %hleg.Title.String = 'Cold electrons';
end
if 0 % Partitioning between drift and thermal energy for cold and hot electrons in df04 run
  hca = h(isub); isub = isub + 1;
  plot(hca,df04.UK(1)/df04.i(1).UB,df04.UT(1)/df04.i(1).UB,'-',...
           df04.UK(46)/df04.i(1).UB,df04.UT(46)/df04.i(1).UB,'--')
  hca.XLabel.String = 'U_{B}/U_B(t=0)';
  hca.YLabel.String = 'U_e^c/U_B(t=0)';
  %hca.XDir = 'reverse';
  %hca.XLim(2) = 1;
  %hleg = legend(hca,{'U_K^{c1}','U_T^{c1}','U_K^{c2}','U_T^{c2}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','off','location','best');
  %hleg.Title.String = 'Cold electrons';
end
if 0 % Partitioning between cold and hot ions and electrons in df04 run
  hca = h(isub); isub = isub + 1;
  uih = df04.UK(1)+df04.UT(1); 
  ueh = df04.UK(2)+df04.UT(2);
  uic = df04.UK(35)+df04.UT(35);
  uec = df04.UK(46)+df04.UT(46);
  utut = uih+ueh+uic+uec;
  plot(hca,df04.UB/df04.i(1).UB,uih/df04.i(1).UB,'-',...
           df04.UB/df04.i(1).UB,ueh/df04.i(1).UB,'--',...
           df04.UB/df04.i(1).UB,uic/df04.i(1).UB,'-',...
           df04.UB/df04.i(1).UB,uec/df04.i(1).UB,'--')
  hold(hca,'on')
  plot(hca,df04.UB/df04.i(1).UB,(df04.UK(135)+df04.UK(246)+df04.UT(135)+df04.UT(246))/df04.i(1).UB,'k-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 'U_{B}/U_B(t=0)';
  hca.YLabel.String = 'U/U_B(t=0)';
  hca.XDir = 'reverse';
  hca.XLim(2) = 1;
  hleg = legend(hca,{'U^{ih}','U^{eh}','U^{ic}','U^{ec}','U^{tot}'},'Box','off','location','best');
  hleg.Title.String = 'Cold electrons';
end
if 0 % Partitioning between cold and hot ions and electrons in df04 run
  hca = h(isub); isub = isub + 1;
  uih = df04.UK(1)+df04.UT(1); uih = uih-uih(2);
  ueh = df04.UK(2)+df04.UT(2); ueh = ueh-ueh(2);
  uic = df04.UK(35)+df04.UT(35); uic = uic-uic(2);
  uec = df04.UK(46)+df04.UT(46); uec = uec-uec(2);
  utot = uih+ueh+uic+uec;
  plot(hca,df04.UB/df04.i(1).UB,uih/df04.i(1).UB,'-',...¨
           df04.UB/df04.i(1).UB,ueh/df04.i(1).UB,'--',...
           df04.UB/df04.i(1).UB,uic/df04.i(1).UB,'-',...
           df04.UB/df04.i(1).UB,uec/df04.i(1).UB,'--')
  hold(hca,'on')
  plot(hca,df04.UB/df04.i(1).UB,utot/df04.i(1).UB,'k-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 'U_{B}/U_B(t=0)';
  hca.YLabel.String = 'U/U_B(t=0)';
  hca.XDir = 'reverse';
  hca.XLim(2) = 1;
  hleg = legend(hca,{'U^{ih}','U^{eh}','U^{ic}','U^{ec}','U^{tot}'},'Box','off','location','best');
  hleg.Title.String = 'Cold electrons';
end
if 0 % Partitioning between cold and hot ions and electrons in df04 run
  hca = h(isub); isub = isub + 1;
  uih = df04.UK(1)+df04.UT(1); uih = uih-uih(2);
  ueh = df04.UK(2)+df04.UT(2); ueh = ueh-ueh(2);
  uic = df04.UK(35)+df04.UT(35); uic = uic-uic(2);
  ieref = 25;
  uec = df04.UK(46)+df04.UT(46); uec = uec-uec(ieref);
  utot = uih+ueh+uic+uec;
  uu = [uih;ueh;uic;uec]'/df04.i(1).UB;
  xx = df04.UB/df04.i(1).UB;
  harea = area(hca,xx,uu);
  c_eval('harea(?).FaceAlpha = 0.5;',1:numel(harea))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  plot(hca,df04.UB/df04.i(1).UB,utot/df04.i(1).UB,'k-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 'U_{B}/U_B(t=0)';
  hca.YLabel.String = '(U-U_{t=0})/U_B_{t=0}';
  hca.XDir = 'reverse';
  hca.XLim(2) = 1;
  hca.YLim(1) = 0;
  hleg = legend(hca,{'U^{ih}','U^{eh}','U^{ic}',sprintf('U^{ec} (t_{ref}w_{ci} = %g)',df04.twci(ieref)),'U^{tot}'},'Box','off','location','best');  
end
if 0 % Partitioning between cold and hot ions and electrons in df04 run / per particle
  hca = h(isub); isub = isub + 1;
  %Nih = sum(sum(df04.i(1).n(1)));
  %Neh = sum(sum(df04.i(1).n(2)));
  %Nic = sum(sum(df04.i(1).n(3)))+sum(sum(df04.i(1).n(5)));
  %Nec = sum(sum(df04.i(1).n(4)))+sum(sum(df04.i(1).n(6)));
  Nih = 0.2;
  Nic = 0.4;
  Neh = 0.2;
  Nec = 0.4;
  Ni = Nih + Nic;
  Ne = Neh + Nec;

  
  uih = df04.UK(1)+df04.UT(1); uih = uih-uih(2);
  ueh = df04.UK(2)+df04.UT(2); ueh = ueh-ueh(2);
  uic = df04.UK(35)+df04.UT(35); uic = uic-uic(2);
  ieref = 25;
  uec = df04.UK(46)+df04.UT(46); uec = uec-uec(ieref);
  utot = uih+ueh+uic+uec;
  uu = [uih*Ni/Nih;ueh*Ne/Neh;uic*Ni/Nic;uec*Ne/Nec]'/df04.i(1).UB;
  xx = df04.UB/df04.i(1).UB;
  harea = area(hca,xx,uu);
  c_eval('harea(?).FaceAlpha = 0.5;',1:numel(harea))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  %plot(hca,df04.UB/df04.i(1).UB,utot/df04.i(1).UB/Ni,'k-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 'U_{B}/U_B(t=0)';
  hca.YLabel.String = '(U-U_{t=0})/U_B_{t=0}';
  hca.XDir = 'reverse';
  hca.XLim(2) = 1;
  hca.YLim(1) = 0;
  hleg = legend(hca,{'U^{ih}','U^{eh}','U^{ic}',sprintf('U^{ec} (t_{ref}w_{ci} = %g)',df04.twci(ieref)),'U^{tot}'},'Box','off','location','best');  
end

if 1 % Partitioning between drift and thermal energy for cold ions in df04 run normalized to UB(1)
  hca = h(isub); isub = isub + 1;
  colors = pic_colors('matlab');
  hold(hca,'on');
  hca.ColorOrder = colors(1:2,:);
  hca.LineStyleOrder = {'-','-.','--'};
  
  plot(hca,df04.twci,df04.UK(4)/df04(1).UB,...
           df04.twci,df04.UT(4)/df04(1).UB,...
           df04.twci,df04.UK(6)/df04(1).UB,...
           df04.twci,df04.UT(6)/df04(1).UB,'linewidth',1.5)
  %hold(hca,'on')
  plot(hca,df04.twci,df04.UK(46)/df04(1).UB,'--',...
           df04.twci,df04.UT(46)/df04(1).UB,'--','linewidth',1.5)  
  hca.ColorOrder = colors(1:2,:);
  hca.LineStyleOrder = {'-','--'};
  hca.XLabel.String = 't\omega_{ci}';  
  hca.YLabel.String = 'U_e^c/U_B(t=0)';
  hleg = legend(hca,{'U_K^{c1}','U_T^{c1}','U_K^{c2}','U_T^{c2}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','on','location','northwest');
  hleg.Title.String = 'n_c = 0.4 n_0';  
  hold(hca,'off')  
end
if 1 % Partitioning between drift and thermal energy for cold ions in df04 run normalized to UB(1)
  hca = h(isub); isub = isub + 1;
  colors = pic_colors('matlab');
  hold(hca,'on');
  hca.ColorOrder = colors(1:2,:);
  hca.LineStyleOrder = {'-','--'};
  
  plot(hca,df04.twci,2*df04.UK(4)/df04(1).UB,...
           df04.twci,2*df04.UT(4)/df04(1).UB,'linewidth',1.5)
  %hold(hca,'on')
  plot(hca,df04.twci,df04.UK(46)/df04(1).UB,'--',...
           df04.twci,df04.UT(46)/df04(1).UB,'--','linewidth',1.5)
  hca.XLabel.String = 't\omega_{ci}';  
  hca.YLabel.String = 'U_e^c/U_B(t=0)';
  hleg = legend(hca,{'2U_K^{c1}','2U_T^{c1}','U_K^{c1+c2}','U_T^{c1+c2}'},'Box','on','location','northwest');
  hleg.Title.String = 'n_c = 0.4 n_0';
  hold(hca,'off')
end
if 0
  hca = h(isub); isub = isub + 1;
end
if 0
  hca = h(isub); isub = isub + 1;
end
if 0
  hca = h(isub); isub = isub + 1;
end

c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on''; h(?).Box = ''on'';',1:npanels)

%% stacked area plot of energy content
h = setup_subplots(2,1);
isub = 1;

if 1 % Partitioning between cold and hot ions and electrons in df04 run
  pic = df04;
  iref = 15;
  hca = h(isub); isub = isub + 1;
  uih = pic.UK(1)+pic.UT(1); uih = uih-uih(iref);
  ueh = pic.UK(2)+pic.UT(2); ueh = ueh-ueh(iref);
  uic = pic.UK(35)+pic.UT(35); uic = uic-uic(iref);
  ieref = 20;
  uec = pic.UK(46)+pic.UT(46); uec = uec-uec(ieref);
  utot = uih+ueh+uic+uec;
  uu = [uih;ueh;uic;uec]'/pic.i(1).UB;
  xx = pic.UB/pic.i(1).UB;
  harea = area(hca,xx,uu);
  c_eval('harea(?).FaceAlpha = 0.5;',1:numel(harea))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  plot(hca,pic.UB/pic.i(1).UB,utot/pic.i(1).UB,'k-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 'U_{B}/U_B(tw_{ci}=4)';
  hca.YLabel.String = sprintf('(U-U_{twci=%g})/U_B_{t=0}',pic.twci(iref));
  hca.XDir = 'reverse';
  hca.XLim(2) = 1;
  %hca.YLim(1) = 0;
  hold(hca,'on')
  plot(hca,hca.XLim,1-hca.XLim,'k-.')
  hold(hca,'off')
  hleg = legend(hca,{'U^{ih}','U^{eh}','U^{ic}',sprintf('U^{ec} (t_{ref}w_{ci} = %g)',pic.twci(ieref)),'U^{tot}','\Delta (U_T+U_K) = -\Delta U_B'},'Box','off','location','best');  
end

if 1 % Partitioning between cold and hot ions and electrons in df08 run
  pic = df08;
  iref = 20;
  hca = h(isub); isub = isub + 1;
  uih = pic.UK(1)+pic.UT(1); uih = uih-uih(iref);
  ueh = pic.UK(2)+pic.UT(2); ueh = ueh-ueh(iref);
  uic = pic.UK(3)+pic.UT(3); uic = uic-uic(iref);
  ieref = 20;
  uec = pic.UK(4)+pic.UT(4); uec = uec-uec(ieref);
  utot = uih+ueh+uic+uec;
  uu = [uih;ueh;uic;uec]'/pic.i(1).UB;
  xx = pic.UB/pic.i(1).UB;
  harea = area(hca,xx,uu);
  c_eval('harea(?).FaceAlpha = 0.5;',1:numel(harea))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  plot(hca,pic.UB/pic.i(1).UB,utot/pic.i(1).UB,'k-','linewidth',1.5)
  hold(hca,'off')
  hca.XLabel.String = 'U_{B}/U_B(tw_{ci}=4)';
  hca.YLabel.String = sprintf('(U-U_{twci=%g})/U_B_{t=0}',pic.twci(iref));
  hca.XDir = 'reverse';
  hca.XLim(2) = 1;
  %hca.YLim(1) = 0;  
  hold(hca,'on')
  plot(hca,hca.XLim,1-hca.XLim,'k-.')
  hold(hca,'off')
  hleg = legend(hca,{'U^{ih}','U^{eh}','U^{ic}',sprintf('U^{ec} (t_{ref}w_{ci} = %g)',pic.twci(ieref)),'U^{tot}','\Delta (U_T+U_K) = -\Delta U_B'},'Box','off','location','best');  
end

%% thermalization
h = setup_subplots(2,1);
isub = 1;

it = 30;

if 1 % 
  pic = df04;  
  hca = h(isub); isub = isub + 1;
  imagesc(pic.xi,pic.zi,pic.i(it).)
  
end
