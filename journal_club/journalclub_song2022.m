h = setup_subplots(3,1);
isub = 1;

pic_comb = no02m;%no02m.twpelim(15000:1000:24000);
Uti_hot = get_timeline_attributes(pic_comb,'Uti_hot');
Ute_hot = get_timeline_attributes(pic_comb,'Ute_hot');
Uki_hot = get_timeline_attributes(pic_comb,'Uki_hot');
Uke_hot = get_timeline_attributes(pic_comb,'Uke_hot');
Uki_cold = get_timeline_attributes(pic_comb,'Uki_cold');
Uke_cold = get_timeline_attributes(pic_comb,'Uke_cold');
Uti_cold = get_timeline_attributes(pic_comb,'Uti_cold');
Ute_cold = get_timeline_attributes(pic_comb,'Ute_cold');

if 1 % RE
  hca = h(isub); isub = isub + 1;
  plot(hca,no02m.twci,no02m.RE,'-')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'R_E';
end
if 0 % UB
  hca = h(isub); isub = isub + 1;
  plot(hca,no02m.twci,no02m.UB,'o')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_B';
end
if 0 % UT
  hca = h(isub); isub = isub + 1;
  plot(hca,no02m.twci,no02m.UT(1:6)-repmat(no02m(1).UT(1:6),no02m.length,1),'*')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_T';
  legend(hca,{'i hot','e hot','i cold top','e cold top','i cold bot','e cold bot'},...
    'location','best')
end
if 0 % UK
  hca = h(isub); isub = isub + 1;
  plot(hca,no02m.twci,no02m.UK(1:6)-repmat(no02m(1).UK(1:6),no02m.length,1),'x')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_K';
  legend(hca,{'i hot','e hot','i cold top','e cold top','i cold bot','e cold bot'},...
    'location','best')
end
if 0 % UT, UK
  hca = h(isub); isub = isub + 1;
  linestyle_cold = '-';
  linestyle_hot = '-';
  plot(hca,no02m.twci,Uti_hot,'-',...
           no02m.twci,Ute_hot,'-',...
           no02m.twci,Uki_hot,'--',...
           no02m.twci,Uke_hot,'--',...
           no02m.twci,Uti_cold,'-',...
           no02m.twci,Ute_cold,'-',...
           no02m.twci,Uki_cold,'--',...
           no02m.twci,Uke_cold,'--')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U';
  legend(hca,{'U_{ti}^{hot}','U_{te}^{hot}','U_{ki}^{hot}','U_{tke}^{hot}',...
              'U_{ti}^{cold}','U_{te}^{cold}','U_{ki}^{cold}','U_{tke}^{cold}'},...
    'location','eastoutside')
end
if 0 % UT+UK
  hca = h(isub); isub = isub + 1;
  linestyle_cold = '-';
  linestyle_hot = '-';
  plot(hca,no02m.twci,Uti_hot+Uki_hot,'-',...
           no02m.twci,Ute_hot+Uke_hot,'-',...
           no02m.twci,Uti_cold+Uki_cold,'-',...
           no02m.twci,Ute_cold+Uke_cold,'-')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U';
  legend(hca,{'U_{i}^{hot}','U_{e}^{hot}',...
              'U_{i}^{cold}','U_{e}^{cold}'},...
    'location','eastoutside')
end
if 1 % dUT, dUK
  hca = h(isub); isub = isub + 1;
  linestyle_cold = '-';
  linestyle_hot = '-';
  plot(hca,no02m.twci,Uti_hot-Uti_hot(1) + Uki_hot-Uki_hot(1),'-',...
           no02m.twci,Ute_hot-Ute_hot(1) + Uke_hot-Uke_hot(1),'-',...
           no02m.twci,Uti_cold-Uti_cold(1) + Uki_cold-Uki_cold(1),'-',...
           no02m.twci,Ute_cold-Ute_cold(1) + Uke_cold-Uke_cold(1),'-')           
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U-U(1)';
  legend(hca,{'U_{i}^{hot}','U_{e}^{hot}',...
              'U_{i}^{cold}','U_{e}^{cold}'},...
    'location','eastoutside')
end
if 1 % dUT, dUK average
  hca = h(isub); isub = isub + 1;
  linestyle_cold = '-';
  linestyle_hot = '-';
  plot(hca,no02m.twci,(Uti_hot-Uti_hot(1) + Uki_hot-Uki_hot(1))/no02m.nx/no02m.nz,'-',...
           no02m.twci,(Ute_hot-Ute_hot(1) + Uke_hot-Uke_hot(1))/no02m.nx/no02m.nz,'-',...
           no02m.twci,(Uti_cold-Uti_cold(1) + Uki_cold-Uki_cold(1))/no02m.nx/no02m.nz,'-',...
           no02m.twci,(Ute_cold-Ute_cold(1) + Uke_cold-Uke_cold(1))/no02m.nx/no02m.nz,'-')           
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = '(U-U(1))/n_xn_z';
  legend(hca,{'U_{i}^{hot}','U_{e}^{hot}',...
              'U_{i}^{cold}','U_{e}^{cold}'},...
    'location','eastoutside')
end
if 0 % UK
  hca = h(isub); isub = isub + 1;
  plot(hca,no02m.twci,no02m.UK(1:6)-repmat(no02m(1).UK(1:6),no02m.length,1),'x')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_K';
  legend(hca,{'i hot','e hot','i cold top','e cold top','i cold bot','e cold bot'},...
    'location','best')
end

if 0 % UB+UK+UT
  hca = h(isub); isub = isub + 1;
  plot(hca,no02m.UB + sum(no02m.UK(1:6),2) + sum(no02m.UT(1:6),2),'x')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_B+U_K+U_T';
end

if 0 % UT, UK cold ions comparison, separate/together
  hca = h(isub); isub = isub + 1;
  iSp = [3 5];
  UT_sep = sum(pic_comb.UT(iSp),2);
  UT_tog = Uti_cold;
  UK_sep = sum(pic_comb.UK(iSp),2);
  UK_tog = Uki_cold;
  plot(hca,pic_comb.twci,UT_sep,'kx',pic_comb.twci,UT_tog,'rx',...
           pic_comb.twci,UK_sep,'ko',pic_comb.twci,UK_tog,'ro')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_T+U_K (cold ions)';
  legend(hca,{'U_{T,sep}','U_{T,comb}','U_{K,sep}','U_{K,comb}'},...
    'location','best')
end

if 0 % UT + UK cold ions comparison, separate/together
  hca = h(isub); isub = isub + 1;
  iSp = [3 5];
  Utot_sep = sum(pic_comb.UT(iSp),2) + sum(pic_comb.UK(iSp),2);
  Utot_tog = UK_cold_ions + UT_cold_ions;c
  plot(hca,pic_comb.twci,Utot_sep,'x',pic_comb.twci,Utot_tog,'o')
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'U_T+U_K (cold ions)';
  legend(hca,{'separate','together'},...
    'location','best')
end

compact_panels(0.04)
hlinks = linkprop(h,{'XLim'});
h(1).XLim = [0 125];
c_eval('h(?).FontSize = 12;',1:numel(h))
c_eval('h(?).Position(3) = 0.7;',1:numel(h))
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
hl = findobj(h,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

%% Edot J