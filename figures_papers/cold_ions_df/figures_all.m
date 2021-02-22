%% Figure 1, change of inflow parameters, compression of current sheet
%% Figure 1, prepare data
pic = no02m;

x0 = no02m.xi(end)/2;
twpe1 = 1000; twpe1 = no02m.twpelim(twpe1).twpe; twci1 = no02m.twpelim(twpe1).twci;
twpe2 = 24000; twpe2 = no02m.twpelim(twpe2).twpe; twci2 = no02m.twpelim(twpe2).twci;

%pic_Bxline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'Babs');
pic_Bxline_z1_ = pic.get_points(pic.x_xline,pic.z_xline+1,pic.twci,[-0.1 0.1],'Babs');
%pic_nxline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'ni');
pic_nxline_z1_ = pic.get_points(pic.x_xline,pic.z_xline+1,pic.twci,[-0.1 0.1],'ni');
%pic_nhot_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'n(1)');
pic_nhot_z1_ = pic.get_points(pic.x_xline,pic.z_xline+1,pic.twci,[-0.1 0.1],'n(1)');
%pic_ncold_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'n(3)');
pic_ncold_z1_ = pic.get_points(pic.x_xline,pic.z_xline+1,pic.twci,[-0.1 0.1],'n(3)');
%pic_tixline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'ti');
%pic_vA_z1 = squeeze(pic_Bxline_z1./sqrt(pic_nxline_z1));
pic_vA_z1_ = squeeze(pic_Bxline_z1_./sqrt(pic_nxline_z1_));



zlim = [-0.5 0.5];
xlim = x0 + 40*[-1 1];
pic = no02m.zlim(zlim).xlim(xlim);
pic_Bz_tx = squeeze(mean(pic.Bz,2));
pic_A_tx = squeeze(mean(pic.A,2));
%% Figure 1, plot
colors = pic_colors('matlab');
Alev = -25:1:25;
doA = 0;
istepA = 3;

nrows = 5;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
hb = gobjects(0);


if 1 % zcut at early times
  hca = h(isub); isub = isub + 1;
  pic = no02m.twpelim(twpe1).xlim(x0+[-1 1]);
  hn1 = plot(hca,pic.zi,mean(pic.n(1),1),'color',colors(2,:));
  hold(hca,'on')
  hn2 = plot(hca,pic.zi,smooth(mean(pic.n([3 5]),1),5),'color',colors(1,:));
  hB = plot(hca,pic.zi,mean(pic.Babs,1),'k');
  hold(hca,'off')
  hca.XLim = [-10 10];
  if 0 % Add axes for A
    %hold(hca,'on')
    ax1 = hca;
    ax1_pos = ax1.Position; % position of first axes
    ax2 = axes('Position',ax1_pos,...
      'XAxisLocation','top',...
      'YAxisLocation','right',...
      'Color','none');
    ax2all(isub-1) = ax2;
    hA = line(ax2,pic.zi,mean(pic.A,1),'color',[0 0 0],'linestyle','--');
    ax2.XLim = hca.XLim;
    ax2.XTick = [];
    ax2.YLabel.String = 'A (B_0d_i)';
  end
  hca.XLabel.String = 'z (d_i)';
  hca.YLabel.String = 'n, B';  
  %legend([hn1,hn2,hB,hA],{'n_{hot}','n_{cold}','|B|','A_y'},'location','northwest','Box','off')
  %pos = hca.Position;
  legend([hn1,hn2,hB],{'n_{hot}','n_{cold}','|B|'},'location','east','Box','off')
  %legend([hn1,hn2,hB],{'n_{hot}','n_{cold}','|B|'},'location','northoutside','Box','off','Orientation','horizontal')
  %hca.Position = pos;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  irf_legend(hca,{['inflow at t\omega_{ci}= ' num2str(twci1)]},[0.02 0.9],'fontsize',14,'color',[0 0 0])
end

if 0 % zcut at, vs A 
  hca = h(isub); isub = isub + 1;
  pic = no02m.twpelim(twpe1).xlim(x0+[-1 1]);
  A = mean(pic.A,1);
  plot(hca,A,mean(pic.n(1),1),'color',colors(2,:))
  hold(hca,'on')
  plot(hca,A,smooth(mean(pic.n([3 5]),1),5),'color',colors(1,:))
  plot(hca,A,mean(pic.Babs,1),'k')
  hold(hca,'off')
  %hca.XLim = [-10 10];  
  hca.XLim = [min(A) max(A)];
  hca.XLabel.String = 'A (B_0d_i)';
  hca.YLabel.String = 'n, B';
  %legend(hca,{'n_{hot}','n_{cold}','|B|'},'location','northwest','box','off')
end
if 0 % xcut vs A
  hca = h(isub); isub = isub + 1;
  pic = no02m.twpelim(twpe2).xgrid(2:no02m.nx-1).zlim(0+[-1 1]);
  A = mean(pic.A,2);
  nh = mean(pic.n(1),2);
  plot(hca,A,nh,'color',colors(2,:))
  hold(hca,'on')
  plot(hca,A,mean(pic.n([3 5]),2),'color',colors(1,:))
  plot(hca,A,mean(pic.Babs,2),'k')
  hold(hca,'off')  
  %hca.XLim = x0 + [-50 50];
  hca.XLim = [min(A) max(A)];
  hca.XLabel.String = 'A (B_0d_i)';
  hca.YLabel.String = 'n, B';
  %legend(hca,{'n_{hot}','n_{cold}','|B|'},'location','northwest','box','off')
end
if 0 % ni, Bx, vA, RE
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic_nxline_z1,pic.twci,pic_Bxline_z1,pic.twci,pic_vA_z1,pic.twci,pic.RE);
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'n, B_x';
  %hca.XLim = [0 0.4];
  legend(hca,{'n(x_x,z_x+1)','B_x(x_x,z_x+1)'},'location','southwest','box','off')
end
if 0 % ni, Bx, vA
  hca = h(isub); isub = isub + 1;  
  [AX,H1,H2] = plotyy(hca,pic.twci,[pic_nxline_z1,pic_Bxline_z1]',pic.twci,smooth(pic_vA_z1,1));  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'n, B_x';
  AX(2).YLabel.String = 'v_A';
  AX(2).YColor = [0 0 0];
  H2.LineStyle = '-.';
  H2.Color = [0 0 0];
  %hca.XLim = [0 0.4];
  legend(hca,{'n(x_X,z_X+1)','B_x(x_X,z_X+1)','v_A(x_X,z_X+1)'},'location','northeast','box','off')
end
if 0 % ni, tot/hot/cold, Bx
  hca = h(isub); isub = isub + 1;  
  H1 = plot(hca,pic.twci,[pic_nxline_z1,pic_nhot_z1,pic_ncold_z1,pic_Bxline_z1]');  
  H1(1).Color = colors(4,:);
  H1(2).Color = colors(2,:);
  H1(3).Color = colors(1,:);
  H1(4).Color = colors(5,:);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'n, |B|';
  AX(2).YLabel.String = 'v_A';
  AX(2).YColor = [0 0 0];
  H2.LineStyle = '-.';
  H2.Color = [0 0 0];
  %hca.XLim = [0 0.4];
  legend(hca,{'n(x_X,z_X+1)','n_{hot}(x_X,z_X+1)','n_{cold}(x_X,z_X+1)','B_x(x_X,z_X+1)'},'location','northeast','box','off')
end
pic = no02m.zlim(zlim).xlim(xlim);
if 1 % ni, tot/hot/cold, Bx, vA
  hca = h(isub); isub = isub + 1; 
  
  [AX,H1,H2] = plotyy(hca,pic.twci,[pic_nxline_z1_,pic_nhot_z1_,pic_ncold_z1_,pic_Bxline_z1_]',pic.twci,smooth(pic_vA_z1_,1));  
  [AX,H1,H2] = plotyy(hca,pic.twci,[pic_nxline_z1,pic_nhot_z1,pic_ncold_z1,pic_Bxline_z1]',pic.twci,smooth(pic_vA_z1_,1));  
  H1(1).Color = colors(4,:);
  H1(2).Color = colors(2,:);
  H1(3).Color = colors(1,:);
  H1(4).Color = [0 0 0];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'n, |B|';
  AX(2).YLabel.String = 'v_A';
  AX(2).YColor = [0 0 0];
  H2.LineStyle = '-.';
  H2.Color = colors(5,:);
  AX(1).YLim = [0 0.75];
  AX(1).YTick = [0 0.25 0.5 0.75];
  AX(2).YLim = [0 1.5];
  AX(2).YTick = [0 0.25 0.5 0.75]*2;
  AX(1).XLim = [5 125];
  AX(2).XLim = [5 125];
  %hca.XLim = [0 0.4];
  %legend(hca,{'n(x_X,z_X+1)','B_x(x_X,z_X+1)','v_A(x_X,z_X+1)'},'location','northeast','box','off')
  %hleg = legend(hca,{'n(x_X,z_X+1)','n_{hot}(x_X,z_X+1)','n_{cold}(x_X,z_X+1)','B_x(x_X,z_X+1)','v_A(x_X,z_X+1)'},'location','northwest','box','off','orientation','horizontal');
  %hleg = legend(hca,{'n','n_{hot}','n_{cold}','|B|','v_A'},'location','northwest','box','off','orientation','horizontal');
  hleg = legend([H1(1) H2],{'n','v_A'},'location','southwest','box','off');
  irf_legend(hca,{'inflow values 1d_i above X line'},[0.02 0.98],'fontsize',14,'color',[0 0 0])
end
if 0 % ni, Bx
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic_nxline_z1,pic.twci,pic_Bxline_z1);
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'n, B_x';
  %hca.XLim = [0 0.4];
  legend(hca,{'n(x_x,z_x+1)','B_x(x_x,z_x+1)'},'location','southwest','box','off')
end
if 0 % vA
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic_vA_z1);
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'v_x';
  %hca.XLim = [0 0.4];
  legend(hca,{'v_A(x_x,z_x+1)'},'location','southwest','box','off')
end
if 1 % RE
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic.RE);
  h_.Color = colors(3,:);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'E_R';
  %hca.XLim = [0 0.4];
  %legend(hca,{'v_A(x_x,z_x+1)'},'location','southwest','box','off')
  hca.XLim = [5 125];
  hca.YLim = [0 0.17];
end
if 1 % Bz(x,t)
  hca = h(isub); isub = isub + 1;  
  pic = no02m.zlim(zlim).xlim(xlim);
  h_ = pcolor(hca,pic.twci,pic.xi,pic_Bz_tx);
  shading(hca,'flat')
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x (d_i)';
  colormap(hca,pic_colors('blue_red'))
  clim = hca.CLim;
  pos = hca.Position;
  hcb = colorbar('peer',hca);
  hcb.Title.String = 'B_z';
  if 1 % A
    hold(hca,'on')
    contour(hca,pic.twci,pic.xi,pic_A_tx,0:0.5:25,'color',[0 0 0]);
    plot(hca,pic.twci,pic.x_xline,'linewidth',2,'color',[0 0 0])
    hold(hca,'off')
  end
  colormap(hca,pic_colors('blue_red'))    
  hca.CLim = clim;
  hca.Position = pos;
  %legend(hca,{'v_A(x_x,z_x+1)'},'location','southwest','box','off')
  hca.XLim = [5 125];
  %hca.YLim = [0 0.17];
end
if 1 % xcut at late times
  hca = h(isub); isub = isub + 1;
  pic = no02m.twpelim(twpe2).xgrid(2:no02m.nx-1).zlim(0+[-1 1]);
  hn1 = plot(hca,pic.xi,mean(pic.n(1),2),'color',colors(2,:));
  hold(hca,'on')
  hn2 = plot(hca,pic.xi,mean(pic.n([3 5]),2),'color',colors(1,:));
  hB = plot(hca,pic.xi,mean(pic.Babs,2),'k');
  hold(hca,'off')  
  hca.XLim = x0 + [-50 50];
  if 0 % Add axes for A
    %hold(hca,'on')
    ax1 = hca;
    ax1_pos = ax1.Position; % position of first axes
    ax2 = axes('Position',ax1_pos,...
      'XAxisLocation','top',...
      'YAxisLocation','right',...
      'Color','none');
    ax2all(isub-1) = ax2;
    hA = line(ax2,pic.xi,mean(pic.A,2),'color',[0 0 0],'linestyle','--');
    ax2.XLim = hca.XLim;
    ax2.XTick = [];
    ax2.YLabel.String = 'A (B_0d_i)';
    legend([hn1,hn2,hB],{'n_{hot}','n_{cold}','|B|'},'location','west','Box','off')
  end
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'n, B';
  %legend([hn1,hn2,hB,hA],{'n_{hot}','n_{cold}','|B|','A_y'},'location','north','Box','off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  irf_legend(hca,{['outflow at t\omega_{ci}= ' num2str(twci2)]},[0.20 0.9],'fontsize',14,'color',[0 0 0])
end

legends = {'a)','b)','c)','d)','e)','f)'};
for ip = 1:numel(h)
  irf_legend(h(ip),legends{ip},[-0.1 1.00],'fontsize',14)  
end

hall = findobj(gcf,'type','axes'); hall = hall(end:-1:1);
for ip = 1:numel(hall)
  hall(ip).FontSize = 12;
  hall(ip).XGrid = 'on';
  hall(ip).YGrid = 'on';
end

compact_panels(h([2 3 4]))

% Set up panels
if 1
  %%
  h(1).Position([2 4]) = [0.83 0.15];
  h(2).Position([2 4]) = [0.61 0.15];
  h(3).Position([2 4]) = [0.505 0.10];
  h(4).Position([2 4]) = [0.3 0.20];
  h(5).Position([2 4]) = [0.08 0.15];
else
  h(2).Position(4) = h(1).Position(4);
  h(2).Position(2) = h(2).Position(2)+0.015;
  h(1).Position(2) = h(1).Position(2)+0.030;
end

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

%% Figure 2, overview at t=120
% what do we want to show here?
varstrs = {'n([3 5])','n(3)';'t([3 5])','t(3)'}';
%varstrs = {'n([3 5])','n(3)','t([3 5])','t(3)'}';
clims = {[0 0.5],[0 0.5],[0 0.5],[0 0.5]};

xlim = [40 165];
zlim = 0.99*[-10 10];
cmapth = pic_colors('thermal');
cmaps = {cmapth,cmapth,cmapth,cmapth};
cbarlabels = {'n_{i,cold}','n_{i,cold}^{top}','T_{i,cold}','T_{i,cold}^{top}'};
pic = no02m.twpelim(24000).xlim(xlim).zlim(zlim);
h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps,'cbarlabels',cbarlabels);

for ip = 1:numel(h)
  h(ip).FontSize = 12;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end
hl = findobj(gcf,'Type','Contour');
c_eval('hl(?).Color = 0.5 + [0 0 0];',1:numel(hl))

legends = {'a)','b)','c)','d)','e)','f)'};
if 1 % 2x2
  %%
  compact_panels(h,0.002,0.002)
  for ip = 1:numel(h)
    irf_legend(h(ip),{[legends{ip} ' ' cbarlabels{ip}]},[0.02 0.98],'color',[1 1 1],'fontweight','bold','fontsize',12)
  end
  hcbar = findobj(gcf,'Type','ColorBar');
  delete(hcbar(2:end))
  hcbar(1).Position = [0.85 0.11 0.02 0.815];
  hcbar(1).YLabel.String = 'n, T';
  h(3).YTickLabels = [];
  h(3).YLabel.String = [];
  h(4).YTickLabels = [];
  h(4).YLabel.String = [];  
  for ip = 1:numel(h)
    h(ip).Position(2) = h(ip).Position(2) + 0.05;
  end
  hcbar(1).Position(2) = hcbar(1).Position(2) + 0.05;
  hcbar(1).FontSize = 12;
  h(1).Title.String = '';
end

%% Figure 3, reduced distributions
%ds100 = PICDist('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');
%% Figure 3, prepare data
twpe = 24000; xlim = [50 155]; zlim = [-15 15];
sep = no02m.twpelim(twpe).separatrix_location;
for zpick = [0 2 4]
  ds = ds100.twpelim(twpe).zfind(zpick).xlim(xlim).findtag({'line horizontal'});
  
  xdist = (ds.xi1{1}+ds.xi2{1})/2;
  zdist = (ds.zi1{1}+ds.zi2{1})/2;
  tdist = repmat(twpe,size(xdist));
  vExBx_tmp = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,[-0.25 0.25],'vExBx'); eval(sprintf('vExBx_z%g = vExBx_tmp;',zpick))
  vExBy_tmp = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,[-0.25 0.25],'vExBy'); eval(sprintf('vExBy_z%g = vExBy_tmp;',zpick))
  vExBz_tmp = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,[-0.25 0.25],'vExBz'); eval(sprintf('vExBz_z%g = vExBz_tmp;',zpick))
  pic_lim = no02m.xlim(xlim).zlim(zlim).twpelim(twpe);
  pic = no02m.twpelim(twpe);
  Bx_ = pic.Bx;
  By_ = pic.By;
  Bz_ = pic.Bz;
  Bx = interpfield(pic.xi,pic.zi,Bx_,xdist,zdist); 
  By = interpfield(pic.xi,pic.zi,By_,xdist,zdist); 
  Bz = interpfield(pic.xi,pic.zi,Bz_,xdist,zdist); 
  %fred5_tmp = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred5_z%g = fred5_tmp;',zpick))
  %fred3_tmp = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred3_z%g = fred3_tmp;',zpick))
  %fred35_tmp = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred35_z%g = fred35_tmp;',zpick))
  %fred46_tmp = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred46_z%g = fred46_tmp;',zpick))    
end
%% Figure 3, plot
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 3;
ncols = 3;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 0; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.0;
doPhi = 0; colorPhi = [0.5 0.5 0];
doSep = 1;

cmap_dist = pic_colors('waterfall');

freds = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
%freds = {fred35_z2,fred35_z2,fred35_z2};
%freds = {fred5_z4,fred5_z4,fred5_z4;fred5_z2,fred5_z2,fred5_z2;fred5_z0,fred5_z0,fred5_z0}';
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)'}; 

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.fv' labstr]);
    pcolor(hca,fred.x,fred.v,log10(fredplot)')
    shading(hca,'flat')
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy4')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %irf_legend(hca,{sprintf('%s z = %g',legends{ifred},unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('%s f(v_%s,z=%g)',legends{ifred},labstr,unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 12;
    if 0*doE
      hold(hca,'on')
      plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
      hold(hca,'off')
    end
    if 0*doV
      hold(hca,'on')
      plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
      hold(hca,'off')
    end
    if doExB
      hold(hca,'on')
      xx = eval(['x_z' num2str(unique(fred.z))]);
      vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
      hExB = plot(hca,xx,vv,'color',colorExB,'linewidth',0.5);
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end    
    if doSep
      hold(hca,'on')   
      if 1           
        zz = unique(fred.z);
        [PKS,LOCS] = findpeaks(-abs(sep.z-zz),'sort','descend');
        [~,iz] = min(abs(sep.z-zz));
        if zz == 0
          xx = sep.x(LOCS(1));
          hSep = plot(hca,xx(1)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
        else
          xx = sep.x(LOCS(1:2));
          hSep = plot(hca,xx(1)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
          hSep = plot(hca,xx(2)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');      
        end
      elseif 0
        xx = unique(fred.x);
        [~,ix] = min(abs(sep.x-xx));
        zz = sep.z(ix);
      end      
      hold(hca,'off')      
    end
  end
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
hl = findobj(h(1),'type','line');
legend(hl([3 2]),{'v_{ExB}','separatrix'},'box','off','location','north')
%legend([hExB,hSep],{'v_{ExB}','separatrix'})
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
hl = findobj(h(8),'type','line'); delete(hl(2));
hl = findobj(h(9),'type','line'); delete(hl(2));

hcb = colorbar('peer',h(9));
hcb.Position = [0.91 0.13 0.015 0.795];
hcb.YLabel.String = 'log_{10} f(x,v_{x,y,z})';
c_eval('h(?).YLabel.String = []; h(?).YTickLabel = [];',[2 3 5 6 8 9])
c_eval('h(?).YLabel.String = ''v'';',[1 4 7])
compact_panels(h,0.002,0.002)
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
h(1).CLim = 0.99*[-4 2];
h(1).YLim = 0.99*4*[-1 1];

%% Figure 4, reduced vertical distribution 
%ds100 = PICDist('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');
%% Figure 4, prepare data
for xpick = 75:5:95  
  ds = ds100.twpelim(twpe).xfind(xpick).findtag({'line vertical'});
  xdist = (ds.xi1{1}+ds.xi2{1})/2;
  zdist = (ds.zi1{1}+ds.zi2{1})/2;
  
  zlim_fred = [min(zdist) max(zdist)];
  xlim_fred = xdist(1) + 0.25*[-1 1];
  eval(sprintf('z_x%g = pic.xlim(xlim_fred).zlim(zlim_fred).zi;',xpick))
  %fred3_tmp = ds.reduce_1d_new('x',[3],[]); eval(sprintf('fred3_x%g = fred3_tmp;',xpick))
  %fred5_tmp = ds.reduce_1d_new('x',[5],[]); eval(sprintf('fred5_x%g = fred5_tmp;',xpick))
  fred35_tmp = ds.reduce_1d_new('x',[3 5],[]); eval(sprintf('fred35_x%g = fred35_tmp;',xpick))

  
  eval(sprintf('vExBx_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBx,1));',xpick))
  eval(sprintf('vExBy_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBy,1));',xpick))
  eval(sprintf('vExBz_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vExBz,1));',xpick))
%   eval(sprintf('vix_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vx([3 5]),1));',xpick))
%   eval(sprintf('viy_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vy([3 5]),1));',xpick))
%   eval(sprintf('viz_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vz([3 5]),1));',xpick))
%   eval(sprintf('vex_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vx([4 6]),1));',xpick))
%   eval(sprintf('vey_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vy([4 6]),1));',xpick))
%   eval(sprintf('vez_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).vz([4 6]),1));',xpick))
%   eval(sprintf('Ex_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Ex,1));',xpick))
%   eval(sprintf('Ey_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Ey,1));',xpick))
  eval(sprintf('Ez_x%g = squeeze(mean(pic.xlim(xlim_fred).zlim(zlim_fred).Ez,1));',xpick))
  eval(sprintf('phiz_x%g = -cumtrapz(z_x%g,Ez_x%g);',xpick,xpick,xpick))
end
%sep = no02m.twpelim(twpe).separatrix_location;
disp('Done.')
%% Figure 4, plot
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 3;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1; 
doE = 0; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doVe = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];
doSep = 1;
% if doSep
%   xline_pos = no02m.twpelim(twpe).xline_position;
%   sep = no02m.twpelim(twpe).separatrix_location;
% end
hleg = gobjects(0);

cmap_dist = pic_colors('waterfall');

%freds = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
%freds = {fred35_x75,fred35_x80,fred35_x85,fred35_x90,fred35_x95};
freds = {fred35_x75,fred35_x85,fred35_x95,fred35_x75,fred35_x85,fred35_x95,fred35_x75,fred35_x85,fred35_x95};
freds = {fred35_x75,fred35_x80,fred35_x85,fred35_x90,fred35_x95,...
         fred35_x75,fred35_x80,fred35_x85,fred35_x90,fred35_x95,...
         fred35_x75,fred35_x80,fred35_x85,fred35_x90,fred35_x95};
freds = {fred35_x75,fred35_x85,...
         fred35_x75,fred35_x85,...
         fred35_x75,fred35_x85};
xpicks = [75 85];
%freds = {fred3_x75,fred3_x85,fred3_x95};
%freds = {fred3_x85};
%freds = {fred5_x75,fred5_x85,fred5_x95};
%freds = {fred46_x75,fred46_x85,fred46_x95};
%freds = {fred46_x75,fred46_x85,fred46_x95,fred46_x75,fred46_x85,fred46_x95,fred46_x75,fred46_x85,fred46_x95};
%freds = {fred46_x75,fred46_x80,fred46_x85,fred46_x90,fred46_x95};
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','x','x','x','x','x','x','x','x'};
labstrs = {'y','y','y','y','y','y','y','y','y'};
labstrs = {'z','z','z','z','z','z','z','z','z'};
labstrs = {'z','z','z','y','y','y','x','x','x'};
labstrs = {'z','z','z','z','z','y','y','y','y','y','x','x','x','x','x'};
labstrs = {'z','z','y','y','x','x'};
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)'};

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.fv' labstr]);
    pcolor(hca,fred.v,fred.z,log10(fredplot))
    shading(hca,'flat')
    hca.YLabel.String = 'z (d_i)';
    hca.XLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy4')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %irf_legend(hca,{sprintf('%s x = %g',legends{isub-1},unique(fred.x))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('%s',legends{isub-1})},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 12;
    if 1*doE
      hold(hca,'on')
      %plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['E' labstr '_x' num2str(unique(fred.x))]);
      hE = plot(hca,smooth(vv,50),zz,'color',0*colorE,'linewidth',1,'linestyle',':');      
      hold(hca,'off')
    end
    if doSep
      hold(hca,'on')      
      xx = unique(fred.x);
      [~,ix] = min(abs(sep.x-xx));
      zz = sep.z(ix);
      hSep = plot(hca,hca.XLim,zz*[1 1],'color',0*colorExB,'linewidth',1,'linestyle',':');      
      hold(hca,'off')      
    end
    if doV
      hold(hca,'on')
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['vi' labstr '_x' num2str(unique(fred.x))]);      
      hv = plot(hca,vv,zz,'color',colorV,'linewidth',1.5);
      hold(hca,'off')
      hleg()
    end
    if doVe
      hold(hca,'on')
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['ve' labstr '_x' num2str(unique(fred.x))]);      
      hv = plot(hca,vv,zz,'color',colorV,'linewidth',1.5);
      hold(hca,'off')
      hleg()
    end
    if doExB
      hold(hca,'on')
      unique(fred.x);
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['vExB' labstr '_x' num2str(unique(fred.x))]);
      hExB = plot(hca,smooth(vv,50),zz,'color',0*colorExB,'linewidth',1,'linestyle','-');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end  
    if doPhi && ifred == 2
      hold(hca,'on')
      unique(fred.x);
      %if 1
      try
      zz = eval(['z_x' num2str(unique(fred.x))]);
      pp = eval(['phi' labstr '_x' num2str(unique(fred.x))]);
      vv0 = eval(['vExB' labstr '_x' num2str(unique(fred.x))]);
      vv0 = vv0(end);
      vv_phi = sqrt(2*(abs(pp-pp(end-150)))).*sign(pp-pp(end-150));
      %(smooth(vv,10)+vv0-vv(end));
      hphi = plot(hca,vv_phi,zz,'color',0*[1 1 1],'linewidth',1,'linestyle','--');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      end
      hold(hca,'off')
    end   
  end
end
drawnow
compact_panels(h(1:end),0.002,0.002)

%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
% colorbars
if 1 % 3 colorbars
  ip = ncols;
  position_hcb_peer = h(ip).Position;
  hcb = colorbar('peer',h(ip));
  h(ip).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})','z');
  hcb.YLabel.FontSize = 14;
  hcb.Position(3)=0.02;
  
  ip = 2*ncols;
  position_hcb_peer = h(ip).Position;
  hcb = colorbar('peer',h(ip));
  h(ip).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})','y');
  hcb.YLabel.FontSize = 14;
  hcb.Position(3)=0.02;
  
  ip = 3*ncols;
  position_hcb_peer = h(ip).Position;
  hcb = colorbar('peer',h(ip));
  h(ip).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})','x');
  hcb.YLabel.FontSize = 14;
  hcb.Position(3)=0.02;
  
else % 1 colorbar
  position_hcb_peer = h(end).Position;
  hcb = colorbar('peer',h(end));
  h(end).Position = position_hcb_peer;
  hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})',labstr);
  hcb.YLabel.FontSize = 14;
end
hl = findobj(h(ncols),'type','line');
hleg = legend(hl,{'v_{ExB}','separatrix'},'edgecolor',[1 1 1],'location','northeast');
hleg = legend(hl([2 1 3]),{'v_{ExB}','\phi_z','separatrix'},'edgecolor',[1 1 1],'location','northeast');
% if all([doV doExB])
%   legend([hv,hExB],{sprintf('v_{i,%s}',labstr),sprintf('v_{ExB,%s}',labstr)})
% elseif doV
%   legend([hv],{sprintf('v_{i,%s}',labstr)})
% elseif doExB
%   legend([hExB],{sprintf('v_{ExB,%s}',labstr)})
% end
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
%compact_panels(0.005,0.01)
h(1).CLim = 0.99*[-4 2];
h(1).CLim = [-3.5 1.8];
%h(1).CLim = 0.99*[-6 0]; % electrons
h(1).XLim = 0.99*3*[-1 1];
for ip = [2:ncols (ncols+2):(2*ncols) (2*ncols+2):(3*ncols)]
  h(ip).YTickLabels = '';
  h(ip).YLabel.String = '';  
end
for ip = 1:npanels
  h(ip).FontSize = 12;
  %h(ip).Position(2) = 0.17;
end
for ip = [2 4 6]
  h(ip).Position(1) = 0.5;
end
%xpicks = xpicks;
for ip = 1:numel(xpicks)
  h(ip).Title.String = sprintf('x = %g',xpicks(ip));
  h(ip).Title.FontWeight = 'normal';
end
for ip = ((nrows-1)*ncols+1):nrows*ncols
  h(ip).XLabel.String = 'v';  
end
for ip = [2 4 6]
  h(ip).Position(1) = 0.4;
end
compact_panels(h,0.002,0.002)
hleg.FontSize = 11;
hcb_ = findobj(gcf,'Type','ColorBar');

%% Figure 5, test particle trajectories

%% Figure 6, magnetic curvature
%% Curvature plot, combined, also including f(vx,vz) at z=2 to illustrate inward vs outward beam temperature
clear h;
nrows = 4;
ncols = 5;
h(1) = subplot(nrows,ncols,1:ncols);
h1pos = h(1).Position;
h(2) = subplot(nrows,ncols,ncols+(1:ncols));
h2pos = h(2).Position;
c_eval('h(?+2) = subplot(nrows,ncols,?+2*ncols);',1:2*ncols);
twpe = 24000;
xpos = [72:2:80];
fontsize = 14;
legfontsize = 11;
xlim = [60 110];
xpos2 = [69 71 73 75 77 79];
xpos2 = [71 73 75 77 79];
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','k)','l)','m)','n)','o)'};

ds = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([71 80]).zfind(0).xfind(xpos);
ds2 = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([71 80]).zfind(2).xfind(xpos);
ds2 = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([50 80]).zfind(2).xfind(xpos2);
pic = no02m.twpelim(twpe).xlim(xlim).zlim([-7 7]);

hca = h(1);
hmap = pic.plot_map(hca,{'log10(curvbrad)'},'A',1,'cmap',pic_colors('blue_red'),'cbarlabels',{'log_{10}r_B'},'smooth',2);
hmap.CLim = [-1 3];
hca.Position = h1pos;
hold(hca,'on')
ds.plot_boxes(hca,'color',[1 1 1]);
ds2.plot_boxes(hca,'color',[0 0 0]);
hold(hca,'off')

hca = h(2);
comp = 'x';
hline = pic.xlim([63 xlim(2)]).zlim([-0.25 0.25]).plot_line(hca,comp,{{'Babs','curvbrad','curvbrad.*Babs'}},'smooth',10);
hca.Position = h2pos;
hca.YLim = [0 3.99];
%legend(hca,{'|B|','r_B','r_B\omega_{ci}'},'location','north','orientation','horizontal')
legend(hca,{'|B|','r_B','r_B\omega_{ci}'},'location','north','orientation','horizontal','box','off')
hca.Title.String = [];
hca.XLim = xlim;
hca.Position(2) = hca.Position(2) + 0.08;
compact_panels(h(1:2),0)
c_eval('h(?).Position(2) = h(?).Position(2) + 0.03;',1:2)
irf_legend(hca,{'z=0\pm0.25'},[0.98 0.98],'color',[0 0 0],'fontsize',fontsize)

ih0 = 2+ncols;
if 1 % f(v_x,v_y)
  sumdim = 3;
  clim = [-4 1];
  xlim = 3.5*0.99*[-1 1];
  ylim = 3.5*0.99*[-1 1];
  x0 = (ds.xi1{1}+ds.xi2{1})/2;
  hds = ds.plot_map(h(ih0+(1:ncols)),[3 5],sumdim,'v',no02m,'log','curv',{no02m,1},'nolabel'); % 
  hlinks = linkprop(hds.ax,{'XLim','YLim','CLim','XTick','YTick'});
  c_eval('hds.ax(?).Position(2) = hds.ax(?).Position(2) + 0.10;',1:ncols)
  hds.ax(1).Position(1) = h(1).Position(1);
  compact_panels(hds.ax,0.00,0.00)
  c_eval('irf_legend(hds.ax(?),sprintf(''%s x = %g, z = 0'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',legfontsize);',1:ncols);
  %[hax,hlab] = label_panels(hds.ax);
  hds.ax(1).CLim = clim;
  %h.ax(1).CLim = [0 1];
  hds.ax(1).XLim = xlim;
  hds.ax(1).YLim = ylim;
  %c_eval('hds.ax(?).XLabel.String = []; hds.ax(?).XTickLabels = [];',[1 2 3]);
  c_eval('hds.ax(?).YLabel.String = []; hds.ax(?).YTickLabels = [];',2:ncols);
  pos = hds.ax(end).Position;
  hcb = colorbar('peer',hds.ax(end));
  hcb.YLabel.String = 'f(v_x,v_y)';
  hds.ax(end).Position = pos;
end

ih0 = 2;
if 1 % f(v_x,v_z)
  hds2 = ds2.plot_map(h(ih0+(1:ncols)),[3 5],2,'v',no02m,'bline',no02m,'log','nolabel','nan','exb',no02m); % 
  hlinks2 = linkprop(hds2.ax,{'XLim','YLim','CLim','XTick','YTick'});
  x0 = (ds2.xi1{1}+ds2.xi2{1})/2;
  c_eval('irf_legend(hds2.ax(?),sprintf(''%s x = %g, z = 2'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',legfontsize);',1:ncols);
  hds2.ax(1).Position(1) = h(1).Position(1);
  %c_eval('hds2.ax(?).Position(2) = hds2.ax(?).Position(2) + 0.10;',1:3)
  compact_panels(hds2.ax,0.00,0.00)
  c_eval('hds2.ax(?).YLabel.String = []; hds2.ax(?).YTickLabels = [];',2:ncols);
  hds2.ax(1).CLim = clim;
  %h.ax(1).CLim = [0 1];
  hds2.ax(1).XLim = xlim;
  hds2.ax(1).YLim = ylim;
  pos = hds2.ax(end).Position;
  hcb = colorbar('peer',hds2.ax(end));
  hcb.YLabel.String = 'f(v_x,v_z)';
  hds2.ax(end).Position = pos;
end

irf_legend(hmap,{'a)'},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)
irf_legend(h(2),{'b)'},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)


hlinks = linkprop(h(3:end),{'XLim','YLim','CLim'});
h(1).Title.String = '';

c_eval('h(?).FontSize = fontsize;',1:numel(h));
c_eval('h(?).Position(1) = h(?).Position(1)-0.04;',1:numel(h));
c_eval('h(?).Position(4) = h(end).Position(4);',3:numel(h))
drawnow
c_eval('h(?).Position(2) = h(?).Position(2)+0.06;',3:(2+ncols));
%compact_panels(h(3:end),0.0,0)

%irf_legend(h(8),{'B'},[0.92 0.27],'color',0.5+[0 0 0],'fontsize',16)

%c_eval('h().Position')
hl = findobj(h(2+ncols),'type','line');
hs = findobj(h(2+ncols),'type','scatter');
legend([hl;hs],{'B','v_{bulk}','v_{ExB}'},'box','off','location','northoutside','orientation','horizontal')

%% Figure 7, f(vy), py
%% Reduced distributions, plot
twpe = 24000;
pic = no02m.twpelim(twpe);
fred = fred35_z0;
twpe = 24000;
fontsize = 12;
% What to include
% - overview of whxx  ere boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];
zlim_line = [-0.5 0.5];
xlim_line = [min(fred3_z0.x) max(fred3_z0.x)];

py0 = [8 7 6.5 5.8 5.2 4.6];
py0 = [8 6.9 6.3 5.6 5.1 4.6];
py0 = [8 6.9 6.3 5.6 5.1 4.6];
py0 = [8 7 6.5 5.8 5.4 4.6];
xmin = [70 70 71 73 80 82];
xmax = [80 85 87 93 95 100];


py0 = [8.5 8.1 7.7 7.4 6.9 6.3 5.8 5.2 4.6];
xmin = [71 71 71 71 71 71.5 73 80 83];
xmax = [74 77 79 81 84 87 92 96 97]+00;
npy = numel(py0);

nrows = 1;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 1; colorE = [0 0.8 0.8];
doV = 1; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];
doTraj = 1; colorTraj = [0 0 0];
trs = tr100.find([tr100.z0] == 0,[tr100.vy0] > 0.5,[tr100.x0] > 73,[tr100.xstart] < 100);
trs = tr100.find([tr100.z0] == 0,[tr100.vy0] > 0.5).lim('t',[23000 25000]/200);

cmap_dist = pic_colors('waterfall');

%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};

if 0 % A(x),
  hca = h(isub); isub = isub + 1; 
  pic_tmp = pic.xlim(xlim_line).zlim(zlim_line);
  plot(hca,pic_tmp.xi,mean(pic_tmp.A,2))  
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'A_{y}';    
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % fi(v_y)
  hca = h(isub); isub = isub + 1;
  fred = fred35_z0;    
  pcolor(hca,fred.x,fred.v,log10(fred.fvy)')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('candy4'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'f_{i,cold}(x,v_{y})';  
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  hca.CLim = 0.99*[-4 2];
  hca.CLim = 0.99*[-4 2];
  hca.YLim = 0.99*4*[-1 1];
  
  for ipy = 1:npy
    hold(hca,'on')
    pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
    A = squeeze(mean(pic_tmp.A,2));    
    xx = pic_tmp.xi;
    yy = f_vy_A(py0(ipy),A);
    plot(hca,xx,yy,':k','linewidth',0.5)
    hold(hca,'off')
    ht = text(hca,xx(end),yy(end),sprintf(' %.1f',py0(ipy)),'fontsize',fontsize);
    ht.HorizontalAlignment = 'left';
    dy = (yy(end)-yy(end-1))./diff(hca.YLim);
    dx = (xx(end)-xx(end-1))./diff(hca.XLim);
    ht.Rotation = atand(dy/dx);
    
  end
  if 0*doE
    hold(hca,'on')
    plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
    hold(hca,'off')
  end
  if 0*doV
    hold(hca,'on')
    plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if 0 % doExB
    hold(hca,'on')
    xx = eval(['x_z' num2str(unique(fred.z))]);
    vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
    plot(hca,xx,vv,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi(v_y)
  hca = h(isub); isub = isub + 1;
  fred = fred3_z0;    
  pcolor(hca,fred.x,fred.v,log10(fred.fvy)')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('candy4'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'f_{i,cold}(x,v_{y})';  
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  hca.CLim = 0.99*[-4 2];
  hca.CLim = 0.99*[-2 1];
  hca.YLim = 0.99*4*[-1 1];
  
%   for ipy = 1:npy
%     hold(hca,'on')
%     pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
%     A = squeeze(mean(pic_tmp.A,2));
%     plot(hca,pic_tmp.xi,f_vy_A(py0(ipy),A),'--k')
%     hold(hca,'off')
%   end

   for itr = 1:trs.ntr
      hold(hca,'on')
      plot(hca,trs(itr).x,trs(itr).vy,'k')
      hold(hca,'off')
   end
end
if 0 % fi_top(v_y)/fi_tot(v_y)
  hca = h(isub); isub = isub + 1;
  fred = fred3_z0;    
  pcolor(hca,fred.x,fred.v,fred3_z0.fvy'./fred35_z0.fvy')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('waterfall'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'f_{i,top}/f_{i,tot}(x,v_{y})';  
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  hca.CLim = 1*[0 1];  
  hca.YLim = 0.99*4*[-1 1];
  
%   for ipy = 1:npy
%     hold(hca,'on')
%     pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
%     A = squeeze(mean(pic_tmp.A,2));
%     plot(hca,pic_tmp.xi,f_vy_A(py0(ipy),A),'--k')
%     hold(hca,'off')
%   end

   for itr = 1:trs.ntr
      hold(hca,'on')
      plot(hca,trs(itr).x,trs(itr).vy,'k')
      hold(hca,'off')
   end
end
if 0 % A0(v_y,A)
  hca = h(isub); isub = isub + 1; 
  A_ = no02m.interp(fred.x,fred.z,no02m.twpelim(twpe).twci,'A');
  [VY,A] = meshgrid(fred.v,A_);  
  A0map = f_A0_vy_py0(A,VY);
  A0map(fred.fvy==0) = NaN;
  %[Ccont,hcont] = contourf(hca,fred.x,fred.v,A0map',-30:1:30);
  [Ccont,hcont] = contourf(hca,fred.x,fred.v,A0map',[0 sort(py0)]);
  clabel(Ccont,hcont,'LabelSpacing',72,'Color','k','FontWeight','bold');
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('waterfall'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'A_{y0}(A_{y,loc},v_y,v_{y0}=0)';  
  irf_legend(hca,{'A_{y0} = A_{y,loc} - v_{y,loc}';'assuming v_{y0}=0'},[0.98 0.1],'color',[0 0 0])
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  %hca.CLim = 0.99*[-4 2];
  %hca.CLim = 0.99*[-2 1];
  hca.YLim = 0.99*4*[-1 1];
  
  for ipy = 1:npy
    hold(hca,'on')
    pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
    A = squeeze(mean(pic_tmp.A,2));
    xx = pic_tmp.xi;
    yy = f_vy_A(py0(ipy),A);
    plot(hca,xx,yy,'--k')
    hold(hca,'off')
    ht = text(hca,xx(end),yy(end),sprintf('p_y = %.1f',py0(ipy)),'fontsize',fontsize);
  end
  if 0*doE
    hold(hca,'on')
    plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
    hold(hca,'off')
  end
  if 0*doV
    hold(hca,'on')
    plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if 0 % doExB
    hold(hca,'on')
    xx = eval(['x_z' num2str(unique(fred.z))]);
    vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
    plot(hca,xx,vv,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % ExB_x
  hca = h(isub); isub = isub + 1; 
  pic_tmp = pic.xlim(xlim_line).zlim(zlim_line);
  plot(hca,pic_tmp.xi,mean(pic_tmp.vExBx,2),pic_tmp.xi,mean(pic_tmp.vix,2))
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{x}';    
  legend(hca,{'v_{ExB}','v_i'},'location','best')
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 0 % A(x), for different times
  hca = h(isub); isub = isub + 1; 
  pic_tmp = no02m.twpelim(19000:1000:24000,'exact').xlim(xlim_line).zlim(zlim_line);
  plot(hca,pic_tmp.xi,squeeze(mean(pic_tmp.A,2)))  
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'A_{y}';    
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim'});
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
for ip = 1:nrows*ncols
  axwidth(ip) = h(ip).Position(3);
  h(ip).FontSize = 14;
end
for ip = 1:nrows*ncols
  h(ip).Position(3) = 0.7;%min(axwidth);
end
% for ip = 1:nrows*ncols
%   h(ip).Position(2) = h(ip).Position(2)-0.05;
% end

%c_eval('h(?).YTickLabel = []; h(?).YLabel = [];',[2 3 5 6 8 9])
%c_eval('h(?).YTick = -10:1:10;',1)


%% Figure 8, Energy content
nrows = 1;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

Utot = no02m.UB+no02m.Uke+no02m.Uki+no02m.Ute+no02m.Uti;
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,no02m.twci,no02m.UB,...
           no02m.twci,no02m.Uti,...
           no02m.twci,no02m.Ute,...
           no02m.twci,no02m.Uki,...
           no02m.twci,no02m.Uke,...
           no02m.twci,Utot)
legend(hca,{'U_B','U_{ti}','U_{te}','U_{ki}','U_{ke}','U_{tot}'},'location','best')
end

%% Figure 9, streaming instability
%% PICDist.reduce_1d, along a line, A tags
%ds = ds01.zfind(3);
%fred = ds.reduce_1d_new('x',[5],[]);
twpe = 10000; xlim = [130 110]; zlim = [-8 8];
twpe = 24000; xlim = [60 80]; zlim = [-8 8];

%ds = ds100.twpelim(twpe).findtag({'A=-6'});
%ds = ds100.twpelim(twpe).zfind(0).findtag({'line horizontal'});
%ds = ds100.twpelim(twpe).zfind(4).findtag({'line horizontal'});
ds = ds100.twpelim(twpe).findtag({'A=7.5'});

xdist = (ds.xi1{1}+ds.xi2{1})/2;
zdist = (ds.zi1{1}+ds.zi2{1})/2;
% arclength = [0 cumsum(sqrt(diff(xdist).^2 + diff(zdist).^2))];

pic_lim = no02m.xlim(xlim).zlim(zlim).twpelim(twpe);
pic = no02m.twpelim(twpe);
Bx_ = pic.Bx;
By_ = pic.By;
Bz_ = pic.Bz;

Bx = interpfield(pic.xi,pic.zi,Bx_,xdist,zdist); 
By = interpfield(pic.xi,pic.zi,By_,xdist,zdist); 
Bz = interpfield(pic.xi,pic.zi,Bz_,xdist,zdist);
%get_points(obj,x,z,t,range,field,varargin)

% reduced distributions
if 0 % saved reduced distributions
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_2.mat')  
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_3.mat')
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_4.mat')
  %save('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_4.mat','fred35_4')
else % make reduced distributions
  %fred35_A9 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A9 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_A8 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A8 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_A7 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A7 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  fred35_A75 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  fred3_A75 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  fred46_A75 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_A4 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A4 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_A5 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A5 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_A6 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A6 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_z4 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_z4 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  
end
%% Plot
fredi_str = '3'; iSpecies = [3];
frede_str = '3'; eSpecies = [3];
fredi = eval(['fred' fredi_str '_A75']);
frede = eval(['fred' frede_str '_A75']);
fred = fredi;
arclength = [0; cumsum(sqrt(diff(fredi.x).^2 + diff(fredi.z).^2))];
arclength = arclength(end:-1:1);
if 1; arclength = arclength - arclength(find(abs(fredi.z)==min(abs(fredi.z)))); end
darc = arclength(2)-arclength(1);
arcedges = [arclength(1)-0.5*darc; arclength+0.5*darc];
narc = 1800; % 900
arclength_interp = linspace(arclength(1),arclength(end),narc);
xinterp = interp1(arclength,xdist,arclength_interp);
zinterp = interp1(arclength,zdist,arclength_interp);
ni = interpfield(pic.xi,pic.zi,pic.ni,fredi.x,fredi.z); 
ni_ = interpfield(pic.xi,pic.zi,pic.ni,xinterp,zinterp); 
ne = interpfield(pic.xi,pic.zi,pic.ne,fredi.x,fredi.z); 
ne_ = interpfield(pic.xi,pic.zi,pic.ne,xinterp,zinterp); 
tic; vipar = interpfield(pic.xi,pic.zi,pic.vpar(iSpecies),fredi.x,fredi.z); toc
%tic; vipar_ = pic.get_points(xdist,zdist,twpe,[-0.1 0.1],'vpar([3])'); toc
vepar = interpfield(pic.xi,pic.zi,pic.vpar([3 5]),fredi.x,fredi.z); 
%vix = interpfield(pic.xi,pic.zi,pic.vx(iSpecies),fredi.x,fredi.z); 
%viy = interpfield(pic.xi,pic.zi,pic.vy(iSpecies),fredi.x,fredi.z); 
%viz = interpfield(pic.xi,pic.zi,pic.vz(iSpecies),fredi.x,fredi.z); 
%vex = interpfield(pic.xi,pic.zi,pic.vx(eSpecies),fredi.x,fredi.z); 
%vey = interpfield(pic.xi,pic.zi,pic.vy(eSpecies),fredi.x,fredi.z); 
%vez = interpfield(pic.xi,pic.zi,pic.vz(eSpecies),fredi.x,fredi.z); 
%vepar = interpfield(pic.xi,pic.zi,pic.vpar(eSpecies),fredi.x,fredi.z); 
Epar = interpfield(pic.xi,pic.zi,pic.Epar,fredi.x,fredi.z); 
Epar_ = interpfield(pic.xi,pic.zi,pic.Epar,xinterp,zinterp); 
%Eparx = interpfield(pic.xi,pic.zi,pic.Eparx,fredi.x,fredi.z); 
%Epary = interpfield(pic.xi,pic.zi,pic.Epary,fredi.x,fredi.z); 
%Eparz = interpfield(pic.xi,pic.zi,pic.Eparz,fredi.x,fredi.z); 
%Ex = interpfield(pic.xi,pic.zi,pic.Ex,fredi.x,fredi.z); 
%Ey = interpfield(pic.xi,pic.zi,pic.Ey,fredi.x,fredi.z); 
%Ez = interpfield(pic.xi,pic.zi,pic.Ez,fredi.x,fredi.z);
%vExBx = interpfield(pic.xi,pic.zi,pic.vExBx,fredi.x,fredi.z); 
%vExBy = interpfield(pic.xi,pic.zi,pic.vExBy,fredi.x,fredi.z); 
%vExBz = interpfield(pic.xi,pic.zi,pic.vExBz,fredi.x,fredi.z);
%vExBabs = sqrt(vExBx.^2 + vExBy.^2 + vExBz.^2);
%EExB = vExBabs.^2/2;
%Bx = interpfield(pic.xi,pic.zi,pic.Bx,fredi.x,fredi.z); 
%By = interpfield(pic.xi,pic.zi,pic.By,fredi.x,fredi.z); 
%Bz = interpfield(pic.xi,pic.zi,pic.Bz,fredi.x,fredi.z); 

fi_clim = [0 0.0499];
fe_clim = [0 1.3e-2];

nrows = 5;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;
doE = 0; colorE = [0.0 0.0 0.0];
doV = 0; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
isMap = [];

if 0 % line position on map, vy
  isMap(end+1) = isub; 
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic_lim.xi,pic_lim.zi,pic_lim.vy(iSpecies)');
  colormap(hca,pic_colors('blue_red'));
  hcb = colorbar('peer',hca);
  hca.CLim = max(max(get(findobj(hca.Children,'type','Image'),'CData')))*[-1 1];
  hcb.YLabel.String = sprintf('v_{y,%s}',fredi_str);
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  if 1 % plot_boxes
    hold(hca,'on')
    ds.plot_boxes(hca)
    hold(hca,'off')
  end
end
if 0 % line position on map, Ez
  isMap(end+1) = isub; 
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic_lim.xi,pic_lim.zi,pic_lim.Ez');
  colormap(hca,pic_colors('blue_red'));
  hcb = colorbar('peer',hca);
  hca.CLim = max(max(get(findobj(hca.Children,'type','Image'),'CData')))*[-1 1];
  hcb.YLabel.String = 'E_z';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  if 1 % plot_boxes
    hold(hca,'on')
    ds.plot_boxes(hca)
    hold(hca,'off')
  end
end
if 1 % line position on map, Epar
  isMap(end+1) = isub; 
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic_lim.xi,pic_lim.zi,smooth2(pic_lim.Epar,3)');
  colormap(hca,pic_colors('blue_red'));
  hcb = colorbar('peer',hca);
  hca.CLim = max(max(get(findobj(hca.Children,'type','Image'),'CData')))*[-1 1];
  hcb.YLabel.String = 'E_{||}';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  if 1 % plot_boxes
    hold(hca,'on')
    ds.plot_boxes(hca)
    hold(hca,'off')
  end
  if 1 % plot A
    hold(hca,'on')
    A = pic_lim.A;
    clim = hca.CLim;
    contour(hca,pic_lim.xi,pic_lim.zi,A',[0:1:25],'color',0.5*[1 1 1])
    hca.CLim = clim;
    hold(hca,'off')
  end
end
if 0 % line position on map, vepar
  isMap(end+1) = isub; 
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic_lim.xi,pic_lim.zi,abs(pic_lim.vepar)');
  colormap(hca,pic_colors('blue_red'));
  hcb = colorbar('peer',hca);
  hca.CLim = max(max(get(findobj(hca.Children,'type','Image'),'CData')))*[-1 1];
  hcb.YLabel.String = 'v_{e,||}';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  if 0 % plot_boxes
    hold(hca,'on')
    ds.plot_boxes(hca)
    hold(hca,'off')
  end
  if 1 % plot A
    hold(hca,'on')
    A = pic_lim.A;
    clim = hca.CLim;
    contour(hca,pic_lim.xi,pic_lim.zi,A',[0:1:25],'color',0.5*[1 1 1])
    contour(hca,pic_lim.xi,pic_lim.zi,A',[7.5 7.5],'color',0.0*[1 1 1])
    hca.CLim = clim;
    hold(hca,'off')
  end
end
if 0 % line position on map, n
  isMap(end+1) = isub; 
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic_lim.xi,pic_lim.zi,pic_lim.n(3)');
  colormap(hca,pic_colors('thermal'));
  hcb = colorbar('peer',hca);
  hca.CLim = max(max(get(findobj(hca.Children,'type','Image'),'CData')))*[-1 1];
  hcb.YLabel.String = 'n';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.CLim = [0 0.5];
  if 1 % plot_boxes
    hold(hca,'on')
    ds.plot_boxes(hca)
    hold(hca,'off')
  end
end
if 0 % line position
  hca = h(isub); isub = isub + 1;
  [ax,h1,h2] = plotyy(hca,arclength,fredi.x,arclength,fredi.z);
  hca.XLabel.String = 'arclength (d_i)';
  ax(1).YLabel.String = 'x';
  ax(2).YLabel.String = 'z';
  legend(hca,{'x','z'},'location','best')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % B
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,Bx,arclength,By,arclength,Bz)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'B';
  legend(hca,{'B_x','B_y','B_z'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % n
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,ni,arclength,ne)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'n';
  legend(hca,{'n_i','n_e'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % v
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,vix,arclength,viy,arclength,viz,arclength,vex,arclength,vey,arclength,vez)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'B';
  legend(hca,{'v_{ix}','v_{iy}','v_{iz}','v_{ex}','v_{ey}','v_{ez}'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % v
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,vExBx,arclength,vExBy,arclength,vExBz,arclength,vExBabs)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{ExB}';
  legend(hca,{'v_{x}','v_{y}','v_{z}','|v|'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % E
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,Ex,arclength,Ey,arclength,Ez)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'E';
  legend(hca,{'E_{x}','E_{y}','E_{z}'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Epar, int(Epar)dl
  hca = h(isub); isub = isub + 1;
  intE = -cumtrapz(arclength_interp,Epar_);
  plot(hca,arclength_interp,Epar_,arclength_interp,intE)  
  legend(hca,{'E_{||}','-\int E_{||}dl_{||}'},'location','eastoutside') 
  if 0
  hold(hca,'on')
  plot(hca,arclength,Eparx,arclength,Epary,arclength,Eparz)  
  hold(hca,'off')
  legend(hca,{'E_{||}','-\int E_{||}dl_{||}','E_{||,x}','E_{||,y}','E_{||,z}'},'location','eastoutside') 
  end
  hca.YLabel.String = 'E_{||}, \int E_{||}dl_{||}'; 
  
  hca.XLabel.String = 'arclength (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % Epar
  hca = h(isub); isub = isub + 1;  
  %plot(hca,arclength_interp,Epar_,arclength_interp,ni_,arclength_interp,ne_,arclength_interp,(ni_-ne_)*10)  
  %legend(hca,{'E_{||}','n_{i,cold}','n_{e,cold}','ni-ne'},'location','eastoutside')   
  plot(hca,arclength_interp,Epar_,'k')  
  legend(hca,{'E_{||}'},'location','eastoutside')   
  hca.YLabel.String = 'E_{||}';   
  hca.XLabel.String = 'arclength (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % dn
  hca = h(isub); isub = isub + 1;  
  %plot(hca,arclength_interp,Epar_,arclength_interp,ni_,arclength_interp,ne_,arclength_interp,(ni_-ne_)*10)  
  %legend(hca,{'E_{||}','n_{i,cold}','n_{e,cold}','ni-ne'},'location','eastoutside')   
  plot(hca,arclength_interp,(ni_-ne_))  
  legend(hca,{'n_{i,cold}-n_{e,cold}'},'location','eastoutside')   
  hca.YLabel.String = 'n';   
  hca.XLabel.String = 'arclength (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Epar, n
  hca = h(isub); isub = isub + 1;  
  %plot(hca,arclength_interp,Epar_,arclength_interp,ni_,arclength_interp,ne_,arclength_interp,(ni_-ne_)*10)  
  %legend(hca,{'E_{||}','n_{i,cold}','n_{e,cold}','ni-ne'},'location','eastoutside')   
  plot(hca,arclength_interp,Epar_,arclength_interp,(ni_-ne_)*10)  
  legend(hca,{'E_{||}','n_{i,cold}-n_{e,cold}'},'location','eastoutside')   
  hca.YLabel.String = 'E_{||}, n';   
  hca.XLabel.String = 'arclength (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % fi(vx)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.v,fredi.fvx')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{x}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},v_{x})'];
  hca.CLim(2) = prctile(fredi.fvx(:),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ex*max(abs(hca.YLim))/max(abs(Ex)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vix,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBx,'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fe(vx)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,frede.v,frede.fvx')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{x}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{e,' frede_str '}(l_{||},v_{x})'];
  hca.CLim(2) = prctile(frede.fvx(:),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ex*max(abs(hca.YLim))/max(abs(Ex)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vix,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBx,'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi(vy)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.v,fredi.fvy')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},v_{y})'];
  hca.CLim(2) = prctile(fredi.fvy(:),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBy,'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi(vz)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.v,fredi.fvz')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{z}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},v_{z})'];
  hca.CLim(2) = prctile(fredi.fvz(:),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ez*max(abs(hca.YLim))/max(abs(Ez)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,viz,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBz,'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred35_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fred.fvpar')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,cold}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(fred.fvpar(:),99);
  hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi3(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fred.fvpar')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,cold}^{top}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(fred.fvpar(:),99);
  hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred35_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  
  if 0 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.vpar_center,log10(fred.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{i,cold}(l_{||},v_{||})'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = [-6 -1];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % log 10 fi3(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred3_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{i,cold}^{top}(l_{||},v_{||})'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = [-6 -1];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % fi3/fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fredall = fred35_A75;
  fred = fred3_A75;
  fplot = (fred.fvpar./fredall.fvpar);
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fplot')
  view(hca,[0 0 1]);   
  shading(hca,'flat')
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.vpar_center,log10(fredall.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  %colormap(hca,pic_colors('candy')) 
  %colormap(hca,pic_colors('blue_red'))
  colormap(hca,pic_colors('pasteljet'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,cold}^{top}/f_{i,cold}^{tot}(l_{||},v_{||})'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.CLim = [0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fe(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred46_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')  
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',fred.fvpar')
  view(hca,[0 0 1]);  
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{e,cold}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(fred.fvpar(:),99);
  hca.CLim = fe_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-7 7];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vepar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log 10 fe(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred46_A75;
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')  
  surf(hca,arcedges,fred.vpar_edges,zeros(numel(arcedges),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]);  
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{e,cold}(l_{||},v_{||})'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fe_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-7 7];
  hca.CLim = [-4 -1.8];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
  if 1%doV
    hold(hca,'on')
    plot(hca,arclength,vepar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi(vabs)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.vabs_center,fredi.fvabs')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},|v|)'];
  hca.CLim(2) = prctile(fredi.fvabs(:),99);
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if 0%doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if 0%doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % defi35(m*abs^2/2)
  %%
  %isub = 6;
  hca = h(isub); isub = isub + 1;
  darc = arclength(2)-arclength(1);
  arclength_edges = [arclength-darc/2; arclength(end)+darc/2];
  [ARC,EN] = meshgrid(arclength_edges,fredi.vabs_edges.^2/2);
  surf(hca,ARC,EN,ARC*0,log10(fredi.fdefE'))
  shading(hca,'flat')
  view(hca,[0 0 1])
  %pcolor(hca,arclength,fredi.vabs_center.^2/2,log10(fredi.fdefE'))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'mv^2/2';
  colormap(hca,pic_colors('candy'))
  colormap(hca,[ 1 1 1; pic_colors('waterfall')])
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['dpf_{i,' fredi_str '}(l_{||},mv^2/2)'];
  hca.CLim = [prctile(log10(fredi.fdefE(fredi.fdefE>0)),1) prctile(log10(fredi.fdefE(fredi.fdefE>0)),99)];
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YScale = 'log';
  hca.YTick = 10.^(-10:1:10);
  hca.YLim(1) = 10^(-2);
  if doExB
    hold(hca,'on')
    plot(hca,arclength,EExB,'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
  if 0%doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if 0%doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % defe35(m*vabs^2/2)
  %%
  %isub = 6;
  hca = h(isub); isub = isub + 1;
  darc = arclength(2)-arclength(1);
  arclength_edges = [arclength-darc/2; arclength(end)+darc/2];
  [ARC,EN] = meshgrid(arclength_edges,frede.vabs_edges.^2/2/25);
  surf(hca,ARC,EN,ARC*0,log10(frede.fdefE'))
  shading(hca,'flat')
  view(hca,[0 0 1])
  %pcolor(hca,arclength,fredi.vabs_center.^2/2,log10(fredi.fdefE'))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'mv^2/2';
  colormap(hca,pic_colors('candy'))
  colormap(hca,[ 1 1 1; pic_colors('waterfall')])
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['dpf_{e,' frede_str '}(l_{||},mv^2/2)'];
  hca.CLim = [prctile(log10(frede.fdefE(frede.fdefE>0)),1) prctile(log10(frede.fdefE(frede.fdefE>0)),99)];
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YScale = 'log';
  hca.YTick = 10.^(-10:1:10);
  hca.YLim(1) = 10^(-2);
  if doExB
    hold(hca,'on')
    plot(hca,arclength,EExB/25,'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
  if 0%doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if 0%doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % f46(vpar)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,frede.vpar_center,frede.fvpar')
  hca.YLim = [-7 7];
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))
  hcb = colorbar('peer',hca);
  hca.CLim = fe_clim;
  hcb.YLabel.String = ['f_{e,' frede_str '}(l_{||},v_{||})'];  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if 1%doV
    hold(hca,'on')
    plot(hca,arclength,vepar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  
end
if 0 % def35(vabs^2)
  hca = h(isub); isub = isub + 1;
  darc = arclength(2)-arclength(1);
  arclength_edges = [arclength-darc/2 arclength(end)+darc/2];
  surf(hca,arclength_edges,frede.vabs_edges.^2/2/25,log10(frede.fdefE'))
  view([0 1 0])
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'mv^2/2';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['dpf_{i,' frede_str '}(l_{||},mv^2/2)'];
  hca.CLim = [prctile(log10(frede.fdefE(:)),5) prctile(log10(frede.fdefE(:)),95)];
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YScale = 'log';
  hca.YTick = 10.^(-10:1:10);
  hca.YLim(1) = 10^(-2);
  if doExB
    hold(hca,'on')
    plot(hca,arclength,EExB/25,'color',colorExB)
    hold(hca,'off')
  end
  if 0%doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if 0%doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
%
legends = {'a)','b)','c)','d)','e)','f)'};
for ip = 1:numel(h)
  irf_legend(h(ip),legends{ip},[-0.06 0.98],'fontsize',12,'color','k')
end
compact_panels(h(setdiff(1:nrows*ncols,isMap)),0.01)
compact_panels(h(isMap),0.01)
h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
%fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks_dist = linkprop(h(setdiff(1:nrows*ncols,isMap)),{'XLim'});
hlinks = linkprop(h(isMap),{'XLim','YLim'});

%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
for ip = 1:nrows*ncols
  axwidth(ip) = h(ip).Position(3);
end
for ip = 1:nrows*ncols
  h(ip).Position(3) = min(axwidth);
end



%ax(2).YAxisLocation = 'right';
if 0
  %%
  nrows = 3;
  ncols = 1;
  h = setup_subplots(nrows,ncols);
  isub = 1;  
  
  if 1 % line position
    hca = h(isub); isub = isub + 1;
    [ax,h1,h2] = plotyy(hca,arclength,fredi.x,arclength,fredi.z);
    hca.XLabel.String = 'arclength (d_i)';
    ax(1).YLabel.String = 'x';
    ax(2).YLabel.String = 'z';
    legend(hca,{'x','z'},'location','eastoutside')
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 1 % Epar
    hca = h(isub); isub = isub + 1;
    plot(hca,arclength,Eparx,arclength,Epary,arclength,Eparz,arclength,Epar)          
    legend(hca,{'E_{||}','-\int E_{||}dl_{||}','E_{||,x}','E_{||,y}','E_{||,z}'},'location','eastoutside') 
    legend(hca,{'E_{||,x}','E_{||,y}','E_{||,z}','E_{||}'},'location','eastoutside') 
    hca.YLabel.String = 'E_{||}'; 

    hca.XLabel.String = 'arclength (d_i)';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 1 % int Epar
    hca = h(isub); isub = isub + 1;
    plot(hca,arclength,-cumtrapz(arclength,Eparx),arclength,-cumtrapz(arclength,Epary),arclength,-cumtrapz(arclength,Eparz),arclength,-cumtrapz(arclength,Epar))
    %legend(hca,{'E_{||}','-\int E_{||}dl_{||}','E_{||,x}','E_{||,y}','E_{||,z}'},'location','eastoutside') 
    legend(hca,{'E_{||,x}','E_{||,y}','E_{||,z}','E_{||}'},'location','eastoutside') 
    hca.YLabel.String = '-\int E_{||}dl_{||}'; 

    hca.XLabel.String = 'arclength (d_i)';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  compact_panels(0.01)
  h = findobj(get(gcf,'children'),'type','axes'); hall = hall(end:-1:1);
  hlink = linkprop(h,{'XLim'});
  h(1).XLim = arclength([1 end]);
  for ip = 1:nrows*ncols
    axwidth(ip) = h(ip).Position(3);
  end
  for ip = 1:nrows*ncols
    h(ip).Position(3) = min(axwidth);
  end
end


