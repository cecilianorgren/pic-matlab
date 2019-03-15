% Data are here: cell_ts_line_x
% Description of data are here: varstrs_ts_line_x
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/';

nvars = numel(varstrs_ts_line_x);
nzpicks = size(cell_ts_line_x{1},2);
zpicks = zval_collect;

xlim = [-200 200];

for ivar = [10]%:nvars
  subdir = 'lines_x/';
  cmax = max(cell_ts_line_x{ivar}(:));
  cmin = min(cell_ts_line_x{ivar}(:));  
  clim = max(abs([cmin cmax]))*[-1 0];
  
  
  if cmax+cmin < 0.1*max(abs([cmin cmax])) % they are symmetrical    
     clim = max(abs([cmin cmax]))*[-1 1];
  end
  clim = 0.2*[-1 1];
  for  izpick = 1:nzpicks
    savestr = sprintf('%s_z_%s',strrep(varstrs_ts_line_x{ivar},'.','_'),strrep(num2str(zpicks(izpick)),'.','_'));    
    hca = subplot(1,1,1);imagesc(x,times,squeeze(cell_ts_line_x{ivar}(:,izpick,:))');
    hca.XDir = 'reverse';
    hca.YDir = 'normal';
    hca.XLim = xlim;
    hca.CLim = clim;
    hca.YLabel.String = sprintf('time (1/wci)');
    hca.XLabel.String = sprintf('x (di)');
    hca.Title.String = sprintf('z = %g (di)',zpicks(izpick));    
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = varstrs_ts_line_x{ivar};
    hca.FontSize = 12;
    if not(exist([savedir_root subdir],'dir')), mkdir([savedir_root subdir]); end
    colormap(pic_colors('blue_red'))
    print('-dpng','-r200',[savedir_root subdir savestr '.png']);
    %pause(0.1)
  end
end

%% Calculate reconnected flux in a few different ways
dx = x(2)-x(1);
dt = times(2)-times(1);
ix0 = find_closest_ind(x,0);

R.times = times;

% Bz
ivar = find(cellfun(@(x)strcmp(x,'B.z'),varstrs_ts_line_x));
izpick = find_closest_ind(zpicks,0);
Bz = squeeze(cell_ts_line_x{ivar}(:,izpick,:)); % Bz(x,z=0,t)
flux = sum(Bz(ix0:end,:),1)*dx;
R.flux = flux;
R.Bz = interp1(times(2:end)-0.5*dt,diff(flux)/dt,times);
  
% Ey
ivar = find(cellfun(@(x)strcmp(x,'E.y'),varstrs_ts_line_x));
izpick = find_closest_ind(zpicks,0);
Ey = squeeze(cell_ts_line_x{ivar}(ix0,izpick,:)); % Ey(x,z=0,t)
R.Ey = reshape(Ey,1,numel(Ey));

%% Plot reconnected flux as a function of time
h = setup_subplots(4,1);
isub = 1;
if 1 % time map of Ey at z = 0
  hca = h(isub); isub = isub + 1;
  ivar = find(cellfun(@(x)strcmp(x,'E.y'),varstrs_ts_line_x));
  izpick = find_closest_ind(zpicks,0);
  imagesc(hca,times,x,squeeze(cell_ts_line_x{ivar}(:,izpick,:)));  
  hca.YDir = 'normal';
  hca.CLim = max(max(abs(cell_ts_line_x{ivar}(:,izpick,:))))*[-1 1];  
  hca.XLabel.String = sprintf('time (1/wci)');
  hca.YLabel.String = sprintf('x (di)');
  %hca.Title.String = sprintf('z = %g (di)',zpicks(izpick));    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstrs_ts_line_x{ivar};  
  colormap(hca,pic_colors('blue_red'))
  hca.YLim = [-100 100];
  htext = text(hca,hca.XLim(1),hca.YLim(2), sprintf('z = %g (di)',zpicks(izpick)),'horizontalalignment','left','verticalalignment','top','fontsize',12);
end
if 1 % time map of Bz at z = 0
  hca = h(isub); isub = isub + 1;
  ivar = find(cellfun(@(x)strcmp(x,'B.z'),varstrs_ts_line_x));
  izpick = find_closest_ind(zpicks,0);
  imagesc(hca,times,x,squeeze(cell_ts_line_x{ivar}(:,izpick,:)));  
  hca.YDir = 'normal';  
  hca.CLim = max(max(abs(cell_ts_line_x{ivar}(:,izpick,:))))*[-1 1];
  hca.XLabel.String = sprintf('time (1/wci)');
  hca.YLabel.String = sprintf('x (di)');
  %hca.Title.String = sprintf('z = %g (di)',zpicks(izpick));    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstrs_ts_line_x{ivar};
  colormap(hca,pic_colors('blue_red'))
  hca.YLim = [-100 100];
  htext = text(hca,hca.XLim(1),hca.YLim(2), sprintf('z = %g (di)',zpicks(izpick)),'horizontalalignment','left','verticalalignment','top','fontsize',12);
end
if 1 % reconnected flux
  hca = h(isub); isub = isub + 1;
  plot(hca,R.times,R.flux)
  hca.YLabel.String = 'Outflow flux (...)';
  hca.XLabel.String = sprintf('time (1/wci)');
  legend(hca,{sprintf('int(Bzdx), x=[%.0f,%.0f]',x(ix0),x(end))},'location','northwest','box','off')
end
if 1 % reconnection rate
  hca = h(isub); isub = isub + 1;
  plot(hca,R.times,[R.Bz; R.Ey])
  hca.YLabel.String = 'Reconnection rate (...)';
  legend(hca,{sprintf('(d/dt)int(Bzdx), x=[%.0f,%.0f]',x(ix0),x(end)),'Ey @ x=z=0'},'location','northwest','box','off')
  hca.XLabel.String = sprintf('time (1/wci)');
end

arrayfun(@(x)set(x,'FontSize',12),h)
arrayfun(@(x)set(x,'XLim',[0 times(end)]),h)
arrayfun(@(x)set(x,'XGrid','on','YGrid','on'),h)
for ipanel = 1:numel(h)
  h(ipanel).Position(3) = 0.7;
end
arrayfun(@(x)set(x,'XTickLabels',[],'XLabel',[]),h(1:end-1))
compact_panels

%% Plot reconnected flux as a function of 
h = setup_subplots(4,1);
isub = 1;
if 1 % time map of Ey at z = 0
  hca = h(isub); isub = isub + 1;
  ivar = find(cellfun(@(x)strcmp(x,'E.y'),varstrs_ts_line_x));
  izpick = find_closest_ind(zpicks,0);
  imagesc(hca,times,x,squeeze(cell_ts_line_x{ivar}(:,izpick,:)));  
  hca.YDir = 'normal';
  hca.CLim = max(max(abs(cell_ts_line_x{ivar}(:,izpick,:))))*[-1 1];  
  hca.XLabel.String = sprintf('time (1/wci)');
  hca.YLabel.String = sprintf('x (di)');
  %hca.Title.String = sprintf('z = %g (di)',zpicks(izpick));    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstrs_ts_line_x{ivar};  
  colormap(hca,pic_colors('blue_red'))
  hca.YLim = [-100 100];
  htext = text(hca,hca.XLim(1),hca.YLim(2), sprintf('z = %g (di)',zpicks(izpick)),'horizontalalignment','left','verticalalignment','top','fontsize',12);
end
if 1 % time map of Bz at z = 0
  hca = h(isub); isub = isub + 1;
  ivar = find(cellfun(@(x)strcmp(x,'B.z'),varstrs_ts_line_x));
  izpick = find_closest_ind(zpicks,0);
  imagesc(hca,times,x,squeeze(cell_ts_line_x{ivar}(:,izpick,:)));  
  hca.YDir = 'normal';  
  hca.CLim = max(max(abs(cell_ts_line_x{ivar}(:,izpick,:))))*[-1 1];
  hca.XLabel.String = sprintf('time (1/wci)');
  hca.YLabel.String = sprintf('x (di)');
  %hca.Title.String = sprintf('z = %g (di)',zpicks(izpick));    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstrs_ts_line_x{ivar};
  colormap(hca,pic_colors('blue_red'))
  hca.YLim = [-100 100];
  htext = text(hca,hca.XLim(1),hca.YLim(2), sprintf('z = %g (di)',zpicks(izpick)),'horizontalalignment','left','verticalalignment','top','fontsize',12);
end
if 1 % reconnected flux
  hca = h(isub); isub = isub + 1;
  plot(hca,R.times,R.flux)
  hca.YLabel.String = 'Outflow flux (...)';
  hca.XLabel.String = sprintf('time (1/wci)');
  legend(hca,{sprintf('int(Bzdx), x=[%.0f,%.0f]',x(ix0),x(end))},'location','northwest','box','off')
end
if 1 % reconnection rate
  hca = h(isub); isub = isub + 1;
  plot(hca,R.times,[R.Bz; R.Ey])
  hca.YLabel.String = 'Reconnection rate (...)';
  legend(hca,{sprintf('(d/dt)int(Bzdx), x=[%.0f,%.0f]',x(ix0),x(end)),'Ey @ x=z=0'},'location','northwest','box','off')
  hca.XLabel.String = sprintf('time (1/wci)');
end

arrayfun(@(x)set(x,'FontSize',12),h)
arrayfun(@(x)set(x,'XLim',[0 times(end)]),h)
arrayfun(@(x)set(x,'XGrid','on','YGrid','on'),h)
for ipanel = 1:numel(h)
  h(ipanel).Position(3) = 0.7;
end
arrayfun(@(x)set(x,'XTickLabels',[],'XLabel',[]),h(1:end-1))
compact_panels
