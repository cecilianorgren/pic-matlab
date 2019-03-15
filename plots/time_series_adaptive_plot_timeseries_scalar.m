% time_series_adaptive_plot_timeseries_scalar

%% Energies
dx_grid = x(2)-x(1);
dz_grid = z(2)-z(1);
dx_box = x(end)-x(1);
dz_box = z(end)-z(1);
dt = times(2)-times(1);
diff_times = times(1:end-1)+0.5*dt;
t_ref_index = 1; % 30

colors = pic_colors('matlab');
colors = [0 0 0; colors([2 1 3 5],:); colors([2 1 3 5],:); 0.8 0.8 0.8];
linestyles = {'-','-','-','-','-','--','--','--','--','-'};

nrows = 3;
ncols = 2;
h = setup_subplots(nrows,ncols,'horizontal');
arrayfun(@(x)set(x,'ColorOrder',colors),h)
arrayfun(@(x)set(x,'LineStyleOrder',linestyles),h)


isub = 1;
if 1 % Mean energy densities
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_mean',...
    'U.Uke1_mean','U.Uke2_mean','U.Uki1_mean','U.Uki2_mean',...
    'U.Ute1_mean','U.Ute2_mean','U.Uti1_mean','U.Uti2_mean'};
  
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));
    variable = cell_ts_scalar_{ivar_super}(1,1:54);
    all_variables(:,ivar) = variable;
    legends{ivar} = varstr;
    hold(hca,'on')
    hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),'-');
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
  hca.YLabel.String = {'Energy density','U(t)'};
end
if 0 % Sum energy densities
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_sum',...
    'U.Uke1_sum','U.Uke2_sum','U.Uki1_sum','U.Uki2_sum',...
    'U.Ute1_sum','U.Ute2_sum','U.Uti1_sum','U.Uti2_sum'};
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));
    variable = cell_ts_scalar_{ivar_super}(1,1:54);
    all_variables(:,ivar) = variable;
    legends{ivar} = varstr;
    hold(hca,'on')
    hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});    
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),linestyles{ivar+1});
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
end
if 1 % Mean energy densities * vol_box
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_mean',...
    'U.Uke1_mean','U.Uke2_mean','U.Uki1_mean','U.Uki2_mean',...
    'U.Ute1_mean','U.Ute2_mean','U.Uti1_mean','U.Uti2_mean'};
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));
    variable = cell_ts_scalar_{ivar_super}(1,1:54)*dx_box*dz_box;
    all_variables(:,ivar) = variable;    
    legends{ivar} = sprintf('%s*dx_box*dz_box',varstr);
    hold(hca,'on')
    hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});    
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),linestyles{ivar+1});
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
  hca.YLabel.String = {'Energy','U(t)'};
end
if 0 % Sum energy densities
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_sum',...
    'U.Uke1_sum','U.Uke2_sum','U.Uki1_sum','U.Uki2_sum',...
    'U.Ute1_sum','U.Ute2_sum','U.Uti1_sum','U.Uti2_sum'};
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));
    variable = cell_ts_scalar_{ivar_super}(1,1:54)*dx_grid*dz_grid;
    all_variables(:,ivar) = variable;
    legends{ivar} = sprintf('%s*dx_grid*dz_grid',varstr);
    hold(hca,'on')
    hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});    
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),linestyles{ivar+1});
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
end
% Relative to given time
if 1 % Mean energy densities 
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_mean',...
    'U.Uke1_mean','U.Uke2_mean','U.Uki1_mean','U.Uki2_mean',...
    'U.Ute1_mean','U.Ute2_mean','U.Uti1_mean','U.Uti2_mean'};
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));
    variable = cell_ts_scalar_{ivar_super}(1,1:54);
    variable = variable - variable(t_ref_index);
    all_variables(:,ivar) = variable;
    legends{ivar} = varstr;
    hold(hca,'on')
    if strcmp(varstr,'U.B_mean')
      hplot = plot(hca,times,-variable);
      legends{ivar} = sprintf('-%s',varstr);
    else
      hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});    
    end
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),linestyles{ivar+1});
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
  hca.YLabel.String = {'Change in energy density',sprintf('U = U(t)-U(%g)',times(t_ref_index))};
  
end
if 0 % Sum energy densities
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_sum',...
    'U.Uke1_sum','U.Uke2_sum','U.Uki1_sum','U.Uki2_sum',...
    'U.Ute1_sum','U.Ute2_sum','U.Uti1_sum','U.Uti2_sum'};
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));
    variable = cell_ts_scalar_{ivar_super}(1,1:54);
    variable = variable - variable(t_ref_index);
    all_variables(:,ivar) = variable;
    legends{ivar} = varstr;
    hold(hca,'on')
    if strcmp(varstr,'U.B_sum')
      hplot = plot(hca,times,-variable);
      legends{ivar} = sprintf('-%s',varstr);
    else
      hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});    
    end
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),linestyles{ivar+1});
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
end
if 1 % Mean energy densities * vol_box
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_mean',...
    'U.Uke1_mean','U.Uke2_mean','U.Uki1_mean','U.Uki2_mean',...
    'U.Ute1_mean','U.Ute2_mean','U.Uti1_mean','U.Uti2_mean'};
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));
    variable = cell_ts_scalar_{ivar_super}(1,1:54)*dx_box*dz_box;
    variable = variable - variable(t_ref_index);
    all_variables(:,ivar) = variable;    
    legends{ivar} = sprintf('%s*dx_box*dz_box',varstr);
    hold(hca,'on')
    if strcmp(varstr,'U.B_mean')
      hplot = plot(hca,times,-variable);
      legends{ivar} = sprintf('-%s',varstr);
    else
      hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});    
    end 
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),linestyles{ivar+1});
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
  hca.YLabel.String = {'Change in energy',sprintf('U = U(t)-U(%g)',times(t_ref_index))};
end
if 0 % Sum energy densities
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_sum',...
    'U.Uke1_sum','U.Uke2_sum','U.Uki1_sum','U.Uki2_sum',...
    'U.Ute1_sum','U.Ute2_sum','U.Uti1_sum','U.Uti2_sum'};
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));
    variable = cell_ts_scalar_{ivar_super}(1,1:54)*dx_grid*dz_grid;
    variable = variable - variable(t_ref_index);
    all_variables(:,ivar) = variable;
    legends{ivar} = sprintf('%s*dx_grid*dz_grid',varstr);
    hold(hca,'on')
    if strcmp(varstr,'U.B_sum')
      hplot = plot(hca,times,-variable);
      legends{ivar} = sprintf('-%s',varstr);
    else
      hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});    
    end
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),linestyles{ivar+1});
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
end
% Rate of changes
if 1 % Mean energy densities
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_mean',...
    'U.Uke1_mean','U.Uke2_mean','U.Uki1_mean','U.Uki2_mean',...
    'U.Ute1_mean','U.Ute2_mean','U.Uti1_mean','U.Uti2_mean'};
  
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));
    variable = cell_ts_scalar_{ivar_super}(1,1:54);    
    variable = interp1(diff_times,diff(variable),times);
    all_variables(:,ivar) = variable;
    legends{ivar} = varstr;
    hold(hca,'on')
    hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),'-');
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
  hca.YLabel.String = {'Rate of change of energy density','dU/dt'};
end
if 0 % Sum energy densities
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_sum',...
    'U.Uke1_sum','U.Uke2_sum','U.Uki1_sum','U.Uki2_sum',...
    'U.Ute1_sum','U.Ute2_sum','U.Uti1_sum','U.Uti2_sum'};
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));    
    variable = cell_ts_scalar_{ivar_super}(1,1:54);
    variable = interp1(diff_times,diff(variable),times);
    all_variables(:,ivar) = variable;
    legends{ivar} = varstr;
    hold(hca,'on')
    hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});    
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),linestyles{ivar+1});
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
end
if 1 % Mean energy densities * vol_box
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_mean',...
    'U.Uke1_mean','U.Uke2_mean','U.Uki1_mean','U.Uki2_mean',...
    'U.Ute1_mean','U.Ute2_mean','U.Uti1_mean','U.Uti2_mean'};
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));
    variable = cell_ts_scalar_{ivar_super}(1,1:54)*dx_box*dz_box;
    variable = interp1(diff_times,diff(variable),times);
    all_variables(:,ivar) = variable;    
    legends{ivar} = sprintf('%s*dx_box*dz_box',varstr);
    hold(hca,'on')
    hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});    
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),linestyles{ivar+1});
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
  hca.YLabel.String = {'Rate of change of energy','dU/dt'};
end
if 0 % Sum energy densities
  hca = h(isub); isub = isub + 1;
  varstrs = {'U.B_sum',...
    'U.Uke1_sum','U.Uke2_sum','U.Uki1_sum','U.Uki2_sum',...
    'U.Ute1_sum','U.Ute2_sum','U.Uti1_sum','U.Uti2_sum'};
  nvars = numel(varstrs);
  legends = cell(ivar,1);
  for ivar = 1:nvars
    varstr = varstrs{ivar};
    ivar_super = find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_scalar));
    variable = cell_ts_scalar_{ivar_super}(1,1:54)*dx_grid*dz_grid;
    variable = interp1(diff_times,diff(variable),times);
    all_variables(:,ivar) = variable;
    legends{ivar} = sprintf('%s*dx_grid*dz_grid',varstr);
    hold(hca,'on')
    hplot = plot(hca,times,variable,'linestyle',linestyles{ivar});    
    hold(hca,'off')
  end 
  hold(hca,'on')
  hplot = plot(hca,times,sum(all_variables,2),linestyles{ivar+1});
  legends{end+1} = 'Sum of all';
  hold(hca,'off')
  legend(hca,legends,'location','eastoutside','interpreter','none')
end


%compact_panels(0.012)
hlink = linkprop(h,{'XLim'});
hlink.Targets(1).XLim = [0 times(end)];
arrayfun(@(x)eval(sprintf('x.XGrid = ''on''; x.YGrid = ''on'';'),x),h)
arrayfun(@(x)eval(sprintf('x.Position(3) = 0.3;'),x),h)
%arrayfun(@(x)set(x,'XLabel','time'),h)
c_eval('h(?).XLabel.String = ''time'';',1:numel(h))
c_eval('h(?).Box = ''on'';',1:numel(h))

