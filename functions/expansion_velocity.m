function [v,x_vals,data_vals] = expansion_velocity(x,t,data,varargin)
% EXPANSION_VELOCITY Calculates expansion velocity
%   [velocity,location,value_at_location] = expansion_velocity(position,time,data)
%     data - time-space series with dimension (n_positions x n_times)

%%
doSingle = 0; % find single maximum
doPlot = 0;
doPos = 0; % find maximum for positive x
doNeg = 0; % find maximum for positive x

nargs = numel(varargin);
have_options = nargs;
args = varargin;
while have_options
  l = 1;
  switch(lower(args{1}))   
    case {'plot'}
      doPlot = 1;
      l = 1;
    case {'sym'}
      doSym = args{l+1};
      if doSym
        doPos = 1;
        doNeg = 1;
      end
      doSingle = 0;
      l = 2;
    case {'pos'}
      doPos = args{l+1};
      doSingle = 0;
      l = 2;
    case {'neg'}
      doNeg = args{l+1};
      doSingle = 0;
      l = 2;
  end
  args = args(l+1:end);  
  if isempty(args), break, end    
end

ntimes = numel(t);
npos = numel(x);
  
x_ind = {};
if doSingle
  x_ind{end+1} = find(x);
end
if doNeg
  x_ind{end+1} = find(x<mean(x));
end
if doPos
  x_ind{end+1} = find(x>mean(x));
end
nx_ind = numel(x_ind);

for itime = 1:ntimes  
  for ix_ind = 1:nx_ind
    [max_val_tmp,max_ind_tmp] = max(abs(data(x_ind{ix_ind},itime)));
    max_val{ix_ind}(itime) = max_val_tmp;      
    max_ind{ix_ind}(itime) = max_ind_tmp + min(x_ind{ix_ind});
    x_vals{ix_ind}(itime) = x(max_ind{ix_ind}(itime));
  end
end
%x_vals = max_val;
data_vals = max_val;
for ix_ind = 1:nx_ind
  v_tmp = diff(x_vals{ix_ind})./diff(t);
  v_tmp_interp = interp1(t(1:end-1)+0.5*(t(2)-t(1)),v_tmp,t);
  v{ix_ind} = v_tmp_interp;
end

if doPlot
  %%
  h(1) = subplot(2,1,1);
  h(2) = subplot(2,1,2);
  
  hca = subplot(2,1,1);
  imagesc(hca,t,x,data);
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'distance';
  for ix_ind = 1:nx_ind
    hold(hca,'on')
    plot(hca,t,x_vals{ix_ind})
    hold(hca,'off')
  end
  
  hca = subplot(2,1,2);
  for ix_ind = 1:nx_ind
    hold(hca,'on')
    plot(hca,t,v{ix_ind});
    hold(hca,'off')
  end
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'velocity';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  hlink = linkprop(h,{'XLim'});
  compact_panels
end