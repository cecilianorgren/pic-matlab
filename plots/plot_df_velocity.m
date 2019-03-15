input_varstr = 'B.z';
input_ivar = find(cellfun(@(x)strcmp(x,input_varstr),varstrs_ts_line_x));
input_data_for_velocity = squeeze(cell_ts_line_x{input_ivar}(:,11,:));
[front_velocity,front_location,value_at_location] = expansion_velocity(x,times,input_data_for_velocity,'pos',1);
n_vel = numel(front_velocity);

% If a particle was reflected at the front at time = t1, where would it be
% at a alter time t2? First we get the velocity as a function of the front
% speed: v2 = -v1+2*v_front, we assume v1 = 0
velocity_before_reflection = 0.5;
velocity_after_reflection = cellfun(@(x)2*x+velocity_before_reflection,front_velocity,'UniformOutput',false);

% n_particles = ntimes;
% particle_velocity = zeros(ntimes,n_particles);
% for itime = 1:ntimes
%   particle_velocity(itime,itime:end) = velocity_after_reflection{1}(itime);
%   particle_location_{iparticle}(1) = front_location{1}(iparticle);
% end

% reflection location is front location
% particle velocity = front_location+(t-t_reflection)*velocity_after_reflection
% calculate this for each time step

% set location prior to/at the time of the reflection = location of front
particle_location = cell(ntimes,1);
for iparticle = 1:n_particles
  particle_location{iparticle}(1) = NaN;%front_location{1}(iparticle);
end
% Integrate particle position 
dt = times(2)-times(1);
for itime = 1:ntimes
  %fprintf('itime = %g, time: %g \n',itime,times(itime))
  for iparticle = 1:n_particles
   % fprintf('Particle: %g \n',iparticle)
    %fprintf('Particle location: %g, Front location: %g \n',particle_location{iparticle}(itime),front_location{1}(itime))
    if iparticle > itime % particle has not yet encountered the front
      particle_location{iparticle}(itime) = NaN;
    elseif iparticle == itime
      particle_location{iparticle}(itime) = front_location{1}(itime);   
    else % particle has encountered the front, advance position            
      particle_location{iparticle}(itime) = particle_location{iparticle}(itime-1) + dt*velocity_after_reflection{1}(iparticle);
    end
  end
end

% Plot
if exist('h','var'); delete(h); end
nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);

isub = 1;
if 0 % Data used to get the expansion velocity
  hca = h(isub); isub = isub + 1;
  himag = imagesc(hca,times,x,input_data_for_velocity);
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = input_varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  for i_vel = 1:n_vel
    hold(hca,'on')
    plot(hca,times,front_location{i_vel},'k')
    hold(hca,'off')
  end  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Location of expansion front
  hca = h(isub); isub = isub + 1;
  prefix = {'1','1'};
  operator = cellfun(@(x)eval(x),prefix);
  for i_vel = 1:n_vel
    %velocity_tmp = eval(sprintf('%s,%s;',velocity{i_vel}));
    if i_vel == 1, hold(hca,'on'); end
    plot(hca,times,operator(i_vel)*front_location{i_vel});
    if i_vel == n_vel, hold(hca,'off'); end
  end
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'front position';
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
end
if 0 % Velocity of expansion
  hca = h(isub); isub = isub + 1;
  prefix = {'1','1'};
  operator = cellfun(@(x)eval(x),prefix);
  for i_vel = 1:n_vel
    %velocity_tmp = eval(sprintf('%s,%s;',velocity{i_vel}));
    if i_vel == 1, hold(hca,'on'); end
    plot(hca,times,operator(i_vel)*front_velocity{i_vel});
    if i_vel == n_vel, hold(hca,'off'); end
  end
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'velocity';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Velocity of particle after having encountered the front at the given time
  hca = h(isub); isub = isub + 1;
  prefix = {'1','1'};
  operator = cellfun(@(x)eval(x),prefix);
  for i_vel = 1:n_vel
    %velocity_tmp = eval(sprintf('%s,%s;',velocity{i_vel}));
    if i_vel == 1, hold(hca,'on'); end
    plot(hca,times,operator(i_vel)*velocity_after_reflection{i_vel});
    if i_vel == n_vel, hold(hca,'off'); end
  end
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = {'particle velocity','after reflection at front'};
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  irf_legend(hca,{'velocity after reflection = - velocity before reflection + 2*vfront';...
    sprintf('velocity before reflection = %g',velocity_before_reflection)},[0.02 0.98],'k')
end
if 1 % Velocity of front and particle after having encountered the front at a given time
  hca = h(isub); isub = isub + 1;
  prefix = {'1','1'};
  operator = cellfun(@(x)eval(x),prefix);
  for i_vel = 1:n_vel
    %velocity_tmp = eval(sprintf('%s,%s;',velocity{i_vel}));
    if i_vel == 1, hold(hca,'on'); end
    plot(hca,times,operator(i_vel)*front_velocity{i_vel});
    plot(hca,times,operator(i_vel)*velocity_after_reflection{i_vel});
    if i_vel == n_vel, hold(hca,'off'); end
  end
  
%   for i_vel = 1:n_vel
%     %velocity_tmp = eval(sprintf('%s,%s;',velocity{i_vel}));
%     if i_vel == 1, hold(hca,'on'); end
%     plot(hca,times,operator(i_vel)*front_velocity{i_vel});
%     if i_vel == n_vel, hold(hca,'off'); end
%   end
  legend(hca,{'velocity of front','particle velocity after reflection at front'},'box','off','location','northwest')
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = {'particle velocity','after reflection at front'};
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hleg = irf_legend(hca,{'velocity after reflection = - velocity before reflection + 2*vfront';...
      sprintf('velocity before reflection = %g',velocity_before_reflection)},[0.1 0.5]);
  arrayfun(@(x)eval(sprintf('x.Color = ''k'';'),x),hleg)
end

if 0 % Location of particle
  hca = h(isub); isub = isub + 1;
  hold(hca,'on');
  for i_particle = 1:n_particles    
    plot(hca,times,particle_location{i_particle});
  end
  hold(hca,'off');
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'front position';
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
end
  zzz = 2;
if 1 % Plot particle lines starting from front
  hca = h(isub); isub = isub + 1;
  zzind = find_closest_ind(zpicks,zzz);
  varstr = 'E.z';
  plot_data = squeeze(cell_ts_line_x{find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_line_x))}(:,zzind,:));
  himag = imagesc(hca,times,x,plot_data);
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  for i_vel = 1:n_vel
    hold(hca,'on')
    hfront = plot(hca,times,front_location{i_vel},'color','k','linewidth',1);
    hold(hca,'off')
  end  
  irf_legend(hca,{sprintf('z = %g',zpicks(zzind))},[0.02 0.98],'k')
  
  hold(hca,'on');  
  for i_particle = 21:2:n_particles
    plot(hca,times,particle_location{i_particle},'k','linewidth',0.5);
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.YLim(1) = 0;
  %hca.YLim(2) = 100;
end
if 1 % Plot particle lines starting from front
  hca = h(isub); isub = isub + 1;
  zzind = find_closest_ind(zpicks,zzz);
  varstr = 'vi2.x';
  plot_data = squeeze(cell_ts_line_x{find(cellfun(@(x)strcmp(x,varstr),varstrs_ts_line_x))}(:,zzind,:));
  himag = imagesc(hca,times,x,plot_data);
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  colormap(pic_colors('blue_red'))
  hca.CLim = max(abs(himag.CData(:)))*[-1 1];
  hca.XDir = 'normal';
  hca.YDir = 'normal';
  hca.XLabel.String = 'time';
  hca.YLabel.String = 'x (di)';
  for i_vel = 1:n_vel
    hold(hca,'on')
    hfront = plot(hca,times,front_location{i_vel},'color','k','linewidth',1);
    hold(hca,'off')
  end  
  irf_legend(hca,{sprintf('z = %g',zpicks(zzind))},[0.02 0.98],'k')
  
  hold(hca,'on');  
  for i_particle = 21:2:n_particles
    plot(hca,times,particle_location{i_particle},'k','linewidth',0.5);
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %hca.YLim(1) = 0;
  %hca.YLim(2) = 100;
end

arrayfun(@(x)eval(sprintf('x.Position(3) = 0.7;'),x),h)
arrayfun(@(x)eval(sprintf('x.XLabel.String = [];'),x),h(1:end-1))
arrayfun(@(x)eval(sprintf('x.Box = ''on'';'),x),h)
hlink = linkprop(h,{'XLim'});
hlink.Targets(1).XLim(1) = 80;
hlink.Targets(1).XLim(2) = times(end);
arrayfun(@(x)eval(sprintf('x.YLim = [0 200];'),x),h(2:3))
compact_panels
