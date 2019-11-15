% Run time_series_adaptive to get desired quantities

%%
% Initialize figure 
fig = figure(109);
nrows = 4;
ncols = 6;
npanels = nrows*ncols;
clear h hb;
for ipanel = 1:npanels  
  h(ipanel) = subplot(nrows,ncols,ipanel); 
end

nvars = numel(varstrs_ts_stacked_adapted);

for zind = 11:21
isub = 1;
  for ivar = 1:nvars%[1 2 15]%[3 7 11 16 20]+2%[1 2]%[3:6 7:10 11:14 16:19 20:23] %1:nvars
    hca = h(isub); isub = isub + 1;
    variable_str = varstrs_ts_stacked_adapted{ivar};
    variable = eval(varstrs_ts_stacked_adapted{ivar});

    himag = imagesc(hca,x,timesteps/wpewce/mass(1),squeeze(variable(:,zind,:))');
    if abs(himag.CData(not(isnan(himag.CData)))) % dont do if is zero
      hca.CLim = max(abs(himag.CData(not(isnan(himag.CData)))))*[-1 1];  
    else
      %pause
      hca.CLim = max(abs(himag.CData(not(isnan(himag.CData)))))*[-1 1];  
    end

    %hca.CLim = max(abs(himag.CData(not(isnan(himag.CData)))))*[-1 1];  % dont do if is zero

    hcb = colorbar('peer',hca);
    hca.Title.String = variable_str;
    hca.Title.Interpreter = 'none';
    hca.YLabel.String = 'time (1/wci)';
    hca.XLabel.String = 'x (di)';
    hca.YDir = 'normal';
    colormap(pic_colors('blue_red'))
  end
  for ipanel = 1:npanels    
    %h(ipanel).Position(3) = h(ipanel).Position(3)*1.1;
    %h(ipanel).Position(4) = h(ipanel).Position(4)*1.1;
  end
  print('-dpng','-r200',sprintf('%s/ts_stacked_all/stacked_time_all_z=%.0f_20190296.png',savedir_root,z(ind_z0(zind))));
end

%% Adaptive
isub = 1;
savedir = 'ts_stacked_all';
varstrs_ts_stacked = {'B.z','E.y',...
              've1.x','ve2.x','vi1.x','vi2.x'...
              };
            

varstrs_ts_stacked = {'sep_vepar',...
           'sep_ne',...
           };            
nx = 6400;
nvars_ts_stacked = numel(varstrs_ts_stacked);
for ivar_ts = 1:nvars_ts_stacked
  ivar_ts
  disp([varstrs_ts_stacked{ivar} '_ts_stacked'])
  %variable = eval([varstr_ts_stacked '_ts_stacked']);
end
%%
for ivar_ts = 1:nvars_ts_stacked
    %disp([varstrs_ts_stacked_adapted{ivar_ts} '(:,itime) = ' varstrs_ts_stacked{ivar_ts} '(:,ind_z0);']);    
    %eval([varstrs_ts_stacked_adapted{ivar_ts} '(:,itime) = ' varstrs_ts_stacked{ivar_ts} '(:,ind_z0);']);    
  hca = h(isub); isub = isub + 1;subdir = 'Bz_at_x=0';
%  savedir = [savedir_root,subdir];
 % mkdir(savedir)
 % savestr = sprintf('%s_t%05.0f',subdir,timestep);    
  figure(33)
  hca = subplot(2,1,1);
  plot(hca,x,eval([varstrs_ts_stacked{ivar_ts} '(:,ind_z0);']))      
  hca.XLabel.String = 'x (di)';
  hca.YLabel.String = 'B.z';
  hca.Title.String = 'B.z at z = 0';
  hca = subplot(2,1,2);
  imagesc(hca,timesteps/wpewce/mass(1),x,eval(varstrs_ts_stacked_adapted{ivar_ts}))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'B.z';
  hca.XLabel.String = 'time (1/wci)';
  hca.YLabel.String = 'x (di)';
  print('-dpng','-r200',[savedir '/' savestr '.png']);
    
  himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
  hca.Title.String = sprintf('%s',varstr); 
  hca.Title.Interpreter = 'none';  
  if abs(himag.CData(not(isnan(himag.CData)))) % dont do if is zero
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
  end
  hcb = colorbar('peer',hca);
  hb(isub-1) = hcb;
  %hcb.YLim = hca.CLim(2)*[-1 1];
  colormap(hca,cn.cmap('blue_red'));
end