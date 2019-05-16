% Compare flux transfer of multiple simulations

% Can be different for different simulations
% Sub-directories have similar structure
data_root_dirs = {
  '/Volumes/Fountain/Data/PIC/df_cold_protons_1/data_separated/',...
  '/Volumes/pic/in_progress/df_cold_protons_04/data_separated/'
  };

timesteps = {200:200:10800,200:200:12000}; % should be possible for different timesteps
ntimes = cellfun(@(x)numel(x),timesteps);

nsims = numel(data_root_dirs);
% NB: populations do not have the same names!
varstrs = {
  {'ni1'},...
  {'ni1'}...
  };
nvars = cellfun(@(x)numel(x),varstrs);

%% load data
for isim = 1:nsims
  varstr_reload = sprintf('%s/%s.mat',[data_root_dirs{isim} 'same_for_all_times'],'sim_info');
  siminfo = load(varstr_reload); 
  eval(sprintf('%s_%g = siminfo;','siminfo',isim))
  for ivar = 1:nvars(isim)
    isInitialized = 0;
    varstr = varstrs{isim}{ivar};
    vardir = [data_root_dirs{isim} varstr];
    
    fprintf('Loading: sim %g: %s \n',isim,varstr)
    tic;    
    for itime = 1:ntimes(isim)
      %eval(sprintf('%s'),varstrs);
      timestep = timesteps{isim}(itime);    
      varstr_reload = sprintf('%s/%s-%05.0f.mat',vardir,varstr,timestep);
      data_tmp = load(varstr_reload,varstr); 
      data_tmp = eval(sprintf('data_tmp.%s',varstr));    
      if isnumeric(data_tmp)
        if not(isInitialized)
          data = zeros([ntimes(isim) size(data_tmp)]); 
          isInitialized = 1;
        end
        data(itime,:,:,:,:,:) = data_tmp;
      elseif isstruct(data_tmp)
        var_fields = fields(data_tmp);
        nfields = numel(var_fields);
        vec_fields = {'x';'y';'z'};
        tens_fields = {'xx';'xy';'xz';'yy';'yz';'zz'};
        if all(cellfun(@isequal,var_fields,vec_fields))
          datasize = size(data_tmp.x);
          if not(isInitialized)
            data = zeros([ntimes(isim) datasize numel(vec_fields)]); 
            isInitialized = 1;
          end
          for ifield = 1:nfields
            %data_field = eval([]);
            data(itime,:,:,ifield) = eval(['data_tmp.' var_fields{ifield}]);
          end
          %data(itime,:,:,) = data_tmp;
        elseif all(cellfun(@isequal,var_fields,tens_fields))  
          % tensor fields not supported yet.
        end
      end
    end
    eval(sprintf('%s_%g = data;',varstr,isim))
    % whos('*_*')
    toc
  end
end

%% Make gif of flux contour A(x,z)
savedir_gifs = '/Users/cno062/Research/PIC/df_cold_protons_04/gifs/';
savedir_gifs = '/Users/cno062/Research/PIC/combination_of_runs/gifs/';
hca = subplot(1,1,1);
set(gcf,'color',[1 1 1]);

nx = 6400;
nz = 1600;
plx = 1:5:nx;
plz = 1:5:nz;
ntimes = 54;
imovie = 1;
cell_movies = cell(1,2);
for itime = 1:ntimes
 if 1%doMovie(iplot) % Collect frame for movie    
    currentBackgroundColor = get(gcf,'color');
    Alevels = -25:1:0;
    contour(hca,-x(plx),z(plz),squeeze(A_1(itime,plx,plz))',Alevels,'-','displayname','n_c=0.8','linewidth',1.5);
    hold(hca,'on')
    contour(hca,x(plx)*0.95,z(plz)*0.95,squeeze(A_2(itime,plx,plz))',Alevels,'-','displayname','n_c=0.4');
    hold(hca,'off')
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = sprintf('t = %g w_{pe} = %g w_{ci}',timesteps{1}(itime),timesteps{1}(itime)/50);
    drawnow
    legend
    pause(0.1)
    %imovie = sum(doMovie(1:iplot));
    tmp_frame = getframe(gcf);
    %cell_movies{imovie}(itime) = tmp_frame;
    if itime == 1 % initialize animated gif matrix
      [im_tmp,map] = rgb2ind(tmp_frame.cdata,256,'nodither');
      %map(end+1,:) = get(gcf,'color');
      im_tmp(1,1,1,ntimes) = 0;
      cell_movies{imovie,1} = map;
      cell_movies{imovie,2} = im_tmp;
    else
      cell_movies{imovie,2}(:,:,1,itime) = rgb2ind(tmp_frame.cdata,cell_movies{imovie,1},'nodither');
    end       
 end
end
%imwrite(cell_movies{ivar,2},cell_movies{ivar,1},[savedir_gifs 'A_12.gif'],'DelayTime',0.0,'LoopCount',0)

%% Rescale time to different quntity C, for example 
% 1. Spent magnetic energy, dB(t)
% 2. A value at X line, AX(t)
% We need the quantity as a function of t. Unless the relation between the
% quantity C and t is linear, we time frames will not be equally spaced.
% Therefore we need to interpolate the fields to equispaced C's. So that if
% we want to make a movie, each frame is the progress for the same change
% in C.
% Practically we can first define our C. For:
% 2. AX(t) 

dA = 1/4;
AX_1 = ceil(A_1_xline(1)*4)/4:dA:floor(A_1_xline(end)*4)/4; % keep within existing boundaries
AX_2 = ceil(A_2_xline(1)*4)/4:dA:floor(A_2_xline(end)*4)/4; % keep within existing boundaries

tAX_1 = interp1(A_1_xline',timesteps{1}/50,AX_1); % interpolated times
tAX_2 = interp1(A_2_xline',timesteps{2}/50,AX_2); % interpolated times

hca = subplot(1,1,1);
plot(hca,timesteps{1}/50,A_1_xline,'.',tAX_1,AX_1,'o'...
  ,timesteps{2}/50,A_2_xline,'.',tAX_2,AX_2,'o'); % looks fine
legend(hca,{sprintf('(n_c=0.8n_0): A_X(t): dt = %g',(timesteps{1}(2)-timesteps{1}(1))/50),sprintf('(n_c=0.8n_0): t(A_X): dA_X = %g',dA),...
            sprintf('(n_c=0.4n_0): A_X(t): dt = %g',(timesteps{1}(2)-timesteps{1}(1))/50),sprintf('(n_c=0.4n_0): t(A_X): dA_X = %g',dA)},...
       'location','northwest')
hca.XLabel.String = 't (w_{ci})';
hca.YLabel.String = 'A_X = A(t,x_X,z_X) (...)';
hca.YMinorTick = 'on';
hca.YMinorGrid = 'on';

% interpolate the entire A
% I think the most time consuming is doing the for-loop in nx and nz, so
% might as well do
A_1_interp = interpolate_time(timesteps{1}/50,A_1,tAX_1); % takes a long time
A_2_interp = interpolate_time(timesteps{2}/50,A_2,tAX_2); % takes a long time
%save('/Users/cno062/Research/PIC/combination_of_runs/A_and_A_interp','AX_1','AX_2','tAX_1','tAX_2','A_1_xline','A_2_xline','A_1_interp','A_2_interp','A_1','A_2')
%load('/Users/cno062/Research/PIC/combination_of_runs/A_and_A_interp');

%% Interpolate A's to spent magnetic energy (calculated below in script)
dt = (timesteps{1}(2)/50-timesteps{1}(1)/50);
spentUB_1 = interp1(timesteps{1}(1:end-1)/50+0.5*dt,-diff(U_1.B),timesteps{1}(1:end)/50);
spentUB_2 = interp1(timesteps{2}(1:end-1)/50+0.5*dt,-diff(U_2.B),timesteps{2}(1:end)/50);

dSpentUB = 0.05e-3;
rounding = 1e4;
spentUB_1_equi = ceil(spentUB_1(2)*rounding)/rounding:dSpentUB:floor(spentUB_1(end-1)*rounding)/rounding; % keep within existing boundaries
spentUB_2_equi = ceil(spentUB_2(2)*rounding)/rounding:dSpentUB:floor(spentUB_2(end-1)*rounding)/rounding; % keep within existing boundaries
% rounding errors gives different vectors, fix below
spentUB_1_equi = spentUB_2_equi(1:numel(spentUB_1_equi));

tUB_1 = interp1(spentUB_1(2:end-1)',timesteps{1}(2:end-1)/50,spentUB_1_equi); % interpolated times
tUB_2 = interp1(spentUB_2(2:end-1)',timesteps{2}(2:end-1)/50,spentUB_2_equi); % interpolated times


hca = subplot(1,1,1);
plot(hca,timesteps{1}/50,spentUB_1,'.',tUB_1,spentUB_1_equi,'o',...
         timesteps{2}/50,spentUB_2,'.',tUB_2,spentUB_2_equi,'o'); % looks fine
%legend(hca,{sprintf('(n_c=0.8n_0): A_X(t): dt = %g',(timesteps{1}(2)-timesteps{1}(1))/50),sprintf('(n_c=0.8n_0): t(A_X): dA_X = %g',dA),...
%            sprintf('(n_c=0.4n_0): A_X(t): dt = %g',(timesteps{1}(2)-timesteps{1}(1))/50),sprintf('(n_c=0.4n_0): t(A_X): dA_X = %g',dA)},...
%       'location','northwest')
hca.XLabel.String = 't (w_{ci})';
hca.YLabel.String = '-\Delta <U_B> (...)';
hca.YMinorTick = 'on';
hca.YMinorGrid = 'on';

A_1_interp_spentUB = interpolate_time(timesteps{1}/50,A_1,tUB_1); % takes a long time
A_2_interp_spentUB = interpolate_time(timesteps{2}/50,A_2,tUB_2); % takes a long time

%% Make gif of flux contour A(x,z) rescaled to A_X(t)
savedir_gifs = '/Users/cno062/Research/PIC/df_cold_protons_04/gifs/';
savedir_gifs = '/Users/cno062/Research/PIC/combination_of_runs/gifs/';
hca = subplot(1,1,1);
set(gcf,'color',[1 1 1]);

x = siminfo_1.x-mean(siminfo_1.x);
z = siminfo_1.z;
nx = 6400;
nz = 1600;
plx = 1:5:nx;
plz = 1:5:nz;
ntimes = 54;
imovie = 1;
cell_movies = cell(1,2);
A_common = unique([AX_1,AX_2]);
nA_common = numel(A_common);
doholdon = 0;
for iA = 1:nA_common
doholdon = 0;
 if 1%doMovie(iplot) % Collect frame for movie    
    currentBackgroundColor = get(gcf,'color');
    Alevels = -25:0.5:0;
    if 0%intersect(AX_1,A_common(iA))
      displayname = sprintf('n_c=0.8,n_0 tw_{ci} = %.0f',tAX_1(iA));
      contour(hca,-x(plx),z(plz),squeeze(A_1_interp(iA,plx,plz))',Alevels,'-','displayname',displayname,'linewidth',1.5);
      doholdon = 1;
    else % keep last lines plot, to see if they eventually becomes similar again    
      iAplot = numel(AX_2);
      iAplot = 12;  
      displayname = sprintf('n_c=0.8n_0, tw_{ci} = %.0f',tAX_1(iAplot));      
      contour(hca,x(plx),z(plz),squeeze(A_1_interp(iAplot,plx,plz))',Alevels,'-','displayname',displayname,'linewidth',1.5);      
      doholdon = 1;
    end
    if intersect(AX_2,A_common(iA))
      displayname = sprintf('n_c=0.4n_0, tw_{ci} = %.0f',tAX_2(iA));
      if doholdon; hold(hca,'on'); end
      contour(hca,x(plx),z(plz),squeeze(A_2_interp(iA,plx,plz))',Alevels,'-','displayname',displayname);
      if doholdon; hold(hca,'off'); end
    else % keep last lines plot, to see if they eventually becomes similar again
      %iAplot = numel(AX_2);
      %iAplot = 12;
      displayname = sprintf('n_c=0.4n_0, tw_{ci} = %.0f',tAX_2(iAplot));
      if doholdon; hold(hca,'on'); end
      contour(hca,x(plx),z(plz),squeeze(A_2_interp(iAplot,plx,plz))',Alevels,'-','displayname',displayname);
      if doholdon; hold(hca,'off'); end
    end
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = sprintf('A_X = %g ',A_common(iA));
    hca.XLim(1) = -120;
    hca.XLim(2) = 10;
    hca.YLim = [-20 20];
    drawnow
    legend
    pause(0.1)
    %imovie = sum(doMovie(1:iplot));
    tmp_frame = getframe(gcf);
    %cell_movies{imovie}(itime) = tmp_frame;
    if iA == 1 % initialize animated gif matrix
      [im_tmp,map] = rgb2ind(tmp_frame.cdata,256,'nodither');
      %map(end+1,:) = get(gcf,'color');
      im_tmp(1,1,1,nA_common) = 0;
      cell_movies{imovie,1} = map;
      cell_movies{imovie,2} = im_tmp;
    else
      cell_movies{imovie,2}(:,:,1,iA) = rgb2ind(tmp_frame.cdata,cell_movies{imovie,1},'nodither');
    end       
 end
end
%imwrite(cell_movies{ivar,2},cell_movies{ivar,1},[savedir_gifs 'A_12.gif'],'DelayTime',0.0,'LoopCount',0)

%% Make gif of flux contour A(x,z) rescaled to spentUB(t)
savedir_gifs = '/Users/cno062/Research/PIC/df_cold_protons_04/gifs/';
savedir_gifs = '/Users/cno062/Research/PIC/combination_of_runs/gifs/';
h(1) = subplot(1,3,[1 2]);
h(2) = subplot(1,3,[3]);

set(gcf,'color',[1 1 1]);

x = siminfo_1.x-mean(siminfo_1.x);
z = siminfo_1.z;
nx = 6400;
nz = 1600;
plx = 1:5:nx;
plz = 1:5:nz;
ntimes = 54;
imovie = 1;
cell_movies = cell(1,2);
UB_common = unique([spentUB_1_equi,spentUB_2_equi]);
nUB_common = numel(UB_common);
doholdon = 0;
for iUB = 1:nUB_common
doholdon = 0;
 if 1%doMovie(iplot) % Collect frame for movie    
    currentBackgroundColor = get(gcf,'color');
    Alevels = -25:0.5:0; % still draw these
    hca = h(1);
    if intersect(spentUB_1_equi,UB_common(iUB))
      [C,IA,IB] = intersect(spentUB_1_equi,UB_common(iUB));
      iUBplot = IA;
      displayname = sprintf('n_c=0.8,n_0 tw_{ci} = %.0f',tUB_1(iUBplot));
      contour(hca,-x(plx),z(plz),squeeze(A_1_interp_spentUB(iUBplot,plx,plz))',Alevels,'-','displayname',displayname,'linewidth',1.5);
      doholdon = 1;
    elseif 1% % keep last lines plot, to see if they eventually becomes similar again    
      iUBplot = numel(spentUB_1_equi);
      %iAplot = 12;  
      displayname = sprintf('n_c=0.8n_0, tw_{ci} = %.0f',tUB_1(iUBplot));      
      contour(hca,-x(plx),z(plz),squeeze(A_1_interp_spentUB(iUBplot,plx,plz))',Alevels,'-','displayname',displayname,'linewidth',1.5);      
      doholdon = 1;
    end
    if intersect(spentUB_2_equi,UB_common(iUB))
      [C,IA,IB] = intersect(spentUB_2_equi,UB_common(iUB));
      iUBplot = IA;
      displayname = sprintf('n_c=0.4n_0, tw_{ci} = %.0f',tUB_2(iUB));
      if doholdon; hold(hca,'on'); end
      contour(hca,x(plx),z(plz),squeeze(A_2_interp_spentUB(iUBplot,plx,plz))',Alevels,'-','displayname',displayname);
      if doholdon; hold(hca,'off'); end
    elseif 0 % keep last lines plot, to see if they eventually becomes similar again
      iUBplot = numel(spentUB_2_equi);
      %iAplot = 12;
      displayname = sprintf('n_c=0.4n_0, tw_{ci} = %.0f',tUB_2(iUBplot));
      if doholdon; hold(hca,'on'); end
      contour(hca,x(plx),z(plz),squeeze(A_2_interp_spentUB(iUBplot,plx,plz))',Alevels,'-','displayname',displayname);
      if doholdon; hold(hca,'off'); end
    end
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = sprintf('-dU_{B} = %g ',UB_common(iUB));
    hca.XLim(1) = -100;
    hca.XLim(2) = 100;
    hca.YLim = [-20 20];
    
    if 1 % plot energy vs time plot
      hca = h(2);
      doholdon = 0;
      hlines = plot(hca,timesteps{1}/50,spentUB_1,'-',tUB_1,spentUB_1_equi,'.',...
                          timesteps{2}/50,spentUB_2,'-',tUB_2,spentUB_2_equi,'.'); % looks fine
      hold(hca,'on')
      if intersect(spentUB_1_equi,UB_common(iUB))
        hlines = plot(hca,tUB_1(iUB),spentUB_1_equi(iUB),'+','displayname',sprintf('(n_c=0.4n_0): t(U_B): dU_B = %g',dSpentUB)); % looks fine       
      end
      if intersect(spentUB_2_equi,UB_common(iUB))
        hlines = plot(hca,tUB_2(iUB),spentUB_2_equi(iUB),'+','displayname',sprintf('(n_c=0.8n_0): t(U_B): dU_B = %g',dSpentUB)); % looks fine       
      end
      hold(hca,'off')
      
      hlines(1).LineWidth = 1.5;
      legend(hca,{sprintf('(n_c=0.8n_0): -dU_B(t): dt = %g',(timesteps{1}(2)-timesteps{1}(1))/50),...
                  sprintf('(n_c=0.8n_0): t(U_B): dU_B = %g',dSpentUB),...
                  sprintf('(n_c=0.4n_0): -dU_B(t): dt = %g',(timesteps{1}(2)-timesteps{1}(1))/50),...
                  sprintf('(n_c=0.4n_0): t(U_B): dU_B = %g',dSpentUB)},...
             'location','northwest','box','off');
      hca.XLabel.String = 't (w_{ci})';
      hca.YLabel.String = '-d<U_B> (...)';
      hca.YMinorTick = 'on';
      hca.YMinorGrid = 'on';
    end
    
    drawnow
    legend;
    pause(0.1)
    %imovie = sum(doMovie(1:iplot));
    tmp_frame = getframe(gcf);
    %cell_movies{imovie}(itime) = tmp_frame;
    if iUB == 1 % initialize animated gif matrix
      [im_tmp,map] = rgb2ind(tmp_frame.cdata,256,'nodither');
      %map(end+1,:) = get(gcf,'color');
      im_tmp(1,1,1,nUB_common) = 0;
      cell_movies{imovie,1} = map;
      cell_movies{imovie,2} = im_tmp;
    else
      cell_movies{imovie,2}(:,:,1,iUB) = rgb2ind(tmp_frame.cdata,cell_movies{imovie,1},'nodither');
    end       
 end
end
%imwrite(cell_movies{ivar,2},cell_movies{ivar,1},[savedir_gifs 'A_12.gif'],'DelayTime',0.0,'LoopCount',0)

%% Calculate reconnection rate, 
% Magnetic energy, for plotting against, as complement to time scaling
for it = 1:numel(timesteps{1})
  PB_tmp = sqrt(B_1(it,:,:,1).^2 + B_1(it,:,:,2).^2 + B_1(it,:,:,3).^2)/2;
  U_1.B(it) = nanmean(PB_tmp(:));
end
for it = 1:numel(timesteps{2})
  PB_tmp = sqrt(B_2(it,:,:,1).^2 + B_2(it,:,:,2).^2 + B_2(it,:,:,3).^2)/2;
  U_2.B(it) = nanmean(PB_tmp(:));
end

% A value at X line
[R_1_A,A_1_xline] = dAdt(timesteps{1}/50,A_1,'xline'); 
R_1.A = R_1_A; 
R_1.Ax = A_1_xline;
R_1.UB = U_1.B;
R_1.time = timesteps{1}/50;
[R_2_A,A_2_xline] = dAdt(timesteps{2}/50,A_2,'xline'); 
R_2.A = R_2_A; 
R_2.Ax = A_2_xline;
R_2.UB = U_2.B;
R_2.time = timesteps{2}/50;

if 1 % plot
  %%
  h = setup_subplots(2,1);
  isub = 1;
  
  if 1 % reconnection rate based on A
    hca = h(isub); isub = isub + 1;
    plot(hca,timesteps{1}/50,R_1.A,timesteps{2}/50,R_2.A)
    hca.XLabel.String = 't (w_{ci})';
    hca.YLabel.String = 'R (v_{A0}B_0)';
    legend(hca,{'n_c = 0.8 n_0','n_c = 0.4 n_0'},'location','best')
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.YLim(1) = 0;
  end
  if 1 % A_xline
    hca = h(isub); isub = isub + 1;
    plot(hca,timesteps{1}/50,A_1_xline,timesteps{2}/50,A_2_xline)
    hca.XLabel.String = 't (w_{ci})';
    hca.YLabel.String = 'A^X (d_iB_0)';
    legend(hca,{'n_c = 0.8 n_0','n_c = 0.4 n_0'},'location','best')
    hca.XGrid = 'on';
    hca.YGrid = 'on';    
  end
  if 0 % Magnetic energy
    hca = h(isub); isub = isub + 1;
    plot(hca,timesteps{1}/50,U_1.B,timesteps{2}/50,U_2.B)
    hca.XLabel.String = 't (w_{ci})';
    hca.YLabel.String = '<U_B> (...)';
    legend(hca,{'n_c = 0.8 n_0','n_c = 0.4 n_0'},'location','best')
  end
  if 0 % Magnetic energy
    hca = h(isub); isub = isub + 1;
    plot(hca,timesteps{1}(2:end)/50/1,-1.0*diff(U_1.B),timesteps{2}(2:end)/50,-diff(U_2.B))
    hca.XLabel.String = 't (w_{ci})';
    hca.YLabel.String = '-\Delta <U_B> (...)';
    legend(hca,{'n_c = 0.8 n_0','n_c = 0.4 n_0'},'location','best')
  end
  if 0 % Reconnectionr ate vs magnetic energy in system
    hca = h(isub); isub = isub + 1;
    plot(hca,U_1.B,R_1.A,U_2.B,R_2.A)
    hca.YLabel.String = 'R (v_{A0}B_0)';
    hca.XLabel.String = 'U_B (...)';
    legend(hca,{'n_c = 0.8 n_0','n_c = 0.4 n_0'},'location','best')
  end
  if 0 % Reconnectionr ate vs magnetic energy in system
    hca = h(isub); isub = isub + 1;
    plot(hca,-diff(U_1.B),R_1.A(2:end),-diff(U_2.B),R_2.A(2:end))
    hca.YLabel.String = 'R (v_{A0}B_0)';
    hca.XLabel.String = '-\Delta U_B (...)';
    legend(hca,{'n_c = 0.8 n_0','n_c = 0.4 n_0'},'location','best')
  end
  if 0 % 
    hca = h(isub); isub = isub + 1;
    plot(hca,-diff(U_1.B),R_1.Ax(2:end),-diff(U_2.B),R_2.Ax(2:end))
    hca.YLabel.String = 'A^x (d_{i}B_0)';
    hca.XLabel.String = '-\Delta U_B (...)';
    legend(hca,{'n_c = 0.8 n_0','n_c = 0.4 n_0'},'location','best')
  end
  if 0 % 
    hca = h(isub); isub = isub + 1;
    plot(hca,(R_1.Ax(2:end)-R_1.Ax(1))*1,-diff(U_1.B)*1,R_2.Ax(2:end)-R_2.Ax(1),-diff(U_2.B))
    hca.XLabel.String = 'A^x-A^{x}(1) (d_{i}B_0)';
    hca.YLabel.String = '-\Delta U_B (...)';
    legend(hca,{'n_c = 0.8 n_0','n_c = 0.4 n_0'},'location','best')
  end
  
end

%% Ey at X line
R_1.Ey = zeros(54,1);
R_2.Ey = zeros(54,1);
for it = 1:54    
  [saddle_locations,saddle_values] = saddle(squeeze(A_1(it,:,:)),'sort');
  saddle_1_x = saddle_locations(1,1); % indices
  saddle_1_y = saddle_locations(1,2);
  R_1.Ey(it) = E_1(it,saddle_1_x,saddle_1_y,2);
  
%   [saddle_locations,saddle_values] = saddle(squeeze(A_2(it,:,:)),'sort');
%   saddle_2_x = saddle_locations(1,1);
%   saddle_2_y = saddle_locations(1,2);
%   R_2.Ey = E_2(it,saddle_1_x,saddle_2_y,2);
  
  
end
%%
% partial time derivative
dAdt(1,:,:) = (A(2,:,:)-A(1,:,:))/dt;
dAdt(2:end-1,:,:) = (A(3:end,:,:)-A(1:end-2,:,:))/(2*dt);
dAdt(end,:,:) = (A(end,:,:)-A(end-1,:,:))/dt;

%% Calculate flux speed
x = siminfo_1.x - mean(siminfo_1.x);
z = siminfo_1.z;
tic; [vAx_1,vAz_1] = flux_speed(timesteps{1}/50,siminfo_1.x,siminfo_1.z,A_1); vA_1.x = vAx_1; vA_1.z = vAz_1; toc
tic; [vAx_2,vAz_2] = flux_speed(timesteps{2}/50,siminfo_2.x,siminfo_2.z,A_2); vA_2.x = vAx_2; vA_2.z = vAz_2; toc
tic; [aAx_1,aAz_1] = flux_acc(timesteps{1}/50,siminfo_1.x,siminfo_1.z,vA_1); aA_1.x = aAx_1; aA_1.z = aAz_1; toc
tic; [aAx_2,aAz_2] = flux_acc(timesteps{2}/50,siminfo_2.x,siminfo_2.z,vA_2); aA_2.x = aAx_2; aA_2.z = aAz_2; toc

if 1 % plot flux velocity at neutral plane as a function of time, i.e. vA_1,2(x,t)
  %%
  npanels = 3;
  h = setup_subplots(npanels,1);
  isub = 1;
  zind = 801+10*[-1,1];
  if 1
    hca = h(isub); isub = isub + 1;
    imagesc(hca,x,timesteps{1}/50,squeeze(mean(vA_1.x(:,:,zind),3))); 
    hcb = colorbar('peer',hca); 
    hca.CLim = [-2 2];
    colormap(pic_colors('blue_red'))
    hca.YDir = 'normal';
    %hold(hca,'on')
    %contour(hca,x,timesteps{2}/50,squeeze(mean(vA_2.x(:,:,zind),3))')
    %hold(hca,'off')
    hca.Title.String = 'n_c = 0.8n_0';
    hcb.YLabel.String = 'v_{A,x}(x,z~0,t)';
    hold(hca,'on')
    [C,h_] = contour(hca,x,timesteps{1}/50,squeeze(mean(A_1(:,:,zind),3)),-25:1:0,'k'); 
    clabel(C,h_)
    hold(hca,'off')
  end
  if 1
    hca = h(isub); isub = isub + 1;
    imagesc(hca,x,timesteps{2}/50,squeeze(mean(vA_2.x(:,:,zind),3))); 
    hcb = colorbar('peer',hca); 
    hca.CLim = [-2 2];
    colormap(pic_colors('blue_red'))
    hca.YDir = 'normal';
    hca.Title.String = 'n_c = 0.4n_0';
    hcb.YLabel.String = 'v_{A,x}(x,z~0,t)';
    hold(hca,'on')
    [C,h_] = contour(hca,x,timesteps{2}/50,squeeze(mean(A_2(:,:,zind),3)),-25:1:0,'k'); 
    clabel(C,h_)
    hold(hca,'off')
  end
  if 1
    hca = h(isub); isub = isub + 1;
    %vAx_diff = squeeze(mean(vA_2.x(1:54,:,zind),3))-flipdim(-squeeze(mean(vA_1.x(1:54,:,zind),3)),2);
    vscale2 = 0.95;
    vAx_diff = vscale2*squeeze(mean(vA_2.x(1:54,:,zind),3))-squeeze(mean(vA_1.x(1:54,:,zind),3));
    imagesc(hca,x,timesteps{2}(1:54)/50,vAx_diff); 
    hcb = colorbar('peer',hca); hca.CLim = [-2 2];
    colormap(pic_colors('blue_red'))
    hca.YDir = 'normal';
    hca.Title.String = sprintf('%.2fv_A(n_c = 0.4n_0)-v_A(n_c = 0.8n_0)',vscale2);
    hcb.YLabel.String = 'v_{A,x}(x,z~0,t)';
  end
  
  for ipanel = 1:npanels
    h(ipanel).CLim = [-1 1];
    h(ipanel).XLabel.String = 'x (d_i)';
    h(ipanel).YLabel.String = 't (w_{ci})';
    h(ipanel).Title.String = {h(ipanel).Title.String,sprintf('z = [%.2f,%.2f]',z(zind(1)),z(zind(end)))};
  end
  
end

if 1 % plot flux acceleration at neutral plane as a function of time, i.e. vA_1,2(x,t)
  %%
  npanels = 2;
  h = setup_subplots(npanels,1);
  isub = 1;
  zind = 801+15*[-1,1];
  if 1
    hca = h(isub); isub = isub + 1;
    tt_diff = timesteps{1}(2:end)/50-(timesteps{1}(2)-timesteps{1}(1))/50;
    imagesc(hca,x,tt_diff,squeeze(mean(diff(vA_1.x(:,:,zind),1),3))); 
    %imagesc(hca,x,tt_diff,squeeze(diff(mean(vA_1.x(:,:,zind),3),1))); 
    hcb = colorbar('peer',hca); 
    hca.CLim = [-2 2];
    colormap(pic_colors('blue_red'))
    hca.YDir = 'normal';
    %hold(hca,'on')
    %contour(hca,x,timesteps{2}/50,squeeze(mean(vA_2.x(:,:,zind),3))')
    %hold(hca,'off')
    hca.Title.String = 'n_c = 0.8n_0';
    hcb.YLabel.String = 'dt*d(v_{A,x}(x,z~0,t))/dt';
    hold(hca,'on')
    [C,h_] = contour(hca,x,timesteps{1}/50,squeeze(mean(A_1(:,:,zind),3)),-25:1:0,'k'); 
    clabel(C,h_)
    hold(hca,'off')
  end
  if 1
    hca = h(isub); isub = isub + 1;
    tt_diff = timesteps{2}(2:end)/50-(timesteps{2}(2)-timesteps{2}(1))/50;
    imagesc(hca,x,tt_diff,squeeze(mean(diff(vA_2.x(:,:,zind),1),3))); 
    hcb = colorbar('peer',hca); 
    hca.CLim = [-2 2];
    colormap(pic_colors('blue_red'))
    hca.YDir = 'normal';
    hca.Title.String = 'n_c = 0.4n_0';
    hcb.YLabel.String = 'dt*d(v_{A,x}(x,z~0,t))/dt';
    hold(hca,'on')
    [C,h_] = contour(hca,x,timesteps{2}/50,squeeze(mean(A_2(:,:,zind),3)),-25:1:0,'k'); 
    clabel(C,h_)
    hold(hca,'off')
  end
  if 0
    hca = h(isub); isub = isub + 1;
    %vAx_diff = squeeze(mean(vA_2.x(1:54,:,zind),3))-flipdim(-squeeze(mean(vA_1.x(1:54,:,zind),3)),2);
    vscale2 = 0.95;
    vAx_diff = vscale2*squeeze(mean(vA_2.x(1:54,:,zind),3))-squeeze(mean(vA_1.x(1:54,:,zind),3));
    imagesc(hca,x,timesteps{2}(1:54)/50,vAx_diff); 
    hcb = colorbar('peer',hca); hca.CLim = [-2 2];
    colormap(pic_colors('blue_red'))
    hca.YDir = 'normal';
    hca.Title.String = sprintf('%.2fv_A(n_c = 0.4n_0)-v_A(n_c = 0.8n_0)',vscale2);
    hcb.YLabel.String = 'dt*d(v_{A,x}(x,z~0,t))/dt';
  end
  
  for ipanel = 1:npanels
    h(ipanel).CLim = 0.5*[-1 1];
    h(ipanel).XLim = 100*[-1 1];
    h(ipanel).XLabel.String = 'x (d_i)';
    h(ipanel).YLabel.String = 't (w_{ci})';
    h(ipanel).Title.String = {h(ipanel).Title.String,sprintf('z = [%.2f,%.2f]',z(zind(1)),z(zind(end)))};
  end
  hlink = linkprop(h,{'XLim','YLim'});
  
end
%%
if 1 % plot arrows vs x,z
  %%
  h = setup_subplots(1,1);
  isub = 1;
  
  x = siminfo_2.x - mean(siminfo_2.x);
  z = siminfo_2.z;

  colors = pic_colors('matlab');

  xlim = [x(1) x(end)]; xlim = [-20 20];
  zlim = [z(1) z(end)]; zlim = 5*[-1 1];
  ipx1 = find(x>xlim(1),1,'first'); ipx2 = find(x<xlim(2),1,'last');
  ipz1 = find(z>zlim(1),1,'first'); ipz2 = find(z<zlim(2),1,'last');
  ipx = ipx1:3:ipx2;
  ipz = ipz1:3:ipz2;

  % Quivers
  doQ = 0;
  nQx = 50;
  nQz = 50;
  [Z,X] = meshgrid(z,x);
  ipxQ = fix(linspace(ipx1,ipx2,nQx));
  ipzQ = fix(linspace(ipz1,ipz2,nQz));
  [dataQx,dataQz] = meshgrid(ipxQ,ipzQ);
  ipXQ = dataQx; ipZQ = dataQz;
  % dataQ.x = E.perp.x;
  % dataQ.y = E.perp.y;
  % dataQ.z = E.perp.z;
  dataQ = vA_2;
  maxQ = 2.5;
  dataQ.abs = sqrt(dataQ.x.^2 + dataQ.z.^2);
  dataQ.x(dataQ.abs>maxQ) = NaN;
  dataQ.y(dataQ.abs>maxQ) = NaN;
  dataQ.z(dataQ.abs>maxQ) = NaN;

  tind = 25;
  %xind = 1:3:6400;
  %zind = 1:3:1600;
  [X,Z] = ndgrid(x,z);
  if 1
    hca = h(isub); isub = isub + 1;
    doQ = 1;
    Alevels = -25:0.1:0;
    contour(hca,x(ipx),z(ipz),squeeze(A_2(tind,ipx,ipz))',Alevels,'-','color',[0.7 0.7 0.7],'displayname','A')    
    hca.XLabel.String = 'x (d_{pi})';
    hca.YLabel.String = 'z (d_{pi})';
    
    vquiverscaling = 3;
    if doQ
      hold(hca,'on')
      hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),squeeze(dataQ.x(tind,ipxQ,ipzQ))*vquiverscaling/10,squeeze(dataQ.z(tind,ipxQ,ipzQ))*vquiverscaling/10,0,'displayname','0.1v_{A}');
      hold(hca,'off')  
    end
    if 1 % doQ1
      hold(hca,'on')
      dataQ = ve2;
      hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),squeeze(dataQ.x(ipxQ,ipzQ))*vquiverscaling,squeeze(dataQ.z(ipxQ,ipzQ))*vquiverscaling,0,'displayname','v_{e2}');
      dataQ = ve3;
      hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),squeeze(dataQ.x(ipxQ,ipzQ))*vquiverscaling,squeeze(dataQ.z(ipxQ,ipzQ))*vquiverscaling,0,'displayname','v_{e3}');
      dataQ.x = (ne2.*ve2.x + ne3.*ve3.x)./(ne2+ne3);
      dataQ.z = (ne2.*ve2.z + ne3.*ve3.z)./(ne2+ne3);
      dataQ.x((ne2+ne3)<0.2) = NaN;
      dataQ.z((ne2+ne3)<0.2) = NaN;
      hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),squeeze(dataQ.x(ipxQ,ipzQ))*vquiverscaling,squeeze(dataQ.z(ipxQ,ipzQ))*vquiverscaling,0,'displayname','v_{e23}');
      hold(hca,'off')  
    end
    axis(hca,'equal')
    legend(hca)
  end
  xind = 1:1:6400;
  zind = 800;
  if 0
    hca = h(isub); isub = isub + 1;
    Alevels = -25:1:0;
    tscaling = 1;
    xscaling = 1;
    contour(hca,x(xind)*xscaling,timesteps{1}/50*tscaling,squeeze(A_1(:,xind,zind)),Alevels,'--')
    hold(hca,'on')
    contour(hca,x(xind),timesteps{2}/50,squeeze(A_2(:,xind,zind)),Alevels,'-')
    hold(hca,'off')
    hca.YLabel.String = 't (\omega_{ci})';
    hca.XLabel.String = 'z (d_{pi})';
    hca.XLim = [-100 100];
  end
end

%%

h = setup_subplots(1,2);
isub = 1;

x = siminfo_1.x - mean(siminfo_1.x);
z = siminfo_1.z;

xind = 3200;
zind = 1:1:1600;
if 1
  hca = h(isub); isub = isub + 1;
  Alevels = -25:1:0;
  contour(hca,timesteps{1}/50,z(zind),squeeze(A_1(:,xind,zind))',Alevels,'--','displayname','nc=0.8')
  hold(hca,'on')
  contour(hca,timesteps{2}/50,z(zind),squeeze(A_2(:,xind,zind))',Alevels,'-','displayname','nc=0.4')
  hold(hca,'off')
  hca.XLabel.String = 't (\omega_{ci})';
  hca.YLabel.String = 'z (d_{pi})';
  legend(hca)
end

xind = 1:1:6400;
zind = 800;
if 1
  hca = h(isub); isub = isub + 1;
  Alevels = -25:1:0;
  tscaling =1.1;
  xscaling = 1.5;
  contour(hca,x(xind)*xscaling,timesteps{1}/50*tscaling,squeeze(A_1(:,xind,zind)),Alevels,'--','displayname','nc=0.8')
  hold(hca,'on')
  contour(hca,x(xind),timesteps{2}/50,squeeze(A_2(:,xind,zind)),Alevels,'-','displayname','nc=0.4')
  hold(hca,'off')
  hca.YLabel.String = 't (\omega_{ci})';
  hca.XLabel.String = 'z (d_{pi})';
  hca.XLim = [-100 100];
  legend(hca)
end
if 0
  hca = h(isub); isub = isub + 1;
  Alevels = -25:1:0;
  contour(hca,timesteps{1}/50,x(xind),squeeze(A_1(:,xind,zind))',Alevels,'--')
  hold(hca,'on')
  contour(hca,timesteps{2}/50,x(xind),squeeze(A_2(:,xind,zind))',Alevels,'-')
  hold(hca,'off')
  hca.XLabel.String = 't (\omega_{ci})';
  hca.YLabel.String = 'z (d_{pi})';
end

