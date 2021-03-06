%% Load PIC object
no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
ox = no02m;

%% Make/load trajectories and interpolate data
traj = 1;
switch traj
  case 1 % % top left quadrant only, single time
    
    twpe = 20000;
    twpe1 = twpe;
    twpe2 = twpe;
    
    i_traj_old = [1 2 3];
    x_traj_old = [80 90 95];
    z_traj_old = [8 2 6];
    t_traj_old = ones(1,numel(x_traj));
    
    % interpolate to more trajectory points    
    n_traj_new = 200;
    i_traj_new = linspace(i_traj_old(1),i_traj_old(end),n_traj_new);        
    x_traj = interp1(i_traj_old,x_traj_old,i_traj_new,'spline');
    z_traj = interp1(i_traj_old,z_traj_old,i_traj_new,'spline');
    
  case 6
    %%
    twpe1 = 8000;
    twpe2 = 10000;    
    t_traj = [twpe1:200:twpe2];
    x_traj = [220 225 235 240 242 243 244 245 246 247 248];
    z_traj = [6 3 3 3 3 3 3.5 4 4.5 5.0 5.5];
    
    % interpolate to more trajectory points below plot
       
    if 1 % Plot density for each time so I can select points from it
      %%
      twpelim = [twpe1 twpe2];
      zlim = [-1 8];
      xlim = [190 280];
      pic_ = ox.xlim(xlim).zlim(zlim).twpelim(twpelim);
      nO_tmp = pic_.n(3);
      
      figure(21)
      ncols = 2;
      nrows = ceil(pic_.nt/ncols);
      npanels = nrows*ncols;

      h = setup_subplots(nrows,ncols,'vertical');
      isub = 1;

      for it = 1:pic_.nt
        pic_tmp = pic_(it);
        pic_tmp.twpe;
        if 1 % n
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(nO_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,sprintf('tw_{pe}=%g',pic_tmp.twpe),[0.02 0.98])
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
      end
      %hlinks = linkprop(h,{'XLim','YLim','CLim'});
      %hlinks.Targets(1)
      compact_panels(h,0,0.02)
    end    
        
    t_traj_old = t_traj;
    t_traj = twpe1:1:twpe2;
    x_traj = interp1(t_traj_old,x_traj,t_traj,'spline');
    z_traj = interp1(t_traj_old,z_traj,t_traj,'spline');
  case 2 % a bit furhter out
    %%
    twpe1 = 8000;
    twpe2 = 10000;    
    t_traj = [twpe1:200:twpe2];
    x_traj = [215 225 240 241 241 241 242 243 244 245 246];
    z_traj = [6 3 3 4 4 4 4 4 6 8.0 10];
    
    % interpolate to more trajectory points below plot
       
    if 1 % Plot density for each time so I can select points from it
      %%
      twpelim = [twpe1 twpe2];
      zlim = [-1 9];
      xlim = [190 280];
      pic_ = ox.xlim(xlim).zlim(zlim).twpelim(twpelim);
      nO_tmp = pic_.n(3);
      
      figure(21)
      ncols = 2;
      %nrows = ceil(pic_.nt/ncols);
      nrows = pic_.nt; % plot second column in spacecraft frame.
      npanels = nrows*ncols;

      %h = setup_subplots(nrows,ncols,'vertical');
      h = setup_subplots(nrows,ncols,'horizontal');
      isub = 1;

      for it = 1:pic_.nt
        pic_tmp = pic_(it);
        pic_tmp.twpe;
        if 1 % n
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(nO_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 1 % n, spacecraft frame
          hca = h(isub); isub = isub + 1;
          pcolor(hca,pic_tmp.xi-0*x_traj(it),pic_tmp.zi-0*z_traj(it),squeeze(nO_tmp(:,:,pic_tmp.it))');
          shading(hca,'flat')
          hca.XLim = x_traj(it) + 0.5*diff(xlim)*[-1 1];
          hca.YLim = z_traj(it) + 0.5*diff(zlim)*[-1 1];
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
      end
      %hlinks = linkprop(h,{'XLim','YLim','CLim'});
      %hlinks.Targets(1)
      compact_panels(h,0,0.02)
    end    
        
    t_traj_old = t_traj;
    dt = 1;
    t_traj = twpe1:dt:twpe2;
    
    x_traj = interp1(t_traj_old,x_traj,t_traj,'spline');
    z_traj = interp1(t_traj_old,z_traj,t_traj,'spline');
    vx_traj = diff(x_traj)/(dt/50);
    vz_traj = diff(z_traj)/(dt/50);
  case 3 % down-up-down-up
    %%
    twpe1 = 8000;
    twpe2 = 10000;    
    t_traj = [twpe1:200:twpe2];
    x_traj = [220 230 240 241 241 241 242 243 244 245 246];
    z_traj = [6 2 3 6 5 4 4 4 6 8.0 10];
    
    % shifted a bit more to later times
    x_traj = [220 225 230 240 243 245 246 247 248 249 250 245 246];
    z_traj = [6 5 3 2 4 5 6 3 6 8 9 10];
    
    x_traj = x_traj(1:numel(t_traj));
    z_traj = z_traj(1:numel(t_traj));
    
    % interpolate to more trajectory points below plot
       
    if 1 % Plot density for each time so I can select points from it
      %%
      twpelim = [twpe1 twpe2];
      zlim = [-1 9];
      xlim = [190 280];
      pic_ = ox.xlim(xlim).zlim(zlim).twpelim(twpelim);
      nO_tmp = pic_.n(3);
      vzO_tmp = pic_.vz(3);
      vex_tmp = pic_.vx([2 4]);
      
      figure(21)
      ncols = 3;
      %nrows = ceil(pic_.nt/ncols);
      nrows = pic_.nt; % plot second column in spacecraft frame.
      npanels = nrows*ncols;

      %h = setup_subplots(nrows,ncols,'vertical');
      h = setup_subplots(nrows,ncols,'horizontal');
      isub = 1;

      for it = 1:pic_.nt
        pic_tmp = pic_(it);
        pic_tmp.twpe;
        if 1 % vex
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(vex_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [-1 1];
          colormap(hca,pic_colors('blue_red'))
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 1 % vz
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(vzO_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [-0.4 0.4];
          colormap(hca,pic_colors('blue_red'))
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 1 % n
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(nO_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 0 % n, spacecraft frame
          hca = h(isub); isub = isub + 1;
          pcolor(hca,pic_tmp.xi-0*x_traj(it),pic_tmp.zi-0*z_traj(it),squeeze(nO_tmp(:,:,pic_tmp.it))');
          shading(hca,'flat')
          hca.XLim = x_traj(it) + 0.5*diff(xlim)*[-1 1];
          hca.YLim = z_traj(it) + 0.5*diff(zlim)*[-1 1];
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')            
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
      end
      %hlinks = linkprop(h,{'XLim','YLim','CLim'});
      %hlinks.Targets(1)
      compact_panels(h,0,0.02)
    end    
        
    t_traj_old = t_traj;
    dt = 1;
    t_traj = twpe1:dt:twpe2;
    
    x_traj = interp1(t_traj_old,x_traj,t_traj,'spline');
    z_traj = interp1(t_traj_old,z_traj,t_traj,'spline');
    vx_traj = diff(x_traj)/(dt/50);
    vz_traj = diff(z_traj)/(dt/50);
    gongStruct = load('chirp.mat');
    %sound(gongStruct.y,gongStruct.Fs)
  case 66 % down-up-down-up
    %%
    twpe1 = 8000;
    twpe2 = 10000;    
    t_traj = [twpe1:200:twpe2];
    x_traj = [220 230 240 241 241 241 242 243 244 245 246];
    z_traj = [6 2 3 6 5 4 4 4 6 8.0 10];
    
    % shifted a bit more to later times
    x_traj = [220 225 230 240 243 245 246 247 248 249 250 245 246];
    x_traj = x_traj(end:-1:1);
    x_traj = [240 240 240 240 240 240 240 240 240 240 240 240 240];
    x_traj = x_traj + linspace(0,1,numel(x_traj));
    z_traj = [6 5 3 2 2 2 3 4 5 6 7 8];
    z_traj = [8 6 5 4 3 2 3 4 5 6 7 9];
    z_traj = [6 5.5 5 4.5 4 4 4 4.5 5 5.5 6 6.5];
    
    x_traj = x_traj(1:numel(t_traj));
    z_traj = z_traj(1:numel(t_traj));
    
    % interpolate to more trajectory points below plot
       
    if 1 % Plot density for each time so I can select points from it
      %%
      twpelim = [twpe1 twpe2];
      zlim = [-1 9];
      xlim = [190 280];
      pic_ = ox.xlim(xlim).zlim(zlim).twpelim(twpelim);
      nO_tmp = pic_.n(3);
      vzO_tmp = pic_.vz(3);
      vex_tmp = pic_.vx([2 4]);
      
      figure(21)
      ncols = 3;
      %nrows = ceil(pic_.nt/ncols);
      nrows = pic_.nt; % plot second column in spacecraft frame.
      npanels = nrows*ncols;

      %h = setup_subplots(nrows,ncols,'vertical');
      h = setup_subplots(nrows,ncols,'horizontal');
      isub = 1;

      for it = 1:pic_.nt
        pic_tmp = pic_(it);
        pic_tmp.twpe;
        if 1 % vex
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(vex_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [-1 1];
          colormap(hca,pic_colors('blue_red'))
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 1 % vz
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(vzO_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [-0.4 0.4];
          colormap(hca,pic_colors('blue_red'))
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 1 % n
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(nO_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 0 % n, spacecraft frame
          hca = h(isub); isub = isub + 1;
          pcolor(hca,pic_tmp.xi-0*x_traj(it),pic_tmp.zi-0*z_traj(it),squeeze(nO_tmp(:,:,pic_tmp.it))');
          shading(hca,'flat')
          hca.XLim = x_traj(it) + 0.5*diff(xlim)*[-1 1];
          hca.YLim = z_traj(it) + 0.5*diff(zlim)*[-1 1];
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')            
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
      end
      %hlinks = linkprop(h,{'XLim','YLim','CLim'});
      %hlinks.Targets(1)
      compact_panels(h,0,0.02)
    end    
        
    t_traj_old = t_traj;
    dt = 1;
    t_traj = twpe1:dt:twpe2;
    
    x_traj = interp1(t_traj_old,x_traj,t_traj,'spline');
    z_traj = interp1(t_traj_old,z_traj,t_traj,'spline');
    vx_traj = diff(x_traj)/(dt/50);
    vz_traj = diff(z_traj)/(dt/50);
    %gongStruct = load('chirp.mat');
    %sound(gongStruct.y,gongStruct.Fs)
  case 77 % down-up-down-up
    %%
    twpe1 = 8000;
    twpe2 = 10000;    
    t_traj = [twpe1:200:twpe2];
    x_traj = 240*ones(size(t_traj));
    z_traj = [8 6 4 3 2 2 3 4 6 7 8 10 11];
    z_traj = [7 6 5 4 3 3 3 4 5 6 6 6 6];  
    z_traj = [5.75 5.50 5.25 5.00 4.75 4.5 4.5 4.75 5.00 5.25 5.50 5.75 6];
    
    
    x_traj = x_traj(1:numel(t_traj));
    z_traj = z_traj(1:numel(t_traj));
    
    diff(z_traj)/(200/50)
    
    % interpolate to more trajectory points below plot
       
    if 1 % Plot density for each time so I can select points from it
      %%
      twpelim = [twpe1 twpe2];
      zlim = [-1 9];
      xlim = [190 280];
      xlim = [210 260];
      pic_ = ox.xlim(xlim).zlim(zlim).twpelim(twpelim);
      nO_tmp = pic_.n(3);
      vzO_tmp = pic_.vz(3);
      vex_tmp = pic_.vx([2 4]);
      
      figure(21)
      ncols = 1;
      %nrows = ceil(pic_.nt/ncols);
      nrows = pic_.nt; % plot second column in spacecraft frame.
      npanels = nrows*ncols;

      %h = setup_subplots(nrows,ncols,'vertical');
      h = setup_subplots(nrows,ncols,'horizontal');
      isub = 1;

      for it = 1:pic_.nt
        pic_tmp = pic_(it);
        pic_tmp.twpe;
        if 0 % vex
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(vex_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [-1 1];
          colormap(hca,pic_colors('blue_red'))
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 0 % vz
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(vzO_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [-0.4 0.4];
          colormap(hca,pic_colors('blue_red'))
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 1 % n
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(nO_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 0 % n, spacecraft frame
          hca = h(isub); isub = isub + 1;
          pcolor(hca,pic_tmp.xi-0*x_traj(it),pic_tmp.zi-0*z_traj(it),squeeze(nO_tmp(:,:,pic_tmp.it))');
          shading(hca,'flat')
          hca.XLim = x_traj(it) + 0.5*diff(xlim)*[-1 1];
          hca.YLim = z_traj(it) + 0.5*diff(zlim)*[-1 1];
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')            
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
      end
      %hlinks = linkprop(h,{'XLim','YLim','CLim'});
      %hlinks.Targets(1)
      compact_panels(h,0,0.02)
    end    
        
    t_traj_old = t_traj;
    dt = 10;
    t_traj = twpe1:dt:twpe2;
    
    x_traj = interp1(t_traj_old,x_traj,t_traj,'spline');
    z_traj = interp1(t_traj_old,z_traj,t_traj,'spline');
    %x_traj = interp1(t_traj_old,x_traj,t_traj,'linear');
    %z_traj = interp1(t_traj_old,z_traj,t_traj,'linear');
    
    vx_traj = diff(x_traj)/(dt/50);
    vz_traj = diff(z_traj)/(dt/50);
    t_traj_middle = t_traj(1:end-1)+0.5*dt;
    vx_traj = interp1(t_traj_middle,vx_traj,t_traj);
    vz_traj = interp1(t_traj_middle,vz_traj,t_traj);
    %gongStruct = load('chirp.mat');
    %sound(gongStruct.y,gongStruct.Fs)
  case 4 % 'almost single time', two fronts, 10000
    %%
    twpe1 = 9600;
    twpe2 = 10000;    
    t_traj = [twpe1:200:twpe2];
    x_traj = [220 240 270];
    z_traj = [6 2 10];
    
    % interpolate to more trajectory points below plot
       
    if 1 % Plot density for each time so I can select points from it
      %%
      twpelim = [twpe1 twpe2];
      zlim = [-1 9];
      xlim = [190 280];
      pic_ = ox.xlim(xlim).zlim(zlim).twpelim(twpelim);
      nO_tmp = pic_.n(3);
      vzO_tmp = pic_.vz(3);
      
      figure(21)
      ncols = 3;
      %nrows = ceil(pic_.nt/ncols);
      nrows = pic_.nt; % plot second column in spacecraft frame.
      npanels = nrows*ncols;

      %h = setup_subplots(nrows,ncols,'vertical');
      h = setup_subplots(nrows,ncols,'horizontal');
      isub = 1;

      for it = 1:pic_.nt
        pic_tmp = pic_(it);
        pic_tmp.twpe;
        if 1 % vz
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(vzO_tmp(:,:,it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [-0.4 0.4];
          colormap(hca,pic_colors('blue_red'))
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 1 % n
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(nO_tmp(:,:,it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 1 % n, spacecraft frame
          hca = h(isub); isub = isub + 1;
          pcolor(hca,pic_tmp.xi,pic_tmp.zi,squeeze(nO_tmp(:,:,it))');
          shading(hca,'flat')
          hca.XLim = x_traj(it) + 0.5*diff(xlim)*[-1 1];
          hca.YLim = z_traj(it) + 0.5*diff(zlim)*[-1 1];
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj,z_traj,'k')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
      end
      %hlinks = linkprop(h,{'XLim','YLim','CLim'});
      %hlinks.Targets(1)
      compact_panels(h,0,0.02)
    end    
        
    t_traj_old = t_traj;
    dt = 1;
    t_traj = twpe1:dt:twpe2;
    
    x_traj = interp1(t_traj_old,x_traj,t_traj,'spline');
    z_traj = interp1(t_traj_old,z_traj,t_traj,'spline');
    vx_traj = diff(x_traj)/(dt/50);
    vz_traj = diff(z_traj)/(dt/50);
    gongStruct = load('chirp.mat');
    %sound(gongStruct.y,gongStruct.Fs)
  case 'hakon2' % a bit further out
    %%
    twci_traj_all = load(['/Users/' localuser ,'/MATLAB/pic-matlab/hakon/tList.txt']);
    x_traj_all = load(['/Users/' localuser ,'/MATLAB/pic-matlab/hakon/Xtraj810.txt']);
    z_traj_all = load(['/Users/' localuser ,'/MATLAB/pic-matlab/hakon/Ztraj810.txt']);

    dt = diff(twci_traj_all);
%     for ii = 2:numel(twci_traj_all)
%       dt = twci_traj_all(ii)-twci_traj_all(ii-1);
%       if dt == 1
%         
%       end
%       
%     end
    %x_traj_all = ox.xi(ix_traj_all);
    %z_traj_all = ox.zi(iz_traj_all);
        
    t_traj = 8000:200:10000;    
    twpe_traj_all = twci_traj_all*50;
    x_traj = interp1(twpe_traj_all,x_traj_all,t_traj);
    z_traj = interp1(twpe_traj_all,z_traj_all,t_traj);
    
    if 1
      %%
      figure(33)
      [AX,H1,H2] = plotyy(twpe_traj_all,x_traj_all,twpe_traj_all,z_traj_all);
      hold(AX(1),'on')
      plot(AX(1),t_traj,x_traj,'*')
      hold(AX(1),'off')      
      hold(AX(2),'on')
      plot(AX(2),t_traj,z_traj,'*')
      hold(AX(2),'off')
    
    
    end
    %x_traj = [250 251 252.3];
    %z_traj = [];
    
    %t_traj = [twpe1:200:twpe2];
    %x_traj = [215 225 240 241 241 241 242 243 244 245 246];
    %z_traj = [6 3 3 4 4 4 4 4 6 8.0 10];
    
    % interpolate to more trajectory points below plot
       
    if 1 % Plot density for each time so I can select points from it
      %%
      twpelim = [twpe1 twpe2];
      zlim = [-1 9];
      xlim = [190 280];
      pic_ = ox.xlim(xlim).zlim(zlim).twpelim(twpelim);
      nO_tmp = pic_.n(3);
      
      figure(21)
      ncols = 2;
      %nrows = ceil(pic_.nt/ncols);
      nrows = pic_.nt; % plot second column in spacecraft frame.
      npanels = nrows*ncols;

      %h = setup_subplots(nrows,ncols,'vertical');
      h = setup_subplots(nrows,ncols,'horizontal');
      isub = 1;

      for it = 1:pic_.nt
        pic_tmp = pic_(it);
        pic_tmp.twpe;
        if 1 % n
          hca = h(isub); isub = isub + 1;
          imagesc(hca,pic_tmp.xi,pic_tmp.zi,squeeze(nO_tmp(:,:,pic_tmp.it))');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
        if 1 % n, spacecraft frame
          hca = h(isub); isub = isub + 1;
          pcolor(hca,pic_tmp.xi-0*x_traj(it),pic_tmp.zi-0*z_traj(it),squeeze(nO_tmp(:,:,pic_tmp.it))');
          shading(hca,'flat')
          hca.XLim = x_traj(it) + 0.5*diff(xlim)*[-1 1];
          hca.YLim = z_traj(it) + 0.5*diff(zlim)*[-1 1];
          hca.YDir = 'normal';
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'z';          
          irf_legend(hca,{sprintf('tw_{pe}=%g',pic_tmp.twpe),sprintf('x=%0.1f, z=%0.1f',x_traj(it),z_traj(it))},[0.02 0.98],'color','k')
          hca.CLim = [0 0.4];
          if 1 % plot trajectory positions for each time
            try
            hold(hca,'on')
            plot(hca,x_traj(it),z_traj(it),'ok')
            hold(hca,'off')
            end
          end
        end
      end
      %hlinks = linkprop(h,{'XLim','YLim','CLim'});
      %hlinks.Targets(1)
      compact_panels(h,0,0.02)
    end    
        
    t_traj_old = t_traj;
    dt = 1;
    t_traj = twpe1:dt:twpe2;
    
    x_traj = interp1(t_traj_old,x_traj,t_traj,'spline');
    z_traj = interp1(t_traj_old,z_traj,t_traj,'spline');
    vx_traj = diff(x_traj)/(dt/50);
    vz_traj = diff(z_traj)/(dt/50);
end
% dt = 10;
%t_traj = linspace(8001:dt:8199;
%nt = numel(t_traj);
% interpolate x and z to t
%x_traj = interp1(1:nx,x_traj,(1:nt)/nx)

% limited pic object for faster loading of data
twpelim = [twpe1 twpe2];
xlim = [min(x_traj) max(x_traj)]+[-0.1 0.1];
zlim = [min(z_traj) max(z_traj)]+[-0.1 0.1];
pic = ox.xlim(xlim).zlim(zlim).twpelim(twpelim);

% 3D matrices of data
if 1 % not necessary to redo every time, only if new times are added, or 
  % x,z expanded beyond xlim,zlim
  Bx = pic.Bx;
  By = pic.By;
  Bz = pic.Bz;
  vOx = pic.vx(3);
  vOy = pic.vy(3);
  vOz = pic.vz(3);
  vPx = pic.vx(1);
  vPy = pic.vy(1);
  vPz = pic.vz(1);
  nO = pic.n(3);
end 
% Make 3D meshgrid, which is reuired input to interp3.
% X,Z,T all have dimensions [nz,nx,nt], thats why I need to
% permute the data below, e.g. Bx have dimension [nx,nz,nt], and the 
% permute(Bx,[2 1 3]) gives it dimensions [nz,nx,nt].
if numel(unique(twpelim)) == 1 % Single time, 2D interpolation
  [X,Z] = meshgrid(pic.xi,pic.zi);

  trBx = interp2(X,Z,permute(Bx,[2 1]),x_traj,z_traj);
  trBy = interp2(X,Z,permute(By,[2 1]),x_traj,z_traj);
  trBz = interp2(X,Z,permute(Bz,[2 1]),x_traj,z_traj);
  trVOx = interp2(X,Z,permute(vOx,[2 1]),x_traj,z_traj);
  trVOy = interp2(X,Z,permute(vOy,[2 1]),x_traj,z_traj);
  trVOz = interp2(X,Z,permute(vOz,[2 1]),x_traj,z_traj);
  trVPx = interp2(X,Z,permute(vPx,[2 1]),x_traj,z_traj);
  trVPy = interp2(X,Z,permute(vPy,[2 1]),x_traj,z_traj);
  trVPz = interp2(X,Z,permute(vPz,[2 1]),x_traj,z_traj);
  trNO = interp2(X,Z,permute(nO,[2 1]),x_traj,z_traj); 
else % Multiple times, 3D interpolation
  [X,Z,T] = meshgrid(pic.xi,pic.zi,pic.twpe);

  trBx = interp3(X,Z,T,permute(Bx,[2 1 3]),x_traj,z_traj,t_traj);
  trBy = interp3(X,Z,T,permute(By,[2 1 3]),x_traj,z_traj,t_traj);
  trBz = interp3(X,Z,T,permute(Bz,[2 1 3]),x_traj,z_traj,t_traj);
  trVOx = interp3(X,Z,T,permute(vOx,[2 1 3]),x_traj,z_traj,t_traj);
  trVOy = interp3(X,Z,T,permute(vOy,[2 1 3]),x_traj,z_traj,t_traj);
  trVOz = interp3(X,Z,T,permute(vOz,[2 1 3]),x_traj,z_traj,t_traj);
  trVPx = interp3(X,Z,T,permute(vPx,[2 1 3]),x_traj,z_traj,t_traj);
  trVPy = interp3(X,Z,T,permute(vPy,[2 1 3]),x_traj,z_traj,t_traj);
  trVPz = interp3(X,Z,T,permute(vPz,[2 1 3]),x_traj,z_traj,t_traj);
  trNO = interp3(X,Z,T,permute(nO,[2 1 3]),x_traj,z_traj,t_traj);
end

%% Make TSeries from data
time = EpochTT(i_traj_new);
tsB = irf.ts_vec_xyz(time,[trBx' trBy' trBz']);

%% Plot figure
figure(22)
nrows = 7;
ncols = 1;
npanels = nrows*ncols;

h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % r(t)
  hca = h(isub); isub = isub + 1;
  [AX,H1,H2] = plotyy(hca,t_traj,x_traj,t_traj,z_traj);
  AX(1).YLabel.String = 'x';
  AX(2).YLabel.String = 'z';
  legend(hca,{'x','z'},'Box','off','orientation','horizontal');
end
if 1 % v(t)
  hca = h(isub); isub = isub + 1;
  %[AX,H1,H2] = plotyy(hca,t_traj(1:end-1),smooth(vx_traj,50),t_traj(1:end-1),vz_traj);
  [AX,H1,H2] = plotyy(hca,t_traj,smooth(vx_traj,50),t_traj,vz_traj);
  AX(1).YLabel.String = 'v_x';
  AX(2).YLabel.String = 'v_z';
  legend(hca,{'v_x','v_z'},'Box','off','orientation','horizontal');
end
if 1 % B(t)
  hca = h(isub); isub = isub + 1;
  plot(hca,t_traj,trBx,t_traj,trBy,t_traj,trBz)
  hca.YLabel.String = 'B';
  legend(hca,{'B_x','B_y','B_z'},'Box','off','orientation','horizontal');
end
if 1 % vP(t)
  hca = h(isub); isub = isub + 1;
  plot(hca,t_traj,trVPx,t_traj,trVPy,t_traj,trVPz)
  hca.YLabel.String = 'V_O';
  legend(hca,{'V_{P,x}','V_{P,y}','V_{P,z}'},'Box','off','orientation','horizontal');
end
if 1 % vO(t)
  hca = h(isub); isub = isub + 1;
  plot(hca,t_traj,trVOx,t_traj,trVOy,t_traj,trVOz)
  hca.YLabel.String = 'V_O';
  legend(hca,{'V_{O,x}','V_{O,y}','V_{O,z}'},'Box','off','orientation','horizontal');
end
if 1 % vO(t) - v
  hca = h(isub); isub = isub + 1;
  plot(hca,t_traj,trVOx-vx_traj,t_traj,trVOy,t_traj,trVOz-vz_traj)
  hca.YLabel.String = 'V_O';
  legend(hca,{'V_{O,x}-v_{x,ref}','V_{O,y}-0','V_{O,z}-v_{z,ref}'},'Box','off','orientation','horizontal');
end
if 1 % n(t)
  hca = h(isub); isub = isub + 1;
  plot(hca,t_traj,trNO)
  hca.YLabel.String = 'n_O';
  hca.XLabel.String = 't\omega_{pe}';
  %legend(hca,{'V_{O,x}','V_{O,y}','V_{O,z}'},'Box','off','orientation','horizontal');
end

compact_panels(0)
for ip = 1:npanels
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end

%% Make corresponding movie
tr.t = t_traj/50; % twci
tr.x = x_traj;
tr.z = z_traj;
tr.ntr = 1;


zlim = [-13 13];
xlim = [100 300];
twpe = [8000 10000];
pic = ox.twpelim(twpe).xlim(xlim).zlim(zlim);

varstrs = {'n(3)'};
clims = {[0 0.4]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmapc2 = pic_colors('candy2');
cmapma = pic_colors('matlab');
cmaps = {cmapbr};
filename = [printpath 'ox_nO_trajectories_twpe08000-10000_test1'];
%pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename,'tr',tr100.z0find(4));
%pic.movie(varstrs,'clim',clims,'filename',filename,'tr',tr,'trajcolordot','t');
pic.movie(varstrs,'clim',clims,'filename',filename,'tr',tr);
