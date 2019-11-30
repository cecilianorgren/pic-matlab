%% Read and plot distributions
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/distributions/';
timestep = 08000;
str_timestep = sprintf('%05.0f',timestep);
txttime = sprintf('timestep = %05.0f',timestep); 
    
idist = 0;
tic
for distnumber = 200%1:139%30:40%39%180:200%:250%:100%:100%:10%40%:40%:4%:40
  read_sub_dir = '/1/';
  txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber); % michael's perturbation
  if not(exist(txtfile,'file'))
    warning(sprintf('File not found: %s', txtfile))
    continue
  end  
    
  idist = idist + 1;
  %idist = distnumber;
  
  % Load data  
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
      = read_distributions(txtfile);
    
  if x(1)<-150
    xlo = xlo-x0;
    xhi = xhi-x0;
    zlo = zlo;
    zhi = zhi;    
  end
  vx = axes;
  vy = axes;
  vz = axes;
  Bloc.x = B.x(find_closest_ind(x,0.5*(xlo+xhi)),find_closest_ind(z,0.5*(zlo+zhi)));
  Bloc.y = B.y(find_closest_ind(x,0.5*(xlo+xhi)),find_closest_ind(z,0.5*(zlo+zhi)));
  Bloc.z = B.z(find_closest_ind(x,0.5*(xlo+xhi)),find_closest_ind(z,0.5*(zlo+zhi)));
  disp(sprintf('%.2f ',xlo,xhi,zlo,zhi))
  ispecies = 3;
  vpeaks = find_dist_peaks(fxyz(:,:,:,ispecies));
  %pause(0.1)  
  
  species_clim_temperature = [0.3 0.18 0.3 0.3];
  species_strs = {'i1','e1','i2','e2'}; % keep same order as in distributions!
  doPatch = 1;
  for ispecies = 3 % Temperature       
    if 1 % t**.scalar, t**.par, t**.perp, anisotropy 3x1 distributions for the different species (**)
      %%
      species_str = species_strs{ispecies};
      print_subdir = sprintf('/t_f_%s/',species_str);

      doContourF = 0;
      doQ = 0;
      doA = 1;
      doPatch = 1;
      doLogF = 0;
      doBDir = 1;

      cPatch = 'k';
      distCmap = flipdim(pic_colors('blue_red'),1); distCmap = distCmap(round(end/2):end,:);
      distCmap = irf_colormap('waterfall');
      distCmap = pic_colors('candy');    
      varstrs = {'B.z','ve1.x','log10(pi2.scalar./pi1.scalar)','log10(ti2.scalar./ti1.scalar)','log10(ni2./ni1)'};

      xlim = torow([-60 0]);
      zlim = torow([0 7]);% zlim = [-10 10];
      ipx1 = find(x>xlim(1),1,'first');
      ipx2 = find(x<xlim(2),1,'last');
      ipz1 = find(z>zlim(1),1,'first');
      ipz2 = find(z<zlim(2),1,'last');
      ipx = ipx1:2:ipx2;
      ipz = ipz1:2:ipz2;

      varstr_Q = {'E'};
      dataQ = eval(varstr_Q{1});
      cQ = [0 0 0];
      scaleQ = 2;
      ipxQ = ipx1:25:ipx2;
      ipzQ = ipz1:15:ipz2; 

      cA = 0*[0.8 0.8 0.8];
      nA = [0:-1:min(A(:))];
      ipxA = ipx1:20:ipx2;
      ipzA = ipz1:20:ipz2;    

      if idist == 1 % Initialize figure, only once
        fig = figure(501);
        screensize = get(groot,'Screensize');
        fig.Position = [screensize(1) screensize(2) screensize(3)*0.5 screensize(4)*0.7];

        nrows_maps = 3;          
        ncols_maps = 2; % which will be merged
        npanels_maps = nrows_maps;

        nrows_dists = 3;
        ncols_dists = 1;    
        npanels_dists = nrows_dists*ncols_dists;      

        ncols = ncols_maps + ncols_dists;
        clear h hb hmap hp;

        isub = 0;
        for imap = 1:npanels_maps
          isub = isub + 1;
          hmap(isub) = subplot(nrows_maps,ncols,ncols*(imap-1)+[1:ncols_maps]);
        end          
        isub = 0;
        for imap = 1:nrows_dists
          isub = isub + 1;
          h(isub) = subplot(nrows_dists,ncols,ncols*(imap-1)+ncols_maps+1);
        end
        linkaxes(hmap)
      end
      if idist == 1 % Plot map panels
        doA = 1;
        if doA    
          cA = [0.8 0.8 0.8];
          cA = [0.7 0.7 0.7];
          nA = 40;
          nA = [0:-1:min(A(:))];
          ipxA = ipx1:5:ipx2;
          ipzA = ipz1:5:ipz2;
        end

        % Panels
        isub = 1;
        if 1 % t**.scalar
          hca = hmap(isub); isub = isub + 1;
          varstr = sprintf('t%s.scalar',species_str);
          variable = eval(varstr);
          himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
          hca.XLabel.String = 'x (d_i)';
          hca.YLabel.String = 'z (d_i)';
          hca.CLim = species_clim_temperature(ispecies)*[0 1];
          hcb = colorbar('peer',hca);
          hcb.YLabel.String = varstr;
          hb(ivar) = hcb;
          colormap(hca,pic_colors('candy'));
          if doA
            hold(hca,'on')
            hcont = contour(hca,x(ipxA),z(ipzA),A(ipxA,ipzA)',nA,'color',cA,'linewidth',0.7);  
            hold(hca,'off')  
          end
        end            
        if 1 % t**.par
          hca = hmap(isub); isub = isub + 1;
          varstr = sprintf('t%s.par',species_str);
          variable = eval(varstr);
          himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
          hca.XLabel.String = 'x (d_i)';
          hca.YLabel.String = 'z (d_i)';
          hca.CLim = species_clim_temperature(ispecies)*[0 1];
          hcb = colorbar('peer',hca);
          hcb.YLabel.String = varstr;
          hb(ivar) = hcb;
          colormap(hca,pic_colors('candy'));
          if doA
            hold(hca,'on')
            hcont = contour(hca,x(ipxA),z(ipzA),A(ipxA,ipzA)',nA,'color',cA,'linewidth',0.7);  
            hold(hca,'off')  
          end
        end          
        if 1 % t**.perp
          hca = hmap(isub); isub = isub + 1;
          varstr = sprintf('t%s.perp',species_str);
          variable = eval(varstr);
          himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
          hca.XLabel.String = 'x (d_i)';
          hca.YLabel.String = 'z (d_i)';
          hca.CLim = species_clim_temperature(ispecies)*[0 1];
          hcb = colorbar('peer',hca);
          hcb.YLabel.String = varstr;      
          hb(ivar) = hcb;      
          colormap(hca,pic_colors('candy'));
          if doA
            hold(hca,'on')
            hcont = contour(hca,x(ipxA),z(ipzA),A(ipxA,ipzA)',nA,'color',cA,'linewidth',0.7);  
            hold(hca,'off')  
          end
        end            
        if 0 % log2(t**.par./t**.perp)
          hca = hmap(isub); isub = isub + 1;
          varstr = sprintf('log2(t%s.par./t%s.perp)',species_str',species_str);
          variable = eval(varstr);
          himag = imagesc(hca,x(ipx),z(ipz),real(variable(ipx,ipz))');
          hca.XLabel.String = 'x (d_i)';
          hca.YLabel.String = 'z (d_i)';
          hca.CLim = 3*[-1 1];
          hcb = colorbar('peer',hca);
          hcb.YLabel.String = varstr;
          hb(ivar) = hcb;
          colormap(hca,pic_colors('blue_red'));
          if doA
            hold(hca,'on')
            hcont = contour(hca,x(ipxA),z(ipzA),A(ipxA,ipzA)',nA,'color',cA,'linewidth',0.7);  
            hold(hca,'off')  
          end
        end
        %compact_panels
        for ipanel = 1:nrows_maps
          hmap(ipanel).YDir = 'normal';
          hmap(ipanel).XDir = 'reverse';
          hmap(ipanel).XLim = xlim;
          hmap(ipanel).YLim = zlim;               
        end
        hmap(1).Title.String = sprintf('time = %g (1/wci) = %g (1/wpe)',time,timestep);
        set(gcf,'color',[1 1 1]);
        c_eval('h(?).XLabel.String = [];',1:nrows_maps-1)        
      end        
      if doLogF % Take log10 of data if specified
        strprint = sprintf('dist_%s_t_%05.0f_x_%g_%g_z_%g_%g_f_map',species_str,timestep,xlo,xhi,zlo,zhi);
      else
        strprint = sprintf('dist_%s_t_%05.0f_x_%g_%g_z_%g_%g_log10f_map',species_str,timestep,xlo,xhi,zlo,zhi);      
  %       strprint = sprintf('dist_%04.0f_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_log10f_map',distnumber,'ie',timestep,xlo,xhi,zlo,zhi);
  %     else
  %       strprint = sprintf('dist_%04.0f_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_map',distnumber,'ie',timestep,xlo,xhi,zlo,zhi);
      end
      if doPatch % Plot/update patch    
        if exist('hpatch','var'); delete(hpatch); end % delete last patch      
        if exist('hp','var'); delete(hp); end % delete last patch      
        for ipanel = 1:nrows_maps
          hca = hmap(ipanel);
          hold(hca,'on')        
          hpatch = patch(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'k');
          hp(ipanel) = hpatch;
          hpatch.FaceAlpha = 0;
          hpatch.LineWidth = 1;
          hold(hca,'off') 
        end                       
      end       

      str_title = {sprintf('timestep = %g, time = %g, x_box = [%g, %g], z_box = [%g, %g], Bloc = [%.2f,%.2f,%.2f]',timestep,time,xlo,xhi,zlo,zhi,Bloc.x,Bloc.y,Bloc.z),'', varstrs{1}}; 
      hmap(1).Title.String = str_title;
      hmap(1).Title.Interpreter = 'none';

      isub = 1;    
      strtitle = sprintf('f.%s',species_str);        
      if 1 % fxy
        hca = h(isub); isub = isub + 1;
        toplot = squeeze(sum(fxy(:,:,ispecies),3));
        if doLogF, toplot = log10(toplot); end
        if doContourF      
          contourf(hca,vx(:,ispecies(1)),vy(:,ispecies(1)),toplot')
        else
          imagesc(hca,vx(:,ispecies(1)),vy(:,ispecies(1)),toplot')
        end
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vy';
        hcb = colorbar('peer',hca);
        hb(isub-1) = hcb;
        if doLogF,hcb.YLabel.String = 'log10(f)'; else, hcb.YLabel.String = 'f'; end
        hca.Title.String = strtitle;
        hca.XLim = axes([1 end],ispecies(1));
        hca.YLim = axes([1 end],ispecies(1));
      end        
      if 1 % fxz
        hca = h(isub); isub = isub + 1;
        toplot = squeeze(sum(fxz(:,:,ispecies),3));
        if doLogF, toplot = log10(toplot); end
        if doContourF      
          contourf(hca,vx(:,ispecies(1)),vx(:,ispecies(1)),toplot')
        else
          imagesc(hca,vx(:,ispecies(1)),vx(:,ispecies(1)),toplot')
        end
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vz';
        hcb = colorbar('peer',hca);
        hb(isub-1) = hcb;
        if doLogF,hcb.YLabel.String = 'log10(f)'; else, hcb.YLabel.String = 'f'; end
        hca.Title.String = strtitle;
        hca.XLim = axes([1 end],ispecies(1));
        hca.YLim = axes([1 end],ispecies(1));
      end        
      if 1 % fzy
        hca = h(isub); isub = isub + 1;
        toplot = squeeze(sum(fyz(:,:,ispecies),3));
        if doLogF, toplot = log10(toplot); end
        if doContourF      
          contourf(hca,vy(:,ispecies(1)),vz(:,ispecies(1)),toplot')
        else
          imagesc(hca,vy(:,ispecies(1)),vz(:,ispecies(1)),toplot')
        end
        hca.XLabel.String = 'vy';
        hca.YLabel.String = 'vz';
        hcb = colorbar('peer',hca);
        hb(isub-1) = hcb;
        if doLogF,hcb.YLabel.String = 'log10(f)'; else, hcb.YLabel.String = 'f'; end
        hca.Title.String = strtitle;
        hca.XLim = axes([1 end],ispecies(1));
        hca.YLim = axes([1 end],ispecies(1));
      end    

      if doBDir
        for ipanel = 1 % xy        
          line_slope = (Bloc.y/Bloc.x);
          xx = min(h(ipanel).XLim(2)*[1 1/abs(line_slope)])*[-1 1];
          hold(h(ipanel),'on')                
          hBline = plot(h(ipanel),xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
          hold(h(ipanel),'off')
        end   
        for ipanel = 2 % xz     
          line_slope = (Bloc.z/Bloc.x);
          xx = min(h(ipanel).XLim(2)*[1 1/abs(line_slope)])*[-1 1];
          hold(h(ipanel),'on')                
          hBline = plot(h(ipanel),xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
          hold(h(ipanel),'off')
        end   
        for ipanel = 3 % yz
          line_slope = (Bloc.z/Bloc.y);
          xx = min(h(ipanel).XLim(2)*[1 1/abs(line_slope)])*[-1 1];
          hold(h(ipanel),'on')                
          hBline = plot(h(ipanel),xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
          hold(h(ipanel),'off')
        end
      end

      for ipanel = 1:npanels_maps  
        hmap(ipanel).FontSize = 12;
        %axis(h(ipanel),'equal')
        colormap(distCmap)
      end    
      all_max_clim_dists = 0;
      for ipanel = 1:npanels_dists
        all_max_clim_dists = min([all_max_clim_dists h(ipanel).CLim]);
        h(ipanel).FontSize = 12;
        h(ipanel).XDir = 'reverse';
        %axis(h(ipanel),'equal')
        switch ispecies
          case 1
            tickstep = 1;
            vlim = 5;
          case 2
            tickstep = 5;
            vlim = 10;
          case 3
            tickstep = 1;
            vlim = 1;
          case 4
            tickstep = 5;
            vlim = 10;              
        end

        h(ipanel).XTick = -15:tickstep:15;
        h(ipanel).YTick = -15:tickstep:15;
        h(ipanel).XLim = vlim*[-1 1];
        h(ipanel).YLim = vlim*[-1 1];

        h(ipanel).XGrid = 'on';
        h(ipanel).YGrid = 'on';
        %h(ipanel).XMinorGrid = 'on';
        %h(ipanel).YMinorGrid = 'on';
        h(ipanel).XDir = 'normal';
        h(ipanel).YDir = 'normal';
        axis(h(ipanel),'square')      
        colormap(distCmap)
      end    
      for ipanel = 1:npanels_dists
        h(ipanel).XDir = 'reverse';
      end
      linkprop(h([1 2 3]),'CLim'); % link CLim's        
      print_path = [savedir_root str_timestep print_subdir];
      if ~exist(print_path,'dir'), mkdir(print_path); end
      drawnow
      print('-dpng','-r200',[print_path strprint '.png']);        
    end
  end
  toc
end
 