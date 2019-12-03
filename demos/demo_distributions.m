timestep = 05000;

%% Read fields file
txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/data/fields-%05.0f.dat',timestep); % michael's perturbation
tic; [x,z,E,B,...
  ni1,ne1,ni2,ne2,...
  vi1,ve1,vi2,ve2,...
  ji1,je1,ji2,je2,...
  pi1,pe1,pi2,pe2,...
  ti1,te1,ti2,te2,...
  dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q] = read_fields(txtfile); toc
A = vector_potential(x,z,B.x,B.z); % vector potential
[saddle_locations,saddle_values] = saddle(A);
ExB = cross_product(E.x,E.y,E.z,B.x,B.y,B.z); % Poynting flux
vExB.x = ExB.x./B.abs./B.abs;
vExB.y = ExB.y./B.abs./B.abs;
vExB.z = ExB.z./B.abs./B.abs;
ti2.scalar(ni2<0.01) = NaN;
vi2.x(ni2<0.01) = NaN;
vi2.y(ni2<0.01) = NaN;
vi2.z(ni2<0.01) = NaN;
ni2(ni2<0.00) = 0;
pi2.scalar(pi2.scalar<0.00) = 0;
ti2.scalar(ti2.scalar<0.00) = 0;
ti2.scalar(ni2<0.01) = NaN;
pi2.scalar(ni2<0.00) = 0;
ti2.scalar(ni2<0.01) = NaN;
jtot.x = ji1.x + ji2.x - je1.x - je2.x;
jtot.y = ji1.y + ji2.y - je1.y - je2.y;
jtot.z = ji1.z + ji2.z - je1.z - je2.z;
je.x = je1.x + je2.x;
je.y = je1.y + je2.y;
je.z = je1.z + je2.z;
ji.x = ji1.x + ji2.x;
ji.y = ji1.y + ji2.y;
ji.z = ji1.z + ji2.z;

%% Rotate into different coordinate system

%% Read and plot distributions
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/distributions/';
timestep = 08000;
str_timestep = sprintf('%05.0f',timestep);
txttime = sprintf('timestep = %05.0f',timestep); 
    
idist = 0;
tic
for distnumber = 200:201%30:40%39%180:200%:250%:100%:100%:10%40%:40%:4%:40
  read_sub_dir = '/1/';
  txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber); % michael's perturbation
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
  pause(0.1)  
  if 0 % 1x3 single species plot
    fig = figure(301);
    %fig.Position = [
    nrows = 1;
    ncols = 3;
    npanels = nrows*ncols;
    clear h hb;
    for ipanel = 1:npanels  
      h(ipanel) = subplot(nrows,ncols,ipanel); 
    end
    for ispecies = 1:4
    isub = 1;
      strtitle = sprintf('species %g\n%s\nx = [%g, %g], z = [%g, %g]',ispecies,txttime,xlo,xhi,zlo,zhi);
      strprint = sprintf('dist_%04.0f_species_%g_t_%05.0f_x_%g_%g_z_%g_%g',ispecies,distnumber,timestep,xlo,xhi,zlo,zhi);
      if 1 % fxy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vx(:,ispecies),vy(:,ispecies),squeeze(fxy(:,:,ispecies))')
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vy';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
      end
      if 1 % fxz
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vx(:,ispecies),vz(:,ispecies),squeeze(fxz(:,:,ispecies))')  
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vz';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
      end
      if 1 % fzy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vy(:,ispecies),vz(:,ispecies),squeeze(fyz(:,:,ispecies))')  
        hca.XLabel.String = 'vz';
        hca.YLabel.String = 'vy';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
      end
      for ipanel = 1:npanels
        h(ipanel).FontSize = 12;  
        h(ipanel).XLim = axes([1 end],ispecies);
        h(ipanel).YLim = axes([1 end],ispecies);
        %axis(h(ipanel),'equal')   
        axis(h(ipanel),'square')
      end
      print('-dpng','-r200',[savedir_root sub_dir '/' strprint '.png']);
      drawnow
      pause(1)
    end
  end
  if 0 % 1x3 single species plot, including peaks, not done
    fig = figure(301);
    distCmap = pic_colors('candy');
    
    %fig.Position = [
    nrows = 1;
    ncols = 3;
    npanels = nrows*ncols;
    clear h hb;
    for ipanel = 1:npanels  
      h(ipanel) = subplot(nrows,ncols,ipanel); 
    end
    for ispecies = 2
    vpeaks = find_dist_peaks(fxyz(:,:,:,ispecies));
    isub = 1;
      strtitle = sprintf('species %g\n%s\nx = [%g, %g], z = [%g, %g]',ispecies,txttime,xlo,xhi,zlo,zhi);
      strprint = sprintf('dist_%04.0f_species_%g_t_%05.0f_x_%g_%g_z_%g_%g',ispecies,distnumber,timestep,xlo,xhi,zlo,zhi);
      if 1 % fxy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vx(:,ispecies),vy(:,ispecies),squeeze(fxy(:,:,ispecies))')
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vy';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
      end
      if 1 % fxz
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vx(:,ispecies),vz(:,ispecies),squeeze(fxz(:,:,ispecies))')  
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vz';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
      end
      if 1 % fzy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vy(:,ispecies),vz(:,ispecies),squeeze(fyz(:,:,ispecies))')  
        hca.XLabel.String = 'vz';
        hca.YLabel.String = 'vy';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
      end
      for ipanel = 1:npanels
        h(ipanel).FontSize = 12;  
        h(ipanel).XLim = axes([1 end],ispecies);
        h(ipanel).YLim = axes([1 end],ispecies);
        %axis(h(ipanel),'equal')   
        axis(h(ipanel),'square')
      end
      %print('-dpng','-r200',[savedir_root sub_dir '/' strprint '.png']);
      %drawnow
      %pause(1)
    end
  end
  if 0 % 4x3 all species plot
    fig = figure(302);
    screensize = get(groot,'Screensize');
    fig.Position = [screensize(1) screensize(2) screensize(3)*0.4 screensize(4)*0.8];
    nrows = 4;
    ncols = 3;
    npanels = nrows*ncols;
    clear h hmap hb;
    for ipanel = 1:npanels  
      h(ipanel) = subplot(nrows,ncols,ipanel); 
    end
    
    strprint = sprintf('dist_%04.0f_species_%s_t_%05.0f_x_%g_%g_z_%g_%g',distnumber,'1234',timestep,xlo,xhi,zlo,zhi);
    
    isub = 1;    
    for ispecies = 1:4
      strtitle = sprintf('species %g\n%s\nx = [%g, %g], z = [%g, %g]',ispecies,txttime,xlo,xhi,zlo,zhi);
      if 1 % fxy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vx(:,ispecies),vy(:,ispecies),squeeze(fxy(:,:,ispecies))')
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vy';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
        hca.XLim = axes([1 end],ispecies);
      	hca.YLim = axes([1 end],ispecies);
      end
      if 1 % fxz
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vx(:,ispecies),vz(:,ispecies),squeeze(fxz(:,:,ispecies))')  
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vz';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
        hca.XLim = axes([1 end],ispecies);
      	hca.YLim = axes([1 end],ispecies);
      end
      if 1 % fzy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vy(:,ispecies),vz(:,ispecies),squeeze(fyz(:,:,ispecies))')  
        hca.XLabel.String = 'vy';
        hca.YLabel.String = 'vz';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
        hca.XLim = axes([1 end],ispecies);
      	hca.YLim = axes([1 end],ispecies);
      end      
    end    
    for ipanel = 1:npanels
      h(ipanel).FontSize = 12;  
      
      %axis(h(ipanel),'equal')   
      axis(h(ipanel),'square')
    end
    %print('-dpng','-r200',[savedir_root sub_dir '/' strprint '.png']);
  end
  if 0 % 4x3 all species plot, plus map to show location
    doA = 0;
    doPatch = 1;
    cPatch = 'k';
    distCmap = flipdim(pic_colors('blue_red'),1); distCmap = distCmap(round(end/2):end,:);
    varstr = 've1.x';
    xlim = torow([80 x([round(end/2)])]);
    zlim = torow(z([1 end])); zlim = [-10 10];
    ipx1 = find(x>xlim(1),1,'first');
    ipx2 = find(x<xlim(2),1,'last');
    ipz1 = find(z>zlim(1),1,'first');
    ipz2 = find(z<zlim(2),1,'last');
    ipx = ipx1:2:ipx2;
    ipz = ipz1:2:ipz2;

    fig = figure(303);
    screensize = get(groot,'Screensize');
    fig.Position = [screensize(1) screensize(2) screensize(3)*0.4 screensize(4)*1];
    nrows = 5;
    ncols = 3;
    npanels = nrows*ncols-ncols;
    clear h hb;
    hmap = subplot(nrows,ncols,1:3);
    isub = 1;
    for ipanel = 4:(npanels+ncols)
      h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;
    end

    strprint = sprintf('dist_%04.0f_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_map',distnumber,'1234',timestep,xlo,xhi,zlo,zhi);
    
    if 1 % Plot map showing location of dists
      hca = hmap;
      himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
      hca.Title.String = sprintf('%s, timestep = %g, x = [%g, %g], z = [%g, %g]',varstr,timestep,xlo,xhi,zlo,zhi);       
      hca.Title.Interpreter = 'none';  
      if any(abs(himag.CData(not(isnan(himag.CData(:)))))) % do if any value is non-zero
        hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      end
      hcb = colorbar('peer',hca);     
      colormap(hca,cn.cmap('blue_red'));
      
      if doPatch
        if exist('hpatch','var'); delete(hpatch); end
        hold(hca,'on')        
        hpatch = patch(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'k');
        hpatch.FaceAlpha = 0;
        hpatch.LineWidth = 1;
        hold(hca,'off')      
      end
      if doA
        hold(hca,'on')
        hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',0.5); 
    %     for ixline = 1:size(saddle_locations,1)
    %       sepA = saddle_values(ixline);
    %       hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',sepA*[1 1],'color',cA.^4,'linewidth',2.0);  
    %     end
        hold(hca,'off')  
      end
      if doQ
        hold(hca,'on')
        hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
        hold(hca,'off')  
      end        
    end
    isub = 1;
    for ispecies = 1:4      
      strtitle = sprintf('species %g',ispecies);
      if 1 % fxy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vx(:,ispecies),vy(:,ispecies),squeeze(fxy(:,:,ispecies))')
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vy';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
        hca.XLim = axes([1 end],ispecies);
      	hca.YLim = axes([1 end],ispecies);
      end
      if 1 % fxz
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vx(:,ispecies),vz(:,ispecies),squeeze(fxz(:,:,ispecies))')  
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vz';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
        hca.XLim = axes([1 end],ispecies);
      	hca.YLim = axes([1 end],ispecies);
      end
      if 1 % fzy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vy(:,ispecies),vz(:,ispecies),squeeze(fyz(:,:,ispecies))')  
        hca.XLabel.String = 'vy';
        hca.YLabel.String = 'vz';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
        hca.XLim = axes([1 end],ispecies);
      	hca.YLim = axes([1 end],ispecies);
      end      
    end    
    for ipanel = 1:npanels
      h(ipanel).FontSize = 12;        
      %axis(h(ipanel),'equal')   
      axis(h(ipanel),'square')
      colormap(distCmap)
    end
    print('-dpng','-r200',[savedir_root sub_dir '/' strprint '.png']);
  end  
  if 0 % 2x3 electrons and ions added together, respectively, plus map to show location
    doQ = 0;
    doA = 1;
    doPatch = 1;
    doLogF = 1;
    
    cPatch = 'k';
    distCmap = flipdim(pic_colors('blue_red'),1); distCmap = distCmap(round(end/2):end,:);
    varstr = 'pi2.scalar';
    varstr = 'ne2.*ve2.x';
    varstr = 'ne1.*ve1.x';
    varstr = 'ni2./ni1';
    varstr = 'B.z';
    %varstr = 'ti2.scalar';
    xlim = torow([100 x([round(end/2)])]);
    zlim = torow(z([1 end])); zlim = [-10 10];
    ipx1 = find(x>xlim(1),1,'first');
    ipx2 = find(x<xlim(2),1,'last');
    ipz1 = find(z>zlim(1),1,'first');
    ipz2 = find(z<zlim(2),1,'last');
    ipx = ipx1:2:ipx2;
    ipz = ipz1:2:ipz2;
        
    cA = 0*[0.8 0.8 0.8];
    nA = [0:-1:min(A(:))];
    ipxA = ipx1:20:ipx2;
    ipzA = ipz1:20:ipz2;    

    fig = figure(303);
    screensize = get(groot,'Screensize');
    fig.Position = [screensize(1) screensize(2) screensize(3)*0.4 screensize(4)*0.6];
    nrows = 3;
    ncols = 3;
    npanels = nrows*ncols-ncols;
    clear h hb;
    hmap = subplot(nrows,ncols,1:3);
    isub = 1;
    for ipanel = 4:(npanels+ncols)
      h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;
    end
    
    if doLogF
      strprint = sprintf('dist_%04.0f_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_log10f_map',distnumber,'ie',timestep,xlo,xhi,zlo,zhi);
    else
      strprint = sprintf('dist_%04.0f_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_map',distnumber,'ie',timestep,xlo,xhi,zlo,zhi);
    end
    
    if 1 % Plot map showing location of dists
      hca = hmap;
      if 1%idist == 1 % only plot once
        variable = eval(varstr);
        himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
        hca.YDir = 'normal';
        hca.XLabel.String = 'x (d_i)';
        hca.YLabel.String = 'z (d_i)';
        %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
        hca.Title.String = sprintf('%s, timestep = %g, time = %g, x = [%g, %g], z = [%g, %g]',varstr,timestep,time,xlo,xhi,zlo,zhi);       
        hca.Title.Interpreter = 'none';  
        if any(abs(himag.CData(not(isnan(himag.CData(:)))))) % do if any value is non-zero
          hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
        end
        hcb = colorbar('peer',hca);     
        colormap(hca,cn.cmap('blue_red'));
        hold(hca,'on')
      end
      if doPatch
        %if exist('hpatch','var'); delete(hpatch); end
        hold(hca,'on')        
        hpatch = patch(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'k');
        hpatch.FaceAlpha = 0;
        hpatch.LineWidth = 1;
        hold(hca,'off')      
      end
      if doA %&& idist == 1
        hold(hca,'on')
        hA = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',0.5); 
    %     for ixline = 1:size(saddle_locations,1)
    %       sepA = saddle_values(ixline);
    %       hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',sepA*[1 1],'color',cA.^4,'linewidth',2.0);  
    %     end
        hold(hca,'off')  
      end
      if doQ
        hold(hca,'on')
        hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ));
        hold(hca,'off')  
      end
    end
    isub = 1;      
    
    % Ions
    ispecies = [1 3];
    strtitle = sprintf('species [%g %g]',ispecies);
    if 1 % fxy, all ions together
      hca = h(isub); isub = isub + 1;
      toplot = squeeze(sum(fxy(:,:,ispecies),3));
      if doLogF, toplot = log10(toplot); end
      imagesc(hca,vx(:,ispecies(1)),vy(:,ispecies(1)),toplot')
      hca.XLabel.String = 'vx';
      hca.YLabel.String = 'vy';
      hcb = colorbar('peer',hca);
      if doLogF,hcb.YLabel.String = 'log10(f)'; else, hcb.YLabel.String = 'f'; end
      hca.Title.String = strtitle;
      hca.XLim = axes([1 end],ispecies(1));
      hca.YLim = axes([1 end],ispecies(1));
    end
    if 1 % fxz, all ions together
      hca = h(isub); isub = isub + 1;
      toplot = squeeze(sum(fxz(:,:,ispecies),3));
      if doLogF, toplot = log10(toplot); end
      imagesc(hca,vx(:,ispecies(1)),vx(:,ispecies(1)),toplot')
      hca.XLabel.String = 'vx';
      hca.YLabel.String = 'vz';
      hcb = colorbar('peer',hca);
      if doLogF,hcb.YLabel.String = 'log10(f)'; else, hcb.YLabel.String = 'f'; end
      hca.Title.String = strtitle;
      hca.XLim = axes([1 end],ispecies(1));
      hca.YLim = axes([1 end],ispecies(1));
    end
    if 1 % fzy, all ions together
      hca = h(isub); isub = isub + 1;
      toplot = squeeze(sum(fyz(:,:,ispecies),3));
      if doLogF, toplot = log10(toplot); end
      imagesc(hca,vy(:,ispecies(1)),vz(:,ispecies(1)),toplot')
      hca.XLabel.String = 'vy';
      hca.YLabel.String = 'vz';
      hcb = colorbar('peer',hca);
      if doLogF,hcb.YLabel.String = 'log10(f)'; else, hcb.YLabel.String = 'f'; end
      hca.Title.String = strtitle;
      hca.XLim = axes([1 end],ispecies(1));
      hca.YLim = axes([1 end],ispecies(1));
    end   
    % Electrons
    ispecies = [2 4]; 
    strtitle = sprintf('species [%g %g]',ispecies);
    if 1 % fxy, all electrons together
      hca = h(isub); isub = isub + 1;
      toplot = squeeze(sum(fxy(:,:,ispecies),3));
      if doLogF, toplot = log10(toplot); end
      imagesc(hca,vx(:,ispecies(1)),vy(:,ispecies(1)),toplot')
      hca.XLabel.String = 'vx';
      hca.YLabel.String = 'vy';
      hcb = colorbar('peer',hca);
      if doLogF,hcb.YLabel.String = 'log10(f)'; else, hcb.YLabel.String = 'f'; end
      hca.Title.String = strtitle;
      hca.XLim = axes([1 end],ispecies(1));
      hca.YLim = axes([1 end],ispecies(1));
    end
    if 1 % fxz, all electrons together
      hca = h(isub); isub = isub + 1;
      toplot = squeeze(sum(fxz(:,:,ispecies),3));
      if doLogF, toplot = log10(toplot); end
      imagesc(hca,vx(:,ispecies(1)),vx(:,ispecies(1)),toplot')
      hca.XLabel.String = 'vx';
      hca.YLabel.String = 'vz';
      hcb = colorbar('peer',hca);
      if doLogF,hcb.YLabel.String = 'log10(f)'; else, hcb.YLabel.String = 'f'; end
      hca.Title.String = strtitle;
      hca.XLim = axes([1 end],ispecies(1));
      hca.YLim = axes([1 end],ispecies(1));
    end
    if 1 % fzy, all electrons together
      hca = h(isub); isub = isub + 1;
      toplot = squeeze(sum(fyz(:,:,ispecies),3));
      if doLogF, toplot = log10(toplot); end
      imagesc(hca,vy(:,ispecies(1)),vz(:,ispecies(1)),toplot')
      hca.XLabel.String = 'vy';
      hca.YLabel.String = 'vz';
      hcb = colorbar('peer',hca);
      if doLogF,hcb.YLabel.String = 'log10(f)'; else, hcb.YLabel.String = 'f'; end
      hca.Title.String = strtitle;
      hca.XLim = axes([1 end],ispecies(1));
      hca.YLim = axes([1 end],ispecies(1));        
    end     

    for ipanel = 1:npanels  
      h(ipanel).FontSize = 12;
      %axis(h(ipanel),'equal')
      h(ipanel).XTick = -15:5:15;
      h(ipanel).YTick = -15:5:15;
      h(ipanel).XGrid = 'on';
      h(ipanel).YGrid = 'on';
      %h(ipanel).XMinorGrid = 'on';
      %h(ipanel).YMinorGrid = 'on';
      axis(h(ipanel),'square')      
      colormap(distCmap)
    end
    %pause
    %print('-dpng','-r200',[savedir_root sub_dir '/species_ie/' strprint '.png']);
  end  
  if 0 % a few maps with location on the left and 3x2 electrons and ions added together on the right, respectively
    print_subdir = sprintf('/species_ie_with_location_candy_B_Qj/');
    doContourF = 0;
    doQ = 0;
    doA = 1;
    doPatch = 1;
    doLogF = 1;
    doBDir = 1;
    
    cPatch = 'k';
    distCmap = flipdim(pic_colors('blue_red'),1); distCmap = distCmap(round(end/2):end,:);
    distCmap = irf_colormap('waterfall');
    distCmap = pic_colors('candy');
    varstrs = {'B.z','ve1.x','ji1.y-je1.y+ji2.y-je2.y','(ni1.*vi1.y + ni2.*vi2.y)./(ni1+ni2)','(ne1.*ve1.y + ne2.*ve2.y)./(ne1+ne2)'};
    varstrs = {'B.z','ve1.y','ve2.y','vi1.y','vi2.y'};
    %varstrs = {'B.z','vExB.x','vExB.y','vExB.z'};
    %varstrs = {'B.z','ve1.y','ve2.y','vi1.y','vi2.y','vExB.y'};
    varstrs = {'B.z','ve1.x','pi2.scalar./pi1.scalar','ti2.scalar./ti1.scalar','ni2./ni1'};
    varstrs = {'B.z','B.y','vi2.y','E.z','ve1.par','pae1'};
    varstrs = {'pae1','pae2','pai1','pai2','pai2-pae2','pai1-pae1'}; 
    varstrs = {'angle_vi1vi2','angle_ve1ve2','angle_jije','jtot.x','jtot.z','ni1+ni2-ne1-ne2'}; 
    varstrs = {'B.z','ve1.x','log10(pi2.scalar./pi1.scalar)','log10(ti2.scalar./ti1.scalar)','log10(ni2./ni1)'};
       
    xlim = torow([120 x([round(end/2)])]);
    zlim = torow(z([1 end])); zlim = [-10 10];
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
    ipxQ = ix1:25:ix2;
    ipzQ = iz1:15:iz2; 
    
    cA = 0*[0.8 0.8 0.8];
    nA = [0:-1:min(A(:))];
    ipxA = ix1:20:ix2;
    ipzA = iz1:20:iz2;    
    
    if idist == 1 % Initialize figure, only once
      fig = figure(304);
      screensize = get(groot,'Screensize');
      fig.Position = [screensize(1) screensize(2) screensize(3)*0.7 screensize(4)*0.7];

      nrows_maps = max([3 numel(varstrs)]);
      nrows_maps = numel(varstrs);
      ncols_maps = 2; % which will be merged
      npanels_maps = nrows_maps;

      nrows_dists = 3;
      ncols_dists = 2;    
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
        isub = isub + 1;
        h(isub) = subplot(nrows_dists,ncols,ncols*(imap-1)+ncols_maps+2);
      end
      linkaxes(hmap)
    end
    
    if doLogF
      strprint = sprintf('dist_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_f_map','ie',timestep,xlo,xhi,zlo,zhi);
    else
      strprint = sprintf('dist_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_log10f_map','ie',timestep,xlo,xhi,zlo,zhi);      
%       strprint = sprintf('dist_%04.0f_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_log10f_map',distnumber,'ie',timestep,xlo,xhi,zlo,zhi);
%     else
%       strprint = sprintf('dist_%04.0f_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_map',distnumber,'ie',timestep,xlo,xhi,zlo,zhi);
    end
    
    if 1 % Plot map showing location of dists
      if exist('hpatch','var'); delete(hpatch); end
      if exist('hp','var'); delete(hp); end
      for ivar = 1:nrows_maps
        varstr = varstrs{ivar};
        if ivar == 1
          str_title = {sprintf('timestep = %g, time = %g, x_box = [%g, %g], z_box = [%g, %g], Bloc = [%.2f,%.2f,%.2f]',timestep,time,xlo,xhi,zlo,zhi,Bloc.x,Bloc.y,Bloc.z),'',varstr}; 
        else
          str_title = sprintf('%s',varstr); 
        end
        hca = hmap(ivar);
        if idist == 1 % only plot once
          variable = eval(varstr);
          himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
          hca.YDir = 'normal';
          hca.XLabel.String = 'x (d_i)';
          hca.YLabel.String = 'z (d_i)';
          %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
          hca.Title.String = str_title;
          hca.Title.Interpreter = 'none';  
          if any(abs(himag.CData(not(isnan(himag.CData(:)))))) % do if any value is non-zero
            hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
          end
          hcb = colorbar('peer',hca);     
          colormap(hca,cn.cmap('blue_red'));          
        end
        if doPatch
          hold(hca,'on')        
          hpatch = patch(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'k');
          hp(ivar) = hpatch;
          hpatch.FaceAlpha = 0;
          hpatch.LineWidth = 1;
          hold(hca,'off')      
        end
        if doA && idist == 1
          hold(hca,'on')
          hA = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',0.5); 
      %     for ixline = 1:size(saddle_locations,1)
      %       sepA = saddle_values(ixline);
      %       hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',sepA*[1 1],'color',cA.^4,'linewidth',2.0);  
      %     end
          hold(hca,'off')  
        end
        if doQ && idist == 1
          hold(hca,'on')
          hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ),dataQ.z(ipxQ,ipzQ),scaleQ,'color',cQ);
          legend(hquiv,varstr_Q)
          hold(hca,'off')  
        end
      end      
    end
    % Manual CLim's
    c_eval('hmap(?).CLim = 2*[-1 1];',2:4)
    c_eval('hmap(?).CLim = 1*[-1 1];',5)
    
    str_title = {sprintf('timestep = %g, time = %g, x_box = [%g, %g], z_box = [%g, %g], Bloc = [%.2f,%.2f,%.2f]',timestep,time,xlo,xhi,zlo,zhi,Bloc.x,Bloc.y,Bloc.z),'', varstrs{1}}; 
    hmap(1).Title.String = str_title;
    isub = 1;    
    % fxy
    ispecies = [1 3];
    strtitle = sprintf('species [%g %g]',ispecies);
    if 1 % fxy, all ions together
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
    ispecies = [2 4]; 
    strtitle = sprintf('species [%g %g]',ispecies);
    if 1 % fxy, all electrons together
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
    
    % fxz
    ispecies = [1 3];
    strtitle = sprintf('species [%g %g]',ispecies);
    if 1 % fxz, all ions together
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
    ispecies = [2 4]; 
    strtitle = sprintf('species [%g %g]',ispecies);
    if 1 % fxz, all electrons together
      hca = h(isub); isub = isub + 1;
      toplot = squeeze(sum(fxz(:,:,ispecies),3));
      if doLogF, toplot = log10(toplot); end
      if doContourF      
        contourf(hca,vx(:,ispecies(1)),vz(:,ispecies(1)),toplot')
      else
        imagesc(hca,vx(:,ispecies(1)),vz(:,ispecies(1)),toplot')
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
    
    % fzy
    ispecies = [1 3];
    strtitle = sprintf('species [%g %g]',ispecies);
    if 1 % fzy, all ions together
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
    ispecies = [2 4]; 
    strtitle = sprintf('species [%g %g]',ispecies);
    if 1 % fzy, all electrons together
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
      for ipanel = 1:2 % xy        
        line_slope = (Bloc.y/Bloc.x);
        xx = min(h(ipanel).XLim(2)*[1 1/abs(line_slope)])*[-1 1];
        hold(h(ipanel),'on')                
        hBline = plot(h(ipanel),xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
        hold(h(ipanel),'off')
      end   
      for ipanel = 3:4 % xz     
        line_slope = (Bloc.z/Bloc.x);
        xx = min(h(ipanel).XLim(2)*[1 1/abs(line_slope)])*[-1 1];
        hold(h(ipanel),'on')                
        hBline = plot(h(ipanel),xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
        hold(h(ipanel),'off')
      end   
      for ipanel = 5:6 % yz
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
      %axis(h(ipanel),'equal')
      h(ipanel).XTick = -15:5:15;
      h(ipanel).YTick = -15:5:15;
      h(ipanel).XGrid = 'on';
      h(ipanel).YGrid = 'on';
      %h(ipanel).XMinorGrid = 'on';
      %h(ipanel).YMinorGrid = 'on';
      h(ipanel).XDir = 'normal';
      h(ipanel).YDir = 'normal';
      axis(h(ipanel),'square')      
      colormap(distCmap)
    end    
    linkprop(h([1 3 5]),'CLim'); % link CLim's
    linkprop(h([2 4 6]),'CLim'); % link CLim's
    print_path = [savedir_root str_timestep print_subdir];
    if ~exist(print_path,'dir'), mkdir(print_path); end
    print('-dpng','-r200',[print_path strprint '.png']);
  end  
  
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
        doLogF = 1;
        doBDir = 1;

        cPatch = 'k';
        distCmap = flipdim(pic_colors('blue_red'),1); distCmap = distCmap(round(end/2):end,:);
        distCmap = irf_colormap('waterfall');
        distCmap = pic_colors('candy');    
        varstrs = {'B.z','ve1.x','log10(pi2.scalar./pi1.scalar)','log10(ti2.scalar./ti1.scalar)','log10(ni2./ni1)'};

        xlim = torow([-100 0]);
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
          fig.Position = [screensize(1) screensize(2) screensize(3)*0.6 screensize(4)*0.7];

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
              vlim = 3;
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
        print('-dpng','-r200',[print_path strprint '.png']);        
      end
    end
  toc
end
 
%% Read and plot distributions, regional spatial map
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/distributions/';
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/distributions/';
txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber); % michael's perturbation
timestep = 05000;
str_timestep = sprintf('%05.0f',timestep);
txttime = sprintf('timestep = %05.0f',timestep); 
    
fontsize = 7;
doBDir = 1;
ticks = -15:1:15;

xlim = [160 180] + 1*[-1 1]; % 08000
zlim = [0 6] + 1*[-1 1]; % 08000
fig_position =  [1 330 2555 1015];

xlim = [184 206] + 0.5*[-1 1]; % 05000
zlim = [0 3.5] + 0.5*[-1 1]; % 05000
fig_position = [1 712 2555 633];

vlim = [2 5 2 5]; % ion, electron, ion, electron
ncols = diff(xlim);
nrows = diff(zlim);

idist = 0;
%tic
h = [];
for distnumber = 1:150%50%:281%39%30:40%39%180:200%:250%:100%:100%:10%40%:40%:4%:40
  %read_sub_dir = '/1/';
  read_sub_dir = '';
  %txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber);
  txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber);
  if not(exist(txtfile,'file'))
    warning(sprintf('File not found: %s', txtfile))
    continue
  end  
    
  idist = idist + 1;
  %idist = distnumber;
  
  % Load data  
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
      = read_distributions(txtfile);  
  if not(all(xlo>xlim(1) && xhi<xlim(2) && zlo>zlim(1) && zhi<zlim(2)))
    disp([sprintf('%.2f %.2f %.2f %.2f outside of box',xlo,xhi,zlo,zhi)])
    continue
  end
  disp(sprintf('%.2f ',xlo,xhi,zlo,zhi))
  axes_position = [(xlo-xlim(1))/diff(xlim) ...
                   (zlo-zlim(1))/diff(zlim) ...
                   (xhi-xlo)/diff(xlim) ...
                   (zhi-zlo)/diff(zlim)];
  
  vx = axes;
  vy = axes;
  vz = axes;
  xloc = 0.5*(xlo+xhi); 
  zloc = 0.5*(zlo+zhi);
  xlocind = find_closest_ind(x,xloc);
  zlocind = find_closest_ind(z,zloc);
  Bloc.x = B.x(xlocind,zlocind);
  Bloc.y = B.y(xlocind,zlocind);
  Bloc.z = B.z(xlocind,zlocind);
  
  %vpeaks = find_dist_peaks(fxyz(:,:,:,ispecies));
  pause(0.1)  
  
  if 1 % xz plane
    nrows = 1;
    ncols = 3;
    npanels = nrows*ncols;
    
    for ispecies = [2 4]
      nfig = 400 + ispecies;
      fig = figure(nfig);
      fig.Position = fig_position;
      hca = subplot('Position',axes_position);
      h(end+1) = hca;
      
      hca = gca;
   
   %   strtitle = sprintf('species %g\n%s\nx = [%g, %g], z = [%g, %g]',ispecies,txttime,xlo,xhi,zlo,zhi);
   %   strprint = sprintf('dist_%04.0f_species_%g_t_%05.0f_x_%g_%g_z_%g_%g',ispecies,distnumber,timestep,xlo,xhi,zlo,zhi);
      if 0 % fxy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vx(:,ispecies),vy(:,ispecies),squeeze(fxy(:,:,ispecies))')
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vy';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
      end
      if 1 % fxy
        imagesc(hca,vx(:,ispecies),vy(:,ispecies),squeeze(fxy(:,:,ispecies))')  
        hca.XLabel.String = '';
        hca.YLabel.String = '';
        %hcb = colorbar('peer',hca);
        %hcb.YLabel.String = 'f';
        %hca.Title.String = strtitle;
        hca.YDir = 'normal';        
        hca.XLim = vlim(ispecies)*[-1 1];
        hca.YLim = vlim(ispecies)*[-1 1];
        hca.XGrid = 'on';
        hca.YGrid = 'on';
        hca.XTick = ticks;
        hca.YTick = ticks;
        %if 
        hca.XTickLabel = [];
        hca.YTickLabel = [];
        hca.Box = 'on';
        colormap(hca,pic_colors('candy'))
        irf_legend(hca,{sprintf('x=%.1f, z=%.1f',xloc,zloc);sprintf('B=[%.2f,%.2f,%.2f]',Bloc.x,Bloc.y,Bloc.z)},[0.01 0.99],'color',[0 0 0],'fontsize',fontsize)
        
        if doBDir           
          line_slope = (Bloc.y/Bloc.x);
          xx = min(hca.XLim(2)*[1 1/abs(line_slope)])*[-1 1];
          hold(hca,'on')                
          hBline = plot(hca,xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
          hold(hca,'off')
        end
      end
      if 0 % fxz
        imagesc(hca,vx(:,ispecies),vz(:,ispecies),squeeze(fxz(:,:,ispecies))')  
        hca.XLabel.String = '';
        hca.YLabel.String = '';
        %hcb = colorbar('peer',hca);
        %hcb.YLabel.String = 'f';
        %hca.Title.String = strtitle;
        hca.YDir = 'normal';        
        hca.XLim = vlim(ispecies)*[-1 1];
        hca.YLim = vlim(ispecies)*[-1 1];
        hca.XGrid = 'on';
        hca.YGrid = 'on';
        hca.XTick = ticks;
        hca.YTick = ticks;
        %if 
        hca.XTickLabel = [];
        hca.YTickLabel = [];
        hca.Box = 'on';
        colormap(hca,pic_colors('candy'))
        irf_legend(hca,{sprintf('x=%.1f, z=%.1f',xloc,zloc);sprintf('B=[%.2f,%.2f,%.2f]',Bloc.x,Bloc.y,Bloc.z)},[0.01 0.99],'color',[0 0 0],'fontsize',9)
        
        if doBDir           
          line_slope = (Bloc.z/Bloc.x);
          xx = min(hca.XLim(2)*[1 1/abs(line_slope)])*[-1 1];
          hold(hca,'on')                
          hBline = plot(hca,xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
          hold(hca,'off')
        end
      end
      if 0 % fzy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vy(:,ispecies),vz(:,ispecies),squeeze(fyz(:,:,ispecies))')  
        hca.XLabel.String = 'vz';
        hca.YLabel.String = 'vy';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
      end
      %print('-dpng','-r200',[savedir_root sub_dir '/' strprint '.png']);
      drawnow
      %pause(1)
    end
  end
  %toc
end

%% Read and plot distributions, regional spatial map, h5
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/distributions/';
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/distributions/';
%txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber); % michael's perturbation
timestep = 05000;
str_timestep = sprintf('%05.0f',timestep);
txttime = sprintf('timestep = %05.0f',timestep); 
    
fontsize = 7;
doBDir = 1;
ticks = -15:1:15;

xlim = [160 180] + 1*[-1 1]; % 08000
zlim = [0 6] + 1*[-1 1]; % 08000
fig_position =  [1 330 2555 1015];

xlim = [184 206] + 0.5*[-1 1]; % 05000
zlim = [0 3.5] + 0.5*[-1 1]; % 05000
fig_position = [1 712 2555 633];

vlim = [2 5 2 5]; % ion, electron, ion, electron
ncols = diff(xlim);
nrows = diff(zlim);

idist = 0;
%tic
iSpecies = 1;
h = [];
it = 1;
for id = dst.indices{it}
  %id
  %read_sub_dir = '/1/';
  read_sub_dir = '';
  %txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber);
  txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,id);
  %if not(exist(txtfile,'file'))
  %  warning(sprintf('File not found: %s', txtfile))
  %  continue
  %end  
    
  idist = idist + 1;
  %idist = distnumber;
%  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
%      = read_distributions(txtfile);  
    
  % Load data  
  f = dst.fxyz(it,id,iSpecies);
  xlo = f.x(1);
  xhi = f.x(2);
  zlo = f.z(1);
  zhi = f.z(2);
  
  if not(all(xlo>xlim(1) && xhi<xlim(2) && zlo>zlim(1) && zhi<zlim(2)))
    disp([sprintf('%.2f %.2f %.2f %.2f outside of box',xlo,xhi,zlo,zhi)])
    continue
  end
  disp(sprintf('%.2f ',xlo,xhi,zlo,zhi))
  axes_position = [(xlo-xlim(1))/diff(xlim) ...
                   (zlo-zlim(1))/diff(zlim) ...
                   (xhi-xlo)/diff(xlim) ...
                   (zhi-zlo)/diff(zlim)];
  
  vx = f.vx;
  vy = f.vy;
  vz = f.vz;
  xloc = 0.5*(xlo+xhi); 
  zloc = 0.5*(zlo+zhi);
  %xlocind = find_closest_ind(x,xloc);
  %zlocind = find_closest_ind(z,zloc);
  %Bloc.x = B.x(xlocind,zlocind);
  %Bloc.y = B.y(xlocind,zlocind);
  %Bloc.z = B.z(xlocind,zlocind);
  
  %vpeaks = find_dist_peaks(fxyz(:,:,:,ispecies));
  pause(0.1)  
  
  if 1 % xz plane
    nrows = 1;
    ncols = 3;
    npanels = nrows*ncols;
    
    for ispecies = [2 4]
      nfig = 400 + ispecies;
      fig = figure(nfig);
      fig.Position = fig_position;
      hca = subplot('Position',axes_position);
      h(end+1) = hca;
      
      hca = gca;
   
   %   strtitle = sprintf('species %g\n%s\nx = [%g, %g], z = [%g, %g]',ispecies,txttime,xlo,xhi,zlo,zhi);
   %   strprint = sprintf('dist_%04.0f_species_%g_t_%05.0f_x_%g_%g_z_%g_%g',ispecies,distnumber,timestep,xlo,xhi,zlo,zhi);
      if 0 % fxy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vx(:,ispecies),vy(:,ispecies),squeeze(fxy(:,:,ispecies))')
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vy';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
      end
      if 1 % fxy
        imagesc(hca,f.vx,f.vy,squeeze(sum(f.f(:,:,:,ispecies),3))')  
        hca.XLabel.String = '';
        hca.YLabel.String = '';
        %hcb = colorbar('peer',hca);
        %hcb.YLabel.String = 'f';
        %hca.Title.String = strtitle;
        hca.YDir = 'normal';        
        hca.XLim = vlim(ispecies)*[-1 1];
        hca.YLim = vlim(ispecies)*[-1 1];
        hca.XGrid = 'on';
        hca.YGrid = 'on';
        hca.XTick = ticks;
        hca.YTick = ticks;
        %if 
        hca.XTickLabel = [];
        hca.YTickLabel = [];
        hca.Box = 'on';
        colormap(hca,pic_colors('candy'))
        irf_legend(hca,{sprintf('x=%.1f, z=%.1f',xloc,zloc);sprintf('B=[%.2f,%.2f,%.2f]',Bloc.x,Bloc.y,Bloc.z)},[0.01 0.99],'color',[0 0 0],'fontsize',fontsize)
        
        if doBDir           
          line_slope = (Bloc.y/Bloc.x);
          xx = min(hca.XLim(2)*[1 1/abs(line_slope)])*[-1 1];
          hold(hca,'on')                
          hBline = plot(hca,xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
          hold(hca,'off')
        end
      end
      if 0 % fxz
        imagesc(hca,vx(:,ispecies),vz(:,ispecies),squeeze(fxz(:,:,ispecies))')  
        hca.XLabel.String = '';
        hca.YLabel.String = '';
        %hcb = colorbar('peer',hca);
        %hcb.YLabel.String = 'f';
        %hca.Title.String = strtitle;
        hca.YDir = 'normal';        
        hca.XLim = vlim(ispecies)*[-1 1];
        hca.YLim = vlim(ispecies)*[-1 1];
        hca.XGrid = 'on';
        hca.YGrid = 'on';
        hca.XTick = ticks;
        hca.YTick = ticks;
        %if 
        hca.XTickLabel = [];
        hca.YTickLabel = [];
        hca.Box = 'on';
        colormap(hca,pic_colors('candy'))
        irf_legend(hca,{sprintf('x=%.1f, z=%.1f',xloc,zloc);sprintf('B=[%.2f,%.2f,%.2f]',Bloc.x,Bloc.y,Bloc.z)},[0.01 0.99],'color',[0 0 0],'fontsize',9)
        
        if doBDir           
          line_slope = (Bloc.z/Bloc.x);
          xx = min(hca.XLim(2)*[1 1/abs(line_slope)])*[-1 1];
          hold(hca,'on')                
          hBline = plot(hca,xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
          hold(hca,'off')
        end
      end
      if 0 % fzy
        hca = h(isub); isub = isub + 1;
        imagesc(hca,vy(:,ispecies),vz(:,ispecies),squeeze(fyz(:,:,ispecies))')  
        hca.XLabel.String = 'vz';
        hca.YLabel.String = 'vy';
        hcb = colorbar('peer',hca);
        hcb.YLabel.String = 'f';
        hca.Title.String = strtitle;
      end
      %print('-dpng','-r200',[savedir_root sub_dir '/' strprint '.png']);
      drawnow
      %pause(1)
    end
  end
  %toc
end


%% Find separate populations
vlim = [2 5 2 5]; % ion, electron, ion, electron
read_sub_dir = '/1/';
doPlotDists = 1;
for distnumber = 1:147%:140%1:147
  ispecies = 3;

  txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%05.0f/%s/%.0f.dat',timestep,read_sub_dir,distnumber); % michael's perturbation
  if not(exist(txtfile,'file'))
    error(sprintf('File not found: %s', txtfile))
  end  

  % Load data  
  [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
      = read_distributions(txtfile);  
  vx = axes;
  vy = axes;
  vz = axes;
  xloc = 0.5*(xlo+xhi);
  zloc = 0.5*(zlo+zhi);

  [XMAX,IMAX,~,~] = extrema2(fxz(:,:,ispecies));
  flim = max(max(fxz(:,:,ispecies)))*0.5e-1;
  rem_xmax = find(XMAX<flim);
  
  cl_XMAX = XMAX; cl_XMAX(rem_xmax) = [];
  cl_IMAX = IMAX; cl_IMAX(rem_xmax) = [];
  %[XMAX,IMAX] = clean_extrema2([XMAX,IMAX],'fmin',max(fxz(:,:,ispecies))*1e-2);
  [I,J] = ind2sub(size(fxyz(:,:,ispecies)),IMAX);
  [cl_I,cl_J] = ind2sub(size(fxyz(:,:,ispecies)),cl_IMAX);
  nmax = numel(IMAX);

  if doPlotDists
    if 1 % distributions
      hca = subplot(1,3,1);
      imagesc(hca,vx(:,ispecies),vz(:,ispecies),squeeze(fxz(:,:,ispecies))')  
      hca.XLabel.String = '';
      hca.YLabel.String = '';
      colormap(pic_colors('candy'))
      %hcb = colorbar('peer',hca);
      %hcb.YLabel.String = 'f';
      %hca.Title.String = strtitle;
      hca.YDir = 'normal';        
      hca.XLim = vlim(ispecies)*[-1 1];
      hca.YLim = vlim(ispecies)*[-1 1];
      hca.XGrid = 'on';
      hca.YGrid = 'on';
      axis(hca,'square')

      for imax = 1:nmax
        hold(hca,'on')

        % plot two colors/symbols, just to make sure they show on both white and colored background
        if XMAX(imax)<flim
          plot(hca,vx(I(imax),ispecies),vz(J(imax),ispecies),'w+') 
          plot(hca,vx(I(imax),ispecies),vz(J(imax),ispecies),'rx')      
        else
          plot(hca,vx(I(imax),ispecies),vz(J(imax),ispecies),'w+') 
          plot(hca,vx(I(imax),ispecies),vz(J(imax),ispecies),'ko')      
        end
        text(vx(I(imax),ispecies),vz(J(imax),ispecies),sprintf(' f=%g',XMAX(imax)))
        hold(hca,'off')                
      end
    end
    if 1 % vectors corresponding to distributions
      hca = subplot(1,3,[2 3]);  
      if numel(find(XMAX>flim))>5, nmax = 1; end
      for imax = 1:nmax
        if XMAX(imax)>flim
          hold(hca,'on')
          plot(hca,xloc,zloc,'k.')
          quiver(hca,xloc,zloc,vx(I(imax),ispecies),vz(J(imax),ispecies),0.5)
          hold(hca,'off')
        end
      end
    end
  end
          pause
end