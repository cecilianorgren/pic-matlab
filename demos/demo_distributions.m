distnumber = 1;
sub_dir = '/t8000_exhaust_1/';

txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%s/%.0f.dat',sub_dir,distnumber); % michael's perturbation
tic; [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
    = read_distributions(txtfile); toc
vx = axes;
vy = axes;
vz = axes;
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

%% Rotate into different coordinate system

%% Read and plot distributions
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/';
timestep = 08000;
txttime = sprintf('timestep = %05.0f',timestep); % michael's perturbation


for distnumber = 1:40%:40%:4%:40
  % Load data
  sub_dir = '/t8000_exhaust_1/';
  txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%s/%.0f.dat',sub_dir,distnumber); % michael's perturbation
  %txtfile = sprintf('/Volumes/Fountain/Data/PIC/paul_oxygen/%.0f.dat',distnumber); % michael's perturbation
  tic; [axes,xlo,xhi,zlo,zhi,ic,fxyz,fxy,fxz,fyz,vxa,vya,vza] ...
      = read_distributions(txtfile); toc
  vx = axes;
  vy = axes;
  vz = axes;
  
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
      strprint = sprintf('dist_%g_species_%g_t_%05.0f_x_%g_%g_z_%g_%g',ispecies,distnumber,timestep,xlo,xhi,zlo,zhi);
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
    
    strprint = sprintf('dist_%g_species_%s_t_%05.0f_x_%g_%g_z_%g_%g',distnumber,'1234',timestep,xlo,xhi,zlo,zhi);
    
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

    strprint = sprintf('dist_%g_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_map',distnumber,'1234',timestep,xlo,xhi,zlo,zhi);
    
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
  if 1 % 2x3 electrons and ions added together, respectively, plus map to show location
    doA = 0;
    doPatch = 1;
    doLogF = 1;
    
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
      strprint = sprintf('dist_%g_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_log10f_map',distnumber,'ie',timestep,xlo,xhi,zlo,zhi);
    else
      strprint = sprintf('dist_%g_species_%s_t_%05.0f_x_%g_%g_z_%g_%g_map',distnumber,'ie',timestep,xlo,xhi,zlo,zhi);
    end
    
    if 1 % Plot map showing location of dists
      hca = hmap;
      himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
      hca.YDir = 'normal';
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
      axis(h(ipanel),'square')      
      colormap(distCmap)
    end
    print('-dpng','-r200',[savedir_root sub_dir '/species_ie/' strprint '.png']);
  end  
end
 