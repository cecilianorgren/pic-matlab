%% Define times
timesteps = 200:200:5000;
ntimes = numel(timesteps);
savedir_root = '/Users/cno062/Research/PIC/df_cold_protons_1/';
screensize = get( groot, 'Screensize' );

%% Time loop
Uke1_ts = nan(1,ntimes); 
Uke2_ts = nan(1,ntimes); 
Uki1_ts = nan(1,ntimes); 
Uki2_ts = nan(1,ntimes); 
Ute1_ts = nan(1,ntimes); 
Ute2_ts = nan(1,ntimes); 
Uti1_ts = nan(1,ntimes); 
Uti2_ts = nan(1,ntimes); 
UB_ts = nan(1,ntimes); 
doTs = 1;
doPatch = 1;
for itime = 1:ntimes
  % Load data
  timestep = timesteps(itime);
  txtfile = sprintf('/Users/cno062/Data/PIC/df_cold_protons/data/fields-%05.0f.dat',timestep); % michael's perturbation
  tic; [x,z,E,B,...
        ni1,ne1,ni2,ne2,...
        vi1,ve1,vi2,ve2,...
        ji1,je1,ji2,je2,...
        pi1,pe1,pi2,pe2,...
        ti1,te1,ti2,te2,...
        dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q]... 
        = read_fields(txtfile); toc

  % Calculate auxillary quantities
  %tic; A = vector_potential(x,z,B.x,B.z); toc % vector potential
  pb = B.abs.^2/2; % magnetic pressure
  bcurv = magnetic_field_curvature(x,z,B.x,B.y,B.z); % magnetic curvature
  c_eval('ve?xB = cross_product(ve?.x,ve?.y,ve?.z,B.x,B.y,B.z);',1:2) % electron motional electric field
  c_eval('vi?xB = cross_product(vi?.x,vi?.y,vi?.z,B.x,B.y,B.z);',1:2) % ion motional electric field
  ExB = cross_product(E.x,E.y,E.z,B.x,B.y,B.z); % Poynting flux
  c_eval('E_ve?xB.x = E.x + ve?xB.x; E_ve?xB.y = E.y + ve?xB.y; E_ve?xB.z = E.z + ve?xB.z;',1:2) % electron motional electric field
  c_eval('E_vi?xB.x = E.x + vi?xB.x; E_vi?xB.y = E.y + vi?xB.y; E_vi?xB.z = E.z + vi?xB.z;',1:2) % electron motional electric field
  c_eval('je?E = je?.x.*E.x + je?.y.*E.y + je?.y.*E.z;',1:2)
  c_eval('ji?E = ji?.x.*E.x + ji?.y.*E.y + ji?.y.*E.z;',1:2)
  UB.x = 0.5*B.x.^2;
  UB.y = 0.5*B.y.^2;
  UB.z = 0.5*B.z.^2;
  UB.tot = 0.5*B.abs.^2;
  c_eval('Uke? = mass(2)*0.5*ne?.*(ve?.x.^2 + ve?.y.^2 + ve?.z.^2);',1:2)
  c_eval('Uki? = mass(1)*0.5*ni?.*(vi?.x.^2 + vi?.y.^2 + vi?.z.^2);',1:2)
  c_eval('Ute? = pe?.scalar;',1:2)
  c_eval('Uti? = pi?.scalar;',1:2)
  Uke = Uke1 + Uke2;
  Uki = Uki1 + Uki2;
  Uktot = Uke + Uki;
  Ute = Ute1 + Ute2;
  Uti = Uti1 + Uti2;
  Uttot = Ute + Uti;
  jtot.x = ji1.x + ji2.x - je1.x - je2.x;
  jtot.y = ji1.y + ji2.y - je1.y - je2.y;
  jtot.z = ji1.z + ji2.z - je1.z - je2.z;
  
  % Collect time series
  Uke1_ts(1,itime) = sum(Uke1(:));
  Uke2_ts(1,itime) = sum(Uke2(:));
  Uki1_ts(1,itime) = sum(Uki1(:));
  Uki2_ts(1,itime) = sum(Uki2(:));
  Ute1_ts(1,itime) = sum(Ute1(:));
  Ute2_ts(1,itime) = sum(Ute2(:));
  Uti1_ts(1,itime) = sum(Uti1(:));
  Uti2_ts(1,itime) = sum(Uti2(:));
  UBx_ts(1,itime) = sum(UB.x(:));
  UBy_ts(1,itime) = sum(UB.y(:));
  UBz_ts(1,itime) = sum(UB.z(:));
  UB_ts(1,itime) = sum(UB.tot(:));
  
  if 1 % Plot, define variable in cell array
    %%
    % Save and print info
    subdir = 'energy_density';
    savedir = [savedir_root,subdir];
    mkdir(savedir)
    savestr = sprintf('%s_t%05.0f',subdir,timestep);
    % Define what variables to plot
    %varstrs = {'ve1.x','ve2.x','ve1.z','ve2.z','ve1.par','ve2.par','-ve1xB.x','-ve2xB.x','-ve1xB.z','-ve2xB.z','E.x','E.z'};
    varstrs = {'UB','Uke1','Uke2','Uki1','Uki2','Ute1','Ute2','Uti1','Uti2'};
    nvars = numel(varstrs);

    % Initialize figure
    fig = figure(5);
    fig.Position = [screensize(1) screensize(2) screensize(3)*0.4 screensize(4)*0.7];
    npanels = nvars + doTs;
    nrows = 5;
    ncols = ceil(npanels/nrows);
    npanels = nrows*ncols;
    isub = 1; 
    for ipanel = 1:npanels  
      h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
    end

    doA = 0;
    if doA    
      cA = [0.8 0.8 0.8];
      nA = 20;
      nA = [0:-2:min(A(:))];
    end
    
    %% Plot part of data
    xlim = [x(1)+140 x(end)-140];
    zlim = [-5 5];
    ix1 = find(x>xlim(1),1,'first');
    ix2 = find(x<xlim(2),1,'last');
    iz1 = find(z>zlim(1),1,'first');
    iz2 = find(z<zlim(2),1,'last');
    ipx = ix1:2:ix2;
    ipz = iz1:2:iz2;
    
    % Panels
    isub = 1;
    if doTs % ts plot of energy
      hca = h(isub); isub = isub + 1;
      ts_varstrs = {'UB','Uke1','Uke2','Uki1','Uki2','Ute1','Ute2','Uti1','Uti2'};
      variables = nan(numel(ts_varstrs),ntimes);
      for ivar = 1:numel(ts_varstrs)
        variables(ivar,:) = eval([ts_varstrs{ivar} '_ts']);
      end
      plot(hca,timesteps/wpewce/mass(1),[UB_ts; Uke1_ts; Uke2_ts; Uki1_ts; Uki2_ts; Ute1_ts; Ute2_ts; Uti1_ts; Uti2_ts])     
      legend(hca,ts_varstrs,'location','eastoutside')
      labels = arrayfun(@(x,y) {[num2str(x) ' > Q_{||} > ' num2str(y)]}, edgesQ(end:-1:2),edgesQ(end-1:-1:1));
      hca.XLim = [0 (timesteps(end)+200)/wpewce/mass(1)];
      hca.XLabel.String = 'time (omega_{ci})';
      hca.YLabel.String = 'Energy density (...)';
    end
    for ivar = 1:nvars  
      hca = h(isub); isub = isub + 1;
      varstr = varstrs{ivar};
      variable = eval(varstr);  
      himag = imagesc(hca,x(ipx),z(ipz),variable(ipx,ipz)');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
      hca.Title.String = sprintf('%s',varstr); 
      hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      hcb = colorbar('peer',hca);
      %hcb.YLim = hca.CLim(2)*[-1 1];
      colormap(hca,cn.cmap('blue_red'));

      if doA
        hold(hca,'on')
        hcont = contour(hca,x(ipx),z(ipz),A(ipx,ipz)',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end

    for ipanel = 2:npanels
      h(ipanel).YDir = 'normal';
      h(ipanel).XLim = xlim;
      h(ipanel).YLim = zlim;
    end
    
    print('-dpng','-r200',[savedir '/' savestr '.png']);
  end
  if 0 % figure 1, n, vx, separate populations
    %%
    subdir = 'n_vx';
    savedir = [savedir_root,subdir];
    mkdir(savedir)
    savestr = sprintf('n_vx_t%05.0f',timestep);
    
    % Initialize figure
    fig = figure(1);
    fig.Position = screensize*0.8;
    nrows = 2;
    ncols = 1;
    npanels = nrows*ncols;
    isub = 1; 
    for ipanel = 1:npanels  
      h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
    end
    
    % Panels
    doA = 0;
    cA = [0.8 0.8 0.8];
    nA = 20;
    isub = 1;
    if 1 % ne1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ne1');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'n_{e1}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % ne2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ne2');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'n_{e2}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % ni1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ni1');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'n_{i1}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];
      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % ni2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ni2');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'n_{i2}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 1 % vex1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ve1.x');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'v_{e1x}';
      hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % vex2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ve2.x');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'v_{e2x}';
      hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % vix1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,vi1.x');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'v_{i1x}';
      hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % vix2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,vi2.x');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'v_{i2x}';
      hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end  
    
    %set(fig,'InvertHardCopy', 'off');
    %set(fig,'paperpositionmode','auto');
    tic; print('-dpng','-r150','-painters',[savedir '/' savestr '.png']); toc;
    %cn.print(savestr,'path',savedir)
  end
  if 0 % figure 1, n, vx, separate populations
    %%
    subdir = 'n_vx';
    savedir = [savedir_root,subdir];
    mkdir(savedir)
    savestr = sprintf('n_vx_t%05.0f',timestep);
    
    % Initialize figure
    fig = figure(2);
    fig.Position = screensize*0.8;
    nrows = 2;
    ncols = 1;
    npanels = nrows*ncols;
    isub = 1; 
    for ipanel = 1:npanels  
      h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
    end
    
    % Panels
    doA = 0;
    cA = [0.8 0.8 0.8];
    nA = 20;
    isub = 1;
    if 0 % ne1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ne1');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'n_{e1}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % ne2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ne2');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'n_{e2}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 1 % ni1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ni1');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'n_{i1}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];
      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % ni2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ni2');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'n_{i2}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % vex1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ve1.x');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'v_{e1x}';
      hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % vex2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,ve2.x');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'v_{e2x}';
      hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 1 % vix1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,vi1.x');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'v_{i1x}';
      hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % vix2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,vi2.x');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'v_{i2x}';
      hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end  
    
    %set(fig,'InvertHardCopy', 'off');
    %set(fig,'paperpositionmode','auto');
    tic; print('-dpng','-r150','-painters',[savedir '/' savestr '.png']); toc;
    %cn.print(savestr,'path',savedir)
  end
  if 0 % figure 2, energy densities, separate populations and total
    %%
    subdir = 'energy_densities';
    savedir = [savedir_root,subdir];
    mkdir(savedir)
    savestr = sprintf('energy_densities_t%05.0f',timestep);
    
    % Initialize figure
    fig = figure(3);
    fig.Position = screensize;
    nrows = 2;
    ncols = 4;
    npanels = nrows*ncols;
    isub = 1; 
    for ipanel = 1:npanels  
      h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
    end
    
    % Panels
    doA = 0;
    cA = [0.8 0.8 0.8];
    nA = 20;
    isub = 1;
    if 1 % Uke1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,Uke1');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{ek1}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % Uke2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,Uke2');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{ek2}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 1 % Uki1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,Uki1');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{ik1}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % Uki2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,Uki2');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{ik2}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 1 % Ute1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,Ute1');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{et1}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % Ute2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,Ute2');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{et2}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 1 % Uti1
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,Uti1');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{it1}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 0 % Uti2
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,Uti2');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{it2}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 1 % Uktot
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,Uktot');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{ktot}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 1 % Uttot
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,Uttot');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{ttot}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 1 % UB
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,UB');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{B}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end
    if 1 % UB + Uk + Ut
      hca = h(isub); isub = isub + 1;
      himag = imagesc(hca,x,z,UB' + Uttot' + Uktot');
      hca.XLabel.String = 'x (d_i)';
      hca.YLabel.String = 'z (d_i)';
      hca.Title.String = 'U_{B} + U_{ttot} + U_{ktot}';

      hcb = colorbar('peer',hca);
      colormap(hca,cn.cmap('blue_red'));
      hcb.YLim(1) = 0;
      hca.CLim = hca.CLim(2)*[-1 1];

      if doA
        hold(hca,'on')
        hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
        hold(hca,'off')  
      end
    end  
    
    %cn.print(savestr,'path',savedir)
    tic; print('-dpng','-r150','-painters',[savedir '/' savestr '.png']); toc;
  end
end
if 0 % Panels
  %% Plot 1, 4 species plasma properties, 1 species per column
  % Initialize figure
  nrows = 5;
  ncols = 4;
  npanels = nrows*ncols;
  isub = 1; 
  for ipanel = 1:npanels  
    h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1;  
  end

  % Panels
  doA = 0;
  cA = [0.8 0.8 0.8];
  nA = 20;
  isub = 1;
  if 1 % babs
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,B.abs');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '|B|'; 
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    hcb.YLim = hca.CLim(2)*[0 1];
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % bx
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,B.x');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'B_x';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % by
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,B.y');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'B_y';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % bz
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,B.z');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'B_z';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % epar
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,E.par');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'E_{||}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);  
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ex
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,E.x');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'E_x';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ey
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,E.y');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'E_y';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ez
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,E.z');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'E_z';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ne1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ne1');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'n_{e1}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ne2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ne2');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'n_{e2}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ni1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ni1');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'n_{i1}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];
    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ni2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ni2');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'n_{i2}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vepar1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve1.par');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{e||1}';

    hca.CLim = max(abs(E.par(:)))*[-1 1];  

    hcb = colorbar('peer',hca);  
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vepar2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve2.par');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{e||2}';

    hca.CLim = max(abs(ve2.par(:)))*[-1 1];  

    hcb = colorbar('peer',hca);  
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vipar1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi1.par');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{i||1}';

    hca.CLim = max(abs(vi1.par(:)))*[-1 1];  

    hcb = colorbar('peer',hca);  
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vipar2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi.par2');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{e||1}';

    hca.CLim = max(abs(vi.par2(:)))*[-1 1];  

    hcb = colorbar('peer',hca);  
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % vex1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve1.x');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{e1x}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % vex2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve2.x');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{e2x}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % vix1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi1.x');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{i1x}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % vix2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi2.x');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{i2x}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vey1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve1.y');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{e1y}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vey2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve2.y');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{e2y}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % viy1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi1.y');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{i1y}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % viy2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi2.y');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{i2y}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vez1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve1.z');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{e1z}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vez2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve2.z');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{e2z}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % viz1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi1.z');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{i1z}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % viz2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi2.z');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'v_{i2z}';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,xe,ze,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pxx e1 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pe1.xx');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{e1xx}';
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pxx e2 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pe2.xx');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{e2xx}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pxx i1 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pi1.xx');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{i1xx}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pxx i2 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pi2.xx');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{i2xx}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pyy e1 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pe1.yy');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{e1yy}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pyy e2 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pe2.yy');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{e2yy}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pyy i1 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pi2.yy');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{i1yy}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pyy i2 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pi1.yy');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{i2yy}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pzz e1 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pe1.zz');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{e1zz}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pzz e2 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pe2.zz');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{e2zz}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pzz i1 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pi1.zz');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{i1zz}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pzz i2 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pi2.zz');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{izz2}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pscalar e1 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pe1.scalar');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{e1scalar}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pscalar e2 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pe2.scalar');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{e2scalar}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pscalar i1 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pi1.scalar');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{i1scalar}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % pscalar i2 
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,pi2.scalar');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'p_{i2scalar}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb.YLim = [0 hca.CLim(2)];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ve1xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve1xB.x');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{e1}xB)_x';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ve2xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve2xB.x');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{e2}xB)_x';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vi1xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi1xB.x');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{i1}xB)_x';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vi2xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi2xB.x');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{i2}xB)_x';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ve1xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve1xB.y');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{e1}xB)_y';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ve2xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve2xB.y');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{e2}xB)_y';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vi1xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi1xB.y');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{i1}xB)_y';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vi2xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi2xB.y');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{i2}xB)_y';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ve1xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve1xB.z');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{e1}xB)_z';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ve2xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ve2xB.z');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{e2}xB)_z';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vi1xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi1xB.z');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{i1}xB)_z';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % vi2xB_x
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,vi2xB.z');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = '(v_{i2}xB)_z';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % je1E
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,je1E');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'j_{e1}\cdot E';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % je2E
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,je2E');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'je2E';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ji1E
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ji1E');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'ji1E';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 0 % ji2E
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,ji2E');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'ji1E';
    hca.CLim = max(abs(himag.CData(:)))*[-1 1];  
    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % Uke1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,Uke1');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{ek1}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % Uke2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,Uke2');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{ek2}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % Uki1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,Uki1');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{ik1}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % Uki2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,Uki2');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{ik2}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % Ute1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,Ute1');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{et1}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % Ute1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,Ute2');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{et2}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % Uti1
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,Uti1');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{it1}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % Uti2
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,Uti2');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{it2}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % Uktot
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,Uktot');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{ktot}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % Uttot
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,Uttot');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{ttot}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % UB
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,UB');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{B}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end
  if 1 % UB + Uk + Ut
    hca = h(isub); isub = isub + 1;
    himag = imagesc(hca,x,z,UB' + Uttot' + Uktot');
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = 'z (d_i)';
    hca.Title.String = 'U_{B} + U_{ttot} + U_{ktot}';

    hcb = colorbar('peer',hca);
    colormap(hca,cn.cmap('blue_red'));
    hcb.YLim(1) = 0;
    hca.CLim = hca.CLim(2)*[-1 1];

    if doA
      hold(hca,'on')
      hcont = contour(hca,x,z,A',nA,'color',cA,'linewidth',1.0);  
      hold(hca,'off')  
    end
  end


  for ipanel = 1:npanels
    h(ipanel).YDir = 'normal';
    h(ipanel).YLim = [-10 10];
  end
end