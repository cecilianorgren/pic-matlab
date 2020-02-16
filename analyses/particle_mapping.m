localuser = datastore('local','user');
%tr04 = PICTraj('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5');

%% Get peaks of f

it = 2;
if 1
  iSpecies = [3];
  %ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([-0.1 0.6]);
%  ds = ds04(it).xlim([190 203]).zlim([-0.1 0.1]); % top row.
else % electrons
  iSpecies = [4];
  %ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([-0.1 0.6]);
 % ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([0.4 0.6]); % top row.

end
ds = ds04(2).zlim([-0.2 0.2]).dxlim([0 0.2]);
nPeaks = 3;
spacingPeaks = 2; % for ions its 0.2 vA
fpeaks = ds.get_peaks(nPeaks,spacingPeaks,iSpecies);
nDists = ds.nd;
doPlot = 1;
doPrint = 0;
if doPlot
  % plot results
  for id = ds.nd{1}:-1:1
    f = ds.f(1,id,iSpecies);
    
    figure(27)
    h = setup_subplots(1,3);
    
    
    hca = h(1);
    imagesc(hca,f.v,f.v,f.fxy')
    hca.YDir = 'normal';
    colormap(pic_colors('candy'))
    hold(hca,'on')
    plot(hca,[fpeaks(:,id).vx],[fpeaks(:,id).vy],'k.')
    for iPeak = 1:nPeaks
      text(hca,[fpeaks(iPeak,id).vx],[fpeaks(iPeak,id).vy],sprintf('%g',iPeak))
    end
    hold(hca,'off')
    hca.XGrid = 'on'; hca.YGrid = 'on';
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_y';
    hca.Title.String = sprintf('x = [%.1f,%.1f], z = [%.1f,%.1f]',ds.xi1{1}(id),ds.xi2{1}(id),ds.zi1{1}(id),ds.zi2{1}(id));
    

    hca = h(2);
    imagesc(hca,f.v,f.v,f.fxz')
    hca.YDir = 'normal';
    colormap(pic_colors('candy'))
    hold(hca,'on')
    plot(hca,[fpeaks(:,id).vx],[fpeaks(:,id).vz],'k.')
    for iPeak = 1:nPeaks
      text(hca,[fpeaks(iPeak,id).vx],[fpeaks(iPeak,id).vz],sprintf('%g',iPeak))
    end
    hold(hca,'off')
    hca.XGrid = 'on'; hca.YGrid = 'on';
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_z';
    

    hca = h(3);
    imagesc(hca,f.v,f.v,f.fyz')
    hca.YDir = 'normal';
    colormap(pic_colors('candy'))
    hold(hca,'on')
    plot(hca,[fpeaks(:,id).vy],[fpeaks(:,id).vz],'k.')
    for iPeak = 1:nPeaks
      text(hca,[fpeaks(iPeak,id).vy],[fpeaks(iPeak,id).vz],sprintf('%g',iPeak))
    end
    hold(hca,'off')
    hca.XGrid = 'on'; hca.YGrid = 'on';
    hca.XLabel.String = 'v_y';
    hca.YLabel.String = 'v_z';
    
    
    for ip = 1:3
      axis(h(ip),'square')
      h(ip).XTick = h(ip).YTick;
      h(ip).FontSize = 14;
       h(ip).Position(2) = 0.2;
      h(ip).Position(4) = 0.7;
    end
    hb = colorbar('peer',h(3));
    hb.Position(1) = h(3).Position(1) + h(3).Position(3) + 0.01;
    if doPrint
      cn.print(sprintf('fpeaks_z=0_id=%04.0f',id))
    end
    pause(0.2)
  end
end

%% Integrate trajectories based on fpeaks
pic = df04;
tspan = [60,160,210];
m = 1; 
q = 1;
tic
for iTr = 25:numel(fpeaks)  
  fprintf('iTr = %g/%g\n',iTr,numel(fpeaks))
  r0 = [fpeaks(iTr).x, fpeaks(iTr).y, fpeaks(iTr).z];
  v0 = [fpeaks(iTr).vx, fpeaks(iTr).vy, fpeaks(iTr).vz];
  tr_tmp = df04.integrate_trajectory(r0,v0,tspan,m,q);    
  [Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(tr_tmp.x,tr_tmp.z,tr_tmp.t);  % interpolate

  tr_tmp.Ex = Ex;
  tr_tmp.Ey = Ey;
  tr_tmp.Ez = Ez;
  tr_tmp.Bx = Bx;
  tr_tmp.By = By;
  tr_tmp.Bz = Bz;
  h5write_trajs('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr_tmp,fpeaks(iPeak,id))

  %tr(iPeak,id) = tr_tmp;
  toc
  %catch
  %  continue
  %end
end

%%
% Pick out trajectories that initially belonged to the vy<0 population
% needs adaption for when angle spans 180, because atan2d spans [-180,180]
trif = tr04.pass('x',[190 204],'z',[-0.25 0.25],'atan2d(vy,vx)',[-180 -90]); % This one returns no trajectories, which is surprising


tr = trif.lim('x',[190 204],'z',[-0.25 0.25]);

%% Plot trajectories (lots of information)

vlim = [-2 2];

for itr = 1:numel(trif)
  tr = trif(itr);
  h = setup_subplots(4,2,'vertical');
  isub = 1;
  cmap = pic_colors('waterfall');
  t0_msize = 10;
  tr = tr(1);
  legloc = 'eastoutside';
  %cmap = interp1(linspace(1,64,size(cmap,1)),cmap,1:numel(tr.t)); 
  if 0 % xyz(t)
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t,tr.x-200,tr.t,tr.y,tr.t,tr.z)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'r';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    legend(hca,{'x-200','y','z'},'location',legloc)
  end
  if 1 % xz, colorcoded by time
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.x,tr.z)
    hold(hca,'on')
    plot(hca,tr.x0,tr.z0,'ko','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter(hca,tr.x,tr.z,1,tr.t)
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'twci';
    colormap(hca,cmap)
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'z';
    
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    
  end
  if 1 % x,vx, colorcoded by time
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.x,tr.vx)
    hold(hca,'on')
    plot(hca,tr.x0,tr.vx0,'ko','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter(hca,tr.x,tr.vx,1,tr.t)
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'twci';
    colormap(hca,cmap)
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'v_x';
    
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.YLim = vlim;
  end
  if 1 % x,vy, colorcoded by time
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.x,tr.vy)
    hold(hca,'on')
    plot(hca,tr.x0,tr.vy0,'ko','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter(hca,tr.x,tr.vy,1,tr.t)
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'twci';
    colormap(hca,cmap)
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'v_y';
    
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.YLim = vlim;
  end
  if 1 % x,vz, colorcoded by time
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.x,tr.vz)
    hold(hca,'on')
    plot(hca,tr.x0,tr.vz0,'ko','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter(hca,tr.x,tr.vz,1,tr.t)
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'twci';
    colormap(hca,cmap)
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'v_z';
    
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.YLim = vlim;
  end
  if 1 % x,vy, colorcoded by z
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.x,tr.vy)
    hold(hca,'on')
    
    doDist = 1;
    if doDist
      pcolor(hca,f35_x_vy.x,f35_x_vy.v(1,:),f35_x_vy.f'); 
      shading(hca,'flat')
      hca.XLabel.String = 'x/d_i';
      hca.YLabel.String = 'v_y/v_A';
      hold(hca,'on')
      hca.CLim = [0 max(f35_x_vy.f(:))];
      hb = colorbar('peer',hca);  
      hb.YLabel.String = 'f(x,vy)';
    end
    plot(hca,tr.x0,tr.vy0,'ko','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter(hca,tr.x,tr.vy,1,tr.z)
    hold(hca,'off')
    
    colormap(hca,cmap)
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'v_y';
    
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    %hca.YLim = vlim;
    hold(hca,'off')
  end
  if 0 % xyz, colorcoded by time
    hca = h(isub); isub = isub + 1;
    plot3(hca,tr.x,tr.y,tr.z)
    hold(hca,'on')
    plot3(hca,tr.x0,tr.y0,tr.z0,'ko','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter3(hca,tr.x,tr.y,tr.z,1,tr.t)
    hold(hca,'off')
    hb = colorbar('peer',hca);
    colormap(hca,cmap)
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'y';
    hca.ZLabel.String = 'z';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.ZGrid = 'on';
  end
  if 1 % (vx,vy), colorcoded by time
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.vx,tr.vy,'k')
    hold(hca,'on')
    plot(hca,tr.vx0,tr.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter(hca,tr.vx,tr.vy,20,tr.t,'filled')
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'twci';
    colormap(hca,cmap)
    hca.XLabel.String = 'vx';
    hca.YLabel.String = 'vy';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XLim = vlim;
    hca.YLim = vlim;
    hca.PlotBoxAspectRatio = [1 1 1];
  end
  if 1 % (vx,vz), colorcoded by time
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.vx,tr.vz,'k')
    hold(hca,'on')
    plot(hca,tr.vx0,tr.vz0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter(hca,tr.vx,tr.vz,20,tr.t,'filled')
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'twci';
    colormap(hca,cmap)
    hca.XLabel.String = 'vx';
    hca.YLabel.String = 'vz';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XLim = vlim;
    hca.YLim = vlim;
    hca.PlotBoxAspectRatio = [1 1 1];
  end
  if 1 % (vy,vz), colorcoded by time
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.vy,tr.vz,'k')
    hold(hca,'on')
    plot(hca,tr.vy0,tr.vz0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
    hs = scatter(hca,tr.vy,tr.vz,20,tr.t,'filled');
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'twci';
    colormap(hca,cmap)
    hca.XLabel.String = 'vy';
    hca.YLabel.String = 'vz';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XLim = vlim;
    hca.YLim = vlim;
    hca.PlotBoxAspectRatio = [1 1 1];
  end
  if 0 % (vx,vy), colorcoded by x
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.vx,tr.vy,'k')
    hold(hca,'on')
    plot(hca,tr.vx0,tr.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter(hca,tr.vx,tr.vy,20,tr.x,'filled')
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'x/di';
    colormap(hca,cmap)
    hca.XLabel.String = 'vx';
    hca.YLabel.String = 'vy';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 0 % (vx,vz), colorcoded by x
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.vx,tr.vz,'k')
    hold(hca,'on')
    plot(hca,tr.vx0,tr.vz0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter(hca,tr.vx,tr.vz,20,tr.x,'filled')
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'x/di';
    colormap(hca,cmap)
    hca.XLabel.String = 'vx';
    hca.YLabel.String = 'vz';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 0 % (vy,vz), colorcoded by x
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.vy,tr.vz,'k')
    hold(hca,'on')
    plot(hca,tr.vy0,tr.vz0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
    hs = scatter(hca,tr.vy,tr.vz,20,tr.x,'filled');
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'x/di';
    colormap(hca,cmap)
    hca.XLabel.String = 'vy';
    hca.YLabel.String = 'vz';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 0 % (vx,vy), colorcoded by z
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.vx,tr.vy,'k')
    hold(hca,'on')
    plot(hca,tr.vx0,tr.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter(hca,tr.vx,tr.vy,20,tr.z,'filled')
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'z/di';
    colormap(hca,cmap)
    hca.XLabel.String = 'vx';
    hca.YLabel.String = 'vy';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 0 % (vx,vz), colorcoded by z
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.vx,tr.vz,'k')
    hold(hca,'on')
    plot(hca,tr.vx0,tr.vz0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
    scatter(hca,tr.vx,tr.vz,20,tr.z,'filled')
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'z/di';
    colormap(hca,cmap)
    hca.XLabel.String = 'vx';
    hca.YLabel.String = 'vz';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 0 % (vy,vz), colorcoded by z
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.vy,tr.vz,'k')
    hold(hca,'on')
    plot(hca,tr.vy0,tr.vz0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
    hs = scatter(hca,tr.vy,tr.vz,20,tr.z,'filled');
    hold(hca,'off')
    hb = colorbar('peer',hca);
    hb.YLabel.String = 'z/di';
    colormap(hca,cmap)
    hca.XLabel.String = 'vy';
    hca.YLabel.String = 'vz';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 0 % v 
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t,tr.vx,tr.t,tr.vy,tr.t,tr.vz)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'v';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    legend(hca,{'x','y','z'},'location',legloc)
  end
  if 0 % B
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t,tr.Bx,tr.t,tr.By,tr.t,tr.Bz)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'B';
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    legend(hca,{'x','y','z'},'location',legloc)
  end
  if 0 % E 
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t,tr.Ex,tr.t,tr.Ey,tr.t,tr.Ez)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'E';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    legend(hca,{'x','y','z'},'location',legloc)
  end
  if 0 % vxB
    hca = h(isub); isub = isub + 1;
    vxB = cross_product(tr.vx,tr.vy,tr.vz,tr.Bx,tr.By,tr.Bz);
    plot(hca,tr.t,vxB.x,tr.t,vxB.y,tr.t,vxB.z)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'vxB';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    legend(hca,{'x','y','z'},'location',legloc)
  end

  %[Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(xvtb.x,xvtb.z,xvtb.t);

  %c_eval('h(?).XLim = [100 210];',1:5)
  %compact_panels(0.01)
  
  %h(6).YLim = [-0.49 0.49];
  for ip = 1:numel(h)
    h(ip).FontSize = 14;
    %h(ip).Position(3) = 0.7;
  end
  for ip = []%1:3
    %h(ip).XLim = [-2 2];
    %h(ip).YLim = [-2 2];
    %axis(h(ip),'square')
    %h(ip).XTick = -15:0.5:15;
    %h(ip).YTick = -15:0.5:15;
  end
  %cn.print(sprintf('trajectories_%03.0f',itr),'path',[savedir '/forces/'])
  cn.print(sprintf('trajectories_overlaid_fxvy_i%03.0f',itr))
end

%% Calculate phase space 
zlim = 0+[-0.2 0.2];
xlim = [160 205];
its = 2;
twci = 160;
iss = [2];
icold = [3 5];

rdim = 1; % x
vdim = 2; % y
f1_x_vx = make_space_time_1d(ds04,its,xlim,zlim,rdim,1,1);
f1_x_vy = make_space_time_1d(ds04,its,xlim,zlim,rdim,2,1);
f1_x_vz = make_space_time_1d(ds04,its,xlim,zlim,rdim,3,1);
f35_x_vx = make_space_time_1d(ds04,its,xlim,zlim,rdim,1,icold);
f35_x_vy = make_space_time_1d(ds04,its,xlim,zlim,rdim,2,icold);
f35_x_vz = make_space_time_1d(ds04,its,xlim,zlim,rdim,3,icold);

%% Make phase space plots for several times
xvals = 150:0.2:200;
ds = ds04.dxlim([0 0.21]).zfind(0);
%f35 = ds.reduce_1d(xvals,0,linspace(-2,2,101),[3 5]);
f1 = ds.reduce_1d(xvals,0,linspace(-5,5,101),[1]);

ff = f1;
h = setup_subplots(5,3,'horizontal');
isub = 1;

for it = 1:5  
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,ff(it).x,vv,ff(it).fvx'); shading(hca,'flat'); 
  
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,ff(it).x,vv,ff(it).fvy'); shading(hca,'flat'); 

  hca = h(isub); isub = isub + 1;  
  pcolor(hca,ff(it).x,vv,ff(it).fvz'); shading(hca,'flat'); 
end

colormap(pic_colors('candy'))
hlinks = linkprop(h,{'CLim','XLim','YLim'});
hlinks.Targets(1).XLim = xvals([1 end]);

h(1).Title.String = 'f(x,z=0,vx)';
h(2).Title.String = 'f(x,z=0,vy)';
h(3).Title.String = 'f(x,z=0,vz)';

for ip = 13:15
  h(ip).XLabel.String = 'x/d_i';
end
compact_panels(0.01,0.02)


