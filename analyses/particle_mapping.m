localuser = datastore('local','user');
%tr04 = PICTraj('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5');

%% Get r0,v0 from peaks of distributions

it = 2;
if 1
  iSpecies = [3 5];
  %ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([-0.1 0.6]);
%  ds = ds04(it).xlim([190 203]).zlim([-0.1 0.1]); % top row.
else % electrons
  iSpecies = [4];
  %ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([-0.1 0.6]);
 % ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([0.4 0.6]); % top row.

end
ds = ds04.twcilim(120).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([180:1:210]);
ds = ds04.twcilim(140).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([177:1:210]);
ds = ds04.twcilim(160).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([166:1:205]);
ds = ds04.twcilim(160).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([167:0.2:175]);
ds = ds04.twcilim(160).zlim(2+[-0.2 0.2]).dxlim([0 0.25]).xfind([165:0.2:175]);
ds = ds04.twcilim(160).zfind(2).dxlim([0 0.25]).xfind([169.8:0.2:170.2]);

nPeaks = 7;
spacingPeaks = 5; % for ions its 2 equals 0.2 vA
fpeaks = ds.get_peaks(nPeaks,spacingPeaks,iSpecies); % ,'vz',[-0.19 0.19]

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
      h(ip).XLim = [-1.5 1];
      h(ip).YLim = [-1 1];
    end
    hb = colorbar('peer',h(3));
    hb.Position(1) = h(3).Position(1) + h(3).Position(3) + 0.01;
    if doPrint
      cn.print(sprintf('fpeaks_z=0_id=%04.0f',id))
    end
    pause
  end
end

%% Get r0,v0 from moments, suitable for inflow at early times
pic = no02m.twpelim(16000);
A = squeeze(pic.A);
x_center = [202.2:0.2:202.9 203.2:0.2:203.9]-100;
z_center = (3:1:7);
x_center = [90:115];
z_center = (0:0.2:4);
[ZC,XC] = meshgrid(z_center,x_center);
XC = reshape(XC,prod(size(XC)),1);
ZC = reshape(ZC,prod(size(ZC)),1);
nboxes = numel(XC);
clear fpeaks
if 0 % A
  Alim = [-24 0];
  ind_keep = zeros(nboxes,1);
  for ibox = 1:nboxes
    xind = find(abs(pic.xi-XC(ibox))==min(abs(pic.xi-XC(ibox))));
    zind = find(abs(pic.zi-ZC(ibox))==min(abs(pic.zi-ZC(ibox))));

    if A(xind,zind)>Alim(1) && A(xind,zind)<Alim(2) % keep
      ind_keep(ibox) = 1;
    else
      ind_keep(ibox) = 0;
    end
  end
  keep_boxes = all_boxes(find(ind_keep==1),:);
  n_boxes = size(keep_boxes,1);
  x0 = XC(find(ind_keep));
  z0 = ZC(find(ind_keep));
else % n
  Alim = [0 8.2];
  vlim = [0 1];
  tlim = [0 0.0005];
  A = pic.A;
  t = pic.t(3);
  n = pic.n(3);
  vx = pic.vx(3);
  vy = pic.vy(3);
  vz = pic.vz(3);
  vabs = sqrt(vx.^2+vy.^2+vz.^2);
  nlim = 0.02;
  ind_keep = zeros(nboxes,1);
  for ibox = 1:nboxes
    xind = find(abs(pic.xi-XC(ibox))==min(abs(pic.xi-XC(ibox))));
    zind = find(abs(pic.zi-ZC(ibox))==min(abs(pic.zi-ZC(ibox))));

    x0(ibox) = XC(ibox);
    y0(ibox) = 0;
    z0(ibox) = ZC(ibox);
    vx0(ibox) = vx(xind,zind);
    vy0(ibox) = vy(xind,zind);
    vz0(ibox) = vz(xind,zind);
    
%     if n(xind,zind)>nlim(1) % keep
%       ind_keep(ibox) = 1;
%     else
%       ind_keep(ibox) = 0;
%     end
    if A(xind,zind)>Alim(1) && A(xind,zind)<Alim(2) % keep
      ind_keep(ibox) = 1;
%       if t(xind,zind)>tlim(2) % too fast discard
%         ind_keep(ibox) = 0;
%       else
%         ind_keep(ibox) = 1;
%       end
    else
      ind_keep(ibox) = 0;
    end
  end
  %x0 = XC(find(ind_keep));
  %z0 = ZC(find(ind_keep));
  x0 = x0(find(ind_keep));
  y0 = y0(find(ind_keep));
  z0 = z0(find(ind_keep));
  vx0 = vx0(find(ind_keep));
  vy0 = vy0(find(ind_keep));
  vz0 = vz0(find(ind_keep));
  for ii = 1:numel(x0)
    fpeaks(ii).x = x0(ii);
    fpeaks(ii).y = y0(ii);
    fpeaks(ii).z = z0(ii);
    fpeaks(ii).vx = vx0(ii);
    fpeaks(ii).vy = vy0(ii);
    fpeaks(ii).vz = vz0(ii);
  end
  hca = subplot(1,1,1);
  imagesc(hca,pic.xi,pic.zi,vx');
  colormap(hca,pic_colors('blue_red'))
  colorbar('peer',hca)
  hca.YDir = 'normal';
  hca.Title.String = sprintf('number of position = %g',numel(fpeaks));
  hold(hca,'on')
  plot(hca,[fpeaks.x],[fpeaks.z],'.k')
  hold(hca,'off')
end

%% Integrate trajectories based on fpeaks (which may be based or moments)
pic = df04;
t0 = 160; 
tspan = [50,t0,239];
%t0 = 60;
%tspan = [t0,239];
m = 1; 
q = 1;
istart = 1;
ntr_pre = trs.ntr-(istart-1); % the +(istart-1) is there becuse i already did the 1 and added it to the datafile

for iTr = istart:numel(fpeaks)
  tic
  fprintf('iTr = %g/%g\n',iTr,numel(fpeaks))
  r0 = [fpeaks(iTr).x, fpeaks(iTr).y, fpeaks(iTr).z];
  v0 = [fpeaks(iTr).vx, fpeaks(iTr).vy, fpeaks(iTr).vz];
  tr_tmp = df04.integrate_trajectory(r0,v0,tspan,m,q);
  
  hca = subplot(2,1,1);
  plot(hca,tr_tmp.x,tr_tmp.z,tr_tmp.x0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca = subplot(2,1,2);
  plot(hca,tr_tmp.y,tr_tmp.z,tr_tmp.y0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  drawnow
%   [Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(tr_tmp.x,tr_tmp.z,tr_tmp.t);  % interpolate
% 
%   tr_tmp.t0 = t0;
%   tr_tmp.Ex = Ex;
%   tr_tmp.Ey = Ey;
%   tr_tmp.Ez = Ez;
%   tr_tmp.Bx = Bx;
%   tr_tmp.By = By;
%   tr_tmp.Bz = Bz;
  %0
  h5write_trajs('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr_tmp,'fpeaks',fpeaks(iTr),'id',ntr_pre+iTr)
  %h5write_trajs('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr_tmp)
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
doLogRED = 0;
doLogDEF = 1;

twci = [100 120 140 160 180];
nt = numel(twci);
% new smaller boxes
xvals = 130:0.2:210;
dxi = [0 0.21]; 
% old larger boxes
%xvals = 150:1:210;
zvals = [0 1 2 3 4 5];
%dxi = [0.3 0.6]; 

ds = ds04.dxlim(dxi); % make so that twcilim work for eaxct values if there are more than 2
nt = numel(ds.twci);
%%
%f35 = ds.reduce_1d(xvals,zvals,linspace(-2,2,101),[3 5],'vabs');
%f1 = ds.reduce_1d(xvals,zvals,linspace(-5,5,101),[1],'vabs');
%f135 = ds.reduce_1d(xvals,zvals,linspace(-4,4,121),[1 3 5],'vabs');

ff = f1;

nrows = nt;
ncols = 5;
npanels = nrows*ncols;
toprow = 1:ncols;
bottomrow = (nrows-1)*ncols+[1:ncols];
leftcolumn = 1:ncols:(nrows*ncols);

isDEF = [];
isRED = [];

h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;

for it = 1:2%nt
  hca = h(isub); isub = isub + 1;  
  if doLogRED
    pcolor(hca,ff(it).x,ff(it).v,log10(ff(it).fvx')); shading(hca,'flat'); 
  else
    pcolor(hca,ff(it).x,ff(it).v,ff(it).fvx'); shading(hca,'flat'); 
  end
  isRED(end+1) = isub - 1;
  
  hca = h(isub); isub = isub + 1;  
  if doLogRED
    pcolor(hca,ff(it).x,ff(it).v,log10(ff(it).fvy')); shading(hca,'flat'); 
  else
    pcolor(hca,ff(it).x,ff(it).v,ff(it).fvy'); shading(hca,'flat'); 
  end
  isRED(end+1) = isub - 1;

  hca = h(isub); isub = isub + 1;  
  if doLogRED
    pcolor(hca,ff(it).x,ff(it).v,log10(ff(it).fvz')); shading(hca,'flat'); 
  else
    pcolor(hca,ff(it).x,ff(it).v,ff(it).fvz'); shading(hca,'flat'); 
  end
  isRED(end+1) = isub - 1;
  
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,ff(it).x,ff(it).v(ff(it).v>=0),log10(ff(it).def')); shading(hca,'flat'); 
  isDEF(end+1) = isub - 1;
  
  hca = h(isub); isub = isub + 1;  
  pcolor(hca,ff(it).x,ff(it).v(ff(it).v>=0).^2,log10(ff(it).def')); shading(hca,'flat'); 
  isDEF(end+1) = isub - 1;
  hca.YScale = 'log';
  hca.YLabel.String = 'log_{10}(v^2)';
  %hca = h(isub); isub = isub + 1;  
  %pcolor(hca,ff(it).x,ff(it).v(ff(it).v>=0),ff(it).fvabssum'); shading(hca,'flat'); 
    
end

colormap(pic_colors('candy'))

hlinksDEF = linkprop(h(isDEF),{'CLim'});
%hlinksDEF.Targets(1).YLim = 2.99*[0 1];

hlinksRED = linkprop(h(isRED),{'CLim','YLim'});
hlinksRED.Targets(1).YLim = 2.99*[-1 1];

hlinksALL = linkprop(h,{'XLim'});
hlinksALL.Targets(1).XLim = xvals([1 end]);


h(isRED(1)).Title.String = 'f(x,z=0,vx)';
h(isRED(2)).Title.String = 'f(x,z=0,vy)';
h(isRED(3)).Title.String = 'f(x,z=0,vz)';
if doLogDEF
  for ip = intersect(isDEF,toprow)
    h(ip).Title.String = 'log_{10}(def)';
  end
else
  for ip = intersect(isDEF,toprow)
    h(isDEF(1)).Title.String = 'def';
  end
end
%h(5).Title.String = 'f(x,z=0,|vz|(mean))';

for ip = bottomrow
  h(ip).XLabel.String = 'x/d_i';
end
for it = 1:nt 
  ip = leftcolumn(it);
  h(ip).YLabel.String = 'v/v_A';
  irf_legend(h(ip),sprintf('twci = %g',ds.twci(it)),[0.02 0.98],'fontsize',14)
end
for ip = 1:npanels
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).Layer = 'top';  
  h(ip).FontSize = 14;
end

compact_panels(0.01)

if not(isempty(isDEF)) % colorbar for def
  pos = h(isDEF(end)).Position;
  hb = colorbar('peer',h(isDEF(end)));
  hb.Position(1) = pos(1)+pos(3);
  hb.Position(4) = pos(4);
  h(isDEF(end)).Position = pos;
end

if not(isempty(isDEF)) % colorbar for red
  pos = h(isRED(end)).Position;
  hb = colorbar('peer',h(isRED(end)));
  hb.Position(1) = pos(1)+pos(3);
  hb.Position(4) = pos(4);
  h(isRED(end)).Position = pos;
end
%cn.print('f135')

%% Make phase space plots for several times, several z
doLogRED = 0;
doLogDEF = 1;

twci = [100 180];
nt = numel(twci);
% new smaller boxes
xvals = 150:0.2:210;
dxi = [0 0.21]; 
% old larger boxes
xvals = 150:1:210;
zvals = [0 1 2 3 4];
dxi = [0.3 0.6]; 

ds = ds04.twcilim(twci).dxlim(dxi).zfind(zvals); % make so that twcilim work for eaxct values if there are more than 2
nt = numel(ds.twci);
%
%f35 = ds.reduce_1d(xvals,zvals,linspace(-2,2,101),[3 5],'vabs');
%f1 = ds.reduce_1d(xvals,zvals,linspace(-5,5,101),[1],'vabs');
%f135 = ds.reduce_1d(xvals,zvals,linspace(-4,4,121),[1 3 5],'vabs');

ff = f5;

if 0 % Get X line location for each time
  xxLine = nan(ds.nt,1);
  for itime = 1:ds.nt
    A = df04.twcilim(ds.twci(itime)).A;
    [A_inds,A_vals] = saddle(A,'sort');  
    xxLine(itime) = df04.xi(A_inds(1,1));
  end
end

% Figure
[nrows,ncols] = size(ff);

npanels = nrows*ncols;
toprow = 1:ncols;
bottomrow = (nrows-1)*ncols+[1:ncols];
leftcolumn = 1:ncols:(nrows*ncols);

isDEF = [];
isRED = [];

h = setup_subplots(nrows,ncols,'horizontal');
h = reshape(h,nrows,ncols);
isub = 1;
h_empty = zeros(nrows,ncols);

for it = 1:nrows
  for iSpace = 1:ncols % ncols:-1:1
    fftmp = ff(it,iSpace);
    if isempty(fftmp.iter)
      h_empty(it,iSpace) = 1;
      continue; 
    end
    z = unique(fftmp.z);
      
    if 0 %f(vx)
      hca = h(isub); isub = isub + 1;  
      if doLogRED
        pcolor(hca,fftmp.x,fftmp.v,log10(fftmp.fvx')); shading(hca,'flat'); 
      else
        pcolor(hca,fftmp.x,fftmp.v,fftmp.fvx'); shading(hca,'flat'); 
      end
      isRED(end+1) = isub - 1;
    end 
    if 0 % f(vy)
      hca = h(isub); isub = isub + 1;  
      if doLogRED
        pcolor(hca,fftmp.x,fftmp.v,log10(fftmp.fvy')); shading(hca,'flat'); 
      else
        pcolor(hca,fftmp.x,fftmp.v,fftmp.fvy'); shading(hca,'flat'); 
      end
      isRED(end+1) = isub - 1;
    end
    if 0 % f(vz)
      hca = h(isub); isub = isub + 1;  
      if doLogRED
        pcolor(hca,fftmp.x,fftmp.v,log10(fftmp.fvz')); shading(hca,'flat'); 
      else
        pcolor(hca,fftmp.x,fftmp.v,fftmp.fvz'); shading(hca,'flat'); 
      end
      isRED(end+1) = isub - 1;
    end
    if 0 % def, on vscale
      hca = h(isub); isub = isub + 1;  
      pcolor(hca,fftmp.x,fftmp.v(fftmp.v>=0),log10(fftmp.def')); shading(hca,'flat'); 
      isDEF(end+1) = isub - 1;
    end
    if 1 % def on log10(v^2) scale
      hca = h(it,iSpace); isub = isub + 1;  
      pcolor(hca,fftmp.x,fftmp.v(fftmp.v>=0).^2,log10(fftmp.def')); shading(hca,'flat'); 
      isDEF(end+1) = isub - 1;
      hca.YScale = 'log';
      hca.YLabel.String = 'log_{10}(v^2)';
      hca.Title.String = sprintf('iteration = %g, z = %g',fftmp.iter,z);
      hb = colorbar('peer',hca);
    end
  %hca = h(isub); isub = isub + 1;  
  %pcolor(hca,ff(it).x,ff(it).v(ff(it).v>=0),ff(it).fvabssum'); shading(hca,'flat'); 
  end
end

h = h(not(h_empty));
colormap(pic_colors('candy'))
hlinksDEF = linkprop(h(:),{'CLim','YLim','YTick'});

hlinksDEF.Targets(1).CLim = [-9 -5];
hlinksDEF.Targets(1).YLim = [1e-2 2e1];
hlinksDEF.Targets(1).YTick = [1e-3 1e-2 1e-1 1e0 1e1 1e2];

%hlinksRED = linkprop(h(isRED),{'CLim','YLim'});
%hlinksRED.Targets(1).YLim = 2.99*[-1 1];

hlinksALL = linkprop(h,{'XLim'});
hlinksALL.Targets(1).XLim = xvals([1 end]);


%h(isRED(1)).Title.String = 'f(x,z=0,vx)';
%h(isRED(2)).Title.String = 'f(x,z=0,vy)';
%h(isRED(3)).Title.String = 'f(x,z=0,vz)';
%if doLogDEF
%  for ip = intersect(isDEF,toprow)
%    h(ip).Title.String = 'log_{10}(def)';
%  end
%else
  %for ip = intersect(isDEF,toprow)
  %  h(isDEF(1)).Title.String = 'def';
  %end
%end
%h(5).Title.String = 'f(x,z=0,|vz|(mean))';
