% Integrate an ensemble of particles starting at a given location, with a
% small variation in velocities. A bulk velocity with a randomly directed
% thermal velocity.

%% Pick pounts based on f.
t0 = 24000;
t0 = 24000;
it = 2;

iSpecies = [3];
% ds = ds04.twcilim(120).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([180:1:210]);
% ds = ds04.twcilim(140).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([177:1:210]);
% ds = ds04.twcilim(160).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([166:1:205]);
% ds = ds04.twcilim(160).zlim([-0.2 0.2]).dxlim([0 0.25]).xfind([167:0.2:175]);
% ds = ds04.twcilim(160).zlim(2+[-0.2 0.2]).dxlim([0 0.25]).xfind([165:0.2:175]);
ds = ds100.twpelim(t0).zfind(0).findtag({'line horizontal'}).xfind([70:0.5:92]);

nPeaks = 3;
spacingPeaks = 0; % for ions its 2 equals 0.2 vA
fpeaks = ds.get_peaks(nPeaks,spacingPeaks,iSpecies,'vz',[-2.0 0.0],'vy',[-2.0 0.7]); % ,'vz',[-0.19 0.19]
%fpeaks = ds.get_peaks(nPeaks,spacingPeaks,iSpecies); % ,'vz',[-0.19 0.19]

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
    imagesc(hca,f.v,f.v,log10(f.fxy)')
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
    imagesc(hca,f.v,f.v,log10(f.fxz)')
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
    imagesc(hca,f.v,f.v,log10(f.fyz)')
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
      h(ip).XLim = [-2.5 2];
      h(ip).YLim = [-2 2];
    end
    hb = colorbar('peer',h(3));
    hb.Position(1) = h(3).Position(1) + h(3).Position(3) + 0.01;
    if doPrint
      cn.print(sprintf('fpeaks_z=0_id=%04.0f',id))
    end
    pause
  end
end

%% Integrate and save

tstop = 125;
tspan = [90 120 tstop];
m = 1;
q = 1;

ntr_pre = tr100.ntr;
hca = subplot(1,1,1); plot(hca,NaN,NaN); hold(hca,'on'); drawnow;
% 1-131 in tr100 is from before

nf = numel(fpeaks);

for iTr = 101:numel(fpeaks)
  tic
  disp(sprintf('iTr = %g',iTr))
%   r0 = [x0,y0,z0];
%   v0 = [vx0,vy0,vz0] + vt0*randn(1,3)*0.05;
  fprintf('iTr = %g/%g\n',iTr,numel(fpeaks))
  r0 = [fpeaks(iTr).x, fpeaks(iTr).y, fpeaks(iTr).z];
  v0 = [fpeaks(iTr).vx, fpeaks(iTr).vy, fpeaks(iTr).vz];
  tr_tmp = no02m.integrate_trajectory(r0,v0,tspan,m,q);  
  Atmp = no02m.interp(tr_tmp.x,tr_tmp.z,tr_tmp.t,'A');
  tr_tmp.Ay = Atmp;
  h5write_trajs('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories.h5',tr_tmp,'id',ntr_pre+iTr)
  plot(hca,tr_tmp.x,tr_tmp.z)
  drawnow;
  toc
end


%% Pick points
x0 = 80;
y0 = 0;
z0 = 6; % Maybe this is too far up? They diverge quite a lot already before 
% hitting separatrix it seems. Try further down and a bit later maybe?
t0 = 75;
t0_ = t0;
%h = no02m.twpelim(14000).xlim([50 150]).zlim([0 15]).plot_map({'vz(3)','vt(3)'}','A',1,'sep','clim',{[-1 1],[-1 1]},'cmap',{cmapbr,cmapbr});

vx0 = no02m.interp(x0,z0,t0,'vx(3)');
vy0 = no02m.interp(x0,z0,t0,'vy(3)');
vz0 = no02m.interp(x0,z0,t0,'vz(3)');
vt0 = no02m.interp(x0,z0,t0,'vt(3)');

% Test how to initialize particles
% This seems ok
vhat = randn(100,3);
v0 = [vx0,vy0,vz0] + vt0*vhat*0.05;
vedges = 0.1*linspace(-1,1,200); 
[N edges mid loc] = histcn(v0,vedges,vedges,vedges);
hca = subplot(1,1,1);
imagesc(hca,mid{[1 3]},squeeze(sum(N(:,:,:),2))')
hca.YDir = 'normal';
hold(hca,'on')
plot(hca,vx0,vz0,'ko')
plot(hca,vx0+vt0*cosd(0:360),vz0+vt0*sind(0:360),'linewidth',1,'color','k')
hold(hca,'off')
axis(hca,'equal')
axis(hca,'square')

%%
% Trajectories

tstop = 125;
tspan = [t0_ tstop];
m = 1;
q = 1;

ntr_pre = tr100.ntr;
hca = subplot(1,1,1); plot(hca,NaN,NaN); hold(hca,'on'); drawnow;
% 1-131 in tr100 is from before

for iTr = 1:50
  disp(sprintf('iTr = %g',iTr))
  r0 = [x0,y0,z0];
  v0 = [vx0,vy0,vz0] + vt0*randn(1,3)*0.05;
  tr_tmp = no02m.integrate_trajectory(r0,v0,tspan,m,q);  
  Atmp = no02m.interp(tr_tmp.x,tr_tmp.z,tr_tmp.t,'A');
  tr_tmp.Ay = Atmp;
  h5write_trajs('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories.h5',tr_tmp,'id',ntr_pre+iTr)
  plot(hca,tr_tmp.x,tr_tmp.z)
  drawnow;
end

%% Pick points, along lines in x
x0 = [60:1:110];
y0 = 0*ones(size(x0));
z0 = 3.01*ones(size(x0)); % Maybe this is too far up? They diverge quite a lot already before 
% hitting separatrix it seems. Try further down and a bit later maybe?
t0_ = 70;
t0 = t0_*ones(size(x0));
%h = no02m.twpelim(14000).xlim([50 150]).zlim([0 15]).plot_map({'vz(3)','vt(3)'}','A',1,'sep','clim',{[-1 1],[-1 1]},'cmap',{cmapbr,cmapbr});

vx0 = no02m.interp(x0,z0,t0,'vx(3)');
vy0 = no02m.interp(x0,z0,t0,'vy(3)');
vz0 = no02m.interp(x0,z0,t0,'vz(3)');
vt0 = no02m.interp(x0,z0,t0,'vt(3)');
 
%%
% Trajectories
tstop = 125;
tspan = [t0_ tstop];
m = 1;
q = 1;

ntr_pre = tr100.ntr;
hca = subplot(1,1,1); plot(hca,NaN,NaN); hold(hca,'on'); drawnow;
% 1-131 in tr100 is from before
for iTr = 1:50  
  r0 = [x0(iTr),y0(iTr),z0(iTr)];
  v0 = [vx0(iTr),vy0(iTr),vz0(iTr)];
  tr_tmp = no02m.integrate_trajectory(r0,v0,tspan,m,q);  
  Atmp = no02m.interp(tr_tmp.x,tr_tmp.z,tr_tmp.t,'A');
  tr_tmp.Ay = Atmp;
  h5write_trajs('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories.h5',tr_tmp,'id',ntr_pre+iTr)
  plot(hca,tr_tmp.x,tr_tmp.z)
  drawnow;
end
%% Plot 
% (for now it's just an exmaple plot, we need also to make the approprtiate
% trajectories)

t0 = 70;
tend = 124;
trs = tr100.find([tr100.tstart] == t0,[tr100.x0] == 100,[tr100.z0] == 5);
trs = tr100.find([tr100.tstart] == t0,[tr100.x0] == 80,[tr100.z0] == 3);
trs = tr100.find([tr100.tstart] == t0,[tr100.x0] == 95,[tr100.z0] == 4);
trs = tr;
% need to make a resample function, to resample quantities to given times,
% so that I can plot them together as a function of time.
nt = 100;
timeline = linspace(t0,tend,nt);
trs_ = trs.resample(timeline); 

zz = [trs_.z];
zmean = mean(zz,2);

EEz = [trs_.Ez]; 
Ezmean = mean(EEz,2);
Ezrms = rms(EEz,2);
Ezprc50 = prctile(EEz,50,2);
Ezprc90 = prctile(EEz,90,2);
Ezprc = prctile(EEz,0:10:100,2);


vv = trs_.vabs; 
vv = [vv.v]; % (nt,ntr)
vmean = mean(vv,2);
vrms = rms(vv,2);
vprc10 = prctile(vv,10,2);
vprc50 = prctile(vv,50,2);
vprc90 = prctile(vv,90,2);
vprc = prctile(vv,0:10:100,2);
vedges = 0:0.05:4;

nrows = 2; ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 0 % (x,z)
  hca = h(isub); isub = isub + 1;
  trs.plot_all_xz(hca)  
  hold(hca,'off')
  hca.XLim = [50 110];
  hca.YLim = [-7 7];
end
if 0 % (x,z), color Ez
  hca = h(isub); isub = isub + 1;
  %trs.plot_all_xz(hca,'color','Ez');
  trs.plot_all('x','z',hca,'color','Ex');
  hold(hca,'off')
  hca.XLim = [50 110];
  hca.YLim = [-7 7];
end
if 0 % (x,z), color Ex
  hca = h(isub); isub = isub + 1;
  %trs.plot_all_xz(hca,'color','Ez');
  trs.plot_all('x','z',hca,'color','Ex');
  hold(hca,'off')
  hca.XLim = [50 110];
  hca.YLim = [-7 7];
end
if 0 % (x,y)
  hca = h(isub); isub = isub + 1;
  trs.plot_all_xy(hca)  
  hca.XLim = [50 110];
end
if 0 % |v|
  hca = h(isub); isub = isub + 1;
  % [N1,EDGES1] = histcounts(vv(1,:),edges);
  % [N2,EDGES2] = histcounts(vv(end,:),edges);
  histogram(hca,vv(1,:),vedges,'displayname','t=70')
  hold(hca,'on')
  histogram(hca,vv(end,:),vedges,'displayname',sprintf('t=%g',tend))
  histogram(hca,vv(fix(end/2),:),vedges,'displayname',sprintf('t=%g',t0+(tend-t0)/2))
  hold(hca,'off')
  legend(hca)
  hca.XLabel.String = '|v|';
  hca.YLabel.String = '# part';
end
if 1 % vmean(t), vprc(t)
  hca = h(isub); isub = isub + 1;
  colors = pic_colors('matlab');
  colors = repmat(colors(1,:),10,1);
  ic = 1;
  plot(hca,timeline,vmean,'linewidth',4,'color',colors(ic,:)); ic = ic + 1;
  hold(hca,'on')
  plot(hca,timeline,vprc,'linewidth',0.5); ic = ic + 1;
  %plot(hca,timeline,vprc50,timeline,vprc50,'linewidth',2,'color',colors(ic,:)); ic = ic + 1;
  %plot(hca,timeline,vprc90,timeline,vprc90,'linewidth',0.5,'color',colors(ic,:)); ic = ic + 1;
  hold(hca,'off')
  hca.XLabel.String = 't';
  hca.YLabel.String = '|v|';
  hca.XLim = [70 125];
end
if 1 % vmean(z), vprc(z)
  hca = h(isub); isub = isub + 1;
  colors = pic_colors('matlab');
  colors = repmat(colors(1,:),10,1);
  ic = 1;
  plot(hca,zmean,vmean,'linewidth',4,'color',colors(ic,:)); ic = ic + 1;
  hold(hca,'on')
  plot(hca,zmean,vprc,'linewidth',0.5); ic = ic + 1;
  %plot(hca,zmean,vprc50,zmean,vprc50,'linewidth',2,'color',colors(ic,:)); ic = ic + 1;
  %plot(hca,zmean,vprc90,zmean,vprc90,'linewidth',0.5,'color',colors(ic,:)); ic = ic + 1;
  hold(hca,'off')
  hca.XLabel.String = 'z';
  hca.YLabel.String = '|v|';
end
if 1 % Ezmean(t), Ezprc(t)
  hca = h(isub); isub = isub + 1;
  colors = pic_colors('matlab');
  colors = repmat(colors(1,:),10,1);
  ic = 1;
  plot(hca,timeline,Ezmean,'linewidth',4,'color',colors(ic,:)); ic = ic + 1;
  hold(hca,'on')
  plot(hca,timeline,Ezprc,'linewidth',0.5);
  %plot(hca,timeline,Ezmean+Ezprc90,timeline,Ezmean-Ezprc90,'linewidth',0.5,'color',colors(ic,:)); ic = ic + 1;
  hold(hca,'off')
  hca.XLabel.String = 't';
  hca.YLabel.String = 'E_z';
  hca.XLim = [70 125];
end
if 1 % Ezmean(z), Ezprc(z)
  hca = h(isub); isub = isub + 1;
  colors = pic_colors('matlab');
  colors = repmat(colors(1,:),10,1);
  ic = 1;
  plotc(hca,zmean,Ezmean,'linewidth',4,'color',colors(ic,:)); ic = ic + 1;
  hold(hca,'on')
  plot(hca,zmean,Ezprc,'linewidth',0.5);
  %plot(hca,timeline,Ezmean+Ezprc90,timeline,Ezmean-Ezprc90,'linewidth',0.5,'color',colors(ic,:)); ic = ic + 1;
  hold(hca,'off')
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'E_z';
  %hca.XLim = [70 125];
end

for ip = 1:npanels
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end
         %%
% test to se ehow messy the plot becomes
hold(hca,'on')
colors = pic_colors('matlab');
colors = repmat(colors(2,:),10,1);
ic = 1;
plot(hca,timeline,vmean+0.5,'linewidth',4,'color',colors(ic,:)); ic = ic + 1;
plot(hca,timeline,vmean+vprc50+0.5,timeline,vmean-vprc50+0.5,'linewidth',2,'color',colors(ic,:)); ic = ic + 1;
plot(hca,timeline,vmean+vprc90+0.5,timeline,vmean-vprc90+0.5,'linewidth',0.5,'color',colors(ic,:)); ic = ic + 1;
hold(hca,'off')

hca.XGrid = 'on';
hca.YGrid = 'on';

%% Movie of ion temperature + trajectories
tr = tr100.find([tr100.x0]==95,[tr100.z0]==4);
tr = tr100.find([tr100.tstart] == 90,[tr100.x0] == 95,[tr100.z0] == 4);
tr = tr100.find([tr100.tstart] == 90,[tr100.x0] == 90,[tr100.z0] == 5);
tr = tr100.find([tr100.tstart] == 70,[tr100.z0] == 3.01);
%tr = tr100;
%tr = tr(1:5);
igroup = 7;
tr = tr100.find([tr100.t0]==unique_starts(igroup,1),[tr100.x0]==unique_starts(igroup,2),[tr100.z0]==unique_starts(igroup,3));
tr = tr100.find([tr100.tstart] == 75,[tr100.x0] == 80,[tr100.z0] == 6);
%tr = tr100.find([tr100.tstart] == 70,[tr100.z0] == 3.01);
tr = tr100(783:917);
tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0);
tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]<85);
tr = tr100(783:917);
tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]>85);
tr = tr100(783:917);
tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]<0);


twpe = [25000];
%twpe = [23000 24000];
xlim = no02m.xi([1 end])+[50 -50]';
zlim = [-8 8];
cmapth = pic_colors('thermal');
cmapbr = pic_colors('blue_red');

pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);

varstrs = {'n(3)','Ez','log10(curvbabs)'}';
clims = {[0 0.5],[-1 1],[-2 1]};
cmaps = {cmapth,cmapbr,cmapbr};
trajargs = {};{'Marker','.'};
cbarlabels = {'Ion density from top','E_z','Magnetic curvature'};
trajcolordot = zeros(tr.ntr,3);
trajcolordot([tr.vy0]>0,:) = 1;
trajcolordot(and([tr.vy0]>0,[tr.x0]>85),:) = 0.5;

filename = [printpath 'no02m_ni_Ez_KB_t0_120_z0=0_manyx_startlim__'];
pic.movie(varstrs,'A',1,'colA',[1 1 1]*0.5,'cmap',cmaps,'clim',clims,...
  'filename',filename,'cbarlabels',cbarlabels,...
  'traj',tr,'trajargs',trajargs,'trajcolordot',trajcolordot);
%pic.movie(varstrs,'A',1,'colA',[1 1 1]*0.5,'cmap',cmaps,'clim',clims,'filename',filename,'traj',tr);
%pic.twpelim([17000 25000]).movie({'Ez'},'A',1,'clim',{[-1 1]},'cmap',{pic_colors('blue_red')},'filename',[printpath 'no02m_Ez']);

%% Make common plot with all bunches that can be compared
[unique_starts,counts] = tr100.unique_starts('sorted');
unique_starts = unique_starts(find(counts>25),:);
counts = counts(find(counts>25),:);
ngroups = size(unique_starts,1);
colors = [pic_colors('matlab'); 0 0 0];

nrows = 2; ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % (x,z)
  hca = h(isub); isub = isub + 1;
  holdon = 0;
  for igroup = 1:ngroups
    trs = tr100.find([tr100.t0]==unique_starts(igroup,1),...
                     [tr100.x0]==unique_starts(igroup,2),...
                     [tr100.z0]==unique_starts(igroup,3));
    trs.plot_all_xz(hca,'color',colors(igroup,:))
    leg{igroup} = sprintf('t_0 = %g, x_0 = %g, z_0 = %g',unique_starts(igroup,1:3));
    if not(holdon)
      hold(hca,'on')
    end    
  end
  hlines = hca.Children;
  legend(hlines(cumsum(counts)),leg,'location','eastoutside')
  hold(hca,'off')
end
if 1 % (x,y)
  hca = h(isub); isub = isub + 1;  
  holdon = 0;
  for igroup = 1:ngroups
    trs = tr100.find([tr100.t0]==unique_starts(igroup,1),...
                     [tr100.x0]==unique_starts(igroup,2),...
                     [tr100.z0]==unique_starts(igroup,3));    
    trs.plot_all_xy(hca,'color',colors(igroup,:))
    leg{igroup} = sprintf('t_0 = %g, x_0 = %g, z_0 = %g',unique_starts(igroup,1:3));
    if not(holdon)
      hold(hca,'on')
    end    
  end
  hlines = hca.Children;
  hline_ind = cumsum(counts);
  legend(hlines(hline_ind(end:-1:1)),leg,'location','eastoutside')  
  hold(hca,'off')
end

%% Plot forces with division of trajectories into colored groups
tr = tr100(783:917);
tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]<0);
colors = pic_colors('matlab');

tr = tr100(783:917);
tr1 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]>85);
tr = tr100(783:917);
tr2 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]<0);
tr = tr100(783:917);
tr3 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]<=75);
tr = tr100(783:917);
tr4 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]>75,[tr.x0]<=85);

tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0);

trajcolordot = nan(tr100.ntr,3);
trajcolordot([tr1.id],:) = repmat(colors(1,:),tr1.ntr,1);
trajcolordot([tr2.id],:) = repmat(colors(2,:),tr2.ntr,1);
trajcolordot([tr3.id],:) = repmat(colors(3,:),tr3.ntr,1);
trajcolordot([tr4.id],:) = repmat(colors(4,:),tr4.ntr,1);
trajcolordot(isnan(trajcolordot)) = [];
trajcolordot = reshape(trajcolordot,numel(trajcolordot)/3,3);

%% Plot grouped forces;
tr = tr100(783:917);
tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]<0);

tr = tr100(783:917);
tr1 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]>85);
tr = tr100(783:917);
tr2 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]<0);
tr = tr100(783:917);
tr3 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]<=73);
tr = tr100(783:917);
tr4 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]>73,[tr.x0]<=85);

tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0);

ntrs = [tr1.ntr,tr2.ntr,tr3.ntr,tr4.ntr];
ntrs

iGroup = nan(tr100.ntr,1);
iGroup([tr1.id],:) = 1;
iGroup([tr2.id],:) = 2;
iGroup([tr3.id],:) = 3;
iGroup([tr4.id],:) = 4;
iGroup(isnan(iGroup)) = [];

gTr = tr.group(iGroup,'mean',[]);

h = gTr.plot_forces;
h(1).XLim(2) = 125;

%% Reduced distributiin plot, with forces plotted on top of the distributions
fred = fred3_z0;
fred.t = fred.t(:)/4;
 
% x = fred.x;
% twci = repmat(120,size(x));
% z = repmat(unique(fred.z),size(x));

Ex = pic.interp(fred.x,fred.z,fred.t,'Ex'); 
Ey = pic.interp(fred.x,fred.z,fred.t,'Ey'); 
Ez = pic.interp(fred.x,fred.z,fred.t,'Ez'); 
Bx = pic.interp(fred.x,fred.z,fred.t,'Bx'); 
By = pic.interp(fred.x,fred.z,fred.t,'By'); 
Bz = pic.interp(fred.x,fred.z,fred.t,'Bz'); 

fred.Ex = repmat(Ex,1,numel(fred.v));% fred.Ex(fred.fvx==0) = NaN; 
fred.Ey = repmat(Ey,1,numel(fred.v));% fred.Ey(fred.fvy==0) = NaN; 
fred.Ez = repmat(Ez,1,numel(fred.v));% fred.Ez(fred.fvz==0) = NaN; 

fred.vyBz = Bz*fred.v; fred.vyBz(fred.fvy==0) = NaN;
fred.vzBy = By*fred.v; fred.vzBy(fred.fvz==0) = NaN;
fred.vzBx = Bx*fred.v; fred.vzBx(fred.fvz==0) = NaN;
fred.vxBz = Bz*fred.v; fred.vxBz(fred.fvx==0) = NaN;
fred.vxBy = By*fred.v; fred.vxBy(fred.fvx==0) = NaN;
fred.vyBx = Bx*fred.v; fred.vyBx(fred.fvy==0) = NaN;

fred.vBx = fred.vyBz-fred.vzBy;% fred.vBx(fred.fvx==0) = NaN;
fred.vBy = fred.vzBx-fred.vxBz;% fred.vBy(fred.fvy==0) = NaN;
fred.vBz = fred.vxBy-fred.vyBx;% fred.vBz(fred.fvz==0) = NaN;

fred.EvBx = fred.Ex + fred.vBx;
fred.EvBy = fred.Ey + fred.vBy;
fred.EvBz = fred.Ez + fred.vBz;

fred.EyvxBz = fred.Ey - fred.vxBz; fred.Ey_vxBz(fred.fvx==0) = NaN;
fred.ExvyBz = fred.Ex + fred.vyBz; fred.Ex_vyBz(fred.fvy==0) = NaN;

fred.NaN = fred.Ex*NaN;
fred_ = fred;
%%
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

tr = tr100(783:917);
tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0);

nrows = 4;
ncols = 3;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 1; colorE = [0 0.8 0.8];
doV = 1; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];
doTraj = 0;

cmap_dist = pic_colors('waterfall');
cmap_forces = pic_colors('blue_red');

%freds = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
%freds = {fred3_z0,fred3_z0,fred3_z0};
freds = {fred_,fred_,fred_,fred_,fred_,fred_,fred_,fred_,fred_,fred_,fred_,fred_};
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'fvx','fvy','fvz','Ex','Ey','Ez','vBx','vBy','vBz','EvBx','EvBy','EvBz'};
labstrs = {'fvx','fvy','fvz','Ex','Ey','Ez','vxBy','vyBz','vzBx','vxBz','vyBx','vxBy'};
labstrs = {'fvx','fvy','fvz','Ex','Ey','Ez','vxBz','vyBz','NaN','EyvxBz','ExvyBz','NaN'};
%cmaps = {cmap_dist,cmap_dist,cmap_dist,cmap_forces,cmap_forces,cmap_forces,cmap_forces,cmap_forces,cmap_forces,cmap_forces,cmap_forces,cmap_forces};

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.' labstr]);
    switch labstr
      case {'fvx','fvy','fvz','NaN'}
        pcolor(hca,fred.x,fred.v,log10(fredplot)')
        colormap(hca,pic_colors('candy4')) 
        hcb = colorbar('peer',hca);          
        %hca.CLim(2) = prctile(fred.fvx(:),99);
        hca.YLabel.String = sprintf('v_{%s}',labstr(end));
      case {'Ex','Ey','Ez','vBx','vBy','vBz','EvBx','EvBy','EvBz','vxBy','vyBz','vzBx','vxBz','vyBx','vxBy','EyvxBz','ExvyBz'}
        pcolor(hca,fred.x,fred.v,fredplot')
        colormap(hca,pic_colors('blue_red')) 
        hcb = colorbar('peer',hca);          
        %hca.CLim = prctile(fredplot(:),90)*[-1 1];
        hca.CLim = [-1 1];
        hca.Color = [0.8 0.8 0.8];
        hold(hca,'on')
        %fredcont = eval(['fred.fv' labstr(end)]);
        %fredcont = log10(fredcont);
        %contour(hca,fred.x,fred.v,fredcont',0.8:0.1:2,'k')
        %contour(hca,fred.x,fred.v,fredcont',-5:1:4,'k')
        hold(hca,'off')
        hca.YLabel.String = sprintf('v_{%s}',labstr(2));
    end
    hcb.YLabel.String = sprintf('%s',labstr);
    shading(hca,'flat')
    hca.XLabel.String = 'x (d_i)';
    
    
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 12;
    if doTraj
      hold(hca,'on')
      xx = [tr.x0];
      vv = eval(['[tr.v' labstr '0]']);
      scatter(hca,xx,vv,20,0*[1 1 1])
      hold(hca,'off')
    end    
  end
end
drawnow
compact_panels(h(1:end),0.01)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim'});
hlinks_ = linkprop(h(1:3),{'CLim'});
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
h(1).CLim = 0.99*[-4 2];
h(1).YLim = 0.99*4*[-1 1];





