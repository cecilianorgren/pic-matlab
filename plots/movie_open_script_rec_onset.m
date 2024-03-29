%% Load PIC and PICTraj objects
no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
trp = PICTraj('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories_paul.h5');

%% First prepare all the data I want
% Loading t for three species requires loading 3*3 different sets of data.
% Doing this in advance saves a lot of time when I might need to adjust
% movie in some ways. E.g. test different trajectories.
% T - temperature
% A - vector potential
localuser = datastore('local','user');
savePath = ['/Users/' localuser '/Data/PIC/no02m_ti_A_derived/'];

twpe_tr = no02m.twpelim([1100 14900]).twpe;
nt = numel(twpe_tr);
xlim = no02m.xi([1 end]);
zlim = no02m.zi([1 end]);
species = [1 3 5];
for it = 1:nt
  twpe = twpe_tr(it);
  pic = no02m.twpelim(twpe);
  twci = pic.twci;
  xi = pic.xi;
  zi = pic.zi;
  ti = pic.t(species);
  A = pic.A;
  save([savePath sprintf('ti_A_twpe=%.0f.mat',twpe)],'ti','A','xi','zi','twpe','twci')
  disp(['Saved ' savePath sprintf('ti_A_twpe=%.0f.mat',twpe)])
end

%% Interpolate data prior to 15000

twpe_old = no02m.twpelim([0 15000]).twpe;
twpe_all = [twpe_old(1):400:14000 14200:200:twpe_old(end)];
[twpe_new,idiff] = setdiff(twpe_all,twpe_old);
twci_new = twpe_new/200;
twpe_new
nt_new = numel(twpe_new);
%%
nt_old = numel(twpe_old);
ti_old = zeros(nx,nz,nt_old);
A_old = zeros(nx,nz,nt_old);

for it = 1:nt_old 
  twpe = twpe_old(it);
  % Load data
  data = load([savePath sprintf('ti_A_twpe=%05.0f.mat',twpe)]);
  A = data.A;
  t = data.ti;
  xi = data.xi;
  zi = data.zi;  
  twpe = data.twpe;
  twci = data.twci;
  nx = numel(xi);
  nz = numel(zi);
  
  % Save all data in matrix, to be used for interpolation later
  ti_old(:,:,it) = t;
  A_old(:,:,it) = A;
end

%% Interpolate to more frequent cadence
savePathInterp = ['/Users/' localuser '/Data/PIC/no02m_ti_A_interp/'];
tinds = {[11:15],[16:20],[21:25],[26:nt_new]};
tind = 6:10;
for tt = 1:numel(tinds)
tind = tinds{tt};
%%
%tind = 14;
%pause
[X,Y,Z] = meshgrid(xi,zi,twpe_old);
[Xq,Yq,Zq] = meshgrid(xi,zi,twpe_new(tind));

tic
V = permute(ti_old,[2 1 3]);
ti_data_new = interp3(X,Y,Z,V,Xq,Yq,Zq);
ti_data_new = permute(ti_data_new,[2 1 3]); 
toc
tic
V = permute(A_old,[2 1 3]);
A_data_new = interp3(X,Y,Z,V,Xq,Yq,Zq);
A_data_new = permute(A_data_new,[2 1 3]);
toc

for it = 1:numel(tind)
  ti = ti_data_new(:,:,it);
  A = A_data_new(:,:,it);
  twpe = twpe_new(tind(it));
  twci = twci_new(tind(it));  
  imagesc(ti')
  caxis([0 0.6])
  colormap(pic_colors('thermal'))
  drawnow
  filename = sprintf('ti_A_twpe=%05.0f.mat',twpe);
  save([savePathInterp filename],'ti','A','xi','zi','twpe','twci')
  disp(['Saved ' savePathInterp filename])
end
disp('Done.')
end
%% Construct artificial field between "pre-onset" and first time step
loadPath = ['/Users/' localuser '/Data/PIC/no02m_ti_A_derived/'];
savePathPreonset = ['/Users/' localuser '/Data/PIC/no02m_ti_A_preonset/'];
% Construct pre-onset field
twpe1 = 1000;
data = load([loadPath sprintf('ti_A_twpe=%05.0f.mat',twpe1)]);
A = data.A;
ti = data.ti;
xi = data.xi;
zi = data.zi;  
twpe = data.twpe;
twci = data.twci;

nx = numel(xi);
nz = numel(zi);

% Extract some non-perturbation part of the field
ix = 201:300;
iz = 1:nz;
T_onset = 10000;
dt_onset = 400;
t_onset = (-T_onset):dt_onset:twpe1;
%t_onset = 0:dt_onset:(twpe1-dt_onset);
nt_onset = numel(t_onset);

% ion temperature
ti_tmp = ti(ix,iz);
ti_preonset = repmat(ti_tmp,nx/numel(ix),1);
ti_preonset_and_zero(:,:,1) = ti_preonset;
ti_preonset_and_zero(:,:,2) = ti;

% vector potential
A_tmp = A(ix,iz);
A_preonset = repmat(A_tmp,nx/numel(ix),1);
A_preonset_and_zero(:,:,1) = A_preonset;
A_preonset_and_zero(:,:,2) = A;

[X,Y,Z] = meshgrid(xi,zi,[t_onset(1) twpe1]);
[Xq,Yq,Zq] = meshgrid(xi,zi,t_onset);
% simulate current sheet thinning
doCurrentSheetThinning = 1;
if doCurrentSheetThinning % modify Yq (=z) to be a bit inside, it will look like the sheet is thicker
  Yq_mod_range = linspace(0.8,1,nt_onset);
  for it = 1:nt_onset
    Yq(:,:,it) = Yq(:,:,it)*Yq_mod_range(it);
  end  
  
end

V = squeeze(ti_preonset_and_zero);
V = permute(V,[2 1 3]);
tic
ti_data = interp3(X,Y,Z,V,Xq,Yq,Zq);
toc
ti_data = permute(ti_data,[2 1 3]);
if doCurrentSheetThinning 
  % adjust temperature a bit so that it looks like the thinning gives rise
  % to higher temperature
  ti_mod_range = linspace(0.80,1,nt_onset);
  for it = 1:nt_onset
    ti_data(:,:,it) = ti_data(:,:,it)*ti_mod_range(it);
  end    
end

V = squeeze(A_preonset_and_zero);
V = permute(V,[2 1 3]);
tic
A_data = interp3(X,Y,Z,V,Xq,Yq,Zq);
toc
A_data = permute(A_data,[2 1 3]);
disp('Done interpolating data.')

if 1 % plot interpolated data
  %%
  for it = 1:numel(t_onset)
    imagesc(xi,zi,ti_data(:,:,it)')

    hold('on')  
    iAx = 1:5:nx;
    iAz = 1:5:nz;
    contour(xi(iAx),zi(iAz),A_data(iAx,iAz,it)',0:0.5:25,'k')
    hold('off')
    colorbar
    caxis([0 0.6])
    colormap(pic_colors('thermal'))
    title(sprintf('it=%g',it))
    pause(1)
  end
end

% Save artificially constructed data

for it = 1:nt_onset
  twpe = t_onset(it);
  twci = twpe/200;
  A = squeeze(A_data(:,:,it));
  ti = squeeze(ti_data(:,:,it));
  filename = sprintf('ti_A_twpe=%06.0f.mat',twpe);
  save([savePathPreonset filename],'ti','A','xi','zi','twpe','twci')
  disp(['Saved ' savePathPreonset filename])
end

%% Setup and run
fileName = [printpath 'test5'];
twpe_all_pic = [1000:200:15000 15100:100:25000];
twpe_pre = [-10000:400:000 200:200:800];

twpe_all_pic = [1000:400:14999,15000:100:25000,...
                1000,1600,2000,2600,3000:1000:14000];
twpe_all_pic = unique(sort(twpe_all_pic));
twpe_pre = [-10000:400:800];

twpe_pic = no02m.twpe;
twpe_int = setdiff(twpe_all_pic,twpe_pic);
twpe_int(twpe_int==6000) = []; % manual fix, don't know why it bugs?
%twpe_int(twpe_int==8000) = []; % manual fix, don't know why it bugs?
twpe_all = [twpe_pre,twpe_pic,twpe_int];
twci_all = twpe_all/200;
folder = [zeros(1,numel(twpe_pre))+1 zeros(1,numel(twpe_pic))+2 zeros(1,numel(twpe_int))+3];

[twpe_all,sort_ind] = sort(twpe_all);
folder = folder(sort_ind); 
%
%plot(twpe_all,folder)
loadPaths{1} = ['/Users/' localuser '/Data/PIC/no02m_ti_A_preonset/'];
loadPaths{2} = ['/Users/' localuser '/Data/PIC/no02m_ti_A_derived/'];
loadPaths{3} = ['/Users/' localuser '/Data/PIC/no02m_ti_A_interp/'];


%twpe_loop = 15000:100:25000;
%twpe_loop_and_interp = no02m.twpelim([0 14999]).twpe;
doFill = 1;
colorBackground = [33 33 33]/255;
colorAxes = [0.8 0.8 0.8];
cmap = pic_colors('thermal');
x0 = (no02m.xi(1)+no02m.xi(end))/2;
xlim = x0 + [-50 50];
zlim = [-10 10];
doA = 1;
colorA = [0.3 0.3 0.3];
levelA = [-10:0.5:15]; %levA = floor(min(A(:))/stepA)*stepA:stepA:ceil(max(A(:))/stepA)*stepA;
doMarkedA = 1;
colorMarkedA = [0.9 0.9 0.9];
linewidthMarkedA = 2;
levelMarkedA = 7; 
doVideo = 1;
doGif = 1;
clim = [0 0.6];
doSmooth = 1;
npSmooth = 2;
doGifBackLoop = 1;

% Setup for trajectories
doTrajectories = 1;
tr = trp;
%tr = trall; % see below, including fake preponset trajectories
trajcolordot = zeros(numel(tr),3)+0.4; 
trajcolordot([trp.charge]==-1,:) = 0.9;
%trajcolordot(37:end,:) = 0.9;
%tr = trp(1:2:end);
%trajcolordot = trajcolordot(1:2:end,:);
%colorTrajDot = [0 0 0];
doColorTrajDot = 1;
doColorTrajLine = 0;
doTrajTail = 0;
ntTail = 4;

%try
% Setup figure
fig = figure(777);
h = subplot(1,1,1);
hca = h(1);

hca.Position(2) = hca.Position(2)+0.05;
hca.Position(4) = hca.Position(4)-0.05;

% Figure colors
fig.Color = colorBackground;
ax = findobj(fig,'type','axes');
for iax = 1:numel(ax)
  ax(iax).XAxis.Color = colorAxes;
  ax(iax).YAxis.Color = colorAxes;
end

% Allow for adjustments
disp('Adjust figure size, then hit any key to continue.')
pause

% Setup for gif an video
if doVideo
  vidObj = VideoWriter([fileName '.mp4'],'MPEG-4');
  vidObj.FrameRate = 10;
  open(vidObj);        
end
if doGif
  iframe = 0;
end

% Loop through times
t_diff = [round(diff(twpe_all)/10)*10 100];
doPad = 0;
toPad = [400];
nPad = 1;
t_ind = 29:numel(twpe_all);
nt = numel(t_ind);
for it = t_ind%numel(twpe_all)
  disp(sprintf('it = %g/%g',it,t_ind(end)))
  twpe = twpe_all(it);
  
  % Load data
  loadPath = loadPaths{folder(it)};
  if folder(it) == 1
    data = load([loadPath sprintf('ti_A_twpe=%06.0f.mat',twpe)]);
  else
    data = load([loadPath sprintf('ti_A_twpe=%05.0f.mat',twpe)]);
  end
  
  A = data.A;
  t = data.ti;
  xi = data.xi;
  zi = data.zi;  
  twpe = data.twpe;
  twci = data.twci;
  
  nx = numel(xi);
  nz = numel(zi);
  
  var = t;
  ivar = 1;
        
  % Smoothing of data
  if doSmooth
    var = smooth2(var,npSmooth);
  end
  
  % Plot map
  imagesc(hca,xi,zi,var');
  
  if doA % In-plane magnetic field lines
    climtmp = hca.CLim;        
    iAx = 1:5:nx;
    iAz = 1:5:nz;
    hold(hca,'on')
    contour(hca,xi(iAx),zi(iAz),A(iAx,iAz)',levelA,'color',colorA)
    hold(hca,'off')
    hca.CLim = climtmp; 
  end
  if doMarkedA % In-plane magnetic field lines, a special marked line
    climtmp = hca.CLim;        
    iAx = 1:5:nx;
    iAz = 1:5:nz;
    hold(hca,'on')
    contour(hca,xi(iAx),zi(iAz),A(iAx,iAz)',levelMarkedA*[1 1],'color',colorMarkedA,'linewidth',linewidthMarkedA)
    hold(hca,'off')
    hca.CLim = climtmp; 
  end    
  if doTrajectories
    hold(hca,'on')
    for itr = 1:numel(tr)%tr.ntr
      twci_tr = tr(itr).t;
      xx = tr(itr).x;
      zz = tr(itr).z;
      idup = find(diff(tr(itr).t)==0);
      twci_tr(idup) = [];
      xx(idup) = [];
      zz(idup) = [];
      if twci == twci_tr(1)
        xnow = xx(1);
        znow = zz(1); 
        xtail = NaN;
        ztail = NaN;
      else
        % Current position
        if twci == twci_tr(end)
          xnow = xx(end);
          znow = zz(end);
        else
          xnow = interp1(twci_tr,xx,twci);
          znow = interp1(twci_tr,zz,twci);
        end
        
        if doTrajTail
          if it == t_ind(1)
            xtail = NaN;
            ztail = NaN;
          else
            ind_ = it + [-ntTail:0];
            ind_ = ind_(ind_ > 0);
            it1_tail = find(twci_tr>twci_all(ind_(1)),1,'first');
            it2_tail = find(twci_tr<twci_all(ind_(end)),1,'last');
            it_tail = it1_tail:it2_tail;
            xtail = xx(it_tail);                    
            ztail = zz(it_tail);
          end
        end
      end

      colline = trajcolordot(itr,:);
      coldot = trajcolordot(itr,:);
      %plot(hca,xnow,znow,'color',coldot,'markerSize',20,'marker','.','linestyle','none')
      sc = scatter(hca,xnow,znow,30,coldot,'Marker','o','MarkerFaceColor','flat','MarkerEdgeColor',[0 0 0]);
      %hp = plot(hca,tr(itr).x0,tr(itr).z0,'ko');
      if doTrajTail
       plot(hca,xtail,ztail,'color',coldot)  
      end
    end
    hold(hca,'off')
  end  
        
  hca.XLim = xlim;
  hca.YLim = zlim;
  hca.CLim = clim;  
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.FontSize = 14;  
  hca.XAxis.Color = colorAxes;
  hca.YAxis.Color = colorAxes;
  hca.Title.String = sprintf('twpe = %.0f, twci = %.1f',twpe,twci);
  hca.Title.Color = colorAxes;
  
  hcb = colorbar('peer',hca);
  hcb.Color = colorAxes;
  hcb.YLabel.String = 'Ion Temperature';  
  colormap(hca,cmap)
      
  if doFill
    hca.Position = [0 0 1 1];
  end
  if doVideo
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    if doPad
      indPad = find(2==toPad);
      if indPad
        for iPad = 1:nPad(indPad)
          writeVideo(vidObj,currFrame);
        end
      end
    end
  end
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = nt;
      currentBackgroundColor = get(gcf,'color');     
      set(gcf,'color',colorBackground);
      drawnow      
      tmp_frame = getframe(gcf);
      %cell_movies{imovie}(itime) = tmp_frame;
      if iframe == 1 % initialize animated gif matrix
        [im_tmp,map] = rgb2ind(tmp_frame.cdata,256,'nodither');
        %map(end+1,:) = get(gcf,'color');
        im_tmp(1,1,1,nframes) = 0;                                                
        all_im = im_tmp;             
      else
        all_im(:,:,1,iframe) = rgb2ind(tmp_frame.cdata,map,'nodither');
      end       
    end    
  end
end

% catch % If something goes wrong, write what we have
%   % Write gif and video
%   if doVideo
%     close(vidObj)
%   end
%   if doGif
%     imwrite(all_im,map,[fileName,'.gif'],'DelayTime',0,'LoopCount',inf)
%   end
%   if doGif && doGifBackLoop
%     imwrite(cat(4,all_im,all_im(:,:,:,end:-1:1)),map,[fileName,'_loopback.gif'],'DelayTime',0,'LoopCount',inf)              
%   end
% end
% Write gif and video
if doVideo
  close(vidObj)
end
if doGif
  imwrite(all_im,map,[fileName,'.gif'],'DelayTime',0,'LoopCount',inf)
end
if doGif && doGifBackLoop
  imwrite(cat(4,all_im,all_im(:,:,:,end:-1:1)),map,[fileName,'_loopback.gif'],'DelayTime',0,'LoopCount',inf)              
end
      



%% Move only existing pic and pre-onset
twpe_loop = 15000:100:25000;
twpe_loop_and_interp = no02m.twpelim([0 14999]).twpe;

colorBackground = [33 33 33]/255;
colorAxes = [0.8 0.8 0.8];
cmap = pic_colors('thermal');
x0 = (no02m.xi(1)+no02m.xi(end))/2;
xlim = x0 + [-50 50];
zlim = [-10 10];
doA = 1;
colorA = [0 0 0];
levelA = [-10:0.5:15]; %levA = floor(min(A(:))/stepA)*stepA:stepA:ceil(max(A(:))/stepA)*stepA;
doVideo = 1;
doGif = 1;
clim = [0 0.6];
doSmooth = 1;
npSmooth = 2;
doGifBackLoop = 1;

% Setup for trajectories
doTrajectories = 0;
tr = trp;
trajcolordot = zeros(numel(tr),3)+0.4; 
trajcolordot(37:end,:) = 0.9;
tr = trp(1:2:end);
trajcolordot = trajcolordot(1:2:end,:);
%colorTrajDot = [0 0 0];
doColorTrajDot = 1;
doColorTrajLine = 0;
doTrajTail = 1;
ntTail = 4;

try
% Setup figure
fig = figure(777);
h = subplot(1,1,1);
hca = h(1);

hca.Position(2) = hca.Position(2)+0.05;
hca.Position(4) = hca.Position(4)-0.05;

% Figure colors
fig.Color = colorBackground;
ax = findobj(fig,'type','axes');
for iax = 1:numel(ax)
  ax(iax).XAxis.Color = colorAxes;
  ax(iax).YAxis.Color = colorAxes;
end

% Allow for adjustments
disp('Adjust figure size, then hit any key to continue.')
pause

% Setup for gif an video
if doVideo
  vidObj = VideoWriter([fileName '.mp4'],'MPEG-4');
  vidObj.FrameRate = 10;
  open(vidObj);        
end
if doGif
  iframe = 0;
end

% Loop through onset
for it = 1:(nt_onset-1) % no need to duplicate first/last time step
  
  twpe = t_onset(it);
  A = A_data(:,:,it);
  t = ti_data(:,:,it);
  twci = twpe/200;
  var = t;
  ivar = 1;
        
  % Smoothing of data
  if doSmooth
    var = smooth2(var,npSmooth);
  end
  
  % Plot map
  imagesc(hca,xi,zi,var');
  
  if doA % In-plane magnetic field lines
    climtmp = hca.CLim;        
    iAx = 1:5:nx;
    iAz = 1:5:nz;
    hold(hca,'on')
    contour(hca,xi(iAx),zi(iAz),A(iAx,iAz)',levelA,'color',colorA)
    hold(hca,'off')
    hca.CLim = climtmp; 
  end  
  
  hca.XLim = xlim;
  hca.YLim = zlim;
  hca.CLim = clim;  
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.FontSize = 14;  
  hca.XAxis.Color = colorAxes;
  hca.YAxis.Color = colorAxes;
  hca.Title.String = sprintf('twpe = %.0f, twci = %.1f',twpe,twci);
  hca.Title.Color = colorAxes;
  
  hcb = colorbar('peer',hca);
  hcb.Color = colorAxes;
  hcb.YLabel.String = 'Ion Temperature';  
  colormap(hca,cmap)
  drawnow
  
  if doVideo
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = nt;
      currentBackgroundColor = get(gcf,'color');     
      set(gcf,'color',colorBackground);
      drawnow      
      tmp_frame = getframe(gcf);
      %cell_movies{imovie}(itime) = tmp_frame;
      if iframe == 1 % initialize animated gif matrix
        [im_tmp,map] = rgb2ind(tmp_frame.cdata,256,'nodither');
        %map(end+1,:) = get(gcf,'color');
        im_tmp(1,1,1,nframes) = 0;                                                
        all_im = im_tmp;             
      else
        all_im(:,:,1,iframe) = rgb2ind(tmp_frame.cdata,map,'nodither');
      end       
    end    
  end
end

% Loop through time pic, and interpolate where data is missing


% Loop through time pic
for it = 1:numel(twpe_all)
  twpe = twpe_all(it);
  % Load data
  data = load([savePath sprintf('ti_A_twpe=%.0f.mat',twpe)]);
  A = data.A;
  t = data.ti;
  xi = data.xi;
  zi = data.zi;  
  twpe = data.twpe;
  twci = data.twci;
  
  nx = numel(xi);
  nz = numel(zi);
  
  var = t;
  ivar = 1;
        
  % Smoothing of data
  if doSmooth
    var = smooth2(var,npSmooth);
  end
  
  % Plot map
  imagesc(hca,xi,zi,var');
  
  if doA % In-plane magnetic field lines
    climtmp = hca.CLim;        
    iAx = 1:5:nx;
    iAz = 1:5:nz;
    hold(hca,'on')
    contour(hca,xi(iAx),zi(iAz),A(iAx,iAz)',levelA,'color',colorA)
    hold(hca,'off')
    hca.CLim = climtmp; 
  end  
  if doTrajectories
    hold(hca,'on')
    for itr = 1:tr.ntr
      twpe_tr = tr(itr).t;
      xx = tr(itr).x;
      zz = tr(itr).z;
      idup = find(diff(tr(itr).t)==0);
      twpe_tr(idup) = [];
      xx(idup) = [];
      zz(idup) = [];
      if twpe == twpe_tr(1)
        xnow = xx(1);
        znow = zz(1); 
        xtail = NaN;
        ztail = NaN;
      else
        % Current position
        if twpe == twpe_tr(end)
          xnow = xx(end);
          znow = zz(end);
        else
          xnow = interp1(twpe_tr,xx,twci);
          znow = interp1(twpe_tr,zz,twci);
        end
        
        if doTrajTail
          if it == 1
            xtail = NaN;
            ztail = NaN;
          else
            ind_ = it + [-ntTail:0];
            ind_ = ind_(ind_ > 0);
            it1_tail = find(twpe_tr>twci_all(ind_(1)),1,'first');
            it2_tail = find(twpe_tr<twci_all(ind_(end)),1,'last');
            it_tail = it1_tail:it2_tail;
            xtail = xx(it_tail);                    
            ztail = zz(it_tail);
          end
        end
      end

      colline = trajcolordot(itr,:);
      coldot = trajcolordot(itr,:);
      %plot(hca,xnow,znow,'color',coldot,'markerSize',20,'marker','.','linestyle','none')
      sc = scatter(hca,xnow,znow,30,coldot,'Marker','o','MarkerFaceColor','flat','MarkerEdgeColor',[0 0 0]);
      %hp = plot(hca,tr(itr).x0,tr(itr).z0,'ko');
      if doTrajTail
       plot(hca,xtail,ztail,'color',coldot)  
      end
    end
    hold(hca,'off')
  end  
        
  hca.XLim = xlim;
  hca.YLim = zlim;
  hca.CLim = clim;  
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.FontSize = 14;  
  hca.XAxis.Color = colorAxes;
  hca.YAxis.Color = colorAxes;
  hca.Title.String = sprintf('twpe = %.0f, twci = %.1f',twpe,twci);
  hca.Title.Color = colorAxes;
  
  hcb = colorbar('peer',hca);
  hcb.Color = colorAxes;
  hcb.YLabel.String = 'Ion Temperature';  
  colormap(hca,cmap)
      
  if doVideo
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = nt;
      currentBackgroundColor = get(gcf,'color');     
      set(gcf,'color',colorBackground);
      drawnow      
      tmp_frame = getframe(gcf);
      %cell_movies{imovie}(itime) = tmp_frame;
      if iframe == 1 % initialize animated gif matrix
        [im_tmp,map] = rgb2ind(tmp_frame.cdata,256,'nodither');
        %map(end+1,:) = get(gcf,'color');
        im_tmp(1,1,1,nframes) = 0;                                                
        all_im = im_tmp;             
      else
        all_im(:,:,1,iframe) = rgb2ind(tmp_frame.cdata,map,'nodither');
      end       
    end    
  end
end

catch % If something goes wrong, write what we have
  % Write gif and video
  if doVideo
    close(vidObj)
  end
  if doGif
    imwrite(all_im,map,[fileName,'.gif'],'DelayTime',0,'LoopCount',inf)
  end
  if doGif && doGifBackLoop
    imwrite(cat(4,all_im,all_im(:,:,:,end:-1:1)),map,[fileName,'_loopback.gif'],'DelayTime',0,'LoopCount',inf)              
  end
end
% Write gif and video
if doVideo
  close(vidObj)
end
if doGif
  imwrite(all_im,map,[fileName,'.gif'],'DelayTime',0,'LoopCount',inf)
end
if doGif && doGifBackLoop
  imwrite(cat(4,all_im,all_im(:,:,:,end:-1:1)),map,[fileName,'_loopback.gif'],'DelayTime',0,'LoopCount',inf)              
end
      

%% Make fake trajectories prior to twci = 5
localuser = datastore('local','user');
loadPath = ['/Users/' localuser '/Data/PIC/no02m_ti_A_preonset/'];
% keep them on the same magnetic field line
twpe_pre = [-10000:400:800];
twci_pre = twpe_pre/200;
nt = numel(twpe_pre);
x0 = 90:5:115; % particles moves straight down
% possibility to add some dx or so later
vx0 = -0.025+1*0.05*rand(size(x0));
vx0 = [0.0044, 0.0234, -0.0093, -0.0030, -0.0156, 0.0044];
vy0 = 0.1;

z_part = zeros(nt,numel(x0));
x_part = repmat(x0,nt,1);
z_phase = rand(1,numel(x0));
t_part = 0.9+0.4*rand(1,numel(x0));
vx_part = repmat(vx0,nt,1);

x_part = x_part + cumsum(vx_part,1);


fci = 1;
fce = fci*200;
rci = 0.1;
rce = rci/100;
dtwci = unique(diff(twci_pre));

for it = 1:nt
  twpe = twpe_pre(it);
  
  % Load data  
  data = load([loadPath sprintf('ti_A_twpe=%06.0f.mat',twpe)]);
  A = data.A;
  t = data.ti;
  xi = data.xi;
  zi = data.zi;  
  twpe = data.twpe;
  twci = data.twci;  
  nx = numel(xi);
  nz = numel(zi);
  
  % Find the correct field line and assign the z-position
  S = contourcs(xi,zi,A',Alev*[1 1]);
  S = S(1);
  
  for ix = 1:numel(x0)
    
    ii = find(x>x0(ix),1,'first');
    z_part(it,ix) = S.Y(ii)+rci*t_part(ix)*sin(twpe/200+z_phase(ix));
  end
  plot(x_part,z_part,'.')
  drawnow
end

vx_part = [zeros(size(x0)); diff(x_part,1)./dtwci];
vz_part = [zeros(size(z0)); diff(z_part,1)./dtwci];

%% Merge fake and real trajectories
trp = PICTraj('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories_paul_t0=5_test.h5');

ntr = numel(trp);
for itr = 1:ntr
  trall(itr).t = [twci_pre'; trp(itr).t];
  trall(itr).twpe = [twpe_pre'; trp(itr).t*200];
  trall(itr).x = [x_part(:,itr); trp(itr).x];
  trall(itr).z = [z_part(:,itr); trp(itr).z];  
end







