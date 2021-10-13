%% Load PIC and PICTraj objects
no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
trp = PICTraj('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories_paul.h5');

%% First prepare all the data I want
% Loading t for three species requires loading 3*3 different sets of data.
% Doing this in advance saves a lot of time when I might need to adjust
% movie in some ways. E.g. test different trajectories.
% T - temperature
% A - vector potential

% missing: ti_A_twpe=17400
savePath = ['/Users/' localuser '/Data/PIC/no02m_derived/'];

twci_tr = 15000:100:25000;
nt = numel(twci_tr);
xlim = no02m.xi([1 end]);
zlim = no02m.zi([1 end]);
species = [1 3 5];
for it = 25%:nt
  twpe = twci_tr(it);
  pic = no02m.twpelim(twpe);
  twci = pic.twci;
  xi = pic.xi;
  zi = pic.zi;
  ti = pic.t(species);
  A = pic.A;
  save([savePath sprintf('ti_A_twpe=%.0f.mat',twpe)],'ti','A','xi','zi','twpe','twci')
  disp(['Saved ' savePath sprintf('ti_A_twpe=%.0f.mat',twpe)])
end

%% Setup
fileName = [printpath 'test2'];
twpe_all = 15000:100:25000;
twci_all = no02m.twpelim(twpe_all,'exact').twci;

colorBackground = [33 33 33]/255;
colorAxes = [0.8 0.8 0.8];
cmap = pic_colors('thermal');
x0 = (no02m.xi(1)+no02m.xi(end))/2;
xlim = x0 + [-50 50];
zlim = [-10 10];
doA = 1;
colorA = [0 0 0];
levelA = [0:0.5:15]; %levA = floor(min(A(:))/stepA)*stepA:stepA:ceil(max(A(:))/stepA)*stepA;
doVideo = 1;
doGif = 1;
clim = [0 0.6];
doSmooth = 1;
npSmooth = 2;
doGifBackLoop = 1;

% Setup for trajectories
doTrajectories = 1;
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

for it = 1:numel(twpe_all) % loop through time
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
      twci_tr = tr(itr).t;
      xx = tr(itr).x;
      zz = tr(itr).z;
      idup = find(diff(tr(itr).t)==0);
      twci_tr(idup) = [];
      xx(idup) = [];
      zz(idup) = [];
      if twpe == twci_tr(1)
        xnow = xx(1);
        znow = zz(1); 
        xtail = NaN;
        ztail = NaN;
      else
        % Current position
        if twpe == twci_tr(end)
          xnow = xx(end);
          znow = zz(end);
        else
          xnow = interp1(twci_tr,xx,twci);
          znow = interp1(twci_tr,zz,twci);
        end
        
        if doTrajTail
          if it == 1
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
  hca.YDir = 'normal';n
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
      

