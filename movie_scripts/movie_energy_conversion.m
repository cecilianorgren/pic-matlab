fileName  = [printpath 'vid'];
doVideo = 1;
doGif = 1;
doGifBackLoop = 0;
doDark = 0;
colors = pic_colors('matlab');
colors = [0 0 0; pic_colors('matlab'); 1 0 0; 0 1 0; 0 0 1; 1 1 0];

fontsize = 16;

% Energy partition
pic = pic; % no02m

h(1) = subplot(1,3,[1 2]);
h(1).Position(2) = 0.15;
%h(2) = subplot(1,2,2);

h(2) = subplot(1,3,3);
h(2).Position(2) = 0.15;


if doVideo
  vidObj = VideoWriter([fileName '.mp4'],'MPEG-4');
  vidObj.FrameRate = 10;
  open(vidObj);        
end
if doGif
  iframe = 0;
end

disp('Adjust figure size, then hit any key to continue.')
pause
nSpecies = numel(pic.mass);
for it = 1:pic.length
  pic_tmp = pic(1:it);
  UKtot = sum(pic_tmp.UK(1:nSpecies),2);
  UTtot = sum(pic_tmp.UT(1:nSpecies),2);
  UPtot = UKtot + UTtot;  
  UB = pic_tmp.UB;
  Utot = UB + UPtot;
  Unorm = Utot(1)/100;
  t = pic_tmp.twci;


  hca = h(1);
  pic_tmp(pic_tmp.nt).plot_map(hca,{'(3/2)*pi'},'A',1,'cmap',pic_colors('thermal'),'clim',{[0 1.3]},'cbarlabels',{'Ion thermal energy'})

  hca = h(2);
  plot(hca,[0 t],[UB(1); UB]/Unorm,[0 t],[UPtot(1);  UPtot]/Unorm,'linewidth',3)
  hca.XLim = [0 pic.twci(pic.nt)];
  hca.YLim = [0 100];
  hca.XLabel.String = 'Time (\omega_{ci}^{-1})';
  hca.YLabel.String = 'Energy (%)';
  %legend(hca,{'Magnetic energy','Plasma energy'},'box','off')
  hca.FontSize = fontsize;
  irf_legend(hca,{'Magnetic energy'},[0.05 0.72],'color',colors(1,:),'fontsize',fontsize,'fontweight','bold')
  irf_legend(hca,{'Plasma energy'},[0.05 0.3],'color',colors(2,:),'fontsize',fontsize,'fontweight','bold')
  hca.FontWeight = 'bold';

  drawnow
  if doVideo
    if doDark
      set(gcf,'color',darkBackgroundColor);          
    else
      set(gcf,'color','white');
    end
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = pic.nt;
      currentBackgroundColor = get(gcf,'color');
      if doDark
        set(gcf,'color',darkBackgroundColor);          
      else
        set(gcf,'color','white');
      end
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
  close(vidObj);   
end
if doGif
  imwrite(all_im,map,[fileName,'.gif'],'DelayTime',0,'LoopCount',inf);
end
if doGif && doGifBackLoop
  imwrite(cat(4,all_im,all_im(:,:,:,end:-1:1)),map,[fileName,'_loopback.gif'],'DelayTime',0,'LoopCount',inf);           
end


%% Compare two simulations
localuser = datastore('local','user');
% Load PICs
pic1 = PIC(['/Users/' localuser '/Data/PIC/varying_tite/tite_05/fields.h5']);
pic2 = PIC(['/Users/' localuser '/Data/PIC/varying_tite/tite_10/fields.h5']);
pics = {pic1,pic1};


fileName  = [printpath 'vid__'];
doVideo = 1;
doGif = 1;
doGifBackLoop = 0;
doDark = 0;
colors = pic_colors('matlab');
fontsize = 16;

% Energy partition

h(1) = subplot(1,1,1);
h(1).Position(2) = 0.15;
colors = pic_colors('matlab');
linestyles = {'-','--'};


if doVideo
  vidObj = VideoWriter([fileName '.mp4'],'MPEG-4');
  vidObj.FrameRate = 10;
  open(vidObj);        
end
if doGif
  iframe = 0;
end

disp('Adjust figure size, then hit any key to continue.')
pause
iSpecies = 1:numel(pic.mass);
iSpecies = [1:4];
times = pic2.twci;
Unorm = (pic1(1).UB + sum(pic1(1).UT(:)) + sum(pic1(1).UK(:)))/100;
for it = 1:numel(times)
  % Collect data
  sp = struct([]);
  for ipic = 1:numel(pics)
    if it > pics{ipic}.nt
      pic_tmp = pics{ipic}; 
    else
      pic_tmp = pics{ipic}(1:it);
    end
      UKtot = sum(pic_tmp.UK(iSpecies),2);
      UTtot = sum(pic_tmp.UT(iSpecies),2);
      UPtot = UKtot + UTtot;
      UB = pic_tmp.UB;
      Utot = UB + UPtot;
      %Unorm = Utot(1)/100;
      t = pic_tmp.twci;
      sp(ipic).t = pic_tmp.twci;
      sp(ipic).UKtot = UKtot;
      sp(ipic).UTtot = UTtot;
      sp(ipic).UPtot = UPtot;
      sp(ipic).UB = UB;
      %sp(ipic).Unorm = Utot(1)/100;        
  end

  hca = h(1);
  %pic_tmp(pic_tmp.nt).plot_map(hca,{'(3/2)*pi'},'A',1,'cmap',pic_colors('thermal'),'clim',{[0 1.3]},'cbarlabels',{'Ion thermal energy'})

  hca = h(1);
  if 1
    hp = plot(hca,sp(1).t,sp(1).UB/Unorm,...
                  sp(1).t,sp(1).UTtot/Unorm,...
                  sp(1).t,sp(1).UKtot/Unorm,...
                  sp(2).t,sp(2).UB/Unorm,...
                  sp(2).t,sp(2).UTtot/Unorm,...
                  sp(2).t,sp(2).UKtot/Unorm,...
                  'linewidth',3);
    hp(1).Color = colors(1,:).^.25; hp(1).LineStyle = linestyles{1};
    hp(2).Color = colors(2,:).^.25; hp(2).LineStyle = linestyles{1};
    hp(3).Color = colors(3,:).^.25; hp(3).LineStyle = linestyles{1};
    hp(4).Color = colors(1,:);      hp(4).LineStyle = linestyles{2};
    hp(5).Color = colors(2,:);      hp(5).LineStyle = linestyles{2};
    hp(6).Color = colors(3,:);      hp(6).LineStyle = linestyles{2};

    irf_legend(hca,{'Magnetic energy'},[0.05 0.60],'color',colors(1,:),'fontsize',fontsize,'fontweight','bold')
    irf_legend(hca,{'Plasma thermal energy'},[0.05 0.32],'color',colors(2,:),'fontsize',fontsize,'fontweight','bold')  
    irf_legend(hca,{'Plasma kinetic energy'},[0.05 0.07],'color',colors(3,:),'fontsize',fontsize,'fontweight','bold')  
    
  elseif 1
    hp = plot(hca,sp(1).t,sp(1).UB/Unorm,...
                  sp(1).t,sp(1).UPtot/Unorm,...
                  sp(2).t,sp(2).UB/Unorm,...
                  sp(2).t,sp(2).UPtot/Unorm,...
                  'linewidth',3);
    hp(1).Color = colors(1,:).^.25; hp(1).LineStyle = linestyles{1};
    hp(2).Color = colors(2,:).^.25; hp(2).LineStyle = linestyles{1};
    hp(3).Color = colors(1,:);      hp(3).LineStyle = linestyles{2};
    hp(4).Color = colors(2,:);      hp(4).LineStyle = linestyles{2};

    irf_legend(hca,{'Magnetic energy'},[0.05 0.60],'color',colors(1,:),'fontsize',fontsize,'fontweight','bold')
    irf_legend(hca,{'Plasma energy'},[0.05 0.22],'color',colors(2,:),'fontsize',fontsize,'fontweight','bold')  
    
  else
    hp = plot(hca,sp(1).t,sp(1).UB/Unorm,...
                  sp(1).t,sp(1).UTtot/Unorm,...
                  sp(2).t,sp(2).UB/Unorm,...
                  sp(2).t,sp(2).UTtot/Unorm,...
                  'linewidth',3);
    hp(1).Color = colors(1,:).^.25; hp(1).LineStyle = linestyles{1};
    hp(2).Color = colors(2,:).^.25; hp(2).LineStyle = linestyles{1};
    hp(3).Color = colors(1,:);      hp(3).LineStyle = linestyles{2};
    hp(4).Color = colors(2,:);      hp(4).LineStyle = linestyles{2};

    irf_legend(hca,{'Magnetic energy'},[0.05 0.60],'color',colors(1,:),'fontsize',fontsize,'fontweight','bold')
    irf_legend(hca,{'Thermal plasma energy'},[0.05 0.22],'color',colors(2,:),'fontsize',fontsize,'fontweight','bold')  
  end
  
  hca.XLim = [0 times(end)];
  hca.YLim = [0 110];
  hca.XLabel.String = 'Time (\omega_{ci}^{-1})';
  hca.YLabel.String = 'Energy (%)';  
  hca.FontSize = fontsize;
  hca.FontWeight = 'bold';
  hca.LineWidth = 1;

  drawnow
  if doVideo
    if doDark
      set(gcf,'color',darkBackgroundColor);          
    else
      set(gcf,'color','white');
    end
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = pic.nt;
      currentBackgroundColor = get(gcf,'color');
      if doDark
        set(gcf,'color',darkBackgroundColor);          
      else
        set(gcf,'color','white');
      end
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
  close(vidObj);   
end
if doGif
  imwrite(all_im,map,[fileName,'.gif'],'DelayTime',0,'LoopCount',inf);
end
if doGif && doGifBackLoop
  imwrite(cat(4,all_im,all_im(:,:,:,end:-1:1)),map,[fileName,'_loopback.gif'],'DelayTime',0,'LoopCount',inf);           
end

%% Compare two simulations
localuser = datastore('local','user');
% Load PICs
pic1 = PIC(['/Users/' localuser '/Data/PIC/varying_tite/tite_05/fields.h5']);
pic2 = PIC(['/Users/' localuser '/Data/PIC/varying_tite/tite_10/fields.h5']);
pics = {pic1,pic2};


fileName  = [printpath 'vid_UP_part'];
doVideo = 1;
doGif = 1;
doGifBackLoop = 0;
doDark = 0;
colors = pic_colors('matlab');
fontsize = 16;

% Energy partition

h(1) = subplot(1,1,1);
h(1).Position(2) = 0.15;
colors = pic_colors('matlab');
colors = [0 0 0; pic_colors('matlab'); 1 0 0; 0 1 0; 0 0 1; 1 1 0];
linestyles = {'-','--'};


if doVideo
  vidObj = VideoWriter([fileName '.mp4'],'MPEG-4');
  vidObj.FrameRate = 10;
  open(vidObj);        
end
if doGif
  iframe = 0;
end

disp('Adjust figure size, then hit any key to continue.')
pause
iSpecies = 1:numel(pic.mass);
iSpecies = [1:4];
times = pic2.twci;
Unorm = (pic1(1).UB + sum(pic1(1).UT(:)) + sum(pic1(1).UK(:)))/100;
for it = 1:numel(times)
  % Collect data
  sp = struct([]);
  for ipic = 1:numel(pics)
    if it > pics{ipic}.nt
      pic_tmp = pics{ipic}; 
    else
      pic_tmp = pics{ipic}(1:it);
    end
      UKall = pic_tmp.UK(:);
      UTall = pic_tmp.UT(:);
    
      UKtot = sum(pic_tmp.UK(iSpecies),2);
      UTtot = sum(pic_tmp.UT(iSpecies),2);
      UPtot = UKtot + UTtot;
      UB = pic_tmp.UB;
      Utot = UB + UPtot;
      %Unorm = Utot(1)/100;
      t = pic_tmp.twci;
      sp(ipic).t = pic_tmp.twci;
      sp(ipic).UKall = UKall;
      sp(ipic).UTall = UTall;
      sp(ipic).UKtot = UKtot;
      sp(ipic).UTtot = UTtot;
      sp(ipic).UPtot = UPtot;
      sp(ipic).UB = UB;
      %sp(ipic).Unorm = Utot(1)/100;        
  end

  hca = h(1);
  %pic_tmp(pic_tmp.nt).plot_map(hca,{'(3/2)*pi'},'A',1,'cmap',pic_colors('thermal'),'clim',{[0 1.3]},'cbarlabels',{'Ion thermal energy'})

  hca = h(1);
  if 1
    hp = plot(hca,sp(1).t,sp(1).UTall/Unorm,...
                  sp(1).t,sp(1).UKall/Unorm,...
                  sp(2).t,sp(2).UTall/Unorm,...
                  sp(2).t,sp(2).UKall/Unorm,...
                  'linewidth',3);
    hp(1).Color = colors(1,:).^.25; hp(1).LineStyle = linestyles{1};
    hp(2).Color = colors(2,:).^.25; hp(2).LineStyle = linestyles{1};
    hp(3).Color = colors(3,:).^.25; hp(3).LineStyle = linestyles{1};
    hp(4).Color = colors(4,:).^.25; hp(4).LineStyle = linestyles{1};
    hp(5).Color = colors(5,:).^.25; hp(5).LineStyle = linestyles{1};
    hp(6).Color = colors(6,:).^.25; hp(6).LineStyle = linestyles{1};
    hp(7).Color = colors(7,:).^.25; hp(7).LineStyle = linestyles{1};
    hp(8).Color = colors(8,:).^.25; hp(8).LineStyle = linestyles{1};
    %hp(9).Color = colors(9,:).^.25; hp(9).LineStyle = linestyles{1};
    hp(1+8).Color = colors(1,:).^.99; hp(1+8).LineStyle = linestyles{2};
    hp(2+8).Color = colors(2,:).^.99; hp(2+8).LineStyle = linestyles{2};
    hp(3+8).Color = colors(3,:).^.99; hp(3+8).LineStyle = linestyles{2};
    hp(4+8).Color = colors(4,:).^.99; hp(4+8).LineStyle = linestyles{2};
    hp(5+8).Color = colors(5,:).^.99; hp(5+8).LineStyle = linestyles{2};
    hp(6+8).Color = colors(6,:).^.99; hp(6+8).LineStyle = linestyles{2};
    hp(7+8).Color = colors(7,:).^.99; hp(7+8).LineStyle = linestyles{2};
    hp(8+8).Color = colors(8,:).^.99; hp(8+8).LineStyle = linestyles{2};
    %hp(9+9).Color = colors(9,:).^.99; hp(9).LineStyle = linestyles{2};

    %irf_legend(hca,{'Magnetic energy'},[0.05 0.60],'color',colors(1,:),'fontsize',fontsize,'fontweight','bold')
    %irf_legend(hca,{'Plasma thermal energy'},[0.05 0.32],'color',colors(2,:),'fontsize',fontsize,'fontweight','bold')  
    %irf_legend(hca,{'Plasma kinetic energy'},[0.05 0.07],'color',colors(3,:),'fontsize',fontsize,'fontweight','bold')  

    irf_legend(hp(1:8),{'U_{T1}','U_{T2}','U_{T3}','U_{T4}','U_{K1}','U_{K2}','U_{K3}','U_{K4}'}',...
      [1.01 0.98],'fontsize',fontsize,'fontweight','bold')

  elseif 1


    hp = plot(hca,sp(1).t,sp(1).UB/Unorm,...
                  sp(1).t,sp(1).UTtot/Unorm,...
                  sp(1).t,sp(1).UKtot/Unorm,...
                  sp(2).t,sp(2).UB/Unorm,...
                  sp(2).t,sp(2).UTtot/Unorm,...
                  sp(2).t,sp(2).UKtot/Unorm,...
                  'linewidth',3);
    hp(1).Color = colors(1,:).^.25; hp(1).LineStyle = linestyles{1};
    hp(2).Color = colors(2,:).^.25; hp(2).LineStyle = linestyles{1};
    hp(3).Color = colors(3,:).^.25; hp(3).LineStyle = linestyles{1};
    hp(4).Color = colors(1,:);      hp(4).LineStyle = linestyles{2};
    hp(5).Color = colors(2,:);      hp(5).LineStyle = linestyles{2};
    hp(6).Color = colors(3,:);      hp(6).LineStyle = linestyles{2};

    irf_legend(hca,{'Magnetic energy'},[0.05 0.60],'color',colors(1,:),'fontsize',fontsize,'fontweight','bold')
    irf_legend(hca,{'Plasma thermal energy'},[0.05 0.32],'color',colors(2,:),'fontsize',fontsize,'fontweight','bold')  
    irf_legend(hca,{'Plasma kinetic energy'},[0.05 0.07],'color',colors(3,:),'fontsize',fontsize,'fontweight','bold')  
    
  elseif 1
    hp = plot(hca,sp(1).t,sp(1).UB/Unorm,...
                  sp(1).t,sp(1).UPtot/Unorm,...
                  sp(2).t,sp(2).UB/Unorm,...
                  sp(2).t,sp(2).UPtot/Unorm,...
                  'linewidth',3);
    hp(1).Color = colors(1,:).^.25; hp(1).LineStyle = linestyles{1};
    hp(2).Color = colors(2,:).^.25; hp(2).LineStyle = linestyles{1};
    hp(3).Color = colors(1,:);      hp(3).LineStyle = linestyles{2};
    hp(4).Color = colors(2,:);      hp(4).LineStyle = linestyles{2};

    irf_legend(hca,{'Magnetic energy'},[0.05 0.60],'color',colors(1,:),'fontsize',fontsize,'fontweight','bold')
    irf_legend(hca,{'Plasma energy'},[0.05 0.22],'color',colors(2,:),'fontsize',fontsize,'fontweight','bold')  
    
  else
    hp = plot(hca,sp(1).t,sp(1).UB/Unorm,...
                  sp(1).t,sp(1).UTtot/Unorm,...
                  sp(2).t,sp(2).UB/Unorm,...
                  sp(2).t,sp(2).UTtot/Unorm,...
                  'linewidth',3);
    hp(1).Color = colors(1,:).^.25; hp(1).LineStyle = linestyles{1};
    hp(2).Color = colors(2,:).^.25; hp(2).LineStyle = linestyles{1};
    hp(3).Color = colors(1,:);      hp(3).LineStyle = linestyles{2};
    hp(4).Color = colors(2,:);      hp(4).LineStyle = linestyles{2};

    irf_legend(hca,{'Magnetic energy'},[0.05 0.60],'color',colors(1,:),'fontsize',fontsize,'fontweight','bold')
    irf_legend(hca,{'Thermal plasma energy'},[0.05 0.22],'color',colors(2,:),'fontsize',fontsize,'fontweight','bold')  
  end
  
  hca.XLim = [0 times(end)];
  hca.YLim = [0 35];
  hca.XLabel.String = 'Time (\omega_{ci}^{-1})';
  hca.YLabel.String = 'Energy (%)';  
  hca.FontSize = fontsize;
  hca.FontWeight = 'bold';
  hca.LineWidth = 1;

  drawnow
  if doVideo
    if doDark
      set(gcf,'color',darkBackgroundColor);          
    else
      set(gcf,'color','white');
    end
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = pic.nt;
      currentBackgroundColor = get(gcf,'color');
      if doDark
        set(gcf,'color',darkBackgroundColor);          
      else
        set(gcf,'color','white');
      end
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
  close(vidObj);   
end
if doGif
  imwrite(all_im,map,[fileName,'.gif'],'DelayTime',0,'LoopCount',inf);
end
if doGif && doGifBackLoop
  imwrite(cat(4,all_im,all_im(:,:,:,end:-1:1)),map,[fileName,'_loopback.gif'],'DelayTime',0,'LoopCount',inf);           
end

%% Compare evolution of reconnection
clear h;
nrows = 2; ncols = 1; ipanel = 0;
for irow = 1:nrows; for icol = 1:ncols; ipanel = ipanel + 1; h(ipanel) = subplot(nrows,ncols,ipanel); end; end 
isub = 1;

fonstsize = 14;
varstrs = {'Ey'};
clims = {0.1*[-1 1]};
cbarlabels = {'E_y'};
cmaps = {pic_colors('blue_red')};
zlim = 0.5*[-1 1];

hca = h(isub); isub = isub + 1;
hs1 = pic1.zlim(zlim).plot_timemap(hca,'xt',varstrs,'clim',clims,'cmap',cmaps,'cbarlabel',cbarlabels,'smooth',1);

hca = h(isub); isub = isub + 1;
hs2 = pic2.zlim(zlim).plot_timemap(hca,'xt',varstrs,'clim',clims,'cmap',cmaps,'cbarlabel',cbarlabels,'smooth',1);


hlinks = linkprop(h,{'XLim','YLim','CLim'});

c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on''; h(?).Layer = ''top''; h(?).GridAlpha = 0.1;',1:numel(h))
c_eval('h(?).Color = [0.9 0.9 0.9]; h(?).FontSize = fontsize; h(?).FontWeight = ''bold''; h(?).LineWidth = 1;',1:numel(h))
h(1).YLim = [0 120];


compact_panels(h,0.04)
h(2).Title.String = '';
%hca.XLabel.String = 'Time (\omega_{ci}^{-1})';


%% Energy partition for no02m
fileName  = [printpath 'vid_pic2'];
doVideo = 1;
doGif = 1;
doGifBackLoop = 0;
doDark = 0;
colors = pic_colors('matlab');
fontsize = 16;

% Energy partition
pic = pic2; % no02m

h(1) = subplot(1,1,1);


if doVideo
  vidObj = VideoWriter([fileName '.mp4'],'MPEG-4');
  vidObj.FrameRate = 10;
  open(vidObj);        
end
if doGif
  iframe = 0;
end

disp('Adjust figure size, then hit any key to continue.')
pause
iSpecies = 1:numel(pic.mass);
iSpeciesHot = [1 2];
iSpeciesCold = [3 4];
for it = 1:pic.length
  it
  pic_tmp = pic(1:it);
  UB = pic_tmp.UB;
  
  UKtot1 = sum(pic_tmp.UK(iSpeciesHot),2);
  UTtot1 = sum(pic_tmp.UT(iSpeciesHot),2);
  UPtot1 = UKtot1 + UTtot1;    
  Utot1 = UB + UPtot1;
  
  UKtot2 = sum(pic_tmp.UK(iSpeciesCold),2);
  UTtot2 = sum(pic_tmp.UT(iSpeciesCold),2);
  UPtot2 = UKtot2 + UTtot2;    
  Utot2 = UB + UPtot2;
  
  Unorm = (UB(1) + UPtot1(1) + UPtot2(1))/100;
  
  t = pic_tmp.twci;

  %hca = h(1);
  %pic_tmp(pic_tmp.nt).plot_map(hca,{'(3/2)*pi'},'A',1,'cmap',pic_colors('thermal'),'clim',{[0 1.3]},'cbarlabels',{'Ion thermal energy'})

  hca = h(1);
  if 1 % change in energy
    
    UB_ = [UB(1); UB]/Unorm;
    UP1_ = [UPtot1(1);  UPtot1]/Unorm;
    UP2_ = [UPtot2(1);  UPtot2]/Unorm;
    t_ = [0 t];
    %plot(hca,t_,[0; diff(UB_)./diff(t_)'],t_,[0; diff(UP1_)./diff(t_)'],t_,[0; diff(UP2_)./diff(t_)'],'linewidth',3)
    plot(hca,t_,UB_-UB_(1),t_,UP1_-UP1_(1),t_,UP2_-UP2_(1),'linewidth',3)
    hca.YLabel.String = '\Delta Energy (%)';
    %irf_legend(hca,{'Magnetic energy'},[0.8 0.20],'color',colors(1,:),'fontsize',fontsize,'fontweight','bold')
    %irf_legend(hca,{{'Plasma energy:','pre-existing plasma sheet'}},[0.5 0.85],'color',colors(2,:),'fontsize',fontsize,'fontweight','bold','horizontalalignment','left')
    %irf_legend(hca,{'Plasma energy: cold inflow'},[0.05 0.6],'color',colors(3,:),'fontsize',fontsize,'fontweight','bold')
    irf_legend(hca,{'Magnetic energy'},[0.8 0.20],'color',colors(1,:),'fontsize',fontsize,'fontweight','bold')
    irf_legend(hca,{{'Plasma energy:','pre-existing plasma sheet'}},[0.5 0.85],'color',colors(2,:),'fontsize',fontsize,'fontweight','bold','horizontalalignment','left')
    irf_legend(hca,{'Plasma energy: cold inflow'},[0.05 0.6],'color',colors(3,:),'fontsize',fontsize,'fontweight','bold')
    hca.YLim = [-25 25]; % no02m
    hca.YLim = [0 30];
  else
    plot(hca,[0 t],[UB(1); UB]/Unorm,[0 t],[UPtot1(1);  UPtot1]/Unorm,[0 t],[UPtot2(1);  UPtot2]/Unorm,'linewidth',3)
    hca.YLabel.String = 'Energy (%)';
    irf_legend(hca,{'Magnetic energy'},[0.05 0.70],'color',colors(1,:),'fontsize',fontsize,'fontweight','bold')
    irf_legend(hca,{'Plasma energy: pre-existing plasma sheet'},[0.05 0.35],'color',colors(2,:),'fontsize',fontsize,'fontweight','bold')
    irf_legend(hca,{'Plasma energy: cold inflow'},[0.05 0.07],'color',colors(3,:),'fontsize',fontsize,'fontweight','bold')
    hca.YLim = [0 100];
    
  end
  hca.XLim = [0 pic.twci(pic.nt)];
  hca.XLim = [0 120];
  hca.XLabel.String = 'Time (\omega_{ci}^{-1})';
  %legend(hca,{'Magnetic energy','Plasma energy'},'box','off')
  hca.FontSize = fontsize;
  hca.FontWeight = 'bold';

  drawnow
  if doVideo
    if doDark
      set(gcf,'color',darkBackgroundColor);          
    else
      set(gcf,'color','white');
    end
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = pic.nt;
      currentBackgroundColor = get(gcf,'color');
      if doDark
        set(gcf,'color',darkBackgroundColor);          
      else
        set(gcf,'color','white');
      end
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
  close(vidObj);   
end
if doGif
  imwrite(all_im,map,[fileName,'.gif'],'DelayTime',0,'LoopCount',inf);
end
if doGif && doGifBackLoop
  imwrite(cat(4,all_im,all_im(:,:,:,end:-1:1)),map,[fileName,'_loopback.gif'],'DelayTime',0,'LoopCount',inf);           
end






