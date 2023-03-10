fileName  = [printpath 'vid'];
doVideo = 1;
doGif = 1;
doGifBackLoop = 0;
doDark = 0;
colors = pic_colors('matlab');
fontsize = 16;

% Energy partition
pic = no02m;

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

for it = 1:pic.length
  pic_tmp = pic(1:it);
  UKtot = sum(pic_tmp.UK(1:6),2);
  UTtot = sum(pic_tmp.UT(1:6),2);
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
