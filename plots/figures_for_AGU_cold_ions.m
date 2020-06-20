df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
df04n = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/data_h5/fields.h5');
df08 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data_h5/fields.h5');
ds04 = PICDist('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/dists.h5');
tr04 = PICTraj('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5');
savedir  = ['/Users/' localuser '/GoogleDrive/DF_Cold_ions_figures_for_AGU/'];

%% Reconnection rate, two panels, absolute and massloading compensated
h = setup_subplots(1,2);
isub = 1;
if 1 % R(twci)
  hca = h(isub); isub = isub + 1;
  t04 = df04.twci; R04 = df04.RA; 
  t08 = df08.twci; R08 = df08.RA; 
  
  plot(hca,t04(5:end),R04(5:end),t08(4:end),R08(4:end))
  hca.YLabel.String = 'R   (B_0v_{A0})';
  hca.XLabel.String = 't (\omega_{ci}^{-1})';
  legend(hca,{'n_c = 0.4n_0','n_c = 0.8n_0'},'location','northwest','box','off')
end
if 1 % R(twci), scaled to density
  hca = h(isub); isub = isub + 1;
  fun_massloading = @(nh,nc) sqrt(1+nc/nh);
  t04 = df04.twci; R04 = df04.RA*fun_massloading(0.2,0.4); 
  t08 = df08.twci; R08 = df08.RA*fun_massloading(0.2,0.8);
  
  plot(hca,t04(5:end),R04(5:end),t08(4:end),R08(4:end))
  hca.YLabel.String = 'R (1+n_c/n_h)^{1/2}   (B_0v_A{0})';
  hca.XLabel.String = 't (\omega_{ci}^{-1})';
  legend(hca,{'n_c = 0.4n_0','n_c = 0.8n_0'},'location','northwest','box','off')
end
h(1).Title.String = 'Reconnection rate';
h(2).Title.String = {'Reconnection rate','rescaled for mass loading'};
c_eval('h(?).Position(2) = 0.20; h(?).Position(4) = 0.65;',1:2)
c_eval('h(?).YLim(2) = 0.19;',1:2)
c_eval('h(?).FontSize = 14;',1:2)
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:2)

%% Reconnection rate, one panel, absolute and massloading compensated
h = setup_subplots(1,1);
isub = 1;
if 1 % R(twci)
  hca = h(isub); isub = isub + 1;
  t04 = df04.twci; R04 = df04.RA; 
  t08 = df08.twci; R08 = df08.RA; 
  fun_massloading = @(nh,nc) sqrt(1+nc/nh);
  t04 = df04.twci; R04_ml = df04.RA*fun_massloading(0.2,0.4); 
  t08 = df08.twci; R08_ml = df08.RA*fun_massloading(0.2,0.8);
  
  colors = pic_colors('matlab');
  hlines = plot(hca,t04(5:end),R04(5:end),t08(4:end),R08(4:end),...
                    t04(5:end),R04_ml(5:end),t08(4:end),R08_ml(4:end));
  hlines(1).Color = colors(1,:); hlines(1).LineStyle = '-';
  hlines(3).Color = colors(1,:); hlines(3).LineStyle = '--';
  hlines(2).Color = colors(2,:); hlines(2).LineStyle = '-';
  hlines(4).Color = colors(2,:); hlines(4).LineStyle = '--';
  hca.YLabel.String = 'R   (B_0v_{A0})';
  hca.XLabel.String = 't (\omega_{ci}^{-1})';
  legend(hca,{'n_c = 0.4n_0','n_c = 0.8n_0','n_c = 0.4n_0: R(1+n_c/n_h)','n_c = 0.8n_0: R(1+n_c/n_h)'},'location','eastoutside','box','off')
end
h(1).Title.String = 'Reconnection rate';

c_eval('h(?).Position(2) = 0.20; h(?).Position(4) = 0.65;',1)
c_eval('h(?).YLim(2) = 0.19;',1)
c_eval('h(?).FontSize = 14;',1)
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1)

%% Partitioning between drift and thermal energy for cold ions in df04 run normalized to UB(1)

nrows = 1;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % Partitioning between drift and thermal energy for cold ions in df04 run normalized to UB(1)
  hca = h(isub); isub = isub + 1;
  colors = pic_colors('matlab');
  hold(hca,'on');
  hca.ColorOrder = colors(1:2,:);
  hca.LineStyleOrder = {'-','--'};
  
  plot(hca,df04.twci,2*df04.UK(3)/df04(1).UB,...
           df04.twci,2*df04.UT(3)/df04(1).UB,'linewidth',1.5)
  %hold(hca,'on')
  plot(hca,df04.twci,df04.UK(35)/df04(1).UB,'--',...
           df04.twci,df04.UT(35)/df04(1).UB,'--','linewidth',1.5)
  hca.XLabel.String = 't\omega_{ci}';  
  hca.YLabel.String = 'U_i^c/U_B(t=0)';
  hleg = legend(hca,{'2U_K^{c1} - drift separate','2U_T^{c1} - thermal separate','U_K^{c1+c2} - drift combined','U_T^{c1+c2} - thermal combined'},'Box','on','location','northwest');
  hleg.Title.String = 'n_c = 0.4 n_0';
  hold(hca,'off')
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
end

c_eval('hca.FontSize = 13;',1:npanels)
c_eval('h(?).Position(2) = 0.20; h(?).Position(4) = 0.65;',1)

%% Movie of ni evolution df04, gif
% [all_im, map] = MAKE_GIF(obj,fields,nrows,ncols)      
% make gif
% imwrite(im,map,'delme.gif','DelayTime',0.0,'LoopCount',0)  

% Default options, values
doAdjustCLim = 1;
clim = [0 2];
cmap = pic_colors('candy');
doA = 1;
levA = -25:1:0;

pic = df04;


% setup figure
fig = figure;
fig.Position = [100 100 920 360];
set(fig,'color','white');
nrows = 1;
ncols = 1;
h = setup_subplots(nrows,ncols); % external function, must include in SMILEI.m
isub = 1;
disp('Adjust figure size, then hit any key to continue.')
pause

times = 2:pic.nt;
ntimes = numel(times);
iframe = 0;
for itime = times
  pic = df04(itime).zlim([-15 15]);
  hca = h(isub); 
  data = pic.ni;  
  imagesc(hca,pic.xi,pic.zi,data')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_i/n_0';
  hca.Title.String = sprintf('twci = %3.0f',pic.twci);
  if doAdjustCLim
    hca.CLim = clim;
    %colormap(hca,cmap)
  end
  if doA
    hold(hca,'on')
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;
    A = pic.A;
    contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(hca,'off')
  end
    hca.YDir = 'normal';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
    %toc  
  if iframe == 0
    colormap(cmap)
  end
  pause(0.1)
  if 1 % collect frames, for making gif
    iframe = iframe + 1;    
    nframes = ntimes;
    currentBackgroundColor = get(gcf,'color');
    set(gcf,'color',[1 1 1]);
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
%out = {all_im,map};

% collect frames

%% Movie of ni evolution df04, try with VideoWriter
% [all_im, map] = MAKE_GIF(obj,fields,nrows,ncols)      
% make gif
% imwrite(im,map,'delme.gif','DelayTime',0.0,'LoopCount',0)  

vidObj = VideoWriter([savedir 'ni_df04.mp4'],'MPEG-4');
open(vidObj);

% Default options, values
doAdjustCLim = 1;
clim = [0 2];
cmap = pic_colors('candy');
doA = 1;
levA = -25:1:0;

pic = df04;


% setup figure
fig = figure;
fig.Position = [100 100 920 360];
nrows = 1;
ncols = 1;
h = setup_subplots(nrows,ncols); % external function, must include in SMILEI.m
isub = 1;
disp('Adjust figure size, then hit any key to continue.')
pause

times = 2:pic.nt;
ntimes = numel(times);
iframe = 0;
for itime = times
  pic = df04(itime).zlim([-15 15]);
  hca = h(isub); 
  data = pic.ni;  
  imagesc(hca,pic.xi,pic.zi,data')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_i';
  hca.Title.String = sprintf('twci = %3.0f',pic.twci);
  hca.FontSize = 14;
  if doAdjustCLim
    hca.CLim = clim;
    %colormap(hca,cmap)
  end
  if doA
    hold(hca,'on')
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;
    A = pic.A;
    contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(hca,'off')
  end
    hca.YDir = 'normal';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
    %toc  
  if iframe == 0
    colormap(cmap)
  end
  pause(0.1)
  currFrame = getframe(gcf);
  writeVideo(vidObj,currFrame);
end
%out = {all_im,map};
close(vidObj);
% collect frames

%% Movie of ni evolution df04, gif
% [all_im, map] = MAKE_GIF(obj,fields,nrows,ncols)      
% make gif
% imwrite(im,map,'delme.gif','DelayTime',0.0,'LoopCount',0)  

% Default options, values
doAdjustCLim = 1;
clim = [-0.7 0.7];
cmap = pic_colors('blue_red');
doA = 1;
levA = -25:1:0;

pic = df04;


% setup figure
fig = figure;
fig.Position = [100 100 920 360];
set(fig,'color','white');
nrows = 1;
ncols = 1;
h = setup_subplots(nrows,ncols); % external function, must include in SMILEI.m
isub = 1;
disp('Adjust figure size, then hit any key to continue.')
pause

times = 2:pic.nt;
ntimes = numel(times);
iframe = 0;
for itime = times
  pic = df04(itime).zlim([-15 15]);
  hca = h(isub); 
  data = pic.Ey;  
  imagesc(hca,pic.xi,pic.zi,data')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_y';
  hca.Title.String = sprintf('twci = %3.0f',pic.twci);
  if doAdjustCLim
    hca.CLim = clim;
    %colormap(hca,cmap)
  end
  if doA
    hold(hca,'on')
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;
    A = pic.A;
    contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(hca,'off')
  end
    hca.YDir = 'normal';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
    %toc  
  if iframe == 0
    colormap(cmap)
  end
  pause(0.1)
  if 1 % collect frames, for making gif
    iframe = iframe + 1;    
    nframes = ntimes;
    currentBackgroundColor = get(gcf,'color');
    set(gcf,'color',[1 1 1]);
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
%out = {all_im,map};

% collect frames

%% Movie of Ey evolution df04, try with VideoWriter
% [all_im, map] = MAKE_GIF(obj,fields,nrows,ncols)      
% make gif
% imwrite(im,map,'delme.gif','DelayTime',0.0,'LoopCount',0)  

vidObj = VideoWriter([savedir 'Ey_df04.mp4'],'MPEG-4');
open(vidObj);

% Default options, values
doAdjustCLim = 1;
clim = [-0.7 0.7];
cmap = pic_colors('blue_red');
doA = 1;
levA = -25:1:0;

pic = df04;


% setup figure
fig = figure;
fig.Position = [100 100 920 360];
nrows = 1;
ncols = 1;
h = setup_subplots(nrows,ncols); % external function, must include in SMILEI.m
isub = 1;
disp('Adjust figure size, then hit any key to continue.')
pause

times = 2:pic.nt;
ntimes = numel(times);
iframe = 0;
for itime = times
  pic = df04(itime).zlim([-15 15]);
  hca = h(isub); 
  data = pic.Ey;  
  imagesc(hca,pic.xi,pic.zi,data')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'E_y';
  hca.Title.String = sprintf('twci = %3.0f',pic.twci);
  hca.FontSize = 14;
  if doAdjustCLim
    hca.CLim = clim;
    %colormap(hca,cmap)
  end
  if doA
    hold(hca,'on')
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;
    A = pic.A;
    contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(hca,'off')
  end
    hca.YDir = 'normal';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
    %toc  
  if iframe == 0
    colormap(cmap)
  end
  pause(0.1)
  currFrame = getframe(gcf);
  writeVideo(vidObj,currFrame);
end
%out = {all_im,map};
close(vidObj);
% collect frames

%% Movie of temperature evolution df04, try with VideoWriter
% [all_im, map] = MAKE_GIF(obj,fields,nrows,ncols)      
% make gif
% imwrite(im,map,'delme.gif','DelayTime',0.0,'LoopCount',0)  

vidObj = VideoWriter([savedir 'vy5_df04_viquivers.mp4'],'MPEG-4');
open(vidObj);


% Default options, values
doAdjustCLim = 1;
clim = [-1 1];
cmap = pic_colors('candy2');
cmap = pic_colors('blue_red');
doA = 1;
levA = -25:1:0;

doQ = 1;
maxQ = 5;

pic = df04;


% setup figure
fig = figure;
nrows = 1;
ncols = 1;
h = setup_subplots(nrows,ncols); % external function, must include in SMILEI.m
isub = 1;
disp('Adjust figure size, then hit any key to continue.')
%pause
hcf = gcf;
hcf.Position = [301   190   896   500];
set(hcf,'color','white');

iSpecies = 5;
times = 2:pic.nt;
ntimes = numel(times);
iframe = 0;
for itime = times
  itime
  pic = df04(itime).zlim([-15 15]);
  hca = h(isub); 
  %[n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = pic.njp(iSpecies);
  %p = (pxx+pyy+pzz)/3; 
  %t = p./n;
  
   if doQ % Quivers        
    colorQ = pic_colors('matlab');
    colorQ = [0 0 0; 0.5 0.5 0.5];
    maxQs = [2,2];
    scalesQ_a = 10*[1 1 1]; % rescaling to input data
    scalesQ_b = 0*[1 1 1]; % rescaling option for the quiver function
        
    nQx = 50;
    nQz = 25;
    [X,Z] = ndgrid(pic.xi,pic.zi);
    ipxQ = fix(linspace(1,pic.nx,nQx));
    ipzQ = fix(linspace(1,pic.nz,nQz));
    %[X,Z] = meshgrid(ipxQ,ipzQ);
    %ipXQ = dataQx; ipZQ = dataQz;      
  end 

  data = pic.vy(5);
  imagesc(hca,pic.xi,pic.zi,data')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_y (cold protons from the south)';
  hca.Title.String = sprintf('twci = %3.0f',pic.twci);
  hca.FontSize = 14;
  if doAdjustCLim
    hca.CLim = clim;
    %colormap(hca,cmap)
  end
  if doA
    hold(hca,'on')
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;
    A = pic.A;
    contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(hca,'off')
  end
  if doQ    
      hold(hca,'on')
      dataQx = pic.vx(5); dataQx(n<0.01) = NaN;
      dataQz = pic.vz(5); dataQz(n<0.01) = NaN;
      dataQ.x = dataQx;
      dataQ.z = dataQz;
      maxQ = maxQs(1);     
      iQ = 1;
      displayname = 'v_{ix}';
      hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ)*scalesQ_a(iQ),dataQ.z(ipxQ,ipzQ)*scalesQ_a(iQ),scalesQ_b(iQ),...
        'color',colorQ(iQ,:),'linewidth',1.0,'displayname',displayname,...
        'ShowArrowHead','off','Marker','o','MarkerSize',2);
      %hquiv = quiver(hca,X(ipxQ,ipzQ),Z(ipxQ,ipzQ),dataQ.x(ipxQ,ipzQ)*scalesQ_a(iQ),dataQ.z(ipxQ,ipzQ)*scalesQ_a(iQ),scalesQ_b(iQ),...
      %  'color',colorQ(iQ,:),'linewidth',1.0,'displayname',displayname,...
      %  'ShowArrowHead','off','Marker','o','MarkerSize',1);
      %return
        %hquiv.ShowArrowHead = 'off';
        %hquiv.Marker = '.';
        hold(hca,'off')  
  end    
    hca.YDir = 'normal';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
    %toc  
  if iframe == 0
    colormap(cmap)
  end
  pause(0.1)
  currFrame = getframe(gcf);
  writeVideo(vidObj,currFrame);
end
%out = {all_im,map};
close(vidObj);
% collect frames

%% Map of distributions, 8000
zlim = [0 6];
xlim = [160 180];
sumdim = 3;
iSpecies = [3 5];
hmap = ds04(2).xlim(xlim).zlim(zlim).plot_map(iSpecies,sumdim);

all_axes = findobj(gcf,'type','axes');
compact_panels(0.00,0.00)
drawnow

cn.print(sprintf('fxy_35_t08000_z%.0f-%.0f_x%.0f-%.0f',xlim(1),xlim(2),zlim(1),zlim(2)),'path',savedir)

doSquare = 1;
for ipanel = 1:numel(all_axes)
  if doSquare, axis(all_axes(ipanel),'square'); end
end
cn.print(sprintf('fxy_35_t08000_z%.0f-%.0f_x%.0f-%.0f_square',xlim(1),xlim(2),zlim(1),zlim(2)),'path',savedir)

set(gcf,'Position',pp)
drawnow
cn.print(sprintf('fxy_35_t08000_z%.0f-%.0f_x%.0f-%.0f_square_tight',xlim(1),xlim(2),zlim(1),zlim(2)),'path',savedir)

%% Map of distributions, 5000
zlim = [0 2];
xlim = [192 205];
sumdim = 3;
iSpecies = [3 5];
hmap = ds04(1).xlim(xlim).zlim(zlim).plot_map(iSpecies,sumdim);

all_axes = findobj(gcf,'type','axes');
compact_panels(0.00,0.00)
drawnow

cn.print(sprintf('fxy_35_t05000_z%.0f-%.0f_x%.0f-%.0f',xlim(1),xlim(2),zlim(1),zlim(2)),'path',savedir)

doSquare = 1;
for ipanel = 1:numel(all_axes)
  if doSquare, axis(all_axes(ipanel),'square'); end
end
cn.print(sprintf('fxy_35_t05000_z%.0f-%.0f_x%.0f-%.0f_square',xlim(1),xlim(2),zlim(1),zlim(2)),'path',savedir)

set(gcf,'Position',pp)
drawnow
cn.print(sprintf('fxy_35_t05000_z%.0f-%.0f_x%.0f-%.0f_square_tight',xlim(1),xlim(2),zlim(1),zlim(2)),'path',savedir)

%% Movie looping through distributions at z = 0, and x, seeing the subsequently entering "trees"

%% Plot smaller subset of distributions, for smaller figures embedded in text, one z row
doVideo = 1;
doGif = 1;
movieName = 'f_vyz_xy_xz_z=0_x=160-200';
zlim = 0+[-0.1 0.1];
xlim = [160 205];
sumdim = 3;
iSpecies = [3];
ds = ds04(2).xlim(xlim).zlim(zlim);
cmap = pic_colors('candy');
ticks = -15:1:15;
doColorbar = 0;
clim = [0 0.5e-2]; doCLim = 1;
it = 1;


if doVideo
  vidObj = VideoWriter([savedir movieName '.mp4'],'MPEG-4');
  open(vidObj);
end
if doGif
  iframe = 0;
end
% Set up figure
h = setup_subplots(1,3);

for id = ds.nd{it}:-1:1
  f = ds.f(1,id,iSpecies);
  isub = 1;
  if 1 % fzy
    hca = h(isub); isub = isub + 1;
    imagesc(hca,f.v,f.v,f.fyz')
    hca.XLabel.String = 'v_y';
    hca.YLabel.String = 'v_z';
  if 1 % fxy
    hca = h(isub); isub = isub + 1;
    imagesc(hca,f.v,f.v,f.fxy')
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_y';
  end
  if 1 % fxz
    hca = h(isub); isub = isub + 1;
    imagesc(hca,f.v,f.v,f.fxz')
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_z';
  end
  end
  colormap(gcf,cmap)
  for ip = 1:3
    h(ip).YDir = 'normal';
    h(ip).XTick = ticks;
    h(ip).YTick = ticks;
    h(ip).XGrid = 'on';
    h(ip).YGrid = 'on';    
    h(ip).FontSize = 14;
    axis(h(ip),'square')
    if doCLim, h(ip).CLim = clim; end
  end
  h(2).Title.String = sprintf('x = %.1f, z = %.1f',(f.x(1)+f.x(2))/2,(f.z(1)+f.z(2))/2);
  
  hlinks = linkprop(h,{'CLim'});
  if doColorbar
    h3pos = h(3).Position;
    hcb = colorbar('peer',h(3));
    h(3).Position = h3pos;
    hcb.Position(1) = h(3).Position(1)+h(3).Position(3)+0.01;
  end
  pause(0.1)
  if doVideo
    set(gcf,'color','white');
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = pic0.nt;
      currentBackgroundColor = get(gcf,'color');
      set(gcf,'color',[1 1 1]);
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
if doVideo, close(vidObj); end
if doGif, imwrite(all_im,map,[savedir movieName '.gif'],'DelayTime',0.0,'LoopCount',inf); end

%% Plot smaller subset of distributions, for smaller figures embedded in text, two z rows
doVideo = 1;
doGif = 1;
movieName = 'f_vyz_xy_xz_z=0-1_x=160-200_is5_vlim2';
zlim1 = 0+[-0.1 0.1];
zlim2 = 1+[-0.1 0.1];
xlim = [166 200];
xs = 205:-1:166;
vlim = [-2 2];
sumdim = 3;
iSpecies = [1];
%ds1 = ds04(2).xlim(xlim).zlim(zlim1);
%ds2 = ds04(2).xlim(xlim).zlim(zlim2);
cmap = pic_colors('candy');
ticks = -15:1:15;
doColorbar = 0;
clim = [0 0.5e-2]; doCLim = 1; % cold ions
clim = [0 0.1e-2]; doCLim = 1; % hot ions, saturates at tend
it = 1;


if doVideo
  vidObj = VideoWriter([savedir movieName '.mp4'],'MPEG-4');
  open(vidObj);
end
if doGif
  iframe = 0;
end
% Set up figure
h = setup_subplots(2,3);

for ix = 1:numel(xs)
  %ix
  id = 1;
  ds1 = ds04(2).xlim(xs(ix)+[-0.1 0.1]).zlim(zlim1);
  ds2 = ds04(2).xlim(xs(ix)+[-0.1 0.1]).zlim(zlim2);
  if isempty(ds1.xi{:}) && isempty(ds2.xi{:})
    continue
  elseif isempty(ds1.xi{:})
    f2 = ds2.f(1,id,iSpecies);
    f1 = f2;
    f1.f = f1.f*0;
    f1.fxy = f1.fxy*0;
    f1.fyz = f1.fyz*0;
    f1.fxz = f1.fxz*0;
  elseif isempty(ds2.xi{:})
    f1 = ds1.f(1,id,iSpecies);
    f2 = f1;
    f2.f = f2.f*0;
    f2.fxy = f2.fxy*0;
    f2.fyz = f2.fyz*0;
    f2.fxz = f2.fxz*0;
  else
    f1 = ds1.f(1,id,iSpecies);
    f2 = ds2.f(1,id,iSpecies);
  end
   
      
  isub = 1;
  f = f2;
  if 1 % fzy
    hca = h(isub); isub = isub + 1;
    imagesc(hca,f.v,f.v,f.fyz')
    hca.XLabel.String = 'v_y';
    hca.YLabel.String = 'v_z';
  end
  if 1 % fxy
    hca = h(isub); isub = isub + 1;
    imagesc(hca,f.v,f.v,f.fxy')
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_y';
  end
  if 1 % fxz
    hca = h(isub); isub = isub + 1;
    imagesc(hca,f.v,f.v,f.fxz')
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_z';
  end 
  f1 = ds1.f(1,id,iSpecies);
  f = f1;
  if 1 % fzy
    hca = h(isub); isub = isub + 1;
    imagesc(hca,f.v,f.v,f.fyz')
    hca.XLabel.String = 'v_y';
    hca.YLabel.String = 'v_z';
  end
  if 1 % fxy
    hca = h(isub); isub = isub + 1;
    imagesc(hca,f.v,f.v,f.fxy')
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_y';
  end
  if 1 % fxz
    hca = h(isub); isub = isub + 1;
    imagesc(hca,f.v,f.v,f.fxz')
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_z';
  end   
  
  colormap(gcf,cmap)
  for ip = 1:6
    h(ip).YDir = 'normal';
    h(ip).XLim = vlim;
    h(ip).YLim = vlim;
    h(ip).XTick = ticks;
    h(ip).YTick = ticks;
    h(ip).XGrid = 'on';
    h(ip).YGrid = 'on';    
    h(ip).FontSize = 14;
    axis(h(ip),'square')
    if doCLim, h(ip).CLim = clim; end
  end
  h(5).Title.String = sprintf('x = %.1f, z = %.1f',(f1.x(1)+f1.x(2))/2,(f1.z(1)+f1.z(2))/2);
  h(2).Title.String = sprintf('x = %.1f, z = %.1f',(f2.x(1)+f2.x(2))/2,(f2.z(1)+f2.z(2))/2);
  
  hlinks = linkprop(h,{'CLim'});
  if doColorbar
    h3pos = h(3).Position;
    hcb = colorbar('peer',h(3));
    h(3).Position = h3pos;
    hcb.Position(1) = h(3).Position(1)+h(3).Position(3)+0.01;
  end
  pause(0.1)
  if doVideo
    set(gcf,'color','white');
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = pic0.nt;
      currentBackgroundColor = get(gcf,'color');
      set(gcf,'color',[1 1 1]);
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
if doVideo, close(vidObj); end
if doGif, imwrite(all_im,map,[savedir movieName '.gif'],'DelayTime',0.0,'LoopCount',inf); end

%% Plot showing select distributions with appearing fingers, vx vy plane
x_sel = [183 187 191 197 201 203];
x_sel = x_sel(end:-1:1);
z_sel = x_sel*0 + 0;
t_sel = 2; % id
iSpecies = [3 5];
clim = [0 0.5e-2]; doCLim = 1;
cmap = pic_colors('candy');
ticks = -15:1:15;
doBdir = 1;

nrows = 2;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;


for idist = 1:numel(x_sel)    
  xx = x_sel(idist);
  zz = z_sel(idist);
  ds = ds04(2).zlim(zz+[-0.1 0.1]).xlim(xx+[-0.1 0.1]);
  f = ds.f(1,1,iSpecies);
  
  Bx_ = df04.twpelim(8000).xlim(xx).zlim(xx).Bx;
  By_ = df04.twpelim(8000).xlim(xx).zlim(xx).By;
  Bz_ = df04.twpelim(8000).xlim(xx).zlim(xx).Bz;
  
  if 1 % fxy
    hca = h(isub); isub = isub + 1;
    imagesc(hca,f.v,f.v,f.fxy')
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_y';
    %hca.Title.String = sprintf('x = %.1f d_i, z = %.1f d_i',(f.x(1)+f.x(2))/2,(f.z(1)+f.z(2))/2);
    %irf_legend(hca,sprintf('x = %.1f d_i, z = %.1f d_i',(f.x(1)+f.x(2))/2,(f.z(1)+f.z(2))/2),[0.02 0.98])
    irf_legend(hca,sprintf('x = %.0f d_i',(f.x(1)+f.x(2))/2),[0.03 0.98],'fontsize',13)
    if doBdir        
      line_slope = (By_/Bx_);
      xx = min(hca.XLim(2)*[1 1/abs(line_slope)])*[-1 1];
      hold(hca,'on')                
      hBline = plot(hca,xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
      hold(hca,'off')        
    end
  end
end


colormap(gcf,cmap)
for ip = 1:npanels
  h(ip).YDir = 'normal';
  h(ip).XTick = ticks;
  h(ip).YTick = ticks;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';    
  h(ip).FontSize = 14;
  %h(ip).Title.FontSize = 13;
  h(ip).Position(2) = h(ip).Position(2)+0.015;
  axis(h(ip),'square')
  if doCLim, h(ip).CLim = clim; end
end
c_eval('h(?).YTickLabel = '''';',[2:3 5:6])
c_eval('h(?).YLabel.String = '''';',[2:3 5:6])
compact_panels(0.01,0.01)

h(2).Title.String = 'Cold ions originating from the north, z = [-0.5 0.5] d_i';
h(2).Title.String = 'Cold ions originating from the south, z = [-0.5 0.5] d_i';
h(2).Title.String = 'Cold ions originating from the north and south, z = [-0.5 0.5] d_i';
%h(2).Title.FontSize = 12;

%% Plot showing select distributions with appearing fingers, vx vy plane, include map showing location of boxes
% not modified yet
x_sel = [183 187 191 197 201 203];
x_sel = x_sel(end:-1:1);
z_sel = x_sel*0 + 0;
t_sel = 2; % id
iSpecies = [3 5];
clim = [0 0.5e-2]; doCLim = 1;
cmap = pic_colors('candy');
ticks = -15:1:15;
doBdir = 1;

nrows = 2;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;


for idist = 1:numel(x_sel)    
  xx = x_sel(idist);
  zz = z_sel(idist);
  ds = ds04(2).zlim(zz+[-0.1 0.1]).xlim(xx+[-0.1 0.1]);
  f = ds.f(1,1,iSpecies);
  
  Bx_ = df04.twpelim(8000).xlim(xx).zlim(xx).Bx;
  By_ = df04.twpelim(8000).xlim(xx).zlim(xx).By;
  Bz_ = df04.twpelim(8000).xlim(xx).zlim(xx).Bz;
  
  if 1 % fxy
    hca = h(isub); isub = isub + 1;
    imagesc(hca,f.v,f.v,f.fxy')
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_y';
    %hca.Title.String = sprintf('x = %.1f d_i, z = %.1f d_i',(f.x(1)+f.x(2))/2,(f.z(1)+f.z(2))/2);
    %irf_legend(hca,sprintf('x = %.1f d_i, z = %.1f d_i',(f.x(1)+f.x(2))/2,(f.z(1)+f.z(2))/2),[0.02 0.98])
    irf_legend(hca,sprintf('x = %.0f d_i',(f.x(1)+f.x(2))/2),[0.03 0.98],'fontsize',13)
    if doBdir        
      line_slope = (By_/Bx_);
      xx = min(hca.XLim(2)*[1 1/abs(line_slope)])*[-1 1];
      hold(hca,'on')                
      hBline = plot(hca,xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
      hold(hca,'off')        
    end
  end
end


colormap(gcf,cmap)
for ip = 1:npanels
  h(ip).YDir = 'normal';
  h(ip).XTick = ticks;
  h(ip).YTick = ticks;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';    
  h(ip).FontSize = 14;
  %h(ip).Title.FontSize = 13;
  h(ip).Position(2) = h(ip).Position(2)+0.015;
  axis(h(ip),'square')
  if doCLim, h(ip).CLim = clim; end
end
c_eval('h(?).YTickLabel = '''';',[2:3 5:6])
c_eval('h(?).YLabel.String = '''';',[2:3 5:6])
compact_panels(0.01,0.01)

h(2).Title.String = 'Cold ions originating from the north, z = [-0.5 0.5] d_i';
h(2).Title.String = 'Cold ions originating from the south, z = [-0.5 0.5] d_i';
h(2).Title.String = 'Cold ions originating from the north and south, z = [-0.5 0.5] d_i';
%h(2).Title.FontSize = 12;

%% Time-x maps of Bz, Ey and vExBx
zlim = [-0.2 0.2];
pic = df04.zlim(zlim).twcilim([2 260]); % just remove first time index
nc = 0.8;
firstfract = fix(0.4*6400);
Alev = -25:1:0;
doA = 1;
doShowMaxVal = 1;

nrows = 1;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

hb = gobjects(0);

if 1 % Bz
  hca = h(isub); isub = isub + 1;
  %variable = mean(pic.Bz,2);
  Bz_mean = squeeze(mean(pic.Bz,2));
  imagesc(hca,pic.xi,pic.twci,Bz_mean')
  hcb = colorbar('peer',hca);
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'B_z';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = [-1.2 1.2];
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';
  if doShowMaxVal
    Bz_mean = squeeze(mean(Bz,2));
    Bz_mean_left = Bz_mean(1:firstfract,:); % right side has has island
    Bz_meanmax = max(abs(Bz_mean_left(:)));    
    irf_legend(hca,sprintf('|B_z|^{max} = %.1f B_0',Bz_meanmax),[0.98 0.98],'color','k','fontsize',12)
  end
end
if 1 % Ey
  hca = h(isub); isub = isub + 1;  
  Ey_mean = squeeze(mean(pic.Ey,2));
  imagesc(hca,pic.xi,pic.twci,Ey_mean')
  hcb = colorbar('peer',hca);
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'E_y';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = [-0.9 0.9];
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';  
  if doShowMaxVal
    Ey_mean = squeeze(mean(Ey,2));
    Ey_mean_left = Ey_mean(1:firstfract,:); % right side has has island
    Ey_meanmax = max(abs(Ey_mean_left(:)));   
    irf_legend(hca,sprintf('|E_y|^{max} = %.1f v_AB_0',Ey_meanmax),[0.98 0.98],'color','k','fontsize',12)
  end
end
if 1 % ExB
  if 1
    Ex = pic.Ez;
    Ey = pic.Ey;
    Ez = pic.Ez;
    Bx = pic.Bx;
    By = pic.By;
    Bz = pic.Bz;  
    Babs = sqrt(Bx.^2+By.^2+Bz.^2);
  end
  
  ExB = cross_product(Ex,Ey,Ez,Bx./(Babs.^2),By./(Babs.^2),Bz./(Babs.^2));
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.twci,smooth2(squeeze(mean(ExB.x,2)),2,2)')
  hcb = colorbar('peer',hca);
  hb(isub-1) = hcb;
  hcb.YLabel.String = '(ExB/B^2)_x';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = [-1 1];
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';   
  if doShowMaxVal
    ExBxmean = smooth2(squeeze(mean(ExB.x,2)),2,2);
    ExBxmean_left = ExBxmean(1:firstfract,:); % right side has has island
    ExB_meanmax = max(abs(ExBxmean_left(:)));
    %irf_legend(hca,sprintf('|ExB/B^2|_x^{max} = %.1f v_A',ExB_meanmax),[0.98 0.98],'color','k','fontsize',12)
  end 
end

irf_legend(h(1),sprintf('n_c = %g n_0',nc),[0.02 0.98],'color','k','fontsize',12)

for ip = 1:npanels
  hb(ip).Location = 'northoutside';
  h(ip).Position(2) = 0.21;
  h(ip).Position(4) = 0.5;
  h(ip).XTick = 0:100:500;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).FontSize = 14;
  hb(ip).FontSize = 14;
  
end
h(2).XTick = [100:100:400];
h(3).XTick = [100:100:400];
compact_panels(0.01,0.01)
h(2).YLabel.String = '';
h(3).YLabel.String = '';
h(2).YTickLabel = [];
h(3).YTickLabel = [];

%% Front speed
zlim = [-0.5 0.5];
xlim = [50 200]; % just check left side
tlim = [90 220];
pic = df04.xlim(xlim).zlim(zlim).twcilim(tlim);
Bz_mean = squeeze(mean(pic.Bz,2));
[Bz_peak,ind_Bz_peak] = max(abs(Bz_mean));
xDF = pic.xi(ind_Bz_peak);
dt = diff(pic.twci)';
tcentered = tocolumn(pic.twci(1:pic.nt-1))+tocolumn(dt);
dxDF = diff(xDF);
vDF = dxDF./dt;
tcentered_2 = tocolumn(pic.twci(2:pic.nt-1));
aDF = diff(vDF)./dt(2:end);
RA = pic.RA;


nrows = 5;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % Bz vz x
  hca = h(isub); isub = isub + 1;  
  plot(hca,pic.xi,Bz_mean)
  hold(hca,'on')
  plot(hca,pic.xi(ind_Bz_peak),-Bz_peak,'*')
  hold(hca,'off')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'B_z/B_0';
  %hca.Xlim = [];
end
if 1 % DF Bz strength
  hca = h(isub); isub = isub + 1;  
  plot(hca,pic.twci,Bz_peak)
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'B_z/B_0';
end
if 1 % DF location
  hca = h(isub); isub = isub + 1;  
  plot(hca,pic.twci,xDF)
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
end
if 1 % DF speed
  hca = h(isub); isub = isub + 1;  
  plot(hca,tcentered,vDF)
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'v_x/d_i';
end
if 1 % DF acceleration
  hca = h(isub); isub = isub + 1;  
  plot(hca,pic.twci(2:pic.nt-1),aDF)
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'a_x/d_i';
end

for ip = 1:npanels
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end

%% Mass loading effect on front speed
tlim = [90 220];
pic04 = df04.twcilim(tlim);
pic04n = df04n.twcilim(tlim).twcilim(4:4:200,'exact');
pic08 = df08.twcilim(tlim);
[xDF04,vDF04,aDF04,BDF04] = pic04.xva_df;
[xDF04n,vDF04n,aDF04n,BDF04n] = pic04n.xva_df;
[xDF08,vDF08,aDF08,BDF08] = pic08.xva_df;

%%
% Get ExB velocity at the same location;
for iDF = 1:pic04.length
  %pc04 = pic04.xlim(xDF04(1,iDF)+[-0.1 0.1]).zlim([-0.5 0.5]).twcilim(pic04.twci(iDF));
  %Ex04_(iDF) = mean(mean(pc04.Ex));
  %Ey04_(iDF) = mean(mean(pc04.Ey));  
end
%vExB04 = 

%%
  Ex04 = pic04.interp(xDF04(1,:),xDF04(1,:)*0,pic04.twci,'Ex');
  Ey04 = pic04.interp(xDF04(1,:),xDF04(1,:)*0,pic04.twci,'Ey');
  Ez04 = pic04.interp(xDF04(1,:),xDF04(1,:)*0,pic04.twci,'Ez');
  Bx04 = pic04.interp(xDF04(1,:),xDF04(1,:)*0,pic04.twci,'Bx');
  By04 = pic04.interp(xDF04(1,:),xDF04(1,:)*0,pic04.twci,'By');
  Bz04 = pic04.interp(xDF04(1,:),xDF04(1,:)*0,pic04.twci,'Bz');
  Ex08 = pic08.interp(xDF08(1,:),xDF08(1,:)*0,pic04.twci,'Ex');
  Ey08 = pic08.interp(xDF08(1,:),xDF08(1,:)*0,pic04.twci,'Ey');
  Ez08 = pic08.interp(xDF08(1,:),xDF08(1,:)*0,pic04.twci,'Ez');
  Bx08 = pic08.interp(xDF08(1,:),xDF08(1,:)*0,pic04.twci,'Bx');
  By08 = pic08.interp(xDF08(1,:),xDF08(1,:)*0,pic04.twci,'By');
  Bz08 = pic08.interp(xDF08(1,:),xDF08(1,:)*0,pic04.twci,'Bz');
%%
Babs04 = sqrt(Bx04.^2+By04.^2+Bz04.^2);
Babs08 = sqrt(Bx08.^2+By08.^2+Bz08.^2);
ExB04 = cross_product(Ex04,Ey04,Ez04,Bx04./Babs04,By04./Babs04,Bz04./Babs04);
ExB08 = cross_product(Ex08,Ey08,Ez08,Bx08./Babs08,By08./Babs08,Bz08./Babs08);
%%
%rem_ind04 = find(abs(diff(vDF04))>0.2); vDF04(rem_ind04) = NaN;
rem_ind08 = find(abs(aDF08)>0.05); vDF08(rem_ind08) = NaN;

R0 = @(Rc,nc,nb) Rc*sqrt(1+nc./nb);

Rc = @(R0,nc,nb) R0sqrt(1+nc./nb);

R0(0.1,0.4,0.2);

aDF0 = @(ac,nc,nb) ac*(1+nc./nb)^0.5;

doPlotExB = 1;
doPlotFit = 1;

nrows = 1;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 0 % vDF_x vz t
  hca = h(isub); isub = isub + 1;  
  hl = plot(hca,pic04.twci,vDF04(1,:),'.',...
                pic08.twci,vDF08(1,:),'.','MarkerSize',10);
  hold(hca,'on')
  %plot(hca,pic.xi(ind_Bz_peak),-Bz_peak,'*')
  hold(hca,'off')
  hca.XLabel.String = 't\omega_{wci}';
  hca.YLabel.String = 'v_{DF,x}/v_A';
  %hca.Xlim = [];
  if doPlotFit
    p04 = polyfit(pic04.twci(~isnan(vDF04(1,:))),vDF04(1,~isnan(vDF04(1,:))),1);
    p08 = polyfit(pic08.twci(~isnan(vDF08(1,:))),vDF08(1,~isnan(vDF08(1,:))),1);
 
    % Evaluate the fitted polynomial p and plot:
    f04 = polyval(p04,pic04.twci);
    f08 = polyval(p08,pic08.twci);
    hold(hca,'on')
    hlines = plot(hca,pic04.twci,f04,'-',...
                      pic08.twci,f08,'-');
    hlines(1).Color = hl(1).Color;
    hlines(2).Color = hl(2).Color;
    hold(hca,'off')
    %legend('data','linear fit')
    fit_str04 = sprintf('n_c= 0.4n_0: v_{DF}/v_A = %.4f tw_{ci}',p04(1));
    fit_str08 = sprintf('n_c= 0.8n_0: v_{DF}/v_A = %.4f tw_{ci}',p08(1));
    irf_legend(hca,{fit_str04;fit_str08},[0.02 0.98],'fontsize',12)
  end
  if doPlotExB
    hold(hca,'on')
    hlines = plot(hca,pic04.twci,ExB04.x,'-',...
                      pic08.twci,ExB08.x,'-');
    hlines(1).Color = hl(1).Color;
    hlines(2).Color = hl(2).Color;
    hold(hca,'off')
  end
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % vDF_x vz t, mass loading scaled
  hca = h(isub); isub = isub + 1;  
  hl = plot(hca,pic04.twci,vDF04(1,:),'.',...
                pic08.twci,vDF08(1,:),'.','MarkerSize',10);
  hold(hca,'on')
  %plot(hca,pic.xi(ind_Bz_peak),-Bz_peak,'*')
  hold(hca,'off')
  hca.XLabel.String = 't\omega_{wci}';
  hca.YLabel.String = 'v_{DF,x}/v_A';
  %hca.Xlim = [];
  if doPlotFit
    p04 = polyfit(pic04.twci(~isnan(vDF04(1,:))),vDF04(1,~isnan(vDF04(1,:))),1);
    p08 = polyfit(pic08.twci(~isnan(vDF08(1,:))),vDF08(1,~isnan(vDF08(1,:))),1);
 
    % Evaluate the fitted polynomial p and plot:
    f04 = polyval(p04,pic04.twci);
    f08 = polyval(p08,pic08.twci);
    hold(hca,'on')
    hlines = plot(hca,pic04.twci,f04,'-',...
                      pic08.twci,f08,'-');
    hlines(1).Color = hl(1).Color;
    hlines(2).Color = hl(2).Color;
    hold(hca,'off')
    %legend('data','linear fit')
    fit_str04 = sprintf('n_c= 0.4n_0: v_{DF}/v_A = %.4f tw_{ci}',p04(1));
    fit_str08 = sprintf('n_c= 0.8n_0: v_{DF}/v_A = %.4f tw_{ci}',p08(1));
    irf_legend(hca,{fit_str04;fit_str08},[0.02 0.98],'fontsize',12)
  end
  if doPlotExB
    hold(hca,'on')
    hExB = plot(hca,pic04.twci,abs(ExB04.x)+0.2,'--',...
                      pic08.twci,abs(ExB08.x)+0.2,'--');
    hExB(1).Color = hl(1).Color;
    hExB(2).Color = hl(2).Color;
    hold(hca,'off')
  end
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % aDF_x vz t
  hca = h(isub); isub = isub + 1;  
  plot(hca,pic.twci,aDF04(1,:),'.-',...
           pic.twci,aDF08(1,:),'.-')
  hold(hca,'on')
  %plot(hca,pic.xi(ind_Bz_peak),-Bz_peak,'*')
  hold(hca,'off')
  hca.XLabel.String = 't\omega_{wci}';
  hca.YLabel.String = 'a_{DF,x}';
  %hca.Xlim = [];
end


for ip = 1:npanels  
  h(ip).Position(2) = 0.21;
  h(ip).Position(4) = 0.7;
  %h(ip).XTick = 0:100:500;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).FontSize = 14;
  %hb(ip).FontSize = 14;
  
end

%% Particle trajectories corresponding to cold ion fingers, interpolated EB
% find spots with increased phase space density

it = 4;
if 1
iSpecies = [3];
%ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([-0.1 0.6]);
ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([0.4 0.6]); % top row.
else % electrons
  iSpecies = [4];
  %ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([-0.1 0.6]);
  ds = ds04(it).xlim([166 169]+[-0.1 0.1]).zlim([0.4 0.6]); % top row.

end
nPeaks = 10;
spacingPeaks = 4; % for ions its 0.2 vA
fpeaks = ds.get_peaks(nPeaks,spacingPeaks,iSpecies);
nDists = ds.nd;
doPlot = 1;
if doPlot
  % plot results
  for id = 1:ds.nd{1}
  f = ds.f(1,id,iSpecies);
  figure(27)
  h = setup_subplots(3,1);
  hca = h(1);
  imagesc(hca,f.v,f.v,f.fxy')
  hca.YDir = 'normal';
  colormap(pic_colors('candy'))
  hold(hca,'on')
  plot(hca,[fpeaks(:,id).vx],[fpeaks(:,id).vy],'k.')
  hold(hca,'off')

  hca = h(2);
  imagesc(hca,f.v,f.v,f.fxz')
  hca.YDir = 'normal';
  colormap(pic_colors('candy'))
  hold(hca,'on')
  plot(hca,[fpeaks(:,id).vx],[fpeaks(:,id).vz],'k.')
  hold(hca,'off')

  hca = h(3);
  imagesc(hca,f.v,f.v,f.fyz')
  hca.YDir = 'normal';
  colormap(pic_colors('candy'))
  hold(hca,'on')
  plot(hca,[fpeaks(:,id).vy],[fpeaks(:,id).vz],'k.')
  hold(hca,'off')
  pause(0.5)
  end
end

%% Loop through points, integrate trajectories
pic = df04;
tspan = [60,160,210];
m = 1; 
q = 1;
tic
for id = 1:ds.nd{1}
  for iPeak = 1:nPeaks
    fprintf('(id/nd,ipeak/npeaks) = (%g/%g,%g/%g)\n',id,ds.nd{1},iPeak,nPeaks)
    r0 = [fpeaks(iPeak,id).x, fpeaks(iPeak,id).y, fpeaks(iPeak,id).z];
    v0 = [fpeaks(iPeak,id).vx, fpeaks(iPeak,id).vy, fpeaks(iPeak,id).vz];
    tr_tmp = df04.integrate_trajectory(r0,v0,tspan,m,q);    
    [Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(tr_tmp.x,tr_tmp.z,tr_tmp.t);  % interpolate
    
    tr_tmp.Ex = Ex;
    tr_tmp.Ey = Ey;
    tr_tmp.Ez = Ez;
    tr_tmp.Bx = Bx;
    tr_tmp.By = By;
    tr_tmp.Bz = Bz;
    tr(iPeak,id) = tr_tmp;
    toc
    %catch
    %  continue
    %end
  end
end
% subset saved (20,5)
% save('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/trajectories','tr')
% (1,25) save('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/trajectories2','tr','fpeaks','ds')
%% Plot what we have, all trajectories
hca = subplot(1,1,1);
hold(hca,'on')
for id = 1:size(tr,2)
  for iPeak = 1:nPeaks
    plot3(tr(iPeak,id).x,tr(iPeak,id).y,tr(iPeak,id).z)
    %plot3(tr_pass(iPeak,id).x,tr_pass(iPeak,id).y,tr_pass(iPeak,id).z)
    %plot(tr(iPeak,id).t,tr(iPeak,id).vz)
  end
end
hold(hca,'off')

%% Plot what we have, all trajectories from Traj
hca = subplot(1,1,1);
hold(hca,'on')
for iTr = 1:tr04.ntr  
  %plot3(tr04(iTr).x,tr04(iTr).y,tr04(iTr).z)
  plot(tr04(iTr).x,tr04(iTr).z)
end
hold(hca,'off')

%% Plot what we have, v binning
clim = [130 190];
vlim = [-2 2];
xx = 180:203;
xx = 188;

for ix = 1:numel(xx)
nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

hb = gobjects(0);


vx_all = [];
vy_all = [];
vz_all = [];
vx0_all = [];
vy0_all = [];
vz0_all = [];
t_all = [];
if 1 % 3 panels
  
  
  xlim = xx(ix) + 1*[-0.25 0.25];
  zlim = 0 + 2*[-0.25 0.25];
  for id = 1:size(tr,2)
    for iPeak = 1:size(tr,1)
      if 0%isempty(intersect(iPeakDist,[iPeak,id],'rows'))
        disp(sprintf('Skipping (iPeak,id) = (%g,%g)',iPeak,id))
        continue; 
      end
      % Complete trajectories
      ix = intersect(find(tr(iPeak,id).x<xlim(2)),find(tr(iPeak,id).x>xlim(1)));
      iz = intersect(find(tr(iPeak,id).z<zlim(2)),find(tr(iPeak,id).z>zlim(1)));
      ixz = intersect(ix,iz);
      vx_all = [vx_all; tr(iPeak,id).vx(ixz)];
      vy_all = [vy_all; tr(iPeak,id).vy(ixz)];
      vz_all = [vz_all; tr(iPeak,id).vz(ixz)];
      t_all = [t_all; tr(iPeak,id).t(ixz)];
      % Starting points
      ix0 = intersect(find(tr(iPeak,id).x0<xlim(2)),find(tr(iPeak,id).x0>xlim(1)));
      iz0 = intersect(find(tr(iPeak,id).z0<zlim(2)),find(tr(iPeak,id).z0>zlim(1)));
      ixz0 = intersect(ix0,iz0);
      if not(isempty(ixz0))
        vx0_all = [vx0_all; tr(iPeak,id).vx0];
        vy0_all = [vy0_all; tr(iPeak,id).vy0];
        vz0_all = [vz0_all; tr(iPeak,id).vz0];
      end
    end
  end
  
  hca = h(isub); isub = isub + 1;
  scatter(hca,vx_all,vy_all,5,t_all)
  colormap(hca,pic_colors('waterfall'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'time of passing';
  hca.CLim = clim;
  hold(hca,'on')
  plot(hca,vx0_all,vy0_all,'k+','MarkerSize',5)
  hold(hca,'off')  
  hca.XLim = vlim;
  hca.YLim = vlim;
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  hca.Title.String = sprintf('x = [%g,%g], z = [%g,%g]',xlim(1),xlim(2),zlim(1),zlim(2));
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  hca = h(isub); isub = isub + 1;
  scatter(hca,vx_all,vz_all,5,t_all)
  colormap(hca,pic_colors('waterfall'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'time of passing';
  hca.CLim = clim;
  hold(hca,'on')
  plot(hca,vx0_all,vz0_all,'k+','MarkerSize',5)
  hold(hca,'off')  
  hca.XLim = vlim;
  hca.YLim = vlim;
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
  hca.Title.String = sprintf('x = [%g,%g], z = [%g,%g]',xlim(1),xlim(2),zlim(1),zlim(2));
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YGrid = 'on';
  
  hca = h(isub); isub = isub + 1;
  scatter(hca,vy_all,vz_all,5,t_all)
  colormap(hca,pic_colors('waterfall'))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'time of passing';
  hca.CLim = clim;
  hold(hca,'on')
  plot(hca,vy0_all,vz0_all,'k+','MarkerSize',5)
  hold(hca,'off')  
  hca.XLim = vlim;
  hca.YLim = vlim;
  hca.XLabel.String = 'v_y';
  hca.YLabel.String = 'v_z';
  hca.Title.String = sprintf('x = [%g,%g], z = [%g,%g]',xlim(1),xlim(2),zlim(1),zlim(2));
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  for ip = 1:npanels
    h(ip).XTick = -5:0.5:5;
    h(ip).YTick = -5:0.5:5;
    h(ip).Box = 'on';
    axis(h(ip),'square')
  end
end
pause(0.02)


end
  
%% Find trajectories that meet certain criteria
xpass = 188 + [-0.5 0.5];
zpass = 0 + [-0.5 0.5];
vxpass = [-1 0];
vypass = [0 1];

tr_fields = fieldnames(tr)';
tr_fields{2,1} = {};
tr_pass = struct(tr_fields{:});
iPeakDist = [];
%tr_pass = struct();
itr = 1;
for id = 1:size(tr,2)
  for iPeak = 1:size(tr,1)
    % Check criteria
    % Position
    ix = intersect(find(tr(iPeak,id).x<xpass(2)),find(tr(iPeak,id).x>xpass(1))); if isempty(ix); continue; end    
    iz = intersect(find(tr(iPeak,id).z<zpass(2)),find(tr(iPeak,id).z>zpass(1))); if isempty(iz); continue; end
    ixz = intersect(ix,iz); if isempty(ixz); continue; end
    % Velocities
    ivx = intersect(find(tr(iPeak,id).vx<vxpass(2)),find(tr(iPeak,id).vx>vxpass(1))); if isempty(ivx); continue; end    
    ivy = intersect(find(tr(iPeak,id).vy<vypass(2)),find(tr(iPeak,id).vy>vypass(1))); if isempty(ivy); continue; end
    ivxy = intersect(ivx,ivy); if isempty(ivxy); continue; end
    
    tr_pass(itr) = tr(iPeak,id); itr = itr + 1;
    iPeakDist(end+1,1:2) = [iPeak,id];
  end

  
  if 0 % plot
  hca = h(isub); isub = isub + 1;
  plot(hca,vx_all,vy_all,'.',vx0_all,vy0_all,'*')
  hca.XLim = vlim;
  hca.YLim = vlim;
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  hca.Title.String = sprintf('x = [%g,%g], z = [%g,%g]',xlim(1),xlim(2),zlim(1),zlim(2));
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  hca = h(isub); isub = isub + 1;
  plot(hca,vx_all,vz_all,'.',vx0_all,vz0_all,'*')
  hca.XLim = vlim;
  hca.YLim = vlim;
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
  hca.Title.String = sprintf('x = [%g,%g], z = [%g,%g]',xlim(1),xlim(2),zlim(1),zlim(2));
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  hca = h(isub); isub = isub + 1;
  plot(hca,vy_all,vz_all,'.',vy0_all,vz0_all,'*')
  hca.XLim = vlim;
  hca.YLim = vlim;
  hca.XLabel.String = 'v_y';
  hca.YLabel.String = 'v_z';
  hca.Title.String = sprintf('x = [%g,%g], z = [%g,%g]',xlim(1),xlim(2),zlim(1),zlim(2));
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  end
end

%%
hca = subplot(1,1,1);
hold(hca,'on')
for itr = 1:numel(tr_pass)
  plot3(tr(itr).x,tr(itr).y,tr(itr).z)    
end
hold(hca,'off')

%% Plot forces , trajectories, and v binning for individual trajectories, based on array tr
traj = tr_pass;
nTr = numel(traj);
for iTr = 141:nTr  
  tr = traj(iTr);
  tr.t0 = 160;
  tic
  its = 1:5:numel(tr.t);
  [Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(tr.x(its),tr.z(its),tr.t(its)); 
  vxB = cross_product(tr.vx(its),tr.vy(its),tr.vz(its),Bx,By,Bz,'components',1);
  toc
  
  h = setup_subplots(3,5,'vertical');
  isub = 1;
  cmap = pic_colors('waterfall');
  t0_msize = 10;

  %cmap = interp1(linspace(1,64,size(cmap,1)),cmap,1:numel(tr.t)); 
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
  end
  if 1 % (vx,vy), colorcoded by x
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
  if 1 % (vx,vz), colorcoded by x
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
  if 1 % (vy,vz), colorcoded by x
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
  if 1 % (vx,vy), colorcoded by z
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
  if 1 % (vx,vz), colorcoded by z
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
  if 1 % (vy,vz), colorcoded by z
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
  if 0 % xyz
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
  if 1 % xyz(t)
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t,tr.x-200,tr.t,tr.y,tr.t,tr.z)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'r';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    hold(hca,'on')
    plot(hca,tr.t0*[1 1],hca.YLim,'k--')
    hold(hca,'off')
    legend(hca,{'x-200','y','z','t0'},'location','best')
  end
  if 1 % v 
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t,tr.vx,tr.t,tr.vy,tr.t,tr.vz)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'v';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    hold(hca,'on')
    plot(hca,tr.t0*[1 1],hca.YLim,'k--')
    hold(hca,'off')
    legend(hca,{'x','y','z','t0'},'location','best')
  end
  if 1 % B
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t(its),Bx,tr.t(its),By,tr.t(its),Bz)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'B';
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    hold(hca,'on')
    plot(hca,tr.t0*[1 1],hca.YLim,'k--')
    hold(hca,'off')
    legend(hca,{'x','y','z','t0'},'location','best')
  end
  if 1 % E 
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t(its),Ex,tr.t(its),Ey,tr.t(its),Ez)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'E';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    hold(hca,'on')
    plot(hca,tr.t0*[1 1],hca.YLim,'k--')
    hold(hca,'off')
    legend(hca,{'x','y','z','t0'},'location','best')
  end
  if 1 % vxB
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t(its),vxB.x,tr.t(its),vxB.y,tr.t(its),vxB.z)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'vxB';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    hold(hca,'on')
    plot(hca,tr.t0*[1 1],hca.YLim,'k--')
    hold(hca,'off')
    legend(hca,{'x','y','z','t0'},'location','best')
  end
  if 1 % vxB - components
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t(its),vxB.x_yz,tr.t(its),vxB.x_zy,...
             tr.t(its),vxB.y_zx,tr.t(its),vxB.y_xz,...
             tr.t(its),vxB.z_xy,tr.t(its),vxB.z_yx)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'vxB';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    hold(hca,'on')
    plot(hca,tr.t0*[1 1],hca.YLim,'k--')
    hold(hca,'off')
    pos = hca.Position;
    legend(hca,{'v_yxB_z','v_zxB_y','v_zxB_x','v_xxB_z','v_xxB_y','v_yxB_z','t0'},'location','eastoutside')
    hca.Position = pos;
  end

  h(1).Title.String = {sprintf('(iDist,iPeak) = (%g,%g)',iPeakDist(iTr,2),iPeakDist(iTr,1)),...
    sprintf('[x_0,z_0] = [%g,%g], [vx_0,vy_0,vz_0] = [%.2f,%.2f,%.2f]',tr.x0,tr.z0,tr.vx0,tr.vy0,tr.vz0)};
  %[Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(xvtb.x,xvtb.z,xvtb.t);

  c_eval('h(?).CLim = [100 210];',1:3)
  c_eval('h(?).CLim = [130 220];',4:6)
  c_eval('h(?).CLim = [-3 3];',7:9)
  h(10).YLim = [-80 30]; % x
  h(11).YLim = [-1.8 1.8]; % v
  h(12).YLim = [-1 1]; % B
  h(13).YLim = [-0.5 0.5];
  h(14).YLim = [-0.5 0.5];
  h(15).YLim = [-0.5 0.5];
  for ip = 1:9
    h(ip).XLim = [-2 2];
    h(ip).YLim = [-2 2];
    axis(h(ip),'square')
    h(ip).XTick = -15:0.5:15;
    h(ip).YTick = -15:0.5:15;
  end
  for ip = 10:15
    h(ip).XLim = [100 210];
  end
  for ip = 13:15
    ylim_tmp = h(ip).YLim;
    h(ip).YTick = -2:0.2:2;
    h(ip).YLim = ylim_tmp;
  end
  print_str =  sprintf('traj_iD%g_iP%g_x0%g_z0%g_vx0_%.3f_vy0%.3f_vz0%.3f',iPeakDist(iTr,2),iPeakDist(iTr,1),tr.x0,tr.z0,tr.vx0,tr.vy0,tr.vz0);
  cn.print(print_str,'path',[savedir '/traj/'])
end


%% Plot on top of field, option to make movie
doVideo = 1;
xlim = [120 240];
zlim = [-10 10];
twci = [100 210];
pic = df04.twcilim(twci).xlim(xlim).zlim(zlim);
doA = 1; Alev = -25:1:0;

if doVideo
  vidObj = VideoWriter([savedir 'trajectories_vy_dfions.mp4'],'MPEG-4');
  open(vidObj);
end

nrows = 2;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);

for it = 1:pic.nt
  pc = pic(it);
  isub = 1;
  if 0 % Ez
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.Ez')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'E_z';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % By
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.By')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'B_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % Ey
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.Ey')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'E_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % Bz
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.Bz')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'B_z';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 1 % v35y
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.vy([3 5])')
    hca.CLim = [-.5 0.5];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'v_{y,i,cold}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % vxB_y
    hca = h(isub); isub = isub + 1;
    vx = pc.vx([3 5]);
    vy = pc.vy([3 5]);
    vz = pc.vz([3 5]);
    Bx = pc.Bx;
    By = pc.By;
    Bz = pc.Bz;
    vxB = cross_product(vx,vy,vz,Bx,By,Bz);
    imagesc(hca,pc.xi,pc.zi,vxB.y')
    hca.CLim = [-.5 0.5];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '(vxB)_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % n vy T
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = njp(obj,iSpecies);
    hca = h(isub); isub = isub + 1;    
    
    imagesc(hca,pc.xi,pc.zi,vxB.y')
    hca.CLim = [-.5 0.5];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '(vxB)_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end 
  
  h(1).Title.String = sprintf('twci = %g',pc.twci);
  drawnow
  compact_panels(0.02)
  for ip = 1:npanels
    hca = h(ip);
    hold(hca,'on')
    hca.FontSize = 14;
    hca.YDir ='normal';
    if doA
      iAx = 1:4:pic.nx;
      iAz = 1:4:pic.nz;
      A = pc.A;
      contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',Alev,'k');        
    end
    if 1 % all trajectories
    for itr = 1:numel(tr)      
      plot(hca,tr(itr).x,tr(itr).z)
      ii = find(abs(tr(itr).t-pc.twci)==min(abs(tr(itr).t-pc.twci)));
      plot(hca,tr(itr).x(ii),tr(itr).z(ii),'ko')
    end
    else %subset
      for itr = 1:size(iPeakDist,1)
        %plot(hca,tr(iPeakDist(itr,1),iPeakDist(itr,2)).x,tr(iPeakDist(itr,1),iPeakDist(itr,2)).z)
        ii = find(abs(tr(iPeakDist(itr,1),iPeakDist(itr,2)).t-pc.twci)==min(abs(tr(iPeakDist(itr,1),iPeakDist(itr,2)).t-pc.twci)));
        plot(hca,tr(iPeakDist(itr,1),iPeakDist(itr,2)).x(ii),tr(iPeakDist(itr,1),iPeakDist(itr,2)).z(ii),'ko')
      end
    end
    hold(hca,'off')
  end
  pause(1)
  if doVideo
    set(gcf,'color','white');
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
end
if doVideo, close(vidObj); end
  
%% Few manually defined orbits
m = 1; q = 1;
tr(1) = df04.integrate_trajectory([203,0,0],[-0.0,0.4,0.5],[120,160,200],m,q);
tr(2) = df04.integrate_trajectory([203,0,0],[-0.0,0.2,0.5],[120,160,200],m,q);
tr(3) = df04.integrate_trajectory([197,0,0],[-0.3,-0.25,-0.6],[120,160,200],m,q);
tr(4) = df04.integrate_trajectory([192,0,0],[-1,0.25,0.4],[120,160,200],m,q);
%[Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(tr.x,tr.z,tr.t);
% vxB = cross_product(tr.vx,tr.vy,tr.vz,Bx,By,Bz,'components',1);
%% plot single trajectory
h = setup_subplots(3,5,'vertical');
isub = 1;
cmap = pic_colors('waterfall');
t0_msize = 10;
%cmap = interp1(linspace(1,64,size(cmap,1)),cmap,1:numel(tr.t)); 
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
end
if 1 % (vx,vy), colorcoded by x
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
if 1 % (vx,vz), colorcoded by x
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
if 1 % (vy,vz), colorcoded by x
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
if 1 % (vx,vy), colorcoded by z
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
if 1 % (vx,vz), colorcoded by z
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
if 1 % (vy,vz), colorcoded by z
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
if 0 % xyz
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
if 1 % xyz(t)
  hca = h(isub); isub = isub + 1;
  plot(hca,tr.t,tr.x,tr.t,tr.y,tr.t,tr.z)
  hca.XLabel.String = 't';
  hca.YLabel.String = 'r';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
  legend(hca,{'x','y','z'},'location','best')
end
if 1 % v 
  hca = h(isub); isub = isub + 1;
  plot(hca,tr.t,tr.vx,tr.t,tr.vy,tr.t,tr.vz)
  hca.XLabel.String = 't';
  hca.YLabel.String = 'v';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
  legend(hca,{'x','y','z'},'location','best')
end
if 1 % B
  hca = h(isub); isub = isub + 1;
  plot(hca,tr.t,Bx,tr.t,By,tr.t,Bz)
  hca.XLabel.String = 't';
  hca.YLabel.String = 'B';
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
  legend(hca,{'x','y','z'},'location','best')
end
if 1 % E 
  hca = h(isub); isub = isub + 1;
  plot(hca,tr.t,Ex,tr.t,Ey,tr.t,Ez)
  hca.XLabel.String = 't';
  hca.YLabel.String = 'E';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
  legend(hca,{'x','y','z'},'location','best')
end
if 1 % vxB
  hca = h(isub); isub = isub + 1;
  plot(hca,tr.t,vxB.x,tr.t,vxB.y,tr.t,vxB.z)
  hca.XLabel.String = 't';
  hca.YLabel.String = 'vxB';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';  
  legend(hca,{'x','y','z'},'location','best')
end

%[Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(xvtb.x,xvtb.z,xvtb.t);

for ip = 1:9
  h(ip).XLim = [-2 2];
  h(ip).YLim = [-2 2];
  axis(h(ip),'square')
  h(ip).XTick = -15:0.5:15;
  h(ip).YTick = -15:0.5:15;
end

%% plot multiple trajectories
h = setup_subplots(1,3,'vertical');
isub = 1;
cmap = pic_colors('waterfall');
t0_msize = 10;
t_msize = 4;
if 1 % (vx,vz), colorcoded by time
  hca = h(isub); isub = isub + 1;
  plot(hca,tr1.vx,tr1.vy,'k')
  hold(hca,'on')
  plot(hca,tr2.vx,tr2.vy,'k')
  plot(hca,tr3.vx,tr3.vy,'k')
  plot(hca,tr1.vx0,tr1.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
  plot(hca,tr2.vx0,tr2.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
  plot(hca,tr3.vx0,tr3.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
  scatter(hca,tr1.vx,tr1.vy,t_msize,tr1.t,'filled')
  scatter(hca,tr2.vx,tr2.vy,t_msize,tr2.t,'filled')
  scatter(hca,tr3.vx,tr3.vy,t_msize,tr3.t,'filled')
  hold(hca,'off')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'twci';
  colormap(hca,cmap)
  hca.XLabel.String = 'vx';
  hca.YLabel.String = 'vy';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % (vx,vz), colorcoded by x
  hca = h(isub); isub = isub + 1;
  plot(hca,tr1.vx,tr1.vy,'k')
  hold(hca,'on')
  plot(hca,tr2.vx,tr2.vy,'k')
  plot(hca,tr3.vx,tr3.vy,'k')
  plot(hca,tr1.vx0,tr1.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
  plot(hca,tr2.vx0,tr2.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
  plot(hca,tr3.vx0,tr3.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
  scatter(hca,tr1.vx,tr1.vy,t_msize,tr1.x,'filled')
  scatter(hca,tr2.vx,tr2.vy,t_msize,tr2.x,'filled')
  scatter(hca,tr3.vx,tr3.vy,t_msize,tr3.x,'filled')
  hold(hca,'off')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'x/di';
  colormap(hca,cmap)
  hca.XLabel.String = 'vx';
  hca.YLabel.String = 'vy';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % (vx,vz), colorcoded by z
  hca = h(isub); isub = isub + 1;
  plot(hca,tr1.vx,tr1.vy,'k')
  hold(hca,'on')
  plot(hca,tr2.vx,tr2.vy,'k')
  plot(hca,tr3.vx,tr3.vy,'k')
  plot(hca,tr1.vx0,tr1.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
  plot(hca,tr2.vx0,tr2.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
  plot(hca,tr3.vx0,tr3.vy0,'ks','markersize',t0_msize,'markerfacecolor',[0 0 0])
  scatter(hca,tr1.vx,tr1.vy,t_msize,tr1.z,'filled')
  scatter(hca,tr2.vx,tr2.vy,t_msize,tr2.z,'filled')
  scatter(hca,tr3.vx,tr3.vy,t_msize,tr3.z,'filled')
  hold(hca,'off')
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'z/di';
  colormap(hca,cmap)
  hca.XLabel.String = 'vx';
  hca.YLabel.String = 'vy';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end

for ip = 1:3
  h(ip).XLim = [-2 2];
  h(ip).YLim = [-2 2];
  axis(h(ip),'square')
  h(ip).XTick = -15:0.5:15;
  h(ip).YTick = -15:0.5:15;
end

%% Particle trajectoreis corresponding to cold ion fingers, constant EB (E smoothed)

xvtf = df04.twcilim(160).integrate_trajectory_constant_EB([192,0,0],[-1,0.25,0.4],160,200,25,1);
xvtb = df04.twcilim(160).integrate_trajectory_constant_EB([192,0,0],[-1,0.25,0.4],160,120,25,1);

[Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(xvtb.x,xvtb.z,xvtb.t);

%% Calculate and plot lorentz forces
Alev = -25:1:0;
zlim = [-10 10];
xlim = 200 + 100*[-1 1];
tlim = 8000; % wpe
pic = df04.twpelim(tlim).xlim(xlim).zlim(zlim);

A = pic.A;
Ex = pic.Ex;
Ey = pic.Ey;
Ez = pic.Ez;
Bx = pic.Bx;
By = pic.By;
Bz = pic.Bz;
vx3 = pic.vx([3]);
vy3 = pic.vy([3]);
vz3 = pic.vz([3]);
vx5 = pic.vx([5]);
vy5 = pic.vy([5]);
vz5 = pic.vz([5]);
vx35 = pic.vx([3 5]);
vy35 = pic.vy([3 5]);
vz35 = pic.vz([3 5]);

v3xB = cross_product(vx3,vy3,vz3,Bx,By,Bz,'components',1);
v5xB = cross_product(vx5,vy5,vz5,Bx,By,Bz,'components',1);
v35xB = cross_product(vx35,vy35,vz35,Bx,By,Bz,'components',1);

%% v_xyz, vxB_xyz, E_xyz
nrows = 4; 
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

if 0 % v3x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vx3')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_x 3';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v3y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vy3')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_y 3';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % v3z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vz3')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_z 3';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v5x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vx5')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_x 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v5y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vy5')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_y 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % v5z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vz5')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_z 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v35x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vx35')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_x 35';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v35y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vy35')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_y 35';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % v35z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vz35')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_z 35';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end

if 1 % v3xB.x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v3xB.x')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_x 3';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % v3xB.y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v3xB.y')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_y';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v3xB.z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v3xB.z')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_z';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v3xB.x_yz
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v3xB.x_yz')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_yxB_z';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v3xB.x_zy
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v3xB.x_zy')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_zxB_y';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end

if 1 % v5xB.x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v5xB.x')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_x 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % v5xB.y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v5xB.y')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_y 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v5xB.z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v5xB.z')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_z 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % v35xB.x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v35xB.x')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_x 35';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % v35xB.y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v35xB.y')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_y';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v35xB.z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v35xB.z')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_z';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end

if 0 % Ex
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,Ex')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_x';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % Ey
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,Ey')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_y';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % Ez
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,Ez')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_z';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end

for ip = 1:npanels
  h(ip).YDir = 'normal';
  h(ip).CLim = 0.5*[-1 1];  
end
colormap(gcf,pic_colors('blue_red'))

hlinks = linkprop(h,{'XLim','YLim'});

%% Check forces on cold ions
nrows = 4; 
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

if 0 % v3x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vx3')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_x 3';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v3y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vy3')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_y 3';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v3z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vz3')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_z 3';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v5x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vx5')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_x 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v5y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vy5')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_y 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v5z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vz5')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_z 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v35x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vx35')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_x 35';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % v35y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vy35')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ic,y}';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v35z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,vz35')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_z 35';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end

if 0 % v3xB.x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v3xB.x')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_x 3';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v3xB.y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v3xB.y')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_y';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v3xB.z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v3xB.z')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_z';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v3xB.x_yz
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v3xB.x_yz')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_yxB_z';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v3xB.x_zy
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v3xB.x_zy')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_zxB_y';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v5xB.x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v5xB.x')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_x 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v5xB.y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v5xB.y')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_y 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v5xB.z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v5xB.z')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_z 5';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v35xB.x
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v35xB.x')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_x 35';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % v35xB.y
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v35xB.y')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = '(v_{ic}xB)_y';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % v35xB.z
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,v35xB.z')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vxB_z';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end

if 0 % Ex
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,Ex')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_x';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % Ey
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,Ey')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_y';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % Ez
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,Ez')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_z';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % Ex
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,Ex')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_x';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 1 % Ey + v35xb
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,Ey'+v35xB.y')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_y+(v_{ic}xB)_y';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end
if 0 % Ez
  hca = h(isub); isub = isub + 1;  
  imagesc(hca,pic.xi,pic.zi,Ez')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_z';
  %hold(hca,'on')
  hca.XLabel.String = 'x/d_i';  
  hca.YLabel.String = 'z/d_i';  
  hca.Box = 'on';  
end


for ip = 1:npanels
  h(ip).YDir = 'normal';
  h(ip).CLim = 0.5*[-1 1];  
end
colormap(gcf,pic_colors('blue_red'))

hlinks = linkprop(h,{'XLim','YLim'});



%% Checking canonical momentum and what it might reveal
x_sel = [183 187 191 197 201 203];
z_sel = x_sel*0;
ds = ds04(2);
twpe = 8000;
%A = df04.twpelim(twpe).A;

iSpecies = [3 5];

idist = 6;
xx = x_sel(idist);
zz = z_sel(idist);

A_c = df04.twpelim(twpe).xlim(xx+[-0.5 0.5]).zlim(zz+[-0.5 0.5]).A;
A_c_ = df04.twpelim(twpe).xlim([min(x_sel) max(x_sel)]+[-0.5 0.5]).zlim(zz+[-0.5 0.5]).A;
A_c_z0 = mean(A_c_,2);
xA_c_ = df04.twpelim(twpe).xlim([min(x_sel) max(x_sel)]+[-0.5 0.5]).xi;

ds = ds04(2).zlim(zz+[-0.1 0.1]).xlim(xx+[-0.1 0.1]);
f = ds.f(1,1,iSpecies);
  
imagesc(f.v,f.v,f.fxy'); colorbar
%%
vx0 = 0.0;
vy0 = 0.5;
Ax0 = 0;
Ay0 = mean(A_c(:));
m = 25;
q = 1;
f_p = @(v,A) m*v + q*A;
f_v = @(p,A) (p - q*A)/m;

px0 = f_p(vx0,Ax0);
py0 = f_p(vy0,Ay0);

plotyy(xA_c_,A_c_z0,xA_c_,f_v(py0,A_c_z0))

%% Scatter plot on Traj
xx = cat(1,tr04.x);
zz = cat(1,tr04.z);
vyy = cat(1,tr04.vy);
vxx = cat(1,tr04.vx);
angle = atan2d(vyy,vxx);
%scatter(xx,zz,1,angle)

nrows = 1; 
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;



if 1 % 
  hca = h(isub); isub = isub + 1;  
  i1 = intersect(find(xx>180),find(xx<200));
  i2 = intersect(find(zz>-1),find(zz<1));
  ind = intersect(i1,i2); 
  scatter(hca,vxx(ind),vyy(ind),2,xx(ind))
  hcb = colorbar('peer',hca);    
  hcb.YLabel.String = 'x/d_i';  
  hca.XLabel.String = 'v_x';  
  hca.YLabel.String = 'v_y';  
  hca.XLim = [-2 2];
  hca.YLim = [-2 2];
  hca.Box = 'on';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % 
  hca = h(isub); isub = isub + 1;  
  tr = tr04.lim('x',[189 190],'z',[-0.5 0.5]);
  %tr = tr04.lim('x',[180 200],'z',[-1 1]);
  %tr = tr04.pass('x',[100 190],'z',[-1 1]);
  % scatter(hca,cat(1,tr.vx),cat(1,tr.vy),2,cat(1,tr.x))
  holdon = 0;
  for iTr = 1:tr.ntr
    %if not(isempty(tr(iTr).t))
      if not(holdon), hold(hca,'on'); end
      scatter(hca,tr(iTr).vx,tr(iTr).vy,2,0*tr(iTr).vx+iTr) 
    %end
  end
  hold(hca,'off')
  hcb = colorbar('peer',hca);    
  hcb.YLabel.String = 'trajectory';  
  hca.XLabel.String = 'v_x';  
  hca.YLabel.String = 'v_y';  
  hca.XLim = [-2 2];
  hca.YLim = [-2 2];
  hca.Box = 'on';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end

%% Check if non-smooth changes in A is associated with the fingers
tic;
pic = df04.xlim([160 240]).zlim([-5 5]).twcilim([120 200]);
%A = pic.A;
tt = pic.twci;
tt = tt(2:end) - mean(diff(tt))/2;
toc;
dA = diff(A,1,3);
levA = -0.0003:0.00003:0.0003;
levA = -25:0.5:0;
doA = 1;
h = setup_subplots(2,1);

for it = 1:size(dA,3)
  isub = 1;
  if 0
    hca = h(isub); isub = isub + 1;
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;  
    %contour(hca,pic.xi(iAx),pic.zi(iAz),dA(iAx,iAz,it)',levA,'k'); 
    imagesc(hca,pic.xi,pic.zi,dA(:,:,it)'); 
    hca.CLim = 1.5*[-1 1];
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '\Delta A';
    colormap(hca,pic_colors('blue_red'))
  end
  if 1
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pic.xi,pic.zi,pic.twcilim(tt(it)).n([3 5])')
    hca.CLim = [0 1];
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'n_{cold ions}';
    colormap(hca,pic_colors('candy'))
    if doA
      hold(hca,'on')
      iAx = 1:4:pic.nx;
      iAz = 1:4:pic.nz;  
      contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz,it)',levA,'k');     
      hold(hca,'off')
    end
  end
  if 1
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pic.xi,pic.zi,pic.twcilim(tt(it)).vy([3 5])')
    hca.CLim = [-0.5 0.5];
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'v_{y,cold ions}';
    colormap(hca,pic_colors('blue_red'))
    if doA
      hold(hca,'on')
      iAx = 1:4:pic.nx;
      iAz = 1:4:pic.nz;  
      contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz,it)',levA,'k');     
      hold(hca,'off')
    end
  end
  
  h(1).Title.String = sprintf('twci = %g',tt(it));
  irf_plot_axis_align(h)
  compact_panels(0.01)
  pause(0.1)
end

%% Vide0 of some ion parameters, n vy, t

%% Plot on top of field, option to make movie, frame of the DF.
doVideo = 0;
doGif = 0;
doPlotTrajectoriesParticles = 0; % lines showing the entire trajectory
doPlotTrajectoriesLines = 0; % lines showing the entire trajectory

particleMarker = '.';
particlerMarkerSize = 12;

movieName = 'n_cold_hot_diff_df04_frameoffront_clim';

twci = [100 210];
doA = 1; Alev = -25:1:0;
pic0 = df04.twcilim(twci);
isC = [3 5];

%[xDF,vDF,aDF,BDF] = df04.xva_df;
%xDF = xDF(1,:);
%tr = trif;
tr = tr04.pass('charge',[0 2]); % all ions
tr = tr04.pass('charge',[0 2]).pass('x0',[120 170]); % all ions
if doVideo
  vidObj = VideoWriter([savedir movieName '.mp4'],'MPEG-4');
  open(vidObj);
end
if doGif
  iframe = 0;
end

nrows = 1;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);

for it = 5%1:pic0.nt
  if 1 % frame of df
    x0 = xDF(pic0.it(it));  
    xlim = x0 + [-70 70];
    xlim = x0 + [-50 50];
  else % simulation frame
    xlim = [120 240];
  end
    
  zlim = [-10 10];
  twci = pic0.twci(it);
  pic = df04.twcilim(twci).xlim(xlim).zlim(zlim);

  
  pc = pic;
  
 
    
  isub = 1;
  if 0 % Ez
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.Ez')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'E_z';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % By
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.By')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'B_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % Ey
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.Ey')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'E_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % Bz
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.Bz')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'B_z';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % magnetic field curvature
    hca = h(isub); isub = isub + 1;
    
    Bx = pc.Bx;
    By = pc.By;
    Bz = pc.Bz;   
    bcurv = magnetic_field_curvature(pc.xi,pc.zi,Bx,By,Bz);
    imagesc(hca,pc.xi,pc.zi,bcurv.abs')
    hca.CLim = [0 2];
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '|B_{curv}|';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
    
  end
  if 0 % v35y
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.vy(isC)')
    hca.CLim = [-.5 0.5];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'v_{y,i,cold}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % v35y
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.n(isC)')
    hca.CLim = [0 1.5];
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'n_{i,cold}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % vxB_y
    hca = h(isub); isub = isub + 1;
    vx = pc.vx(isC);
    vy = pc.vy(isC);
    vz = pc.vz(isC);
    Bx = pc.Bx;
    By = pc.By;
    Bz = pc.Bz;
    vxB = cross_product(vx,vy,vz,Bx,By,Bz);
    imagesc(hca,pc.xi,pc.zi,vxB.y')
    hca.CLim = [-.5 0.5];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '(vxB)_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % n vy T 35
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = njp(pic,isC);
    hca = h(isub); isub = isub + 1;   
    t = (pxx+pyy+pzz)./n/3;
    
    imagesc(hca,pc.xi,pc.zi,t')
    hca.CLim = [0 0.2];
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'T_{i,cold}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % n vy T 1
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = njp(pic,[1]);
    hca = h(isub); isub = isub + 1;   
    t = (pxx+pyy+pzz)./n/3;
    
    imagesc(hca,pc.xi,pc.zi,t')
    hca.CLim = [0 0.2];
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'T_{i,hot}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % nP 35
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = njp(pic,isC);
    hca = h(isub); isub = isub + 1;   
    p = (pxx+pyy+pzz)/3;
    
    imagesc(hca,pc.xi,pc.zi,p')
    hca.CLim = [-0.4 0.4];
    colormap(hca,pic_colors('candy2'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'P_{i,cold}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % P 1 - P0
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = njp(pic,[1]);
    hca = h(isub); isub = isub + 1;   
    p = (pxx+pyy+pzz)/3;
    p0 = 0.5*5/6/5;
    imagesc(hca,pc.xi,pc.zi,p'-p0)
    hca.CLim = [-0.4 0.4];
    colormap(hca,pic_colors('candy2'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'P_{i,hot}-P_{i,hot,lobe}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end 
  if 0 % n 35
    n = pic.n(isC);
    hca = h(isub); isub = isub + 1;       
    
    imagesc(hca,pc.xi,pc.zi,n'/0.4)
    hca.CLim = [0 6];
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'n_{i,cold}/n_{i,cold,lobe}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % n
    n = pic.n(1);
    hca = h(isub); isub = isub + 1; 
    imagesc(hca,pc.xi,pc.zi,n'/0.2)
    hca.CLim = [0 6];
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'n_{i,hot}/n_{i,hot,lobe}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end 
  if 1 % nhot-ncold
    nhot = pic.n(1);
    ncold = pic.n(isC);
    hca = h(isub); isub = isub + 1; 
    imagesc(hca,pc.xi,pc.zi,nhot'-ncold')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'n_{i,hot}-n_{i,cold}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end 
  
  
  h(1).Title.String = sprintf('twci = %g',pc.twci);
  drawnow
  compact_panels(0.02)
  for ip = 1:npanels
    hca = h(ip);
    hold(hca,'on')
    hca.FontSize = 14;
    hca.YDir ='normal';
    if doA
      iAx = 1:4:pic.nx;
      iAz = 1:4:pic.nz;
      A = pc.A;
      contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',Alev,'k');        
    end    
    for itr = 1:numel(tr)
      if doPlotTrajectoriesLines
        plot(hca,tr(itr).x,tr(itr).z)
      end
      if doPlotTrajectoriesParticles
      ii = find(abs(tr(itr).t-pc.twci)==min(abs(tr(itr).t-pc.twci)));
      plot(hca,tr(itr).x(ii),tr(itr).z(ii),'color',[0 0 0],'Marker',particleMarker,'MarkerSize',particlerMarkerSize)
      end
    end
    hold(hca,'off')
  end
  pause(1)
  if doVideo
    set(gcf,'color','white');
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = pic0.nt;
      currentBackgroundColor = get(gcf,'color');
      set(gcf,'color',[1 1 1]);
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
if doVideo, close(vidObj); end
if doGif
  imwrite(all_im,map,[savedir movieName '.gif'],'DelayTime',0.0,'LoopCount',inf)
end

%% Figure with forces of example ion that makes a secondary finger
% Select ions that corresponds to certain criteria
% location, speed range, speed angle

%tr1 = tr04.pass('x',[189 190],'z',[-0.25 0.25]);
%tr = tr04.pass('x',[189 190],'z',[-0.25 0.25],'atan2d(vy,vx)',[90 135]);
%tr = tr04.pass('x',[189 190],'z',[-0.25 0.25]).pass('atan2d(vy,vx)',[90 135]);
trif = tr04.pass('x',[189 190],'z',[-0.25 0.25],'atan2d(vy,vx)',[90 135]);
%% % plot of "ifingers"
for itr = 1:numel(trif)
  tr = trif(itr);
  h = setup_subplots(5,1,'vertical');
  isub = 1;
  cmap = pic_colors('waterfall');
  t0_msize = 10;
  tr = tr(1);
  legloc = 'eastoutside';
  %cmap = interp1(linspace(1,64,size(cmap,1)),cmap,1:numel(tr.t)); 
  if 0 % (vx,vy), colorcoded by time
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
  end
  if 0 % (vx,vz), colorcoded by time
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
  end
  if 0 % (vy,vz), colorcoded by time
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
  if 0 % xyz
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
  if 1 % xyz(t)
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t,tr.x-200,tr.t,tr.y,tr.t,tr.z)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'r';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    legend(hca,{'x-200','y','z'},'location',legloc)
  end
  if 1 % v 
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t,tr.vx,tr.t,tr.vy,tr.t,tr.vz)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'v';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    legend(hca,{'x','y','z'},'location',legloc)
  end
  if 1 % B
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t,tr.Bx,tr.t,tr.By,tr.t,tr.Bz)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'B';
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    legend(hca,{'x','y','z'},'location',legloc)
  end
  if 1 % E 
    hca = h(isub); isub = isub + 1;
    plot(hca,tr.t,tr.Ex,tr.t,tr.Ey,tr.t,tr.Ez)
    hca.XLabel.String = 't';
    hca.YLabel.String = 'E';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';  
    legend(hca,{'x','y','z'},'location',legloc)
  end
  if 1 % vxB
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

  c_eval('h(?).XLim = [100 210];',1:5)
  compact_panels(0.01)
  h(2).YLim = [-1.7 1.7];
  h(3).YLim = [-1.02 1.02];
  h(3).YLim = [-0.99 0.99];
  h(4).YLim = [-0.49 0.49];
  h(5).YLim = [-0.49 0.49];
  %h(6).YLim = [-0.49 0.49];
  for ip = 1:5
    h(ip).FontSize = 14;
    h(ip).Position(3) = 0.7;
  end
  for ip = []%1:3
    h(ip).XLim = [-2 2];
    h(ip).YLim = [-2 2];
    axis(h(ip),'square')
    h(ip).XTick = -15:0.5:15;
    h(ip).YTick = -15:0.5:15;
  end
  cn.print(sprintf('ifinger_forces_%g',itr),'path',[savedir '/forces/'])
end

%% Plot distributions that correponds to 'DF-frame' movie
ds = ds04(2).zlim([-0.2 0.2]).xlim([166 169]+[-0.2 0.2]);

nrows = 1;
ncols = 4;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

for id = 1:ds.nd{1}
  hca = h(isub); isub = isub + 1;
  f = ds.f(1,id,[3 5]);
  imagesc(hca,f.v,f.v,f.fxz')
  hca.XLim = [-1.5 1.5];
  hca.YLim = [-1.5 1.5];
  axis(hca,'square')
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
  colormap(hca,pic_colors('candy'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XTick = -5:1:5;
  hca.YTick = -5:1:5;
  hca.Position(2) = 0.15;
  irf_legend(hca,sprintf('x = %.0f, z = %.0f',mean(f.x),mean(f.z)),[0.02 0.98]);
end
compact_panels(0,0.01)
c_eval('h(?).YLabel.String = '''';',2:4)
c_eval('h(?).YTickLabel = '''';',2:4)

%% dEF of ions along a line, foremost between df and xline
twci = 100; it = 1;
zlim = 1 + [-0.24 0.24];
xlim = [155 200];

twci = 160; it = 2;
zlim = 7 + [-0.24 0.24];
xlim = [155 200];
nE = 50;

ds = ds04(it).xlim(xlim).zlim(zlim);
pic = df04.twcilim(twci).xlim(xlim).zlim(zlim);
dEF35 = ds.dEF(1,ds.indices{1},[3 5],nE);
dEF1 = ds.dEF(1,ds.indices{1},[1],nE);
% get ExB velocity to overplot spectrogram
Ex = pic.Ex; 
Ey = pic.Ey;
Ez = pic.Ez;
Bx = pic.Bx;
By = pic.By;
Bz = pic.Bz;
Babs = sqrt(Bx.^2+By.^2+Bz.^2);

vExB = cross_product(Ex,Ey,Ez,Bx./Babs,By./Babs,Bz./Babs);
if 1 % mean first
  vExB.x = mean(vExB.x,2);
  vExB.y = mean(vExB.y,2);
  vExB.z = mean(vExB.z,2);
else
  
end
eExB = (vExB.x.^2+vExB.y.^2+vExB.z.^2)/2;


% plot
doExB = 1;

nrows = 5;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
cmap = pic_colors('candy');
cmap = colormap('jet');
%cmap = irf_colormap('waterfall');


if 1 % B
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.xi,-mean(pic.Bz,2)');
  shading(hca,'flat')
  colormap(cmap)
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = '-B_z/B_0';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % n
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.xi,mean(pic.n([1]),2)',...
    pic.xi,mean(pic.n([3 5]),2)',...
    pic.xi,mean(pic.n([1]),2)'+mean(pic.n([3 5]),2)');
  shading(hca,'flat')
  colormap(cmap)
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'n/n_0';
  irf_legend(hca,{'n_{hot ions}';'n_{cold ions}';'n_{all ions}'},[0.98 0.98])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % eExB
  hca = h(isub); isub = isub + 1;
  semilogy(hca,pic.xi,mean(eExB,2)');
  shading(hca,'flat')
  colormap(cmap)
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'E_{ExB}';
  irf_legend(hca,{'n_{hot ions}';'n_{cold ions}';'n_{all ions}'},[0.98 0.98])
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,dEF1.x,log10(dEF1.energy),log10(dEF1.dEF)'); 
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'diff. energy flux';
  colormap(cmap)
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'log_{10}E_{hot ions}';
end
if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,dEF35.x,log10(dEF35.energy),log10(dEF35.dEF)'); 
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'diff. energy flux';
  colormap(cmap)
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'log_{10}E_{cold ions}';
  if doExB
    hold(hca,'on')
    plot(hca,pic.xi,smooth(log10(mean(eExB,2)),7)','k');
    hold(hca,'on')
    irf_legend(hca,{'E_{ExB}'},[0.02 0.1],'color',[0 0 0])
  end
end
if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,dEF1.x,log10(dEF1.energy),log10(dEF1.dEF+dEF35.dEF)'); 
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'diff. energy flux';
  colormap(cmap)
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 'log_{10}E_{all ions}';
  if doExB
    hold(hca,'on')
    plot(hca,pic.xi,smooth(log10(mean(eExB,2)),7)','k');
    hold(hca,'on')
    irf_legend(hca,{'E_{ExB}'},[0.02 0.1],'color',[0 0 0])
  end
end

h(1).Title.String =sprintf('twci = %g, z = %g d_i',twci,mean(zlim));
hlinks = linkprop(h(3:5),{'CLim'});
hlinks.Targets(1).CLim = [-4 2];
compact_panels(0.01)
irf_plot_axis_align(h)

for ip = 1:npanels
  h(ip).XLim = xlim;
  
end

%% Forces across front option to make video and gif
doVideo = 1;
doGif = 1;
doPlotTrajectories = 0; % lines showing the entire trajectory
particleMarker = '.';
particlerMarkerSize = 15;

movieName = 'forces_at_front_1';
twci = [100 210];
doA = 1; Alev = -25:1:0;
pic0 = df04.twcilim(twci);

%[xDF,vDF,aDF,BDF] = df04.xva_df;
%xDF = xDF(1,:);
%tr = trif;
%tr = tr04.pass('charge',[0 2]); % all ions
if doVideo
  vidObj = VideoWriter([savedir movieName '.mp4'],'MPEG-4');
  open(vidObj);
end
if doGif
  iframe = 0;
end

nrows = 6;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);

for it = 1:pic0.nt
  if 0
    x0 = xDF(pic0.it(it));  
    xlim = x0 + [-70 70];
  else
    xlim = [120 210];
  end
    
  zlim = [-1 1]; % mean over this distance
  twci = pic0.twci(it);
  pic = df04.twcilim(twci).xlim(xlim).zlim(zlim);
  
  pc = pic;
  % load everything thats needed
  [n1,jx1,jy1,jz1,pxx1,pxy1,pxz1,pyy1,pyz1,pzz1] = pic.njp([1]);
  [n35,jx35,jy35,jz35,pxx35,pxy35,pxz35,pyy35,pyz35,pzz35] = pic.njp([3 5]);
  
  n1 = mean(n1,2);
  n35 = mean(n35,2);
  jx1 = mean(jx1,2);
  jy1 = mean(jy1,2);
  jz1 = mean(jz1,2);
  jx35 = mean(jx35,2);
  jy35 = mean(jy35,2);
  jz35 = mean(jz35,2);
  vx1 = jx1./n1;
  vy1 = jy1./n1;
  vz1 = jz1./n1;
  vx35 = jx35./n35;
  vy35 = jy35./n35;
  vz35 = jz35./n35;
  
  ns = 3;
  P1.xx = smooth2(pxx1,ns);
  P1.xy = smooth2(pxy1,ns);
  P1.xz = smooth2(pxz1,ns);
  P1.yy = smooth2(pyy1,ns);
  P1.yz = smooth2(pyz1,ns);
  P1.zz = smooth2(pzz1,ns);
  P35.xx = smooth2(pxx35,ns);
  P35.xy = smooth2(pxy35,ns);
  P35.xz = smooth2(pxz35,ns);
  P35.yy = smooth2(pyy35,ns);
  P35.yz = smooth2(pyz35,ns);
  P35.zz = smooth2(pzz35,ns);
  
  divP1 = div_tensor(pic.xi, pic.zi, P1,'comp',1);
  divP35 = div_tensor(pic.xi, pic.zi, P35,'comp',1);
  divP1xB = cross_product(divP1.x,divP1.y,divP1.z,Bx,By,Bz);
  Ex = mean(pic.Ex,2);
  Ey = mean(pic.Ey,2);
  Ez = mean(pic.Ez,2);
  Bx = mean(pic.Bx,2);
  By = mean(pic.By,2);
  Bz = mean(pic.Bz,2);  
  Babs = sqrt(Bx.^2+By.^2+Bz.^2);    
  vxB1 = cross_product(vx1,vy1,vz1,Bx,By,Bz);
  vxB35 = cross_product(vx35,vy35,vz35,Bx,By,Bz);
 
    
  isub = 1;
  if 1 % B_x,E_y
    hca = h(isub); isub = isub + 1;
    plot(hca,pc.xi,abs(Bz'),pc.xi,Ey')        
    hca.YLabel.String = '|B_z|,Ey';
    hca.XLabel.String = 'x/d_i';
  end
  if 1 % divP
    hca = h(isub); isub = isub + 1;
    colors = pic_colors('matlab');
    plot(hca,pc.xi,-mean(divP1.x,2),pc.xi,-mean(divP1.y,2),pc.xi,-mean(divP1.y,2))
    hold(hca,'on')
    set(hca,'colororder',colors)
    plot(hca,pc.xi,-mean(divP35.x,2),'--',pc.xi,-mean(divP35.y,2),'--',pc.xi,-mean(divP35.y,2),'--')
    hold(hca,'off')
    
    hca.YLabel.String = 'divP';
    hca.XLabel.String = 'x/d_i';      
  end
  if 1 % divP1_x_components
    hca = h(isub); isub = isub + 1;
    colors = pic_colors('matlab');
    plot(hca,pc.xi,-mean(divP1.x_xx,2),pc.xi,-mean(divP1.x_yy,2),pc.xi,-mean(divP1.x_zz,2))    
    hca.YLabel.String = 'divP';
    hca.XLabel.String = 'x/d_i';    
    legend(hca,{'dxPxx','dyPxy','dzPxz'},'location','northwest')
  end
  if 1 % divP1_y_components
    hca = h(isub); isub = isub + 1;
    colors = pic_colors('matlab');
    plot(hca,pc.xi,-mean(divP1.y_xx,2),pc.xi,-mean(divP1.x_yy,2),pc.xi,-mean(divP1.x_zz,2))    
    hca.YLabel.String = 'divP';
    hca.XLabel.String = 'x/d_i';    
    legend(hca,{'dxPxy','dyPyy','dzPyz'},'location','northwest')
  end
  if 1 % E, vxB,divP, y
    hca = h(isub); isub = isub + 1;
    colors = pic_colors('matlab');
    plot(hca,pc.xi,Ey',pc.xi,vxB1.y',pc.xi,-mean(divP1.y,2)',...
             pc.xi,Ey'+vxB1.y'-mean(divP1.y,2)')
    hold(hca,'on')
    set(hca,'colororder',colors)
    %plot(hca,pc.xi,-mean(divP1.x,2),'--',pc.xi,-mean(divP1.y,2),'--',pc.xi,-mean(divP1.y,2),'--')
    hold(hca,'off')
    
    hca.YLabel.String = '...';
    hca.XLabel.String = 'x/d_i';  
    legend(hca,{'E_y','vxB_y','-divP_y','E_y+vxB_y-divP_y'},'location','northwest')
  end
  if 1 % E, vxB,divP, x
    hca = h(isub); isub = isub + 1;
    colors = pic_colors('matlab');
    plot(hca,pc.xi,Ex',pc.xi,vxB1.x',pc.xi,-mean(divP1.x,2)',...
             pc.xi,Ex'+vxB1.x'-mean(divP1.x,2)')
    hold(hca,'on')
    set(hca,'colororder',colors)
    %plot(hca,pc.xi,-mean(divP1.x,2),'--',pc.xi,-mean(divP1.y,2),'--',pc.xi,-mean(divP1.y,2),'--')
    hold(hca,'off')
    
    hca.YLabel.String = '...';
    hca.XLabel.String = 'x/d_i';  
    legend(hca,{'E_x','vxB_x','-divP_x','E_y+vxB_y-divP_x'},'location','northwest')
  end
  if 0 % n(Ey+vxB)_x
    hca = h(isub); isub = isub + 1;
    plot(hca,pc.xi,n1'.*(Ey+vxB1.x)',...
             pc.xi,n35'.*(Ey+vxB35.x)')    
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'n(E+vxB)_x';
    hca.XLabel.String = 'x/d_i';  
    legend(hca,{'hot ions: n(E+vxB)','cold ions: n(E+vxB)'},'location','northwest')
  end
  if 0 % n(Ey+vxB)_y
    hca = h(isub); isub = isub + 1;
    plot(hca,pc.xi,n1'.*(Ey+vxB1.y)',...
             pc.xi,n35'.*(Ey+vxB35.y)')    
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'n(E+vxB)_y';
    hca.XLabel.String = 'x/d_i';  
    legend(hca,{'hot ions: n(E+vxB)','cold ions: n(E+vxB)'},'location','northwest')
  end
  if 0 % -divP.y
    hca = h(isub); isub = isub + 1;
    plot(hca,pc.xi,-smooth(mean(divP1.y,2),10)',...
             pc.xi,-smooth(mean(divP35.y,2),10)')        
    hca.YLabel.String = '-((divP)_y';
    hca.XLabel.String = 'x/d_i';  
    legend(hca,{'hot ions: -divP','cold ions: -divP'},'location','northwest')
  end
  if 0 % n(Ey+vxB)-divP.y
    hca = h(isub); isub = isub + 1;
    plot(hca,pc.xi,n1'.*(Ey+vxB1.y)'-smooth(mean(divP1.y,2),10)',...
             pc.xi,n35'.*(Ey+vxB35.y)'-smooth(mean(divP35.y,2),10)')    
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'n(E+vxB)_y-((divP)_y';
    hca.XLabel.String = 'x/d_i';  
    legend(hca,{'hot ions: n(E+vxB)-divP','cold ions: n(E+vxB)-divP'},'location','northwest')
  end
  if 0 % neEy,divP.y
    hca = h(isub); isub = isub + 1;
    plot(hca,pc.xi,n1'.*Ey',...
             pc.xi,n35'.*Ey',...
             pc.xi,n1'.*vxB1.y',...
             pc.xi,n35'.*vxB35.y',...
             pc.xi,smooth(mean(divP1.y,2),10)',...
             pc.xi,mean(divP35.y,2)')    
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'nE_y';
    hca.XLabel.String = 'x/d_i';    
  end
  if 0 % Ey,divP.y
    hca = h(isub); isub = isub + 1;
    plot(hca,pc.xi,Ey',...
      pc.xi,mean(divP35.y,2)',...
      pc.xi,mean(divP1.y,2)')    
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'E_y';
    hca.XLabel.String = 'x/d_i';    
  end
  if 0 % By
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.By')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'B_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % Ey
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.Ey')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'E_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % Bz
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.Bz')
    hca.CLim = [-1 1];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'B_z';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % magnetic field curvature
    hca = h(isub); isub = isub + 1;
    
    Bx = pc.Bx;
    By = pc.By;
    Bz = pc.Bz;   
    bcurv = magnetic_field_curvature(pc.xi,pc.zi,Bx,By,Bz);
    imagesc(hca,pc.xi,pc.zi,bcurv.abs')
    hca.CLim = [0 2];
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '|B_{curv}|';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
    
  end
  if 0 % v35y
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.vy([3 5])')
    hca.CLim = [-.5 0.5];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'v_{y,i,cold}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % v35y
    hca = h(isub); isub = isub + 1;
    imagesc(hca,pc.xi,pc.zi,pc.n([3 5])')
    hca.CLim = [0 1.5];
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'n_{i,cold}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % vxB_y
    hca = h(isub); isub = isub + 1;
    vx = pc.vx([3 5]);
    vy = pc.vy([3 5]);
    vz = pc.vz([3 5]);
    Bx = pc.Bx;
    By = pc.By;
    Bz = pc.Bz;
    vxB = cross_product(vx,vy,vz,Bx,By,Bz);
    imagesc(hca,pc.xi,pc.zi,vxB.y')
    hca.CLim = [-.5 0.5];
    colormap(hca,pic_colors('blue_red'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '(vxB)_y';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end
  if 0 % n vy T
    [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = njp(pic,[3 5]);
    hca = h(isub); isub = isub + 1;   
    t = (pxx+pyy+pzz)./n/3;
    
    imagesc(hca,pc.xi,pc.zi,t')
    hca.CLim = [0 0.2];
    colormap(hca,pic_colors('candy'))
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'T_{i,cold}';
    hca.XLabel.String = 'x/d_i';
    hca.YLabel.String = 'z/d_i';
  end 
  
  h(1).Title.String = sprintf('twci = %g',pc.twci);
  drawnow
  compact_panels(0.02)
  compact_panels(0.01)
  for ip = 1:npanels
    hca = h(ip);
    hca.FontSize = 14;
    hca.YDir ='normal';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  irf_plot_axis_align
  pause(1)
  if doVideo
    set(gcf,'color','white');
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  if doGif
    if 1 % collect frames, for making gif
      iframe = iframe + 1;    
      nframes = pic0.nt;
      currentBackgroundColor = get(gcf,'color');
      set(gcf,'color',[1 1 1]);
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
if doVideo, close(vidObj); end
if doGif
  imwrite(all_im,map,[savedir movieName '.gif'],'DelayTime',0.0,'LoopCount',inf)
end
