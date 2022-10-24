h5filepath = '/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/fields.h5';
datapath = '/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data/';
pic = PIC(h5filepath); % If you have many times saved, this can take up to a minute


%%



var1 = pic.twpelim(23000).xlim(xlim).zlim(zlim).(varstr);
var2 = pic.twpelim(24000).xlim(xlim).zlim(zlim).(varstr);

%%

hca = subplot(1,1,1);

xlim = [98 115];
zlim = [-3 3];
varstr = 'Epar'; cmap = pic_colors('blue_red'); clim = 1*[-1 1];
%varstr = 'Jy'; cmap = pic_colors('blue_red'); clim = 5*[-1 1];
varstr = 'te'; cmap = flipdim(pic_colors('thermal'),1); clim = [0 0.3];
%varstr = 'log10(tepar./tperp)'; cmap = flipdim(pic_colors('thermal'),1); clim = [-1 1];

times = 23000:200:23600;
var = pic.twpelim(times,'exact').xlim(xlim).zlim(zlim).(varstr);
A = pic.twpelim(times,'exact').xlim(xlim).zlim(zlim).A;

yloc = linspace(0,50,numel(times));
Alev = -25:0.2:25;

holdon = 0;
for iy = 1:numel(yloc)
  pic_tmp = pic.twpelim(times(iy),'exact').xlim(xlim).zlim(zlim);
  [X,Z] = ndgrid(pic_tmp.xi,pic_tmp.zi);
  Y = X*0+yloc(iy);
  VAR = squeeze(var(:,:,iy));
  hs = surf(hca,X,Y,Z,VAR);
  hs.FaceAlpha = 0.9;
  shading(hca,'flat')
  
  
  if not(holdon)
    hold(hca,'on')    
  end
  
  S = contourcs(pic_tmp.xi,pic_tmp.zi,squeeze(A(:,:,iy))',Alev);
  for iline = 1:numel(S)
    xx = S(iline).X;
    yy = S(iline).X*0+yloc(iy);
    zz = S(iline).Y;
    plot3(hca,xx,yy,zz,'k')
  end
  
  %contour3(hca,X,Y,Z,A,Alev)
  
end

colormap(hca,cmap)
hca.CLim = clim;
axis(hca,'equal')

% Add out if plane wave structure 
dy = yloc(1)-yloc(2);
ly = dy;
ky = 2*pi/ly;
yvec = linspace(yloc(1)-10,yloc(end),100);
zwave = sin(ky*yvec);

xvec = yvec*0 + 105;
%plot3(hca,xvec,yvec,zwave,'k')

xpatchdata = [xvec xvec + 3];
ypatchdata = [yvec yvec(end:-1:1)];
zpatchdata = [zwave zwave(end:-1:1)];

plot3(hca,xvec,yvec,zwave,'k')
hpatch = patch(hca,xpatchdata,ypatchdata,zpatchdata,'k');
hpatch.FaceAlpha = 0.4;
hpatch.FaceColor = [0.7 0 0];


view(hca,[-1 -1.5 0.3])
hold(hca,'off')

hca.Visible = 'off';
fig = gcf;
fig.Color = [1 1 1];


%%
var = pic.twpelim([20000:1000:24000],'exact').xlim(xlim).zlim(zlim).(varstr);

hf2 = figure ;
hs = slice(permute(var,[3 2 1]),[],[],1:4);
shading interp
set(hs,'FaceAlpha',0.8);