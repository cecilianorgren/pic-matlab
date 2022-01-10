% AGU2021
%% Load SMILEI object
filepath = '/Users/cecilia/Data/PIC/Smilei/AGU2/Fields0.h5';
namelist = '/Users/cecilia/Data/PIC/Smilei/AGU2/Diamagnetic.py';
particlebinningpath = '/Users/cecilia/Data/PIC/Smilei/AGU2/';

sm = SMILEI(filepath,namelist,particlebinningpath);

%% Get reconnection rate
clear X Y ER AX
pic = sm; %.ylim([10 18])
for it  = 1:pic.nt
  Atmp = pic(it).A;
  [inds,vals] = saddle(-Atmp,'sort','plot');
  drawnow
  
  if isempty(inds)
    ER(it) = NaN;
    continue
  else
    Etmp = pic(it).Ez;
    xx = pic.xi(inds(1,1));
    yy = pic.yi(inds(1,2));
    
    X(it) = xx;
    Y(it) = yy;
    ER(it) = Etmp(inds(1,1),inds(1,2));
    AX(it) = Atmp(inds(1,1),inds(1,2));
    
    if 0
      imagesc(pic.xi,pic.yi,Atmp')
      hold on
      plot(xx,yy,'x')
      hold off
      drawnow
    end
  end
end

%% Get alfven speed one di below the xline
t = pic.twci;
for it  = 1:pic.nt
  xlim = X(it) + [-0.5 0.5];
  ylim_bot = Y(it) - 1 + [-0.5 0.5];
  ylim_top = Y(it) + 1 + [-0.5 0.5];
  Bbot(it) = mean(mean(pic(it).xlim(xlim).ylim(ylim_bot).Babs));
  Btop(it) = mean(mean(pic(it).xlim(xlim).ylim(ylim_top).Babs));
  
  Bxbot(it) = mean(mean(pic(it).xlim(xlim).ylim(ylim_bot).Bx));
  Bxtop(it) = mean(mean(pic(it).xlim(xlim).ylim(ylim_top).Bx));
  
  Bybot(it) = mean(mean(pic(it).xlim(xlim).ylim(ylim_bot).By));
  Bytop(it) = mean(mean(pic(it).xlim(xlim).ylim(ylim_top).By));
  
  Bzbot(it) = mean(mean(pic(it).xlim(xlim).ylim(ylim_bot).Bz));
  Bztop(it) = mean(mean(pic(it).xlim(xlim).ylim(ylim_top).Bz));
  
  nbot(it) = mean(mean(pic(it).xlim(xlim).ylim(ylim_bot).ni));
  ntop(it) = mean(mean(pic(it).xlim(xlim).ylim(ylim_top).ni));
  
  shearbot(it) = atan2d(Bybot(it),Bxbot(it));
  sheartop(it) = atan2d(Bytop(it),Bxtop(it));
  
  vAbot(it) = Bxbot(it)./sqrt(nbot(it));
  vAtop(it) = Bxtop(it)./sqrt(ntop(it));
end
%%
hca = subplot(2,1,1);
RR = -smooth(ER,9);
RR(RR<0) = 0;

hca = subplot(2,1,2);
plot(hca,t,RR)

hca = subplot(2,1,1);
plot(hca,t,Bxbot,t,Bxtop,t,nbot,t,ntop,t,vAbot,t,vAtop)
legend(hca,{'n_{bot}','n_{top}','B_{x,bot}','B_{x,top}','v_{A,bot}','v_{A,top}'},'location','best')
%% Plot reconnection rate
hca = subplot(1,1,1);
plot(hca,pic.twci,-smooth(ER,9))
%plot(hca,pic.twci,-smooth(ER,2))
hca.XLabel.String = 't\omega_{ci}';
hca.YLabel.String = 'Reconnection rate';
hca.YLim(1) = 0;

%% Plot reconnection rate + xline position
XX = X;
XX(75) = 0.5*(X(73)+X(77));
XX(76) = 0.5*(X(73)+X(77));
RR = -smooth(ER,9);
RR(RR<0) = 0;
hca = subplot(1,1,1);
[h,l1,l2] = plotyy(hca,pic.twci,RR,pic.twci,XX);
%plot(hca,pic.twci,-smooth(ER,2))
hca.XLabel.String = 't\omega_{ci}';
hca.YLabel.String = 'Reconnection rate';

h(2).YLabel.String = 'X line location';

%% Make movie of density
pic = sm(1:2:sm.nt).twcilim([20 100]);
h = pic.movie({'ni'}','A',0.5,'cmap',{pic_colors('thermal')},'clim',{[0 2]});


%% Plot of Jx Jy
pic = sm.twcilim(60);

cmaps = {pic_colors('blue_red'),pic_colors('blue_red')};
clims = {[-1 1],[-1 1]};
h = pic.plot_map({'Jx','Jz'}','A',0.5,'cmap',cmaps,'clim',clims);

%% Plot of vex vExBx
pic = sm.twcilim(70);

cmaps = {pic_colors('blue_red'),pic_colors('blue_red')};
clims = {[-2 2],[-2 2]};
h = pic.plot_map({'vex','vExBx'}','A',0.5,'cmap',cmaps,'clim',clims);

%% Make movie of Jx
pic = sm;
filename = [printpath 'diam_Jx'];

h = pic.movie({'Jx'}','A',0.5,'cmap',{pic_colors('blue_red')},'clim',{[-1 1]},'filename',filename);

%% Make movie of vex
pic = sm;
filename = [printpath 'diam_vex'];

h = pic.movie({'vex'}','A',0.5,'cmap',{pic_colors('blue_red')},'clim',{[-2 2]},'filename',filename,'smooth',2);

