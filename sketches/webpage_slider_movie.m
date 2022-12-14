%% Load PIC object
no02m = PIC('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/fields.h5');

%% Make movie


twpe = [15000:200:25000];
xlim = no02m.xi([1 end])+[60 -60]';
zlim = [-10 10];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);

cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmapjet = colormap('jet');
cmapth = pic_colors('thermal');
cmapth_flip = flipdim(pic_colors('thermal'),1);

varstrs = {'te'}';
clims = {[0 0.25],[-1 1]};
cmaps = {cmapth};
cbarlabels = {'Electron temperature'};
color_A = [0.4 0.4 0.4];

filename = [printpath 'slider_Te'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'cbarlabels',cbarlabels,'filename',filename,'fill','colA',color_A,'figpos',[0 0 1200 400]);

%pic.twpelim([17000 25000]).movie({'Ez'},'A',1,'clim',{[-1 1]},'cmap',{pic_colors('blue_red')},'filename',[printpath 'no02m_Ez']);
