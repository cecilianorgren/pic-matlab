filepath = '/Users/cno062/tesla/software/Smilei4.5/Smilei/cold_run/Fields0.h5';
namelist = '/Users/cno062/tesla/software/Smilei4.5/Smilei/cold_run/Harris_new.py';
sm = PIC(filepath,namelist);

%% Figure: Diagnostic 1
pic = sm(sm.nt);
xlim = pic.xi([1 end]);
zlim = pic.xi([1 end]);

varstrs = {'n(1)','n(3)','n(5)','Bx','Bz','Jy','vex','vix'}';
clims = {[0 1.2],[0 1.2],[0 1.2],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1]};
cmap_th = pic_colors('thermal');
cmap_pst = pic_colors('pasteljet');
cmap_br = pic_colors('blue_red');

cmaps = {cmap_th,cmap_th,cmap_th,cmap_br,...
  cmap_br,cmap_br,cmap_br,cmap_br,cmap_br};

pic.xlim(xlim).zlim(zlim).plot_map(varstrs,'clim',clims,'cmap',cmaps)

%% Figure: Diagnostic 2
pic = sm(sm.nt);
xlim = pic.xi([1 end]);
zlim = pic.xi([1 end]);

varstrs = {'n(1)','n(3)','Ex','Ey','Ez'}';
clims = {[0 1.2],[0 1.2],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1]};
cmap_th = pic_colors('thermal');
cmap_pst = pic_colors('pasteljet');
cmap_br = pic_colors('blue_red');

cmaps = {cmap_th,cmap_th,cmap_br,...
  cmap_br,cmap_br,cmap_br,cmap_br,cmap_br};

pic.xlim(xlim).zlim(zlim).plot_map(varstrs,'clim',clims,'cmap',cmaps)

%% Figure: Diagnostic 3, quick
pic = sm(sm.nt);
xlim = pic.xi([1 end]);
zlim = pic.xi([1 end]);

varstrs = {'n(1)','n(3)','Bx','Bz','Jy','vex','vix'}';
clims = {[0 1.2],[0 1.2],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1]};
cmap_th = pic_colors('thermal');
cmap_pst = pic_colors('pasteljet');
cmap_br = pic_colors('blue_red');

cmaps = {cmap_th,cmap_th,cmap_br,...
  cmap_br,cmap_br,cmap_br,cmap_br,cmap_br};

pic.xlim(xlim).zlim(zlim).plot_map(varstrs,'clim',clims,'cmap',cmaps)
