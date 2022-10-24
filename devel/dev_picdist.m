h5filePath = '/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists.h5';
no02m = PIC(h5filepath);
ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists.h5');
%%
ds100.xlim([69 71]).plot_map(2,1,'fieldaligned_cs',no02m);