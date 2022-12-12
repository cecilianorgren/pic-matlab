% 
ons = PIC('/Volumes/DataRaid/Susanne-onset/data_h5/fields.h5');
%%

pic = ons.twcilim([0 65]).xlim([1 25]).zlim([-2 2]);
%pic = ons.twcilim([65]).xlim([1 25]).zlim([-2 2]);

clims = {2*[-1 1],0.2*[0 1]};
varstrs = {'vex','te'}';
cmaps = {pic_colors('blue_red'),flipdim(pic_colors('thermal'),1),pic_colors('blue_red'),pic_colors('blue_red')};

clims = {1*[-1 1]*0.99,3*[-1 1]*0.99};
varstrs = {'vex','Jy'}';
cmaps = {pic_colors('blue_red'),pic_colors('blue_red')};
cbarlabels = {'v_{ex}','J_y'};
pic.movie(varstrs,'clim',clims,'cmap',cmaps,'cbarlabels',cbarlabels,'A',0.2,'smooth',5,'filename',[printpath 'ons_onset2_longer']) 


%%

pic = no02m.twpelim([15000 25000]).xlim(mean(no02m.xi)+50*[-1 1]).zlim([-8 8]);
%pic = ons.twcilim([65]).xlim([1 25]).zlim([-2 2]);


clims = {1*[-1 1]*0.99,3*[-1 1]*0.99};
varstrs = {'vex','Jy'}';
cbarlabels = {'v_{ex}','J_y'};
cmaps = {pic_colors('blue_red'),pic_colors('blue_red')};

clims = {3*[-1 1]*0.99,0.2*[0 1]*0.99};
varstrs = {'vex','te'}';
cbarlabels = {{'Electron outflow','speed v_{ex}'},'Electron temperature'};
cmaps = {pic_colors('blue_red'),flipdim(pic_colors('thermal'),1),pic_colors('blue_red'),pic_colors('blue_red')};



pic.movie(varstrs,'clim',clims,'cmap',cmaps,'cbarlabels',cbarlabels,'A',0.5,'smooth',5,'filename',[printpath 'no02m_2']) 

