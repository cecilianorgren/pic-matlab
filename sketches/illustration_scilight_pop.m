%no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

twpe = 24000;
xlim = [60 140];
zlim = [-10 10];

pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);

cmap = pic_colors('blue_white'); cmap = flipdim(cmap,1);

pic.plot_map({'te'},'cmap',{cmap},'sep')

%% Sphere/ellipsoid with potential

[X,Y,Z] = ellipsoid(0,0,0,2,2,1,100);

surf(X,Y,Z)
axis equal