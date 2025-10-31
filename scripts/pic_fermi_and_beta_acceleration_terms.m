% no02m = PIC('/Users/cecilianorgren/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
twpe = 20000;
pic = no02m.tlim(twpe).xlim([60 145]).zlim([-10 10]);

%curv = pic.magnetic_curvature;


varstrs = {'curvbx','curvby','curvbz','vExBx','vExBy','vExBz'}';
clims = {[-1 1],[-1 1],[-1 1];[-1 1],[-1 1],[-1 1];[-1 1],[-1 1],[-1 1]}';


varstrs = {'curvbx','curvby','curvbz';'vExBx','vExBy','vExBz';'curvbx.*vExBx','curvby.*vExBy','curvbz.*vExBz'}';
clims = {[-1 1],[-1 1],[-1 1];[-1 1],[-1 1],[-1 1];[-1 1],[-1 1],[-1 1]}';

varstrs = {'curvbx.*vExBx','curvby.*vExBy','curvbz.*vExBz';'tepar','vepar','ne'}';
clims = {[-1 1],[-1 1],[-1 1];[-1 1],[-1 1],[-1 1];[-1 1],[-1 1],[-1 1]}';

h = pic.plot_map(varstrs,'clim',clims,'A',1);

colormap(pic_colors('blue_red'))


%%


pic = no02m.xlim([60 145]).zlim([-1 1]);

varstrs = {'curvbx.*vExBx'}';
clims = {[-1 1],[-1 1],[-1 1];[-1 1],[-1 1],[-1 1];[-1 1],[-1 1],[-1 1]}';

h = pic.plot_timemap('xt',varstrs,'clim',clims,'A',1);

colormap(pic_colors('blue_red'))

