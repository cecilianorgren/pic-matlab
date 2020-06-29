%% plotmap
twpe = 10000;
xlim = [140 210];
zlim = [0 15];
pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Ex','Ez';'ni','ne';'vix','vex'}';
clims = {[-1 1],[-1 1];[0 0.5],[0 0.5];[-4 4],[-4 4]}';
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapbr,cmapbr;cmapwa,cmapwa;cmapbr,cmapbr}';

h = pic.plotmap(varstrs,'A',1,'clim',clims,'cmap',cmaps);

%% plotmap
twpe = 10000;
xlim = [140 210];
zlim = [0 15];
pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Ex';'Ez';'Jx';'Jz'};
clims = {[-1 1];[-1 1];0.4*[-1 1];0.4*[-1 1]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapbr;cmapbr;cmapbr;cmapbr};

h = pic.plotmap(varstrs,'A',1,'clim',clims,'cmap',cmaps);

%% plotline, magnetic moment, equatorial plane
comp = 'x';
twpe = 6000:1000:12000;%:1000:10000;
xlim = [100 200];
zlim = 0+0.1*[-1 1];
pic = df04.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bz'};{'tperp([3 5])'};{'magmom([3 5])'}};

h = pic.plotline(comp,varstrs);
h(2).YLim(1) = 0;
h(3).YLim(1) = 0;

%% plotline, for wenya
comp = 'x';
twpe = 7000:1000:10000;
xlim = [100 210];
zlim = 1+0.1*[-1 1];
pic = nobg.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bz'};{'By'};{'Bx'};{'ni'};{'vix'};{'Ey'};{'Ez'};{'t([1 3 5])'}};

h = pic.plotline(comp,varstrs);

%% plotline, typical df
comp = 'x';
twpe = 8000;
xlim = [20 210];
zlim = 2+0.05*[-1 1];
pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {{'Bz'};{'Ey'};{'vix','vex','vExBx'};{'ni','n([1])','n([3 5])'};{'txx([2 4 6])','tyy([2 4 6])','tzz([2 4 6])'};{'txx([1 3 5])','tyy([1 3 5])','tzz([1 3 5])'}};

h = pic.plotline(comp,varstrs);

%% plotline, several timesteps
comp = 'x';
twpe = 5000:1000:10000;
xlim = [20 210];
zlim = 0+0.05*[-1 1];
pic = nobg.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bz'};{'Ey'};{'ni'};{'txx([2 4 6])'};{'tyy([2 4 6])'};{'tzz([2 4 6])'};{'txx([1 3 5])'};{'tyy([1 3 5])'};{'tzz([1 3 5])'}};

h = pic.plotline(comp,varstrs);

%% plotline, vertical, Ez balance
comp = 'z';
twpe = 8000;
xlim = 170+1*[-1 1];
pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);

twpe = 10000;
xlim = 160+1*[-1 1];
pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);
zlim = [0 10];

varstrs = {{'Bx'};{'By'};{'Ez','gradpz([1 3 5])'};{'Jx'};{'vix','vex','vExBx'}};

h = pic.plotline(comp,varstrs);

%% plotline
comp = 'x';
twpe = 10000;
xlim = [120 210];
zlim = 0+0.05*[-1 1];
pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {{'ni','ne'};{'Ex'};{'Ez'};{'txx([4 6])','tyy([4 6])','tzz([4 6])'};{'txx([3 5])','tyy([3 5])','tzz([3 5])'}};
varstrs = {{'ni','ne','n([1])','n([3 5])','n([4 6])','n([4 6])'};{'Ey'};{'Ez'};{'txx([2 4 6])','tyy([2 4 6])','tzz([2 4 6])'};{'txx([1 3 5])','tyy([1 3 5])','tzz([1 3 5])'}};

h = pic.plotline(comp,varstrs);

%% movie
twpe = 3000:1000:10000;
xlim = [140 210];
zlim = [0 15];
pic = df04.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {'Ex';'Ez';'Jx';'Jz'};
clims = {[-1 1];[-1 1];0.4*[-1 1];0.4*[-1 1]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapbr;cmapbr;cmapbr;cmapbr};
filename = [printpath 'df04_ExEzJxJz'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename);