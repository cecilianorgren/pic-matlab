%% plotmap
twpe = 2000;
xlim = [100 300];
xlim = [0 410];
zlim = [-12.5 12.5];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
pic = pic(1);
varstrs = {'Bz';'Jy';'Jx';'Ey';'ni'};
%varstrs = {'Ez';'pxy([2 4 6])';'pyz([2 4 6])'};
h = pic.plot_map(varstrs,'A',1);
%h = pic.plot_map(varstrs,'A',1);

%% plot_map/movie, simulation box sizing
twpe = [7000:200:11000];
twpe = 1000;
xlim = [000 400];
zlim = [0 25];
zlim = [-12.5 12.5];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Ey';'ne';'vx(2)';'vx(4)';'vx(6)';'vz([4 6])'};

clims = {[-1 1];[0 0.5];[-2 2];[-2 2];[-2 2];[-0.5 0.5]}';
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapbr,cmapwa;cmapbr,cmapbr;cmapbr,cmapbr}';

h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps);

filename = [printpath 'nobg_simulation_box_sizing'];
%pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename);

%% plotline, vertical
pic = no02m;
comp = 'z';
twpe = [20];
xlim = 103+1*[-1 1];
zlim = [-15 15];
zlim = pic.zi([1 end]);
pic = pic.xlim(xlim).zlim(zlim);
varstrs = {{'Bx'};{'Jy'};{'n(1)','n([3 5])','n([1 3 5])'};{'Ey','Ez'};{'txx(1)','txx(3)','txx(5)','txx([1 3 5])'}};
varstrs = {{'Bx'};{'Jy'};{'n([1 3 5])'};{'Ey'};{'Ez'}};

h = pic.plot_line(comp,varstrs);

%% plotmap
twpe = 20;
xlim = [190 220];
zlim = [-15 15];
pic = no02.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Ex','Ez';'ni','ne';'vix','vex'}';
clims = {[-1 1],[-1 1];[0 0.5],[0 0.5];[-4 4],[-4 4]}';
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapbr,cmapbr;cmapwa,cmapwa;cmapbr,cmapbr}';

h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps);

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

%% plot_map Epar, 
twpe = 11000;
xlim = [100 200];
zlim = [-5 15];
pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Epar';'ni';'ne'};
clims = {[-1.5 1.5];[0 0.5];[0 0.5]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapbr;cmapwa;cmapwa;cmapbr};

h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps);

%% plotline, magnetic moment, equatorial plane
comp = 'x';
twpe = 6000:1000:12000;%:1000:10000;
xlim = [100 200];
zlim = 0+0.1*[-1 1];
pic = df04.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bz'};{'tperp([3 5])'};{'magmom([3 5])'}};

h = pic.plot_line(comp,varstrs,'smooth',5);
h(2).YLim(1) = 0;
h(3).YLim(1) = 0;

%% Plot PIC.time_map, magnetic moment, equatorial plane
comp = 'xt';
twpe = 6000:1000:12000;%:1000:10000;
twpe = [4000 12000];
xlim = [120 200];
zlim = 0+0.1*[-1 1];
pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Bz';'tperp([3 5])';'magmom([3 5])'};

h = pic.plot_timemap(comp,varstrs);
h(2).CLim = [0 0.15];
h(3).CLim = [0 1.5];

%% Plot PIC.time_map
comp = 'xt';
pic = no02;
twpe = 8000:50:10000;%:1000:10000;
twpe = [3000 4000];pic.twpe;
xlim = diff(pic.xi([1 end]))/2 + [-100 100];
zlim = 4+0.1*[-1 1];
pic = pic.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Bz';'n([1 3])';'t([1 3])'};
varstrs = {'Ey';'ni';'Ex'};

h = pic.plot_timemap(comp,varstrs);

%h(1).CLim = [0 0.5];
%h(2).CLim = [0 1.5];
%h(2).CLim = [0 0.2];

%% plotline, for wenya
comp = 'x';
twpe = 7000:1000:10000;
xlim = [100 210];
xlim = [160 260];
zlim = 0+0.1*[-1 1];
pic = nobg.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bz'};{'By'};{'Bx'};{'ni'};{'vix'};{'Ey'};{'Ez'};{'t([1 3 5])'}};

h = pic.plotline(comp,varstrs);

%% plotline, for wenya
comp = 'x';
twpe = 2000;
xlim = [160 260];
zlim = 1+0.1*[-1 1];
pic = no02.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bx','By','Bz'};{'Ex','Ey','Ez'};{'ni'};{'vex','vix','vExBx'};{'t([1 3 5])'}};

h = pic.plot_line(comp,varstrs);

%% plotline, typical df
comp = 'x';
twpe = 2800;
xlim = [120 210];
zlim = 0+0.05*[-1 1];
pic = no02.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {{'Bz'};{'Ey'};{'vix','vex','vExBx'};{'ni','n([1])','n([3 5])'};{'txx([2 4 6])','tyy([2 4 6])','tzz([2 4 6])'};{'txx([1 3 5])','tyy([1 3 5])','tzz([1 3 5])'}};
%varstrs = {{'Bz'};{'Ey'};{'vix','vex','vExBx'};{'ni','n([1])','n([3])'};{'txx([2 4])','tyy([2 4])','tzz([2 4])'};{'txx([1 3])','tyy([1 3])','tzz([1 3])'}};

h = pic.plot_line(comp,varstrs);

%% plotline, several timesteps
comp = 'x';
twpe = 4000:1000:12000;
xlim = [20 210];
zlim = 0+0.05*[-1 1];
pic = df04.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bz'};{'Ey'};{'ni'};{'txx([2 4 6])'};{'tyy([2 4 6])'};{'tzz([2 4 6])'};{'txx([1 3 5])'};{'tyy([1 3 5])'};{'tzz([1 3 5])'}};
varstrs = {{'t([1 3 5])'};{'t([2 4 6])'}};

h = pic.plot_line(comp,varstrs);

%% plotline, vertical, Ez balance
comp = 'z';
twpe = [9000];
xlim = 204+1*[-1 1];
zlim = [-3 3];
pic = df04.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);

% twpe = 10000;
% xlim = 160+1*[-1 1];

% pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);

varstrs = {{'Bx','Bz','Ey'};{'Ez','-vxBz([1 3 5])','divpz([1 3 5])','dvzdt([1 3 5])'};{'Ez','-vxBz([2 4 6])','-divpz([2 4 6])'};{'Ey','-vxBy([2 4 6])','divpy([2 4 6])','divpy([1 3 5])'};{'pxx([1 3 5])','pyy([1 3 5])','pzz([1 3 5])','pxy([1 3 5])','pxz([1 3 5])','pyz([1 3 5])'}};

h = pic.plot_line(comp,varstrs,'smooth',10);

%% plotline
comp = 'x';
twpe = 8000;
xlim = [160 207];
zlim = 0+0.1*[-1 1];
pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {{'ni','ne'};{'Ex'};{'Ez'};{'txx([4 6])','tyy([4 6])','tzz([4 6])'};{'txx([3 5])','tyy([3 5])','tzz([3 5])'}};
varstrs = {{'ni','ne','n([1])','n([3 5])','n([4 6])','n([4 6])'};{'Ey'};{'Ez'};{'txx([2 4 6])','tyy([2 4 6])','tzz([2 4 6])'};{'txx([1 3 5])','tyy([1 3 5])','tzz([1 3 5])'}};
varstrs = {{'Bz','By'};{'vix','vex','vExBx'};{'Jx','Jz'};{'Ex','Ez'}};
varstrs = {{'Bz','By'};{'vix','vex','vExBx'};{'Jx','Jz'};{'Ex','Ez'}};

h = pic.plot_line(comp,varstrs);

%% plot_line
comp = 'x';
twpe = [8000];
xlim = [130 207];
zlim = 0+0.5*[-1 1];

pic = df04.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bz','Jx'};{'vix','vex','vExBx'};{'Ey','-vexBy','-vixBy'};{'divpy([1 3 5])','divpy([2 4 6])'};{'divpx([1 3 5])','divpx([2 4 6])'};{'PB','pDxx([1 3 5])','pDxx([2 4 6])','pxx([1 3 5])','pxx([2 4 6])'}};
%ylim = {[-1 1]*0.99;[-2 0]*0.99;[-0.1 0.5]*0.99;[-0.1 0.1]*0.99;[-0.1 0.1]*0.99;[0 0.4]*0.99};

varstrs = {{'Bz','Jx'};{'vix','vex','vExBx'};{'Ey','-vixBy','divpy([1 3 5])'};{'Ey','-vexBy','-divpy([2 4 6])'};{'pxy([1 3 5])','pyz([1 3 5])','pxy([2 4 6])','pyz([2 4 6])'}};
%ylim = {[-1 1]*0.99;[-2 0]*0.99;[-0.1 0.4]*0.99;[-0.05 0.4]*0.99;[-0.1 0.1]*0.99};


h = pic.plot_line(comp,varstrs,'smooth',20);

%% plot_map
twpe = 8000;
xlim = [150 205];
zlim = [-10 10];
pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Jx';'Ez';'pxy([3 5])';'pyz([3 5])';'pxy([2 4 6])';'pyz([2 4 6])';'ni'};
clims = {0.7*[-1 1];[-1 1];0.2*[-1 1];0.2*[-1 1];0.1*[-1 1];0.1*[-1 1];[0 2]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapbr;cmapbr;cmapbr;cmapbr;cmapbr;cmapbr;cmapwa};
inds = [3 4];
cmaps = cmaps([3 4]);
clims = clims([3 4]);
varstrs = varstrs([3 4]);

h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps);

%% plottimeseries
comp = 'x';
twpe = df04.twpe;
xlim = 195 + 0.1*[-1 1];
zlim = 0.05 + 0.05*[-1 1];
pic = df04.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bx','By','Bz'};{'Ex','Ey','Ez'}};

h = pic.plottimeseries(varstrs);

%% movie_line
comp = 'x';
twpe = [3000 12000];
xlim = [130 207];
zlim = 0+0.1*[-1 1];

pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {{'ni','ne'};{'Ex'};{'Ez'};{'txx([4 6])','tyy([4 6])','tzz([4 6])'};{'txx([3 5])','tyy([3 5])','tzz([3 5])'}};
varstrs = {{'ni','ne','n([1])','n([3 5])','n([4 6])','n([4 6])'};{'Ey'};{'Ez'};{'txx([2 4 6])','tyy([2 4 6])','tzz([2 4 6])'};{'txx([1 3 5])','tyy([1 3 5])','tzz([1 3 5])'}};
varstrs = {{'Bz','By'};{'vix','vex','vExBx'};{'Jx','Jz'};{'Ex','Ez'}};
varstrs = {{'Bz','By'};{'vix','vex','vExBx'};{'Jx','Jz'};{'viy','vey'}};
ylim = {[-0.5 0.1]*0.99;[-2 0]*0.99;[-1 1]*0.99;[-3 3]*0.99};

h = pic.movie_line(comp,varstrs,'ylim',ylim,'filename',[printpath 'BzBy_vixvexvExBx_JxJz_viyvey_z=0']);

%% movie_line
comp = 'x';
twpe = [4000 12000];
xlim = [130 207];
zlim = 0+1*[-1 1];

pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {{'Bz','Jx'};{'vix','vex','vExBx'};{'Ey','-vexBy','-vixBy'};{'divpy([1 3 5])','divpy([2 4 6])'};{'divpx([1 3 5])','divpx([2 4 6])'};{'PB','pDxx([1 3 5])','pDxx([2 4 6])','pxx([1 3 5])','pxx([2 4 6])'}};
ylim = {[-1 1]*0.99;[-2 0]*0.99;[-0.1 0.5]*0.99;[-0.1 0.1]*0.99;[-0.1 0.1]*0.99;[0 0.4]*0.99};

varstrs = {{'Bz','Jx'};{'vix','vex','vExBx'};{'Ey','-vixBy','divpy([1 3 5])'};{'Ey','-vexBy','-divpy([2 4 6])'};{'pxy([1 3 5])','pyz([1 3 5])','pxy([2 4 6])','pyz([2 4 6])'}};
ylim = {[-1 1]*0.99;[-2 0]*0.99;[-0.1 0.4]*0.99;[-0.1 0.4]*0.99;[-0.1 0.1]*0.99;[0 0.4]*0.99};


h = pic.movie_line(comp,varstrs,'ylim',ylim,'smooth',7,'filename',[printpath 'Ey_z=0_smoothpm7']);

%% movie
twpe = 2000:1000:12000;
xlim = [90 210];
zlim = [0 20];
pic = df04.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {'Ex';'Ez';'Jx';'Jz'};
clims = {[-1 1];[-1 1];0.4*[-1 1];0.4*[-1 1]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapbr;cmapbr;cmapbr;cmapbr};
filename = [printpath 'df04_ExEzJxJz_03000-12000'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename);

%% movie
twpe = [7000 10000];
xlim = [150 210];
zlim = [-10 10];
pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'vx(4)';'vx(6)'};
clims = {10*[-1 1];10*[-1 1]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmapc2 = pic_colors('candy2');
cmapma = pic_colors('matlab');
cmaps = {cmapma;cmapma;cmapbr;cmapbr;cmapbr};
filename = [printpath 'nobg_vx3_vx5_matlab'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename);

%% movie, cold and hot pressure
twpe = 2000:1000:12000;
xlim = [90 210];
zlim = [0 20];
pic = df04.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {'p(1)';'p([3 5])'};
clims = {[0 0.25];[0 0.25];0.4*[-1 1];0.4*[-1 1]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapwa;cmapwa;cmapbr;cmapbr};
filename = [printpath 'df04_p1_p35_02000-12000'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename);

%% movie, ni, Epar
pic = no02;
twpe = [7000 12000];
twpe = [2999 4000];

xlim = pic.xi([1 end])+[100 -100];
xlim = [160 180];
zlim = [0 10];
pic = pic.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'ni';'Epar'};
clims = {[0 0.8];0.5*[-1 1]};
cmapc2 = pic_colors('candy2');
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapwa;cmapbr;cmapbr;cmapbr};
filename = [printpath 'no02_niEpar'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename);

%% movie, Bz
pic = no02.twpelim([1 10000]);
twpe = pic.twpe;
x0 = (pic.xi(end)-pic.xi(1))/2;
xlim = x0 + 100*[-1 1];
%xlim = pic.xi([1 end]);
zlim = [-25 25];
pic = pic.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {'Bz';'Ey'};
clims = {0.5*[-1 1];0.5*[-1 1]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapbr;cmapbr;cmapbr;cmapbr};
filename = [printpath 'no02_BzEy'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename);

%% movie, p, n
twpe = 12000;[5000 13000];
x0 = (nobg.xi(end)-nobg.xi(1))/2;
xlim = x0 + 130*[-1 1];
zlim = [-17 17];
pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'n([1 3 5])';'p([1 3 5])';'t([1 3 5])'};
clims = {[0 0.4],[0 0.2],[0 1]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapwa;cmapwa;cmapwa;cmapbr};
filename = [printpath 'nobg_n135_p135_t135'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename);

%% PICTraj.plot_single
tr = trs.find([trs.t0]==160,[trs.vy0]<-0.25,[trs.x0]>185);
tr([tr.id]==181).tlim([100 240]).plot_single({'xz','ty,tz','tvx,tvy,tvz','tEx,tEy,tEz','tBx,tBy,tBz','tvx,tvy,tvz','tvzBx'}')

tr = trs.find([trs.t0]==160,[trs.vy0]<-0.1,[trs.x0]>195);
tr(1).tlim([100 240]).plot_single({'xz','ty,tz','tvx,tvy,tvz','tEx,tEy,tEz','tBx,tBy,tBz','tvx,tvy,tvz','tvzBx,tEy,tEz'}')

%% PICDist.plot_map
twpe = 10000;
ds = ds04.dxlim([0.3 1]).twpelim(twpe).zfind(0:1:4).xfind(165:176);
twpe = 7000;
ds = ds04.dxlim([0 0.3]).twpelim(twpe).zfind(0:1:4).xfind(173:190);
twpe = 8000;
ds = ds04.dxlim([0.3 1]).twpelim(twpe).zfind(0:2:4).xfind(165:2:176);

%ds = ds04.dxlim([0.3 1]).twpelim([8000]).zfind(2).xfind(170);
clim = [-2 2];
xlim = 0.99*[-1.5 .5];
ylim = 0.99*[-1 1];
sumdim = 2;
fontsie = 16;
h3 = ds.plot_map([3],sumdim,'bline',df04,'v',df04,'log');
h3links = linkprop(h3.ax,{'XLim','YLim','CLim'});
compact_panels(0.00,0.00)
[hax,hlab] = label_panels(h3.ax);
h3.ax(1).CLim = clim;
h3.ax(1).XLim = xlim;
h3.ax(1).YLim = ylim;
c_eval('hlab(?).FontSize = 14;',1:18)
c_eval('h3.ax(?).FontSize = 16;',1:numel(h3.ax))
c_eval('h3.leg(?).FontSize = 16;',1:numel(h3.leg))

%%
h5 = ds.plot_map([5],sumdim,'bline',df04,'v',df04,'log');
h5links = linkprop(h5.ax,{'XLim','YLim','CLim'});
compact_panels(0.00,0.00)
[hax,hlab] = label_panels(h5.ax);
%h5.ax(1).CLim = [0 200];
h5.ax(1).CLim = clim;
h5.ax(1).XLim = xlim;
h5.ax(1).YLim = ylim;
c_eval('hlab(?).FontSize = 16;',1:18)
c_eval('h5.ax(?).FontSize = 16;',1:numel(h5.ax))
c_eval('h5.leg(?).FontSize = 16;',1:numel(h5.leg))
%%
h35 = ds.plot_map([3 5],sumdim,'bline',df04,'v',df04,'log');
h35links = linkprop(h35.ax,{'XLim','YLim','CLim'});
compact_panels(0.00,0.00)
%h35.ax(1).CLim = [0 2];
h35.ax(1).CLim = clim;
h35.ax(1).XLim = xlim;
h35.ax(1).YLim = ylim;
c_eval('hlab(?).FontSize = 16;',1:18)
c_eval('h35.ax(?).FontSize = 16;',1:numel(h35.ax))
c_eval('h35.leg(?).FontSize = 16;',1:numel(h35.leg))
%%
h135 = ds.plot_map([1 3 5],sumdim,'bline',df04,'v',df04,'log');
h135links = linkprop(h135.ax,{'XLim','YLim','CLim'});
compact_panels(0.00,0.00)
%h35.ax(1).CLim = [0 2];
h135.ax(1).CLim = clim;
h135.ax(1).XLim = [-2.5 2.5];
h135.ax(1).YLim = [-2.5 2.5];

%% PICDist.plot_map
twpe = 9000;
ds = ds01.twpelim(twpe).zfind([0 1 2 3 4]).xfind(170:1:185);

sumdim = 2;
%ds = ds04.dxlim([0.3 1]).twpelim([8000]).zfind(2).xfind(170);
clim = [-3 1.5];
xlim = 2*0.99*[-1 1];
ylim = 2*0.99*[-1 1];
fontsie = 16;
h = ds.plot_map([2],sumdim,'bline',nobg,'v',nobg,'log');
hlinks = linkprop(h.ax,{'XLim','YLim','CLim','XTick','YTick'});
compact_panels(0.00,0.00)
[hax,hlab] = label_panels(h.ax);
%h.ax(1).CLim = clim;
%h.ax(1).XLim = xlim;
%h.ax(1).YLim = ylim;
c_eval('hlab(?).FontSize = 12;',1:18)
c_eval('h.ax(?).FontSize = 14;',1:numel(h.ax))
c_eval('h.leg(?).FontSize = 12;',1:numel(h.leg))

%% PICDist.reduce_1d
%ds = ds01.zfind(3);
%fred = ds.reduce_1d_new('x',[5],[]);
twpe = 10000;
xlim = [130 165];
zlim = [-10 10];
ds = ds01.twpelim(twpe).zlim(zlim).xlim(xlim);
id_line1 = 1:274;
id_line2 = 275:462;
ds = ds.update_inds({id_line2});

xdist = (ds.xi1{1}+ds.xi2{1})/2;
zdist = (ds.zi1{1}+ds.zi2{1})/2;
% arclength = [0 cumsum(sqrt(diff(xdist).^2 + diff(zdist).^2))];

pic = nobg.xlim(xlim).zlim(zlim).twpelim(twpe);
Bx_ = pic.Bx;
By_ = pic.By;
Bz_ = pic.Bz;

Bx = interpfield(pic.xi,pic.zi,Bx_,xdist,zdist); 
By = interpfield(pic.xi,pic.zi,By_,xdist,zdist); 
Bz = interpfield(pic.xi,pic.zi,Bz_,xdist,zdist); 

if 1 % saved reduced distributions
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_2.mat')
else % make reduced distributions
  fred1_2 = ds.reduce_1d_new('x',[1],[],'vpar',{Bx,By,Bz});
  fred3_2 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz});
  fred4_2 = ds.reduce_1d_new('x',[4],[],'vpar',{Bx,By,Bz});
  fred5_2 = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz});
  fred6_2 = ds.reduce_1d_new('x',[6],[],'vpar',{Bx,By,Bz});
  fred35_2 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz});
  fred46_2 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz});
end
%%
fredi_str = '35'; iSpecies = [1];
frede_str = '46'; eSpecies = [4 6];
fredi = eval(['fred' fredi_str '_2']);
frede = eval(['fred' frede_str '_2']);
arclength = [0; cumsum(sqrt(diff(fredi.x).^2 + diff(fredi.z).^2))];
if 1; arclength = arclength - arclength(find(abs(fredi.z)==min(abs(fredi.z)))); end
ni = interpfield(pic.xi,pic.zi,pic.ni,fredi.x,fredi.z); 
ne = interpfield(pic.xi,pic.zi,pic.ne,fredi.x,fredi.z); 
vipar = interpfield(pic.xi,pic.zi,pic.vpar(iSpecies),fredi.x,fredi.z); 
vix = interpfield(pic.xi,pic.zi,pic.vx(iSpecies),fredi.x,fredi.z); 
viy = interpfield(pic.xi,pic.zi,pic.vy(iSpecies),fredi.x,fredi.z); 
viz = interpfield(pic.xi,pic.zi,pic.vz(iSpecies),fredi.x,fredi.z); 
vex = interpfield(pic.xi,pic.zi,pic.vx(eSpecies),fredi.x,fredi.z); 
vey = interpfield(pic.xi,pic.zi,pic.vy(eSpecies),fredi.x,fredi.z); 
vez = interpfield(pic.xi,pic.zi,pic.vz(eSpecies),fredi.x,fredi.z); 
vepar = interpfield(pic.xi,pic.zi,pic.vpar(eSpecies),fredi.x,fredi.z); 
Epar = interpfield(pic.xi,pic.zi,pic.Epar,fredi.x,fredi.z); 
vExBx = interpfield(pic.xi,pic.zi,pic.vExBx,fredi.x,fredi.z); 
vExBy = interpfield(pic.xi,pic.zi,pic.vExBy,fredi.x,fredi.z); 
vExBz = interpfield(pic.xi,pic.zi,pic.vExBz,fredi.x,fredi.z);
vExBabs = sqrt(vExBx.^2 + vExBy.^2 + vExBz.^2);
EExB = vExBabs.^2/2;
Bx = interpfield(pic.xi,pic.zi,pic.Bx,fredi.x,fredi.z); 
By = interpfield(pic.xi,pic.zi,pic.By,fredi.x,fredi.z); 
Bz = interpfield(pic.xi,pic.zi,pic.Bz,fredi.x,fredi.z); 

fi_clim = [0 0.013];
fe_clim = [0 7.99e-3];

nrows = 7;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;
doE = 0; colorE = 0*[1 1 1];
doV = 1; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1];

if 1 % line position
  hca = h(isub); isub = isub + 1;
  [ax,h1,h2] = plotyy(hca,arclength,fredi.x,arclength,fredi.z);
  hca.XLabel.String = 'arclength (d_i)';
  ax(1).YLabel.String = 'x';
  ax(2).YLabel.String = 'z';
  legend(hca,{'x','z'},'location','best')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % B
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,Bx,arclength,By,arclength,Bz)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'B';
  legend(hca,{'B_x','B_y','B_z'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % n
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,ni,arclength,ne)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'n';
  legend(hca,{'n_i','n_e'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % v
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,vix,arclength,viy,arclength,viz,arclength,vex,arclength,vey,arclength,vez)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'B';
  legend(hca,{'v_{ix}','v_{iy}','v_{iz}','v_{ex}','v_{ey}','v_{ez}'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % Epar, int(Epar)dl
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,Epar,arclength,-cumtrapz(arclength,Epar))  
  hca.YLabel.String = 'E_{||}, \int E_{||}dl_{||}'; 
  legend(hca,{'E_{||}','-\int E_{||}dl_{||}'},'location','eastoutside') 
  hca.XLabel.String = 'arclength (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % f35(vx)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.v,fredi.fvx')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{x}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},v_{x})'];
  hca.CLim(2) = prctile(fredi.fvx(:),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if 0*doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vExBx,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % f35(vy)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.v,fredi.fvy')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},v_{y})'];
  hca.CLim(2) = prctile(fredi.fvy(:),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if 0*doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vExBy,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % f35(vz)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.v,fredi.fvz')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{z}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},v_{z})'];
  hca.CLim(2) = prctile(fredi.fvz(:),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if 0*doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vExBz,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % f35(vpar)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.vpar_center,fredi.fvpar')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(fredi.fvpar(:),99);
  hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % f35(vabs)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.vabs_center,fredi.fvabs')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},|v|)'];
  hca.CLim(2) = prctile(fredi.fvabs(:),99);
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if 0%doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if 0%doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % def35(vabs^2)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.vabs_center.^2/2,log10(fredi.fdefE'))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'mv^2/2';
  colormap(hca,pic_colors('candy'))
  colormap(hca,[ 1 1 1; pic_colors('waterfall')])
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['dpf_{i,' fredi_str '}(l_{||},mv^2/2)'];
  hca.CLim = [prctile(log10(fredi.fdefE(:)),5) prctile(log10(fredi.fdefE(:)),95)];
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YScale = 'log';
  hca.YTick = 10.^(-10:1:10);
  hca.YLim(1) = 10^(-2);
  if doExB
    hold(hca,'on')
    plot(hca,arclength,EExB,'color',colorExB)
    hold(hca,'off')
  end
  if 0%doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if 0%doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % f46(vpar)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,frede.vpar_center,frede.fvpar')
  hca.YLim = [-7 7];
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))
  hcb = colorbar('peer',hca);
  hca.CLim = fe_clim;
  hcb.YLabel.String = ['f_{e,' frede_str '}(l_{||},v_{||})'];  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vepar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  
end
if 0 % def35(vabs^2)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,frede.vabs_center.^2/2/25,log10(frede.fdefE'))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'mv^2/2';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['dpf_{i,' frede_str '}(l_{||},mv^2/2)'];
  hca.CLim = [prctile(log10(frede.fdefE(:)),5) prctile(log10(frede.fdefE(:)),95)];
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YScale = 'log';
  hca.YTick = 10.^(-10:1:10);
  hca.YLim(1) = 10^(-2);
  if doExB
    hold(hca,'on')
    plot(hca,arclength,EExB/25,'color',colorExB)
    hold(hca,'off')
  end
  if 0%doE
    hold(hca,'on')
    plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    hold(hca,'off')
  end
  if 0%doV
    hold(hca,'on')
    plot(hca,arclength,vipar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
compact_panels(0.01)
h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h_,{'XLim'});
hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
for ip = 1:nrows*ncols
  axwidth(ip) = h(ip).Position(3);
end
for ip = 1:nrows*ncols
  h(ip).Position(3) = min(axwidth);
end



ax(2).YAxisLocation = 'right';

%% Plot quantities along a field line
sep = separatrix_location(nobg);
