%% plotmap
twpe = 2000;
xlim = [100 300];
zlim = [-25 25];
pic = no02.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Bz';'Jy';'Ex';'Ey';'ni';'vix';'vex'};
%varstrs = {'Ez';'pxy([2 4 6])';'pyz([2 4 6])'};
h = pic.plot_map(varstrs,'A',1);

%% plotline, vertical
comp = 'z';
twpe = [1000];
xlim = 205+1*[-1 1];
zlim = [-15 15];
pic = no02(1).xlim(xlim).zlim(zlim);
varstrs = {{'Bx'};{'n(1)','n([3 5])','n([1 3 5])'};{'Ey'};{'txx(1)','txx(3)','txx(5)','txx([1 3 5])'}};

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
pic = bs;
twpe = 7000:1000:10000;%:1000:10000;
twpe = [6000 12000];
twpe = bs.twpe;

xlim = diff(pic.xi([1 end]))/2 + [-100 100];
zlim = 0+0.1*[-1 1];
pic = pic.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {'Bz';'n([1 3])';'t([1 3])'};

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
twpe = 8000;
xlim = [20 210];
zlim = 2+0.05*[-1 1];
pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {{'Bz'};{'Ey'};{'vix','vex','vExBx'};{'ni','n([1])','n([3 5])'};{'txx([2 4 6])','tyy([2 4 6])','tzz([2 4 6])'};{'txx([1 3 5])','tyy([1 3 5])','tzz([1 3 5])'}};

h = pic.plotline(comp,varstrs);

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
twpe = 8000;
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
ds = ds01.twpelim(twpe).zfind([0 1 2 3]).xfind(179:2:210);

sumdim = 2;
%ds = ds04.dxlim([0.3 1]).twpelim([8000]).zfind(2).xfind(170);
clim = [-3 1.5];
xlim = 2*0.99*[-1 1];
ylim = 2*0.99*[-1 1];
fontsie = 16;
h = ds.plot_map([4 6],sumdim,'bline',df04,'v',df04,'log');
hlinks = linkprop(h.ax,{'XLim','YLim','CLim'});
compact_panels(0.00,0.00)
[hax,hlab] = label_panels(h.ax);
%h.ax(1).CLim = clim;
%h.ax(1).XLim = xlim;
%h.ax(1).YLim = ylim;
c_eval('hlab(?).FontSize = 12;',1:18)
c_eval('h.ax(?).FontSize = 14;',1:numel(h.ax))
c_eval('h.leg(?).FontSize = 12;',1:numel(h.leg))

%% PICDist.reduce_1d
ds = ds01.zfind(3);
fred = ds.reduce_1d_new('x',[5],[]);