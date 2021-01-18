%% reconnection rate comparison
if 0
  nobg = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_test/data_h5/fields.h5');
  no02 = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02/data_h5/fields.h5');
  no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
  df04n = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/data_h5/fields.h5');
  bs = PIC('/Volumes/Fountain/Data/PIC/baselinev4/data_h5/fields.h5');
end

%pics = {nobg,no02,no02m,df04,df04n,bs};
pics = {no02m};
mass_cs = [1 1 1 1.2 1.2 1.2];
mass_infl = [0.1 0.2 0.2 0.6 0.6 0.2];
% ER = nin/n0
massload = sqrt(mass_infl)./sqrt(mass_cs);
% fun_massloading = @(nh,nc) sqrt(1+nc/nh);
pic_legs = {'nc=0.1','nc=0.2','nc=0.2,mi=100','nc=0.4,nh=0.2','nc=0.4,nh=0.2,further out','nh=0.2'};
npics = numel(pics);
nrows = 3; ncols = 2; npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
colors = pic_colors('matlab');
if 0 % initial A(x=x_xline,z)
  hca = h(isub); isub = isub + 1;
  for ipic = 1:npics
    pic = pics{ipic};
    x0 = (pic.xi(end)-pic.xi(1))/2;
    A = squeeze(mean(pic(1).xlim(x0+1*[-1 1]).A,1));
    plot(hca,pic.zi,A-min(A))    
    if ipic == 1, hold(hca,'on'); end
    if ipic == npics, hold(hca,'off'); end
  end
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'A(x=x_{line},z)';
  legend(hca,pic_legs,'location','best','box','off')
end
if 0 % initial A(x,z=z_xline)
  hca = h(isub); isub = isub + 1;
  for ipic = 1:npics
    pic = pics{ipic};
    z0 = 0;
    x0 = (pic.xi(end)-pic.xi(1))/2;
    A = squeeze(mean(pic(1).zlim(z0+1*[-1 1]).A,2));
    plot(hca,pic.xi-x0,A-min(A))    
    if ipic == 1, hold(hca,'on'); end
    if ipic == npics, hold(hca,'off'); end
  end
  hca.XLim = [-25 25];
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'A(x,z=z_{line})';
  legend(hca,pic_legs,'location','best','box','off')
end
if 0 % initial A(x,z)
  hca = h(isub); isub = isub + 1;
  for ipic = 1:npics
    pic = pics{ipic};
    A = pic(1).A;
    levA = -25:1:10;
    step = 5;
    x0 = (pic.xi(end)-pic.xi(1))/2;
    contour(hca,pic.xi(1:step:end)-x0,pic.zi(1:step:end),A(1:step:end,1:step:end)','color',colors(ipic,:))
    drawnow
    if ipic == 1, hold(hca,'on'); end
    if ipic == npics, hold(hca,'off'); end
  end
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  %legend(hca,pic_legs,'location','best','box','off')
end
if 1 % ER(twci)
  hca = h(isub); isub = isub + 1;
  for ipic = 1:npics
    pic = pics{ipic};
    plot(hca,pic.twci,pic.RE)
    if ipic == 1, hold(hca,'on'); end
    if ipic == npics, hold(hca,'off'); end
  end
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'E_R';
  legend(hca,pic_legs,'location','best','box','off')
end
if 1 % ER(twpe)
  hca = h(isub); isub = isub + 1;
  for ipic = 1:npics
    pic = pics{ipic};
    plot(hca,pic.twpe,pic.RE)
    if ipic == 1, hold(hca,'on'); end
    if ipic == npics, hold(hca,'off'); end
  end
  hca.XLabel.String = 't\omega_{pe}';
  hca.YLabel.String = 'E_R';
  legend(hca,pic_legs,'location','best','box','off')
end
if 1 % ER(dUB)
  hca = h(isub); isub = isub + 1;
  for ipic = 1:npics
    pic = pics{ipic};
    plot(hca,(pic.UB(1)-pic.UB)/pic.UB(1),pic.RE)
    if ipic == 1, hold(hca,'on'); end
    if ipic == npics, hold(hca,'off'); end
  end
  hca.XLabel.String = '(U_B(t=0)-U_B)/U_B(t=0)';
  hca.YLabel.String = 'E_R';
  legend(hca,pic_legs,'location','best','box','off')
end
if 1 % ER(twci) massloading
  hca = h(isub); isub = isub + 1;
  for ipic = 1:npics
    pic = pics{ipic};
    plot(hca,pic.twci,pic.RE*massload(ipic))
    if ipic == 1, hold(hca,'on'); end
    if ipic == npics, hold(hca,'off'); end
  end
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'E_Rv_{A,in}/v_{A,0}';
  legend(hca,pic_legs,'location','best','box','off')
end
if 1 % ER(dUB)
  hca = h(isub); isub = isub + 1;
  for ipic = 1:npics
    pic = pics{ipic};
    plot(hca,(pic.UB(1)-pic.UB)/pic.UB(1),pic.RE*massload(ipic))
    if ipic == 1, hold(hca,'on'); end
    if ipic == npics, hold(hca,'off'); end
  end
  hca.XLabel.String = '(U_B(t=0)-U_B)/U_B(t=0)';
  hca.YLabel.String = 'E_Rv_{A,in}/v_{A,0}';
  legend(hca,pic_legs,'location','best','box','off')
end

for ip = 1:npanels
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end

%% plotmap
twpe = 24000;

xlim = [200 210];
zlim = [-0 12];
xlim = [00 110];
xlim = [140 260];
xlim = [60 150]+100;
xlim = [50 150];
%pic = no02.twpelim(twpe).xlim(xlim).zlim(zlim);
%pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
%pic = nobg.twpelim(8200).xlim(xlim).zlim(zlim);
%pic = no02.twpelim(twpe).xlim(xlim).zlim(zlim);
%pic = df04n.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Ex+vixBx';'Ey+vixBy';'Ez+vixBz'};
varstrs = {'Ex+vexBx';'Ey+vexBy';'Ez+vexBz'};





varstrs = {'viy';'ti';'te';'log10(ni([3 5]))'};

varstrs = {'By';'Ey';'Jy';'n([1])';'n([3 5])';'A'};
varstrs = {'By';'Jy';'Jx';'Ey';'ni'};
varstrs = {'Epar';'log10(ne)';'ti';'te';'n([3 5])'};
varstrs = {'log10(ne)';'Ez';'viy';'vex';'n([3 5])'};
clims = {[-2 0.2],[-2 2],[-2 2],[-6 6],[0 0.4],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1]};
varstrs = {'Ez','vix.*By','-viy.*Bx','vixBz','Ez+vixBz','-divpiz./ni';'vdviz','-divpez./ni','vdvez','Ez+vixBz-divpiz./ni','JxBz','JxBz./ne'}';
varstrs = {'Ey','viz.*Bx','-vix.*Bz','vixBy','Ey+vixBy','-divpiy./ni';'vdviy','-divpey./ni','vdvey','Ey+vixBy-divpiy./ni','JxBy','JxBy./ne'}';
%varstrs = {'Ez','vix.*By','-viy.*Bx','vixBz','Ez+vixBz';'-divpiz./ni','vex.*By','-vey.*Bx','vexBz','Ez+vexBz','-divpez./ni'}';
clims = {[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1],[-1 1]};
varstrs = {'Epar';'log10(ne)';'ti';'te';'n([3 5])'};
varstrs = {'By';'Jy';'Jx';'Ey';'ni'};
varstrs = {'Ez','vix.*By','-viy.*Bx','vixBz','Ez+vixBz','-divpiz./ni';'vex.*By','-vey.*Bx','vexBz','Ez+vexBz','-divpez./ne','Ez+vexBz+divpez./ne'}';
varstrs = {'Ez','vix.*By','-viy.*Bx','vixBz','Ez+vixBz','-divpiz./ni','Ez+vixBz-divpiz./ni','vdviz','Ez+vixBz-divpiz./ni-vdviz'}';
varstrs = {'Ez','vex.*By','-vey.*Bx','vexBz','Ez+vexBz','+divpez./ne','Ez+vexBz+divpez./ne','-vdvez','Ez+vexBz+divpez./ne+vdvez'}';
varstrs = {'Ez','By','vx([4 6])'}';
clims = {0.5*[-1 1],0.5*[-1 1],2*[-1 1]};
%varstrs = {'Jx';'Epar';'log10(abs(n([1])))';'log10(abs(n([3 5])))'};
%varstrs = {'Ey';'Ey+vexBy';'Ey+vixBy';'vex';'hvey';'vez';'vix'};
%varstrs = {'Babs','By','log10(ni)','log10(magmom([3 5]))','log10(magmom([4 6]))'}';
%varstrs = {'Ez';'pxy([2 4 6])';'pyz([2 4 6])'};
varstrs = {'Epar';'ne';'ni';'log10(ni)'};
clims = {0.5*[-1 1],[0 1],[0 1],[-1.5 1]};
varstrs = {'Ez';'Bx';'By';'Jy'};
clims = {0.5*[-1 1],[-1 1],[-1 1],[-1 1]};

varstrs = {'dvzdt([3 5])';'vdvz([3 5])';'vz([3 5])'};
clims = {0.5*[-1 1],0.5*[-1 1],[-1 1],[-1 1]};

cmapbr = pic_colors('blue_red');
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr}';
%pic = gf05(gf05.nt);
pic = no02m.twpelim(24000).xlim(xlim).zlim(zlim);
h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps);
%h = pic.plot_map(varstrs,'A',0.5);
colormap(pic_colors('blue_red'))

%% plotmap, ion and electron edge
twpe = 24000;

xlim = [50 150];
zlim = [0 15];

varstrs = {'n([5])','n([6])','n([3])','n([4])','n([3])-n([4])','n([5])-n([6])','Ez'}';
clims = {[0 0.3],[0 0.3],[0 0.3],[0 0.3],[-0.3 0.3],[-0.3 0.3],[-1 1]};

cmapcan = pic_colors('candy4');
cmapbr = pic_colors('blue_red');
cmaps = {cmapcan,cmapcan,cmapcan,cmapcan,cmapbr,cmapbr,cmapbr}';

pic = no02m.twpelim(24000).xlim(xlim).zlim(zlim);
h = pic.plot_map(varstrs,'A',1,'sep','clim',clims,'cmap',cmaps);

%h = no02m.twpelim(17000).xlim([50 150]).zlim([0 15]).plot_map({'vz(3)','vt(3)'}','A',1,'sep','clim',{[-1 1],[-1 1]},'cmap',{cmapbr,cmapbr});

%% plot_map/movie, simulation box sizing
twpe = [7000:200:11000];
twpe = [7000:1000:24000];
%twpe = [200:200:6000];
xlim = [0 110];
zlim = [0 15];
%pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);
%pic = nobg.twpelim(10000,'exact').xlim(xlim).zlim(zlim);
pic = no02m.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {' By';'n([4 6])';'vx(2)';'vx(4)';'vx(6)';'vz([4 6])'};
%pic = pic(pic.nt);
clims = {[-1 1];[0 0.5];[-2 2];[-2 2];[-2 2];[-0.5 0.5]}';
clims = {0.2*[-1 1];[0 0.5];[-2 2];[-2 2];[-2 2];[-0.5 0.5]}';


varstrs = {'Ey';'Jy';'n([4 6])';'vx(4)'};
varstrs = {'Epar';'log10(ne)';'ti';'te'};
clims = {[-1 1];[-2.0 0.5];[0 1];[0 0.7]}';
varstrs = {'log10(ne)';'ti';'te';'By'};
clims = {[-2.0 0.5];[0 1];[0 0.7];[-0.5 0.5];}';

varstrs = {'By';'Ez';'vx([4 6])'};
clims = {0.5*[-1 1];0.5*[-1 1];2*[-1 1]}';
%clims = {[-2.0 0.5];[0 1];[0 0.7];[-0.5 0.5];}';


cmapbr = pic_colors('blue_red');
cmapjet = colormap('jet');
cmapwa = pic_colors('waterfall');
cmaps = {cmapbr,cmapwa;cmapbr,cmapbr;cmapbr,cmapbr}';
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr}';

%h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps);

filename = [printpath 'no02m_By_Ez_vx46'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename);

%% plot_line, horizontal
pic = no02(15);
comp = 'x';
twpe = [1000:1000:12000];
twci = [20 100];
twpe = 19000;
xlim = pic.xi([1 end])+[70 -70]';
zlim = 0*[-0.5 0.5];
pic = pic.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');
%pic = pic.xlim(xlim).zlim(zlim).twcilim(twci);
varstrs = {{'Bx'};{'Jy'};{'n(1)','n([3 5])','n([1 3 5])'};{'Ey','Ez'};{'txx(1)','txx(3)','txx(5)','txx([1 3 5])'}};
varstrs = {{'Bx'};{'Jy'};{'n([1 3 5])'};{'Ey'};{'Ez'};{'viz'};{'vez'}};
varstrs = {{'Ey'};{'Jy'};{'ne'};{'vex'};{'vix'}};
varstrs = {{'Ey'};{'Jy'};{'ne'};{'n([4 6])'};{'vex'};{'vix'}};
varstrs = {{'Ex'};{'Ey'};{'vex'};{'vix'}};
%varstrs = {{'Ey'};{'n(4)'};{'tzz(4)'}};
varstrs = {{'Bz','By','Bz'};{'Ex','Ey','Ez'};{'Ey+vexBy'};{'Ey+vixBy'};{'vex';'vey';'vez'};{'vix';'viy';'viz'}};
varstrs = {{'Bz','By','Bz'};{'Ex','Ey','Ez'};{'Ey+vexBy'};{'Ey+vxBy([3 5])'};{'Ey+vxBy([1])'};{'vex';'vey';'vez'};{'vix';'viy';'viz'}};
varstrs = {{'Bz','By','Bz'};{'PB','p(1)','p([3 5])'}};

h = pic.plot_line(comp,varstrs,'smooth',10);

%% plotline, vertical
pic = gf05(gf05.nt);
comp = 'z';
twpe = [100];
twpe = 1000;
xlim = 53+0.5*[-1 1];
zlim = [-15 15];
zlim = pic.zi([1 end]);
zlim = [-15 15];
xlim = 60+0.5*[-1 1];
pic = pic.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');
varstrs = {{'Bx'};{'Jy'};{'n(1)','n([3 5])','n([1 3 5])'};{'Ey','Ez'};{'txx(1)','txx(3)','txx(5)','txx([1 3 5])'}};
varstrs = {{'Bx'};{'Jy'};{'n([1 3 5])'};{'Ey'};{'Ez'};{'viz'};{'vez'}};
varstrs = {{'Ey'};{'Ez'};{'vez'};{'tzz(4)'};{'txx(4)'}};
varstrs = {{'tzz([1])'};{'tzz(2)'};{'ni'}};
%varstrs = {{'Ey'};{'n(4)'};{'tzz(4)'}};
varstrs = {{'n([1])','n([3 5])'};{'Ey'};{'Ez'}};
varstrs = {{'Ez';'Bx';'By';'Jy'}};
h = pic.plot_line(comp,varstrs);

%% plotline, vertical
pic = df04n;
comp = 'z';
twpe = [3000:1000:8000];
%twpe = 3000;
xlim = 205+[-1 1];
zlim = [0 15];
%zlim = pic.zi([1 end]);
zlim = [0 10];
pic = pic.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');
%pic = pic.xlim(xlim).zlim(zlim);
varstrs = {{'Bx'};{'Jy'};{'n(1)','n([3 5])','n([1 3 5])'};{'Ey','Ez'};{'txx(1)','txx(3)','txx(5)','txx([1 3 5])'}};
varstrs = {{'Bx'};{'Jy'};{'n([1 3 5])'};{'Ey'};{'Ez'};{'viz'};{'vez'}};
varstrs = {{'n(2)','n(4)','n([2 4])'};{'tzz(2)','tzz(4)','txx([2 4])'};{'txx(2)','txx(4)','txx([2 4])'}};
varstrs = {{'n(1)','n(4)','n([2 4])'};{'tzz(2)','tzz(4)','txx([2 4])'};{'txx(2)','txx(4)','txx([2 4])'}};
varstrs = {{'PB'};{'pi'};{'pe'}};
varstrs = {{'PB'};{'ni'};{'n(1)'};{'n([3 5])'};{'Jy'}};
%varstrs = {{'n([3 5])'}};

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

%% plotmap, magnetic curvature
twpe = 9000; xlim = [160 210]; zlim = [-10 10];

pic = nobg.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');

varstrs = {'curvbx';'curvby';'curvbz';'curvbabs';'curvbrad';'rc([3 5])'};
clims = {[-1 1];[-1 1];[-1 1];[-1 1];[0 1];[0 1];[0 2]};
cmaps = {cmapbr;cmapbr;cmapbr;cmapbr;cmapbr;cmapwa};

varstrs = {'curvbabs';'curvbrad';'rc([3 5])'};
clims = {[-1 1];[0 5];[0 5];[0 2]};

varstrs = {'vtperp([3 5])';'vtpar([3 5])';'curvbabs';'log10(curvbrad)';'log10(rc([3 5]))';'log10(curvbrad./rc([3 5]))'};
cmaps = {cmapbr;cmapbr;cmapbr;cmapwa;cmapwa;cmapwa};
clims = {[0 2];[0 2];[0 5];[0 5];[-1 1];[-3 3]};
varstrs = {'log10(ni)';'vtperp([3 5])';'vtpar([3 5])';'vtperp([4 6])';'vtpar([4 6])'};
clims = {[-2 1];[0 2];[0 2];[0 7];[0 7]};
cmaps = {cmapbr;cmapbr;cmapbr;cmapbr;cmapbr;cmapwa;cmapwa};


h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps);

%% plot_map Epar, 
twpe = 10000;
xlim = [100 200];
zlim = [-5 15];
pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Epar';'log10(ni)';'ni-ne';'te'};
clims = {0.5*[-1.0 1.0];[-2 0.5];0.02*[-1 1];[0 0.5]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapbr;cmapbr;cmapbr;cmapbr};

h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps,'smooth',40);

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

%% PIC.time_map, magnetic moment, equatorial plane
comp = 'xt';
twpe = 6000:1000:12000;%:1000:10000;
twpe = [4000 12000];
xlim = [120 200];
zlim = 0+0.1*[-1 1];
pic = df04.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Bz';'tperp([3 5])';'magmom([3 5])'};

h = pic.plot_timemap(comp,varstrs);
%h(2).CLim = [0 0.15];
%h(3).CLim = [0 1.5];

%% PIC.time_map, equatorial plane
comp = 'xt';
twpe = 6000:1000:12000;%:1000:10000;
twpe = 5000:1000:24000;%:1000:10000;
twpe = [5000 24000];
xlim = [50 150];
zlim = 0+0.1*[-1 1];
twci = [1 400];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Bz','Ey','ni','ti','te'}';

h = pic.plot_timemap(comp,varstrs,'A',1);
%h(2).CLim = [0 0.15];
%h(3).CLim = [0 1.5];

%% PIC.time_map
comp = 'zt';
pic = no02m;
twpe = 7000:100:10000;%:1000:10000;
twpe = [3000 4000];pic.twpe;
twpe = [1000 24000];
twci = [70:2:150];
twci = [95:1:110];
twci = [120 180];
twci = [1 300];
xlim = diff(pic.xi([1 end]))/2 + 0.5*[-1 1];
zlim = 15*[-1 1]; 
pic = pic.twcilim(twci).xlim(xlim).zlim(zlim); 
pic = pic.twpelim(twpe).xlim(xlim).zlim(zlim); 
varstrs = {'Bz';'n([1 3])';'t([1 3])'};
varstrs = {'Ey';'ni';'vy(1)';'vy([3 5])';'viy'};
varstrs = {'Ey';'Jy'};
varstrs = {'Ey'};

h = pic.plot_timemap(comp,varstrs,'A',1);

%h(1).CLim = [0 0.5];
%h(2).CLim = [0 1.5];
%h(2).CLim = [0 0.2];

%% plotline, for wenya
comp = 'x';
twpe = 7000:1000:10000;
xlim = [100 210];
%xlim = [160 260];
zlim = 0.5+0.1*[-1 1];
pic = nobg.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bz'};{'By'};{'Bx'};{'ni'};{'vix'};{'viy'};{'Ey'};{'Ez'};{'t([1 3 5])'}};

h = pic.plot_line(comp,varstrs);

%% plotline, Ohm's law
comp = 'x';
twpe = [1000:1000:9000];
twpe = 3000;
xlim = [120 200];
zlim = 0.0+0.5*[-1 1];
pic = no02.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bx','By','Bz'};{'Ex','Ey','Ez'};{'ni'};{'divpx([1 3 5])','divpx([1])','divpx([3 5])','divpx(1+[1 3 5])'};{'dvxdt([1])','vdvx(1)'};{'vex','vix','vExBx'};{'t([1 3 5])'}};
varstrs = {{'Bx','By','Bz'};{'Ex','Ey','Ez'};{'divpx([1 3 5])','vxBx([1 3 5])','Ex','-divpx([1 3 5])+vxBx([1 3 5])'};{'dvxdt([1])','vdvx(1)'}};
varstrs = {{'Bx','By','Bz'};{'Ex','Ey','Ez'};{'divpx([1 3 5])','vxBx([1 3 5])','Ex','-divpx([1 3 5])+vxBx([1 3 5])'};{'t([1 3 5])'}};
varstrs = {{'Bz','Ey','Jy'};   {'pi','pe'};{'JxBx','JxBy','JxBz'};{'Ex','vixBx','divpex','JxBx./ni','JxBx./ni-vixBx+divpex'};{'Ex','vixBx','divpix','-vixBx+divpix'};{'Ex','-vexBx','-divpex','-vexBx-divpex'};{'Ex+vixBx','Ex+vxBx(1)','Ex+vxBx([3 5])','Ex+vexBx'}};
varstrs = {{'Bz','Ey','Jy','pi','pe'};{'jiy','jey'};{'Ex','vixBx','divpex','JxBx./ni','JxBx./ni-vixBx+divpex'};{'Ex','vixBx','divpix','-vixBx+divpix'};{'Ex','-vexBx','-divpex','-vexBx-divpex'};{'Ex+vixBx','Ex+vexBx'}};
%varstrs = {{'Bz','Ey','Jy','pi','pe'};{'Ex','vixBx','divpex','JxBx./ni','JxBx./ni-vixBx+divpex'};{'Ex','vixBx','divpix','dvixdt','vdvix','-vixBx+divpix+dvixdt+vdvix'};{'Ex','-vexBx','-divpex','dvexdt','vdvex','-vexBx-divpex'}};     
%varstrs = {{'Bz','Ey','Jy'};{'Jx','Jy','Jz'};{'JxBx','JxBy','JxBz'};{'Ey','vixBy','divpey','JxBy','JxBy-vixBy+divpey'};{'Ey','divpiy','vixBy','-vixBy+divpiy'};{'Ey','-vexBy','-divpey','-vexBy-divpey'};{'Ey+vixBy','Ey+vexBy'}};
%varstrs = {{'Bz','Ey','Jy'};   {'JxBx','JxBy','JxBz'};{'Ey','vixBy','divpey','JxBy','JxBy-vixBy+divpey'};{'Ey','divpiy','vixBy','-vixBy+divpiy'};{'Ey','-vexBy','-divpey','-vexBy-divpey'};{'Ey+vixBy','Ey+vexBy'}};
%varstrs = {{'Bz','Ey'};{'Jz','Jy','Jz'};{'JxBx','JxBy','JxBz'};{'Ex','JxBx','vxBx([1 3])','divpx(1+[1 3])','JxBx-vxBx([1 3])+divpx(1+[1 3])'};{'Ex','divpx([1 3])','vxBx([1 3])','-vxBx([1 3])+divpx([1 3])'};{'Ex','-vxBx(1+[1 3])','-divpx(1+[1 3])','-vxBx(1+[1 3])-divpx(1+[1 3])'}};
%varstrs = {{'Bz'};{'t([1 3 5])'}};

h = pic.plot_line(comp,varstrs,'smooth',20);

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
twpe = [1000];
xlim = 50+0.2*[-1 1];
zlim = [-10 10];
pic = nobg.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);

% twpe = 10000;
% xlim = 160+1*[-1 1];

% pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);

varstrs = {{'Bx','Bz','Ey'};{'Ez','-vxBz([1 3 5])','divpz([1 3 5])','dvzdt([1 3 5])'};{'Ez','-vxBz([2 4 6])','-divpz([2 4 6])'};{'Ey','-vxBy([2 4 6])','divpy([2 4 6])','divpy([1 3 5])'};{'pxx([1 3 5])','pyy([1 3 5])','pzz([1 3 5])','pxy([1 3 5])','pxz([1 3 5])','pyz([1 3 5])'}};
varstrs = {{'n(1)','n([3 5])','n(2)','n([4 6])'};{'Ez','-vxBz([1 3 5])','-vxBz([3 5])','-vxBz([2])','-vxBz([4 6])'};{'Ez','-vxBz([1 3 5])','divpz([1 3 5])','dvzdt([1 3 5])','vdvz([1 3 5])','-vxBz([1 3 5])+divpz([1 3 5])+dvzdt([1 3 5])+vdvz([1 3 5])'};{'Ez','-vxBz([2 4 6])','-divpz([2 4 6])'}};
varstrs = {{'PB','pi','pe','PB+pi+pe'};{'n(1)','n([3 5])','n(2)','n([4 6])'};{'Ez','-vxBz([1 3 5])','-vxBz([3 5])','-vxBz([2])','-vxBz([4 6])'};{'Ez','-vxBz([1])','divpz([1])','dvzdt([1])','vdvz([1])','-vxBz([1])+divpz([1])+dvzdt([1])+vdvz([1])'};{'Ez','-vxBz([3 5])','divpz([3 5])','dvzdt([3 5])','vdvz([3 5])','-vxBz([3 5])+divpz([3 5])+dvzdt([3 5])+vdvz([3 5])'};{'Ez','-vxBz([2 4 6])','-divpz([2 4 6])'}};
varstrs = {{'n(1)','n([3 5])','n([1 3 5])'};{'t(1)','t([3 5])'}};

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
twpe = [24000];
xlim = [50 150];
zlim = 0+0.5*[-1 1];

pic = no02m.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bz','Jx'};{'vix','vex','vExBx'};{'Ey','-vexBy','-vixBy'};{'divpy([1 3 5])','divpy([2 4 6])'};{'divpx([1 3 5])','divpx([2 4 6])'};{'PB','pDxx([1 3 5])','pDxx([2 4 6])','pxx([1 3 5])','pxx([2 4 6])'}};
%ylim = {[-1 1]*0.99;[-2 0]*0.99;[-0.1 0.5]*0.99;[-0.1 0.1]*0.99;[-0.1 0.1]*0.99;[0 0.4]*0.99};

varstrs = {{'Bz','Jx'};{'vix','vex','vExBx'};{'Ey','-vixBy','divpy([1 3 5])'};{'Ey','-vexBy','-divpy([2 4 6])'};{'pxy([1 3 5])','pyz([1 3 5])','pxy([2 4 6])','pyz([2 4 6])'}};
%ylim = {[-1 1]*0.99;[-2 0]*0.99;[-0.1 0.4]*0.99;[-0.05 0.4]*0.99;[-0.1 0.1]*0.99};
varstrs = {{'n(1)','n([3 5])'};{'p(1)','p([3 5])'};{'n(1)./n([1 3 5])'};{'p(1)./p([1 3 5])'};{'t(1)','t([1 3 5])'};{'p(1)+p([3 5])','p([1 3 5])'};{'t(1)+t([3 5])','t([1 3 5])'}};

h = pic.plot_line(comp,varstrs,'smooth',10);

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
twpe = [7000 12000];
xlim = [130 207];
xlim = [140 280];
zlim = 0+0.5*[-1 1];
twpe = 24000;

pic = no02m.twpelim(twpe).zlim(zlim);
%pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {{'ni','ne'};{'Ex'};{'Ez'};{'txx([4 6])','tyy([4 6])','tzz([4 6])'};{'txx([3 5])','tyy([3 5])','tzz([3 5])'}};
varstrs = {{'ni','ne','n([1])','n([3 5])','n([4 6])','n([4 6])'};{'Ey'};{'Ez'};{'txx([2 4 6])','tyy([2 4 6])','tzz([2 4 6])'};{'txx([1 3 5])','tyy([1 3 5])','tzz([1 3 5])'}};
varstrs = {{'Bz','By'};{'vix','vex','vExBx'};{'Jx','Jz'};{'Ex','Ez'}};
varstrs = {{'Bz','By'};{'vix','vex','vExBx'};{'Jx','Jz'};{'viy','vey'}};
varstrs = {{'Bx','By','Bz'};{'vix','viy','viz'};{'vex','vey','vez'}};
varstrs = {{'n(1)','n([3 5])'};{'p(1)','p([3 5])'}};
%ylim = {[-1 1]*0.99;[-2 2]*0.99;[-5 5]*0.99};

h = pic.movie_line(comp,varstrs,'ylim',ylim,'filename',[printpath 'B_vi_ve_z=0']);

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
twpe = 7000:1000:12000;
%twpe = 24000;
twpe = 10000:1000:24000;
twpe = [17000 25000];
xlim = no02m.xi([1 end])+[60 -60]';
zlim = [0 8];
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmapjet = colormap('jet');

pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'Ex';'Ez';'Jx';'Jz'};
clims = {[-1 1];[-1 1];0.4*[-1 1];0.4*[-1 1]};
cmaps = {cmapbr;cmapbr;cmapbr;cmapbr};
varstrs = {'vix';'vex';'vey'};
clims = {2*[-1 1];5*[-1 1];5*[-1 1]};
cmaps = {cmapbr;cmapbr;cmapbr};

varstrs = {'ni','log10(ni)'}';
varstrs = {'jepar','vepar','log10(ne)'}';
clims = {[-2 2],[-7 7],[-2 0.5]};
cmaps = {cmapbr,cmapbr,cmapbr};

varstrs = {'ni','log10(ni)'}';
varstrs = {'jepar','vepar','log10(ne)','Eparx','Epary','Eparz'}';
clims = {[-2 2],[-7 7],[-2 0.5],[-1 1],[-1 1],[-1 1]};
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr};

varstrs = {'Ez'}';
clims = {[-1 1]};
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr};

varstrs = {'Epar'}';
clims = {[-1 1]};
cmaps = {cmapbr};


filename = [printpath 'no02m_Epar'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'sep','clim',clims,'filename',filename);
%pic.twpelim([17000 25000]).movie({'Ez'},'A',1,'clim',{[-1 1]},'cmap',{pic_colors('blue_red')},'filename',[printpath 'no02m_Ez']);

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

%% movie, with trajectories
twpe = [18000 24000];
xlim = [60 100];
zlim = [-10 10];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'ni'};
clims = {[0 1];10*[-1 1]};
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmapc2 = pic_colors('candy2');
cmapma = pic_colors('matlab');
cmaps = {cmapbr};
filename = [printpath 'no02m_ni_trajectories'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'filename',filename,'tr',tr100.z0find(4));

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
ds = ds04.dxlim([0.3 1]).twpelim(twpe).zfind(0:1:4).xfind(165:1:176);
%twpe = 7000;
%ds = ds04.dxlim([0 0.3]).twpelim(twpe).zfind(0:1:4).xfind(173:190);
%twpe = 8000;
%ds = ds04.dxlim([0.3 1]).twpelim(twpe).zfind(0:2:4).xfind(165:2:176);

%ds = ds04.dxlim([0.3 1]).twpelim([8000]).zfind(2).xfind(170);
clim = [-2 2];
xlim = 0.99*[-1.5 .5];
ylim = 0.99*[-1 1];
sumdim = 2;
fontsie = 16;
h3 = ds.plot_map([3],sumdim,'bline',df04,'v',df04,'log');
%h3 = ds.plot_map([3],sumdim,'bline',df04,'v',df04,'diff');
%h3 = ds.plot_map([3],sumdim,'ratio',[3 5]);
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
twpe = 24000;
%ds = ds01.twpelim(twpe).zfind([0 1 2 3]).xfind(140:3:210);
%ds = ds100.twpelim(twpe).xlim([80 90]).zlim([2 7]);
ds = ds100.twpelim(twpe).xlim([60 100]).zlim([0 1]);
ds = ds100.twpelim(twpe).findtag({'A=-6'}).xlim([65 100]);
ds = ds100.twpelim(twpe).findtag({'A=-7'}).xlim([70 90]).zlim([-1 4]);
ds = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([70 80]).zlim([1 3]).xfind(73);
pic = no02m.twpelim(twpe).xlim([40 110]).zlim([-1 10]);

sumdim = 1;
%ds = ds04.dxlim([0.3 1]).twpelim([8000]).zfind(2).xfind(170);
clim = [-4 1];
xlim = 3*0.99*[-1 1];
ylim = 3*0.99*[-1 1];
fontsize = 16;
%h = ds.plot_map([1],sumdim,'bline',no02m,'v',no02m,'frac',[1 3 5]); % ,'log'
%h = ds.plot_map([4],sumdim,'bline',no02m,'v',no02m,'frac',[4 6]); % ,'log'
%h = ds.plot_map([3],sumdim,'bline',nobg,'v',nobg,'diff',[5]); % ,'log'
%h = ds.plot_map([3 5],sumdim,'bline',df04n,'v',df04n); % ,'log'
h = ds.plot_map([3 5],sumdim,'bline',no02m,'v',no02m,'log'); % 
%h = ds.plot_map([3],sumdim,'bline',nobg,'v',nobg,'diff',[5],'log'); %
%h = ds.plot_map([3],sumdim,'bline',nobg,'v',nobg,'frac',[5],'log'); % 
hlinks = linkprop(h.ax,{'XLim','YLim','CLim','XTick','YTick'});
%compact_panels(0.00,0.00)
[hax,hlab] = label_panels(h.ax);
%h.ax(1).CLim = clim;
%h.ax(1).CLim = [0 1];
h.ax(1).XLim = xlim;
h.ax(1).YLim = ylim;
%c_eval('hlab(?).FontSize = 12;',1:18)
%c_eval('h.ax(?).FontSize = 14;',1:numel(h.ax))
%c_eval('h.leg(?).FontSize = 12;',1:numel(h.leg))
%colormap([1 1 1; pic_colors('blue_red')])
%colormap([1 1 1; pic_colors('blue_red'); 1 1 1])
%%
figure;
h2 = pic.plot_map({'Ez','vy([3 5])','vExBy','viz-vExBz'}','A',0.5,'clim',{[-1 1],[-1 1],[-1 1],[-1 1]}');
for ih = 1:numel(h2)
  hold(h2(ih),'on')
  ds.plot_boxes(h2(ih))
  hold(h2(ih),'on')
  colormap(h2(ih),pic_colors('blue_red'))
end
%% PICDist.plot_map
twpe = 21000;
%ds = ds01.twpelim(twpe).zfind([0 1 2 3]).xfind(140:3:210);
ds = ds100.twpelim(twpe);

sumdim = 1;
%ds = ds04.dxlim([0.3 1]).twpelim([8000]).zfind(2).xfind(170);
clim = [-3 1.5];
xlim = 2*0.99*[-1 1];
ylim = 2*0.99*[-1 1];
fontsie = 16;
%h = ds.plot_map([3],sumdim,'bline',nobg,'v',nobg,'frac',[3 5]); % ,'log'
%h = ds.plot_map([3],sumdim,'bline',nobg,'v',nobg,'diff',[5]); % ,'log'
%h = ds.plot_map([3 5],sumdim,'bline',df04n,'v',df04n); % ,'log'
h = ds.plot_map([3],sumdim,'bline',no02m,'v',no02m,'frac',[3 5]); % 
%h = ds.plot_map([3],sumdim,'bline',nobg,'v',nobg,'diff',[5],'log'); %
%h = ds.plot_map([3],sumdim,'bline',nobg,'v',nobg,'frac',[5],'log'); % 
hlinks = linkprop(h.ax,{'XLim','YLim','CLim','XTick','YTick'});
compact_panels(0.00,0.00)
[hax,hlab] = label_panels(h.ax);
%h.ax(1).CLim = clim;
%h.ax(1).XLim = xlim;
%h.ax(1).YLim = ylim;
%c_eval('hlab(?).FontSize = 12;',1:18)
%c_eval('h.ax(?).FontSize = 14;',1:numel(h.ax))
%c_eval('h.leg(?).FontSize = 12;',1:numel(h.leg))
%colormap([1 1 1; pic_colors('blue_red')])
colormap([1 1 1; pic_colors('blue_red'); 1 1 1])

%% PICDist.plot_map
twpe = 09000;
ds = ds01.twpelim(twpe).zfind(0:1:4).xfind(150:3:210);
ds = ds01.twpelim(twpe).zfind(0:1:4).xfind(170:1:182);
ds = ds01.twpelim(twpe).zfind(0:1:4).xfind(170:2:182);
%ds = ds01.twpelim(twpe).zfind(0:1:4).xfind(194:1:206);
%twpe = 7000;
%ds = ds04.dxlim([0 0.3]).twpelim(twpe).zfind(0:1:4).xfind(173:190);
%twpe = 8000;
%ds = ds04.dxlim([0.3 1]).twpelim(twpe).zfind(0:2:4).xfind(165:2:176);

%ds = ds04.dxlim([0.3 1]).twpelim([8000]).zfind(2).xfind(170);
clim = [-2 2];
xlim = 0.99*[-1.5 .5];
ylim = 0.99*[-1 1];
sumdim = 2;
fontsie = 16;
%h = ds.plot_map([3 5],sumdim,'bline',nobg,'v',nobg,'log');
%h = ds.plot_map([3],sumdim,'bline',df04,'v',df04,'diff');
h = ds.plot_map([3],sumdim,'ratio',[3 5]);
%h = ds.plot_map([1],sumdim,'forces','EvBy',nobg,'bline',nobg,'v',nobg,'log');
%h = ds.plot_map([3],sumdim,'bline',nobg,'v',nobg,'diff',[5],'log');
%
hlinks = linkprop(h.ax,{'XLim','YLim','CLim'});
compact_panels(0.00,0.00)
[hax,hlab] = label_panels(h.ax);
%h3.ax(1).CLim = clim;
%h3.ax(1).XLim = xlim;
%h3.ax(1).YLim = ylim;
c_eval('hlab(?).FontSize = 14;',1:18)
c_eval('h.ax(?).FontSize = 16;',1:numel(h.ax))
c_eval('h.leg(?).FontSize = 12;',1:numel(h.leg))

cdata = get(findobj(h.ax(1).Children,'type','Image'),'CData');
%h.ax(1).CLim = prctile(abs(cdata(:)),99)*[-1 1];

%% PICDist.reduce_1d, along a line
%ds = ds01.zfind(3);
%fred = ds.reduce_1d_new('x',[5],[]);
twpe = 10000; xlim = [130 210]; zlim = [-15 15];
twpe = 24000; xlim = [50 210]; zlim = [-15 15];

if 0
ds = ds01.twpelim(twpe).zlim(zlim).xlim(xlim);
ds = ds01.twpelim(twpe);
id_line1 = 1:274;
id_line2 = 275:462;
id_line3 = 463:667;
id_line4 = 668:918;
id_line = id_line3;

id_line_tmp = 150;
ds = ds.update_inds({id_line_tmp});
else
  ds = ds100.findtag({'A=-9.5'});
end


xdist = (ds.xi1{1}+ds.xi2{1})/2;
zdist = (ds.zi1{1}+ds.zi2{1})/2;
% arclength = [0 cumsum(sqrt(diff(xdist).^2 + diff(zdist).^2))];

pic_lim = no02m.xlim(xlim).zlim(zlim).twpelim(twpe);
pic = no02m.twpelim(twpe);
Bx_ = pic.Bx;
By_ = pic.By;
Bz_ = pic.Bz;

Bx = interpfield(pic.xi,pic.zi,Bx_,xdist,zdist); 
By = interpfield(pic.xi,pic.zi,By_,xdist,zdist); 
Bz = interpfield(pic.xi,pic.zi,Bz_,xdist,zdist); 

if 0 % saved reduced distributions
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_2.mat')  
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_3.mat')
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_4.mat')
  %save('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_4.mat','fred35_4')
else % make reduced distributions
  %fred1_2 = ds.reduce_1d_new('x',[1],[],'vpar',{Bx,By,Bz});
  %fred3_2 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz});
  %fred4_2 = ds.reduce_1d_new('x',[4],[],'vpar',{Bx,By,Bz});
  %fred5_2 = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz});
  %fred6_2 = ds.reduce_1d_new('x',[6],[],'vpar',{Bx,By,Bz});
  %fred35_2 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz});
  %fred46_2 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz});
 % fred35_3 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz});
  fred35_A10 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_3 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz});
  %fred1_4 = ds.reduce_1d_new('x',[1],[],'vpar',{Bx,By,Bz});
  %fred3_4 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz});
  %fred4_4 = ds.reduce_1d_new('x',[4],[],'vpar',{Bx,By,Bz});
  %fred5_4 = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz});
  %fred6_4 = ds.reduce_1d_new('x',[6],[],'vpar',{Bx,By,Bz});
  %fred35_4 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz});
  %fred46_4 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz});
end
%% Plot
fredi_str = '35'; iSpecies = [3 5];
frede_str = '46'; eSpecies = [4 6];
fredi = eval(['fred' fredi_str '_3']);
frede = eval(['fred' frede_str '_3']);
fred = fredi;
arclength = [0; cumsum(sqrt(diff(fredi.x).^2 + diff(fredi.z).^2))];
if 1; arclength = arclength - arclength(find(abs(fredi.z)==min(abs(fredi.z)))); end
ni = interpfield(pic.xi,pic.zi,pic.ni,fredi.x,fredi.z); 
ne = interpfield(pic.xi,pic.zi,pic.ne,fredi.x,fredi.z); 
vipar = interpfield(pic.xi,pic.zi,pic.vpar(iSpecies),fredi.x,fredi.z); 
vepar = interpfield(pic.xi,pic.zi,pic.vpar(eSpecies),fredi.x,fredi.z); 
vix = interpfield(pic.xi,pic.zi,pic.vx(iSpecies),fredi.x,fredi.z); 
viy = interpfield(pic.xi,pic.zi,pic.vy(iSpecies),fredi.x,fredi.z); 
viz = interpfield(pic.xi,pic.zi,pic.vz(iSpecies),fredi.x,fredi.z); 
vex = interpfield(pic.xi,pic.zi,pic.vx(eSpecies),fredi.x,fredi.z); 
vey = interpfield(pic.xi,pic.zi,pic.vy(eSpecies),fredi.x,fredi.z); 
vez = interpfield(pic.xi,pic.zi,pic.vz(eSpecies),fredi.x,fredi.z); 
vepar = interpfield(pic.xi,pic.zi,pic.vpar(eSpecies),fredi.x,fredi.z); 
Epar = interpfield(pic.xi,pic.zi,pic.Epar,fredi.x,fredi.z); 
Eparx = interpfield(pic.xi,pic.zi,pic.Eparx,fredi.x,fredi.z); 
Epary = interpfield(pic.xi,pic.zi,pic.Epary,fredi.x,fredi.z); 
Eparz = interpfield(pic.xi,pic.zi,pic.Eparz,fredi.x,fredi.z); 
Ex = interpfield(pic.xi,pic.zi,pic.Ex,fredi.x,fredi.z); 
Ey = interpfield(pic.xi,pic.zi,pic.Ey,fredi.x,fredi.z); 
Ez = interpfield(pic.xi,pic.zi,pic.Ez,fredi.x,fredi.z);
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

nrows = 6;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;
doE = 1; colorE = [0 0.8 0.8];
doV = 1; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
isMap = [];

if 1 % line position on map, vy
  isMap(end+1) = isub; 
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic_lim.xi,pic_lim.zi,pic_lim.vy(iSpecies)');
  colormap(hca,pic_colors('blue_red'));
  hcb = colorbar('peer',hca);
  hca.CLim = max(max(get(findobj(hca.Children,'type','Image'),'CData')))*[-1 1];
  hcb.YLabel.String = sprintf('v_{y,%s}',fredi_str);
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  if 1 % plot_boxes
    hold(hca,'on')
    ds.plot_boxes(hca)
    hold(hca,'off')
  end
end
if 1 % line position on map, Ez
  isMap(end+1) = isub; 
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic_lim.xi,pic_lim.zi,pic_lim.Ez');
  colormap(hca,pic_colors('blue_red'));
  hcb = colorbar('peer',hca);
  hca.CLim = max(max(get(findobj(hca.Children,'type','Image'),'CData')))*[-1 1];
  hcb.YLabel.String = 'E_z';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hca.YDir = 'normal';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  if 1 % plot_boxes
    hold(hca,'on')
    ds.plot_boxes(hca)
    hold(hca,'off')
  end
end
if 0 % line position
  hca = h(isub); isub = isub + 1;
  [ax,h1,h2] = plotyy(hca,arclength,fredi.x,arclength,fredi.z);
  hca.XLabel.String = 'arclength (d_i)';
  ax(1).YLabel.String = 'x';
  ax(2).YLabel.String = 'z';
  legend(hca,{'x','z'},'location','best')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % B
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,Bx,arclength,By,arclength,Bz)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'B';
  legend(hca,{'B_x','B_y','B_z'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % n
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,ni,arclength,ne)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'n';
  legend(hca,{'n_i','n_e'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % v
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,vix,arclength,viy,arclength,viz,arclength,vex,arclength,vey,arclength,vez)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'B';
  legend(hca,{'v_{ix}','v_{iy}','v_{iz}','v_{ex}','v_{ey}','v_{ez}'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % v
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,vExBx,arclength,vExBy,arclength,vExBz,arclength,vExBabs)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{ExB}';
  legend(hca,{'v_{x}','v_{y}','v_{z}','|v|'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % E
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,Ex,arclength,Ey,arclength,Ez)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'E';
  legend(hca,{'E_{x}','E_{y}','E_{z}'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Epar, int(Epar)dl
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,Epar,arclength,-cumtrapz(arclength,Epar))  
  legend(hca,{'E_{||}','-\int E_{||}dl_{||}'},'location','eastoutside') 
  if 0
  hold(hca,'on')
  plot(hca,arclength,Eparx,arclength,Epary,arclength,Eparz)  
  hold(hca,'off')
  legend(hca,{'E_{||}','-\int E_{||}dl_{||}','E_{||,x}','E_{||,y}','E_{||,z}'},'location','eastoutside') 
  end
  hca.YLabel.String = 'E_{||}, \int E_{||}dl_{||}'; 
  
  hca.XLabel.String = 'arclength (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % fi(vx)
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
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ex*max(abs(hca.YLim))/max(abs(Ex)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vix,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBx,'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fe(vx)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,frede.v,frede.fvx')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{x}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{e,' frede_str '}(l_{||},v_{x})'];
  hca.CLim(2) = prctile(frede.fvx(:),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ex*max(abs(hca.YLim))/max(abs(Ex)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vix,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBx,'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % fi(vy)
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
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBy,'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % fi(vz)
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
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ez*max(abs(hca.YLim))/max(abs(Ez)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,viz,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBz,'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi(vpar)
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
if 0 % fe(vpar)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,frede.vpar_center,frede.fvpar')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{e,' frede_str '}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(frede.fvpar(:),99);
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
    plot(hca,arclength,vepar,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi(vabs)
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
if 1 % defi35(m*abs^2/2)
  %%
  %isub = 6;
  hca = h(isub); isub = isub + 1;
  darc = arclength(2)-arclength(1);
  arclength_edges = [arclength-darc/2; arclength(end)+darc/2];
  [ARC,EN] = meshgrid(arclength_edges,fredi.vabs_edges.^2/2);
  surf(hca,ARC,EN,ARC*0,log10(fredi.fdefE'))
  shading(hca,'flat')
  view(hca,[0 0 1])
  %pcolor(hca,arclength,fredi.vabs_center.^2/2,log10(fredi.fdefE'))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'mv^2/2';
  colormap(hca,pic_colors('candy'))
  colormap(hca,[ 1 1 1; pic_colors('waterfall')])
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['dpf_{i,' fredi_str '}(l_{||},mv^2/2)'];
  hca.CLim = [prctile(log10(fredi.fdefE(fredi.fdefE>0)),1) prctile(log10(fredi.fdefE(fredi.fdefE>0)),99)];
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YScale = 'log';
  hca.YTick = 10.^(-10:1:10);
  hca.YLim(1) = 10^(-2);
  if doExB
    hold(hca,'on')
    plot(hca,arclength,EExB,'color',colorExB,'linewidth',1.5)
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
if 0 % defe35(m*vabs^2/2)
  %%
  %isub = 6;
  hca = h(isub); isub = isub + 1;
  darc = arclength(2)-arclength(1);
  arclength_edges = [arclength-darc/2; arclength(end)+darc/2];
  [ARC,EN] = meshgrid(arclength_edges,frede.vabs_edges.^2/2/25);
  surf(hca,ARC,EN,ARC*0,log10(frede.fdefE'))
  shading(hca,'flat')
  view(hca,[0 0 1])
  %pcolor(hca,arclength,fredi.vabs_center.^2/2,log10(fredi.fdefE'))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'mv^2/2';
  colormap(hca,pic_colors('candy'))
  colormap(hca,[ 1 1 1; pic_colors('waterfall')])
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['dpf_{e,' frede_str '}(l_{||},mv^2/2)'];
  hca.CLim = [prctile(log10(frede.fdefE(frede.fdefE>0)),1) prctile(log10(frede.fdefE(frede.fdefE>0)),99)];
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YScale = 'log';
  hca.YTick = 10.^(-10:1:10);
  hca.YLim(1) = 10^(-2);
  if doExB
    hold(hca,'on')
    plot(hca,arclength,EExB/25,'color',colorExB,'linewidth',1.5)
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
  darc = arclength(2)-arclength(1);
  arclength_edges = [arclength-darc/2 arclength(end)+darc/2];
  surf(hca,arclength_edges,frede.vabs_edges.^2/2/25,log10(frede.fdefE'))
  view([0 1 0])
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
%
compact_panels(h(setdiff(1:nrows*ncols,isMap)),0.01)
compact_panels(h(isMap),0.01)
h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
%fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(setdiff(1:nrows*ncols,isMap)),{'XLim'});
hlinks = linkprop(h(isMap),{'XLim','YLim'});

%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
for ip = 1:nrows*ncols
  axwidth(ip) = h(ip).Position(3);
end
for ip = 1:nrows*ncols
  h(ip).Position(3) = min(axwidth);
end



%ax(2).YAxisLocation = 'right';
if 0
  %%
  nrows = 3;
  ncols = 1;
  h = setup_subplots(nrows,ncols);
  isub = 1;  
  
  if 1 % line position
    hca = h(isub); isub = isub + 1;
    [ax,h1,h2] = plotyy(hca,arclength,fredi.x,arclength,fredi.z);
    hca.XLabel.String = 'arclength (d_i)';
    ax(1).YLabel.String = 'x';
    ax(2).YLabel.String = 'z';
    legend(hca,{'x','z'},'location','eastoutside')
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 1 % Epar
    hca = h(isub); isub = isub + 1;
    plot(hca,arclength,Eparx,arclength,Epary,arclength,Eparz,arclength,Epar)          
    legend(hca,{'E_{||}','-\int E_{||}dl_{||}','E_{||,x}','E_{||,y}','E_{||,z}'},'location','eastoutside') 
    legend(hca,{'E_{||,x}','E_{||,y}','E_{||,z}','E_{||}'},'location','eastoutside') 
    hca.YLabel.String = 'E_{||}'; 

    hca.XLabel.String = 'arclength (d_i)';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 1 % int Epar
    hca = h(isub); isub = isub + 1;
    plot(hca,arclength,-cumtrapz(arclength,Eparx),arclength,-cumtrapz(arclength,Epary),arclength,-cumtrapz(arclength,Eparz),arclength,-cumtrapz(arclength,Epar))
    %legend(hca,{'E_{||}','-\int E_{||}dl_{||}','E_{||,x}','E_{||,y}','E_{||,z}'},'location','eastoutside') 
    legend(hca,{'E_{||,x}','E_{||,y}','E_{||,z}','E_{||}'},'location','eastoutside') 
    hca.YLabel.String = '-\int E_{||}dl_{||}'; 

    hca.XLabel.String = 'arclength (d_i)';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  compact_panels(0.01)
  h = findobj(get(gcf,'children'),'type','axes'); hall = hall(end:-1:1);
  hlink = linkprop(h,{'XLim'});
  h(1).XLim = arclength([1 end]);
  for ip = 1:nrows*ncols
    axwidth(ip) = h(ip).Position(3);
  end
  for ip = 1:nrows*ncols
    h(ip).Position(3) = min(axwidth);
  end
end
  
%% PICDist.reduce_1d, along a line, A tags
%ds = ds01.zfind(3);
%fred = ds.reduce_1d_new('x',[5],[]);
twpe = 10000; xlim = [130 210]; zlim = [-15 15];
twpe = 24000; xlim = [50 210]; zlim = [-15 15];

ds = ds100.twpelim(twpe).findtag({'A=-6'});
ds = ds100.twpelim(twpe).zfind(0).findtag({'line horizontal'});
ds = ds100.twpelim(twpe).zfind(4).findtag({'line horizontal'});
ds = ds100.twpelim(twpe).findtag({'A=7.5'});

xdist = (ds.xi1{1}+ds.xi2{1})/2;
zdist = (ds.zi1{1}+ds.zi2{1})/2;
% arclength = [0 cumsum(sqrt(diff(xdist).^2 + diff(zdist).^2))];

pic_lim = no02m.xlim(xlim).zlim(zlim).twpelim(twpe);
pic = no02m.twpelim(twpe);
Bx_ = pic.Bx;
By_ = pic.By;
Bz_ = pic.Bz;

Bx = interpfield(pic.xi,pic.zi,Bx_,xdist,zdist); 
By = interpfield(pic.xi,pic.zi,By_,xdist,zdist); 
Bz = interpfield(pic.xi,pic.zi,Bz_,xdist,zdist); 

if 0 % saved reduced distributions
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_2.mat')  
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_3.mat')
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_4.mat')
  %save('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_4.mat','fred35_4')
else % make reduced distributions
  %fred1_2 = ds.reduce_1d_new('x',[1],[],'vpar',{Bx,By,Bz});
  %fred3_2 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz});
  %fred4_2 = ds.reduce_1d_new('x',[4],[],'vpar',{Bx,By,Bz});
  %fred5_2 = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz});
  %fred6_2 = ds.reduce_1d_new('x',[6],[],'vpar',{Bx,By,Bz});
  %fred35_2 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz});
  %fred46_2 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz});
 % fred35_3 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz});
  %fred35_A9 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A9 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_A8 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A8 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_A7 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A7 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_A75 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  fred3_A75 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred5_A75 = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_A4 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A4 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_A5 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A5 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_A6 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_A6 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred35_z4 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_z4 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  
  %fred35_A2 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  %fred46_3 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz});
  %fred1_4 = ds.reduce_1d_new('x',[1],[],'vpar',{Bx,By,Bz});
  %fred3_4 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz});
  %fred4_4 = ds.reduce_1d_new('x',[4],[],'vpar',{Bx,By,Bz});
  %fred5_4 = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz});
  %fred6_4 = ds.reduce_1d_new('x',[6],[],'vpar',{Bx,By,Bz});
  %fred35_4 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz});
  %fred46_4 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz});
end

%% PICDist.reduce_1d, vertical cut
twpe = 5000;
xlim = [190 210];
zlim = [-10 10];
xpick = 204;
ds = ds04n.twpelim(twpe).xfind(xpick);

xdist = (ds.xi1{1}+ds.xi2{1})/2;
zdist = (ds.zi1{1}+ds.zi2{1})/2;
% arclength = [0 cumsum(sqrt(diff(xdist).^2 + diff(zdist).^2))];

pic = df04n.xlim(xlim).zlim(zlim).twpelim(twpe);
pic = df04n.xlim(xlim).zlim(zlim).twpelim(5000);
pic_sm = df04n.xlim([min(ds.xi1{1}) max(ds.xi2{1})]).zlim([min(ds.zi1{1}) max(ds.zi2{1})]).twpelim(5000);
%pic = df04n.twpelim(twpe);
Bx_ = pic.Bx;
By_ = pic.By;
Bz_ = pic.Bz;

Bx = interpfield(pic.xi,pic.zi,Bx_,xdist,zdist); 
By = interpfield(pic.xi,pic.zi,By_,xdist,zdist); 
Bz = interpfield(pic.xi,pic.zi,Bz_,xdist,zdist); 

if 0 % saved reduced distributions
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_2.mat')  
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_3.mat')
  load('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_4.mat')
  %save('/Volumes/Fountain/Data/PIC/no_hot_bg_test/matlab/fred_4.mat','fred35_4')
else % make reduced distributions
  %fred1_2 = ds.reduce_1d_new('x',[1],[],'vpar',{Bx,By,Bz});
  %fred3_2 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz});
  %fred4_2 = ds.reduce_1d_new('x',[4],[],'vpar',{Bx,By,Bz});
  %fred5_2 = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz});
  %fred6_2 = ds.reduce_1d_new('x',[6],[],'vpar',{Bx,By,Bz});
  %fred35_2 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz});
  %fred46_2 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz});
  fred35_204 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz});
  fred46_204 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz});
  fred1_204 = ds.reduce_1d_new('x',[1],[],'vpar',{Bx,By,Bz});
  fred2_204 = ds.reduce_1d_new('x',[2],[],'vpar',{Bx,By,Bz});
%   fred35_196 = ds.reduce_1d_new('z',[3 5],[],'vpar',{Bx,By,Bz});
%   fred46_196 = ds.reduce_1d_new('z',[4 6],[],'vpar',{Bx,By,Bz});
%   fred1_196 = ds.reduce_1d_new('z',[1],[],'vpar',{Bx,By,Bz});
%   fred2_196 = ds.reduce_1d_new('z',[2],[],'vpar',{Bx,By,Bz});
%   fred35_200 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz});
%   fred46_200 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz});
%   fred1_200 = ds.reduce_1d_new('x',[1],[],'vpar',{Bx,By,Bz});
%   fred2_200 = ds.reduce_1d_new('x',[2],[],'vpar',{Bx,By,Bz});
  %fred1_4 = ds.reduce_1d_new('x',[1],[],'vpar',{Bx,By,Bz});
  %fred3_4 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz});
  %fred4_4 = ds.reduce_1d_new('x',[4],[],'vpar',{Bx,By,Bz});
  %fred5_4 = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz});
  %fred6_4 = ds.reduce_1d_new('x',[6],[],'vpar',{Bx,By,Bz});
  %fred35_4 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz});
  %fred46_4 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz});
end
%% Plot
fredi_str = '35'; iSpecies = [3 5];
frede_str = '46'; eSpecies = [4 6];
fredi = eval(['fred' fredi_str '_A75']);
frede = eval(['fred' frede_str '_A75']);
fred = frede;
arclength = [0; cumsum(sqrt(diff(fredi.x).^2 + diff(fredi.z).^2))];
arclength_interp = min(arclength):0.03:max(arclength);
%arclength = arclength(end:-1:1);
%arclength_interp = arclength_interp(end:-1:1);
x_interp = interp1(arclength,fredi.x,arclength_interp);
z_interp = interp1(arclength,fredi.z,arclength_interp);
%if 1; arclength = arclength - arclength(find(abs(fredi.z)==min(abs(fredi.z)))); end
ni = interpfield(pic.xi,pic.zi,pic.n(iSpecies),fredi.x,fredi.z); 
ne = interpfield(pic.xi,pic.zi,pic.n(eSpecies),fredi.x,fredi.z); 
vipar = interpfield(pic.xi,pic.zi,pic.vpar(iSpecies),fredi.x,fredi.z); 
vepar = interpfield(pic.xi,pic.zi,pic.vpar(eSpecies),fredi.x,fredi.z); 
vix = interpfield(pic.xi,pic.zi,pic.vx(iSpecies),fredi.x,fredi.z); 
viy = interpfield(pic.xi,pic.zi,pic.vy(iSpecies),fredi.x,fredi.z); 
viz = interpfield(pic.xi,pic.zi,pic.vz(iSpecies),fredi.x,fredi.z); 
vex = interpfield(pic.xi,pic.zi,pic.vx(eSpecies),fredi.x,fredi.z); 
vey = interpfield(pic.xi,pic.zi,pic.vy(eSpecies),fredi.x,fredi.z); 
vez = interpfield(pic.xi,pic.zi,pic.vz(eSpecies),fredi.x,fredi.z); 
vepar = interpfield(pic.xi,pic.zi,pic.vpar(eSpecies),fredi.x,fredi.z); 
Epar = interpfield(pic.xi,pic.zi,pic.Epar,fredi.x,fredi.z); 
Epar_interp = interpfield(pic.xi,pic.zi,pic.Epar,x_interp,z_interp); 
Eparx = interpfield(pic.xi,pic.zi,pic.Eparx,fredi.x,fredi.z); 
Epary = interpfield(pic.xi,pic.zi,pic.Epary,fredi.x,fredi.z); 
Eparz = interpfield(pic.xi,pic.zi,pic.Eparz,fredi.x,fredi.z); 
Ex = interpfield(pic.xi,pic.zi,pic.Ex,fredi.x,fredi.z); 
Ey = interpfield(pic.xi,pic.zi,pic.Ey,fredi.x,fredi.z); 
Ez = interpfield(pic.xi,pic.zi,pic.Ez,fredi.x,fredi.z);
vExBx = interpfield(pic.xi,pic.zi,pic.vExBx,fredi.x,fredi.z); 
vExBy = interpfield(pic.xi,pic.zi,pic.vExBy,fredi.x,fredi.z); 
vExBz = interpfield(pic.xi,pic.zi,pic.vExBz,fredi.x,fredi.z);
vExBabs = sqrt(vExBx.^2 + vExBy.^2 + vExBz.^2);
EExB = vExBabs.^2/2;
Bx = interpfield(pic.xi,pic.zi,pic.Bx,fredi.x,fredi.z); 
By = interpfield(pic.xi,pic.zi,pic.By,fredi.x,fredi.z); 
Bz = interpfield(pic.xi,pic.zi,pic.Bz,fredi.x,fredi.z); 

%%
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

arclength = arclength-mean(arclength);

nrows = 4;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;
doE = 0; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 0; colorExB = 0*[1 1 1]+0.5;
doPhi = 0; colorPhi = [0.5 0.5 0];

cmap_dist = pic_colors('candy4');

if 0 % overview with box locations, vepar
  hca = h(isub); isub = isub + 1;
  pic.zlim([0 10]).plot_map(hca,{'vepar'});
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  ds.plot_boxes(hca);
  hold(hca,'off')
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = 10*[-1 1];
end
if 1 % overview with box locations, ni
  hca = h(isub); isub = isub + 1;
  pic.zlim([-10 10]).plot_map(hca,{'ni'});
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  ds.plot_boxes(hca);
  hold(hca,'off')
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = [0 1];
end
if 0 % line position
  hca = h(isub); isub = isub + 1;
  [ax,h1,h2] = plotyy(hca,arclength,fredi.x,arclength,fredi.z);
  hca.XLabel.String = 'arclength (d_i)';
  ax(1).YLabel.String = 'x';
  ax(2).YLabel.String = 'z';
  legend(hca,{'x','z'},'location','best')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % B
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,Bx,arclength,By,arclength,Bz)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'B';
  legend(hca,{'B_x','B_y','B_z'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % n
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,ni,arclength,ne)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'n';
  legend(hca,{'n_i','n_e'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % v
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,vix,arclength,viy,arclength,viz,arclength,vex,arclength,vey,arclength,vez)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'B';
  legend(hca,{'v_{ix}','v_{iy}','v_{iz}','v_{ex}','v_{ey}','v_{ez}'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % v
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,vExBx,arclength,vExBy,arclength,vExBz,arclength,vExBabs)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{ExB}';
  legend(hca,{'v_{x}','v_{y}','v_{z}','|v|'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % E
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength,Ex,arclength,Ey,arclength,Ez)
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'E';
  legend(hca,{'E_{x}','E_{y}','E_{z}'},'location','eastoutside')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % Epar, int(Epar)dl
  hca = h(isub); isub = isub + 1;
  plot(hca,arclength_interp,Epar_interp,arclength_interp,-cumtrapz(arclength_interp,Epar_interp))  
  legend(hca,{'E_{||}','-\int E_{||}dl_{||}'},'location','eastoutside') 
  if 0
  hold(hca,'on')
  plot(hca,arclength,Eparx,arclength,Epary,arclength,Eparz)  
  hold(hca,'off')
  legend(hca,{'E_{||}','-\int E_{||}dl_{||}','E_{||,x}','E_{||,y}','E_{||,z}'},'location','eastoutside') 
  end
  hca.YLabel.String = 'E_{||}, \int E_{||}dl_{||}'; 
  
  hca.XLabel.String = 'arclength (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 0 % fi(vx)
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
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ex*max(abs(hca.YLim))/max(abs(Ex)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vix,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBx,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBx,1),'color',colorExB,'linewidth',1.5)
    
    hold(hca,'off')
  end
end
if 0 % fe(vx)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,frede.v,frede.fvx')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{x}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{e,' frede_str '}(l_{||},v_{x})'];
  hca.CLim(2) = prctile(frede.fvx(:),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ex*max(abs(hca.YLim))/max(abs(Ex)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vix,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBx,'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi(vy)
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
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBy,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi(vz)
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
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ez*max(abs(hca.YLim))/max(abs(Ez)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,viz,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBz,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBz,1),'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log10 fi(vx)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.v,log10(fredi.fvx)')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{x}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},v_{x})'];
  hca.CLim(2) = prctile(log10(fredi.fvx(:)),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ex*max(abs(hca.YLim))/max(abs(Ex)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,vix,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBx,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBx,1),'color',colorExB,'linewidth',1.5)
    
    hold(hca,'off')
  end
end
if 1 % log10 fi(vy)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.v,log10(fredi.fvy'))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},v_{y})'];
  hca.CLim(2) = prctile(log10(fredi.fvy(:)),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBy,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 1 % log10 fi(vz)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.v,log10(fredi.fvz'))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{z}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{i,' fredi_str '}(l_{||},v_{z})'];
  hca.CLim(2) = prctile(log10(fredi.fvz(:)),99);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  if doE
    hold(hca,'on')
    plot(hca,arclength,Ez*max(abs(hca.YLim))/max(abs(Ez)),'color',colorE)
    hold(hca,'off')
  end
  if doV
    hold(hca,'on')
    plot(hca,arclength,viz,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if doExB
    hold(hca,'on')
    plot(hca,arclength,vExBz,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBz,1),'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % log 10 fi(vpar)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,fredi.vpar_center,log10(fredi.fvpar)')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{i,' fredi_str '}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(log10(fredi.fvpar(:)),99);
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
if 0 % fi(vpar)
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
if 0 % fe(vpar)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,arclength,frede.vpar_center,log10(frede.fvpar)')
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['f_{e,' frede_str '}(l_{||},v_{||})'];
  hca.CLim(2) = prctile(frede.fvpar(:),99);
  hca.CLim = fe_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = 10*[-1 1];
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
  if doPhi
    intEdl = -cumtrapz(arclength_interp,Epar_interp);
    vshift = 1;
    vphi = 1*sqrt(2*intEdl/(1/pic.mime))-vshift; % 1 is ion mass
    hold(hca,'on')
    plot(hca,arclength_interp,vphi,'color',colorPhi,'linewidth',1.5)
    hold(hca,'off')    
  end  
end
if 0 % fi(vabs)
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
if 0 % defi35(m*abs^2/2)
  %%
  %isub = 6;
  hca = h(isub); isub = isub + 1;
  darc = arclength(2)-arclength(1);
  arclength_edges = [arclength-darc/2; arclength(end)+darc/2];
  [ARC,EN] = meshgrid(arclength_edges,fredi.vabs_edges.^2/2);
  surf(hca,ARC,EN,ARC*0,log10(fredi.fdefE'))
  shading(hca,'flat')
  view(hca,[0 0 1])
  %pcolor(hca,arclength,fredi.vabs_center.^2/2,log10(fredi.fdefE'))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'mv^2/2';
  colormap(hca,pic_colors('candy'))
  colormap(hca,[ 1 1 1; pic_colors('waterfall')])
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['def_{i,' fredi_str '}(l_{||},mv^2/2)'];
  hca.CLim = [prctile(log10(fredi.fdefE(fredi.fdefE>0)),1) prctile(log10(fredi.fdefE(fredi.fdefE>0)),99)];
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YScale = 'log';
  hca.YTick = 10.^(-10:1:10);
  hca.YLim(1) = 10^(-2);
  if doExB
    hold(hca,'on')
    plot(hca,arclength,EExB,'color',colorExB,'linewidth',1.5)
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
if 0 % defe35(m*vabs^2/2)
  %%
  %isub = 6;
  hca = h(isub); isub = isub + 1;
  darc = arclength(2)-arclength(1);
  arclength_edges = [arclength-darc/2; arclength(end)+darc/2];
  [ARC,EN] = meshgrid(arclength_edges,frede.vabs_edges.^2/2/25);
  surf(hca,ARC,EN,ARC*0,log10(frede.fdefE'))
  shading(hca,'flat')
  view(hca,[0 0 1])
  %pcolor(hca,arclength,fredi.vabs_center.^2/2,log10(fredi.fdefE'))
  shading(hca,'flat')
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'mv^2/2';
  colormap(hca,pic_colors('candy'))
  colormap(hca,[ 1 1 1; pic_colors('waterfall')])
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['dpf_{e,' frede_str '}(l_{||},mv^2/2)'];
  hca.CLim = [prctile(log10(frede.fdefE(frede.fdefE>0)),1) prctile(log10(frede.fdefE(frede.fdefE>0)),99)];
  %hca.CLim = fi_clim;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YScale = 'log';
  hca.YTick = 10.^(-10:1:10);
  hca.YLim(1) = 10^(-2);
  if doExB
    hold(hca,'on')
    plot(hca,arclength,EExB/25,'color',colorExB,'linewidth',1.5)
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
  darc = arclength(2)-arclength(1);
  arclength_edges = [arclength-darc/2 arclength(end)+darc/2];
  surf(hca,arclength_edges,frede.vabs_edges.^2/2/25,log10(frede.fdefE'))
  view([0 1 0])
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
%compact_panels(h(2:end),0.01)
h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(2:end),{'XLim'});
hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
for ip = 1:nrows*ncols
  axwidth(ip) = h(ip).Position(3);
end
for ip = 1:nrows*ncols
  h(ip).Position(3) = min(axwidth);
end

linkprop(h(2:4),{'YLim','CLim','XLim'})
linkprop(h(2:4),{'CLim','XLim'})
h(2).CLim = [-4 2];
%h(2).YLim = [-2 2];
%h(3).CLim = [-8 -0.5];
%linkprop(h,{'YLim','XLim'})
%ax(2).YAxisLocation = 'right';
if 0
  %%
  nrows = 3;
  ncols = 1;
  h = setup_subplots(nrows,ncols);
  isub = 1;  
  
  if 1 % line position
    hca = h(isub); isub = isub + 1;
    [ax,h1,h2] = plotyy(hca,arclength,fredi.x,arclength,fredi.z);
    hca.XLabel.String = 'arclength (d_i)';
    ax(1).YLabel.String = 'x';
    ax(2).YLabel.String = 'z';
    legend(hca,{'x','z'},'location','eastoutside')
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 1 % Epar
    hca = h(isub); isub = isub + 1;
    plot(hca,arclength,Eparx,arclength,Epary,arclength,Eparz,arclength,Epar)          
    legend(hca,{'E_{||}','-\int E_{||}dl_{||}','E_{||,x}','E_{||,y}','E_{||,z}'},'location','eastoutside') 
    legend(hca,{'E_{||,x}','E_{||,y}','E_{||,z}','E_{||}'},'location','eastoutside') 
    hca.YLabel.String = 'E_{||}'; 

    hca.XLabel.String = 'arclength (d_i)';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  if 1 % int Epar
    hca = h(isub); isub = isub + 1;
    plot(hca,arclength,-cumtrapz(arclength,Eparx),arclength,-cumtrapz(arclength,Epary),arclength,-cumtrapz(arclength,Eparz),arclength,-cumtrapz(arclength,Epar))
    %legend(hca,{'E_{||}','-\int E_{||}dl_{||}','E_{||,x}','E_{||,y}','E_{||,z}'},'location','eastoutside') 
    legend(hca,{'E_{||,x}','E_{||,y}','E_{||,z}','E_{||}'},'location','eastoutside') 
    hca.YLabel.String = '-\int E_{||}dl_{||}'; 

    hca.XLabel.String = 'arclength (d_i)';  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  compact_panels(0.01)
  h = findobj(get(gcf,'children'),'type','axes'); hall = hall(end:-1:1);
  hlink = linkprop(h,{'XLim'});
  h(1).XLim = arclength([1 end]);
  for ip = 1:nrows*ncols
    axwidth(ip) = h(ip).Position(3);
  end
  for ip = 1:nrows*ncols
    h(ip).Position(3) = min(axwidth);
  end
end
  

%% Plot quantities along a field line
%sep = separatrix_location(nobg);
pic = nobg.twpelim([7500 10400]).xlim([150 210]);
%sep = separatrix_location(pic);
sepEx  = pic.interpfield(sep.x,sep.z,sep.twci,'Ex');
sepEy  = pic.interpfield(sep.x,sep.z,sep.twci,'Ey');
sepEz  = pic.interpfield(sep.x,sep.z,sep.twci,'Ez');
sepvex = pic.interpfield(sep.x,sep.z,sep.twci,'vex');
sepvey = pic.interpfield(sep.x,sep.z,sep.twci,'vey');
sepvez = pic.interpfield(sep.x,sep.z,sep.twci,'vez');
sepvix = pic.interpfield(sep.x,sep.z,sep.twci,'vix');
sepviy = pic.interpfield(sep.x,sep.z,sep.twci,'viy');
sepviz = pic.interpfield(sep.x,sep.z,sep.twci,'viz');
sepni  = pic.interpfield(sep.x,sep.z,sep.twci,'ni');
%%
h = setup_subplots(5,2,'vertical');
isub = 1;
if 0 % z(x)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep.x(:,1),sep.twci,sep.z')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'z';
end
if 1 % Ex
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep.x(:,1),sep.twci,sepEx')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_x';
end
if 1 % Ey
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep.x(:,1),sep.twci,sepEy')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_y';
end
if 1 % Ey
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep.x(:,1),sep.twci,sepEz')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_z';
end
if 1 % vix
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep.x(:,1),sep.twci,sepvix')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ix}';
end
if 1 % viy
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep.x(:,1),sep.twci,sepviy')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{iy}';
end
if 1 % viz
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep.x(:,1),sep.twci,sepviz')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{iz}';
end
if 1 % vex
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep.x(:,1),sep.twci,sepvex')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ex}';
end
if 1 % vey
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep.x(:,1),sep.twci,sepvey')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ey}';
end
if 1 % vez
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep.x(:,1),sep.twci,sepvez')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ez}';
end
if 1 % ni
  hca = h(isub); isub = isub + 1;
  pcolor(hca,sep.x(:,1),sep.twci,sepni')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'n_i';
end

for ip = 1:numel(h)
  hca = h(ip);
  if 1 % z(x) contours
    hold(hca,'on')
    clim = hca.CLim;
    %contour(hca,sep.x(:,1),sep.twci,sep.z',0:1:20,'k')
    contour(hca,sep.x(:,1),sep.twci,smooth2(sepviy,2)',-3:0.1:3,'k')
    hca.CLim = clim;
    hold(hca,'off')    
  end
end

hlinks = linkprop(h,{'XLim','YLim'});
compact_panels(0.01)
h(1).CLim = [-1 1];
h(1).CLim = 0.5*[-1 1];
h(2).CLim = 0.5*[-1 1];
h(3).CLim = 0.5*[-1 1];
h(4).CLim = 0.5*[-1 1];
h(4).CLim = 0.7*[-1 1];
h(4).CLim = 0.5*[-1 1];
h(5).CLim = 0.5*[-1 1];
h(5).CLim = 1*[-1 1];
h(6).CLim = 0.5*[-1 1];
h(7).CLim = 5*[-1 1];
h(8).CLim = 5*[-1 1];
h(9).CLim = 5*[-1 1];
h(10).CLim = [0 0.1];
%% PICDist.plot_map, with something with curvature
twpe = 24000;
%ds = ds01.twpelim(twpe).zfind([0 1 2 3]).xfind(140:3:210);
%ds = ds100.twpelim(twpe).xlim([80 90]).zlim([2 7]);
xpos = [70:2:80];
%Bz = no02m.twpelim(twpe).zlim(0).xlim(xpos).Bz;
mi = 1; qi = 1;
%fciz = Bz*qi/mi;
vvec = 1:1:3;
%rciz = vvec/fciz;

ds = ds100.twpelim(twpe).xlim([60 100]).zlim([0 1]);
ds = ds100.twpelim(twpe).findtag({'A=-6'}).xlim([65 100]);
ds = ds100.twpelim(twpe).findtag({'A=-7'}).xlim([70 90]).zlim([-1 4]);
ds = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([70 80]).zfind(0).xfind(xpos);
pic = no02m.twpelim(twpe).xlim([40 110]).zlim([-1 10]);

sumdim = 3;
%ds = ds04.dxlim([0.3 1]).twpelim([8000]).zfind(2).xfind(170);
clim = [-4 1];
xlim = 4*0.99*[-1 1];
ylim = 4*0.99*[-1 1];
fontsize = 16;
%h = ds.plot_map([1],sumdim,'bline',no02m,'v',no02m,'frac',[1 3 5]); % ,'log'
%h = ds.plot_map([4],sumdim,'bline',no02m,'v',no02m,'frac',[4 6]); % ,'log'
%h = ds.plot_map([3],sumdim,'bline',nobg,'v',nobg,'diff',[5]); % ,'log'
%h = ds.plot_map([3 5],sumdim,'bline',df04n,'v',df04n); % ,'log'
h = ds.plot_map([3 5],sumdim,'v',no02m,'log','curv',{no02m,1}); % 
%h = ds.plot_map([3],sumdim,'bline',nobg,'v',nobg,'diff',[5],'log'); %
%h = ds.plot_map([3],sumdim,'bline',nobg,'v',nobg,'frac',[5],'log'); % 
hlinks = linkprop(h.ax,{'XLim','YLim','CLim','XTick','YTick'});
compact_panels(0.00,0.00)
[hax,hlab] = label_panels(h.ax);
h.ax(1).CLim = clim;
%h.ax(1).CLim = [0 1];
h.ax(1).XLim = xlim;
h.ax(1).YLim = ylim;
%c_eval('hlab(?).FontSize = 12;',1:18)
%c_eval('h.ax(?).FontSize = 14;',1:numel(h.ax))
%c_eval('h.leg(?).FontSize = 12;',1:numel(h.leg))
%colormap([1 1 1; pic_colors('blue_red')])
%colormap([1 1 1; pic_colors('blue_red'); 1 1 1])

%% plotline, curvature
comp = 'x';
twpe = 24000;
xlim = [10 190];
%xlim = [160 260];
zlim = 0.0+0.25*[-1 1];
pic = no02m.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
varstrs = {{'Bz','By','Bx','Babs'};{'curvbx','curvby','curvbz','curvbabs'};{'curvbrad'};{'1./(curvbrad.*Babs)'};{'curvbrad.*Babs'}};

h = pic.plot_line(comp,varstrs);

%% Curvature plot, combined
clear h;
nrows = 4;
ncols = 3;
h(1) = subplot(nrows,ncols,1:3);
h1pos = h(1).Position;
h(2) = subplot(nrows,ncols,4:6);
h2pos = h(2).Position;
c_eval('h(?+2) = subplot(nrows,ncols,?+6);',1:6);
twpe = 24000;
xpos = [70:2:80];
fontsize = 14;
xlim = [40 110];

ds = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([70 80]).zfind(0).xfind(xpos);
pic = no02m.twpelim(twpe).xlim(xlim).zlim([-7 7]);

hca = h(1);
hmap = pic.plot_map(hca,{'log10(curvbrad)'},'A',1,'cmap',pic_colors('blue_red'),'cbarlabels',{'log_{10}r_B'},'smooth',2);
hmap.CLim = [-1 3];
hca.Position = h1pos;
hold(hca,'on')
ds.plot_boxes(hca);
hold(hca,'off')

hca = h(2);
comp = 'x';
hline = pic.xlim([63 xlim(2)]).zlim([-0.25 0.25]).plot_line(hca,comp,{{'Babs','curvbrad','curvbrad.*Babs'}},'smooth',10);
hca.Position = h2pos;
hca.YLim = [0 5];
legend(hca,{'|B|','r_B','r_B\omega_{ci}'},'location','west')
hca.Title.String = [];
hca.XLim = xlim;
hca.Position(2) = hca.Position(2) + 0.08;
compact_panels(h(1:2),0)
c_eval('h(?).Position(2) = h(?).Position(2) + 0.03;',1:2)
irf_legend(hca,{'z=0\pm0.25'},[0.98 0.98],'color',[0 0 0],'fontsize',fontsize)

ih0 = 2;
sumdim = 3;
clim = [-4 1];
xlim = 4*0.99*[-1 1];
ylim = 4*0.99*[-1 1];
hds = ds.plot_map(h(ih0+(1:6)),[3 5],sumdim,'v',no02m,'log','curv',{no02m,1},'nolabel'); % 
hlinks = linkprop(hds.ax,{'XLim','YLim','CLim','XTick','YTick'});
c_eval('hds.ax(?).Position(2) = hds.ax(?).Position(2) + 0.10;',1:3)
compact_panels(hds.ax,0.00,0.00)
%[hax,hlab] = label_panels(hds.ax);
hds.ax(1).CLim = clim;
%h.ax(1).CLim = [0 1];
hds.ax(1).XLim = xlim;
hds.ax(1).YLim = ylim;
c_eval('hds.ax(?).XLabel.String = []; hds.ax(?).XTickLabels = [];',[1 2 3]);
c_eval('hds.ax(?).YLabel.String = []; hds.ax(?).YTickLabels = [];',[2 3 5 6]);

legends = {'a)','b)','c)','d)','e)','f)','g)','h)'};
x0 = (ds.xi1{1}+ds.xi2{1})/2;
c_eval('irf_legend(hds.ax(?),sprintf(''%s x = %g'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',fontsize);',1:6);
%c_eval('hlab(?).FontSize = 12;',1:18)
%c_eval('h.ax(?).FontSize = 14;',1:numel(h.ax))
%c_eval('h.leg(?).FontSize = 12;',1:numel(h.leg))
%colormap([1 1 1; pic_colors('blue_red')])
%colormap([1 1 1; pic_colors('blue_red'); 1 1 1])

irf_legend(hmap,{'a)'},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)
irf_legend(h(2),{'b)'},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)

c_eval('h(?).FontSize = fontsize;',1:numel(h));
c_eval('h(?).Position(1) = h(?).Position(1)-0.04;',1:numel(h));

hlines = findobj(hds.ax(6).Children,'type','line');
hscat = findobj(hds.ax(6).Children,'type','scatter');
hleg=legend([hlines;hscat],{'r_B','v_{ExB}','v_{bulk}'},'location','southeast');

hpos = hds.ax(3).Position;
hcb = colorbar('peer',hds.ax(3));
hds.ax(3).Position = hpos;
hcb.Position(2) = hds.ax(6).Position(2);
hcb.Position(4) = 2*hcb.Position(4);
hcb.YLabel.String = 'log_{10}f(v_x,v_y)';
%% Curvature plot, combined, also including f(vx,vz) at z=2 to illustrate inward vs outward beam temperature
clear h;
nrows = 4;
ncols = 6;
h(1) = subplot(nrows,ncols,1:ncols);
h1pos = h(1).Position;
h(2) = subplot(nrows,ncols,ncols+(1:6));
h2pos = h(2).Position;
c_eval('h(?+2) = subplot(nrows,ncols,?+2*ncols);',1:12);
twpe = 24000;
xpos = [70:2:80];
fontsize = 14;
xlim = [40 110];

ds = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([70 80]).zfind(0).xfind(xpos);
ds2 = ds100.twpelim(twpe).findtag({'line horizontal'}).xlim([70 80]).zfind(2).xfind(xpos);
pic = no02m.twpelim(twpe).xlim(xlim).zlim([-7 7]);

hca = h(1);
hmap = pic.plot_map(hca,{'log10(curvbrad)'},'A',1,'cmap',pic_colors('blue_red'),'cbarlabels',{'log_{10}r_B'},'smooth',2);
hmap.CLim = [-1 3];
hca.Position = h1pos;
hold(hca,'on')
ds.plot_boxes(hca,'color',[1 1 1]);
ds2.plot_boxes(hca,'color',[0 0 0]);
hold(hca,'off')

hca = h(2);
comp = 'x';
hline = pic.xlim([63 xlim(2)]).zlim([-0.25 0.25]).plot_line(hca,comp,{{'Babs','curvbrad','curvbrad.*Babs'}},'smooth',10);
hca.Position = h2pos;
hca.YLim = [0 3.99];
legend(hca,{'|B|','r_B','r_B\omega_{ci}'},'location','west')
hca.Title.String = [];
hca.XLim = xlim;
hca.Position(2) = hca.Position(2) + 0.08;
compact_panels(h(1:2),0)
c_eval('h(?).Position(2) = h(?).Position(2) + 0.03;',1:2)
irf_legend(hca,{'z=0\pm0.25'},[0.98 0.98],'color',[0 0 0],'fontsize',fontsize)

ih0 = 2;
sumdim = 3;
clim = [-4 1];
xlim = 3.5*0.99*[-1 1];
ylim = 3.5*0.99*[-1 1];
hds = ds.plot_map(h(ih0+(1:6)),[3 5],sumdim,'v',no02m,'log','curv',{no02m,1},'nolabel'); % 
hlinks = linkprop(hds.ax,{'XLim','YLim','CLim','XTick','YTick'});
c_eval('hds.ax(?).Position(2) = hds.ax(?).Position(2) + 0.10;',1:6)
hds.ax(1).Position(1) = h(1).Position(1);
compact_panels(hds.ax,0.00,0.00)
c_eval('irf_legend(hds.ax(?),sprintf(''%s x = %g, z = 0'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',fontsize);',1:6);
%[hax,hlab] = label_panels(hds.ax);
hds.ax(1).CLim = clim;
%h.ax(1).CLim = [0 1];
hds.ax(1).XLim = xlim;
hds.ax(1).YLim = ylim;
%c_eval('hds.ax(?).XLabel.String = []; hds.ax(?).XTickLabels = [];',[1 2 3]);
c_eval('hds.ax(?).YLabel.String = []; hds.ax(?).YTickLabels = [];',2:6);
%

ih0 = 8;
hds2 = ds2.plot_map(h(ih0+(1:6)),[3 5],2,'v',no02m,'bline',no02m,'log','nolabel'); % 
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','k)','l)','m)','n)','o)'};
x0 = (ds2.xi1{1}+ds2.xi2{1})/2;
c_eval('irf_legend(hds2.ax(?),sprintf(''%s x = %g, z = 2'',legends{?+ih0},x0(?)),[0.02 0.98],''color'',[0 0 0],''fontsize'',fontsize);',1:6);
hds2.ax(1).Position(1) = h(1).Position(1);
%c_eval('hds2.ax(?).Position(2) = hds2.ax(?).Position(2) + 0.10;',1:3)
compact_panels(hds2.ax,0.00,0.00)
c_eval('hds2.ax(?).YLabel.String = []; hds2.ax(?).YTickLabels = [];',2:6);
hds2.ax(1).CLim = clim;
%h.ax(1).CLim = [0 1];
hds2.ax(1).XLim = xlim;
hds2.ax(1).YLim = ylim;

irf_legend(hmap,{'a)'},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)
irf_legend(h(2),{'b)'},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)

c_eval('h(?).FontSize = fontsize;',1:numel(h));
c_eval('h(?).Position(1) = h(?).Position(1)-0.04;',1:numel(h));
%%
hlines = findobj(hds.ax(6).Children,'type','line');
hscat = findobj(hds.ax(6).Children,'type','scatter');
hleg=legend([hlines;hscat],{'r_B','v_{ExB}','v_{bulk}'},'location','southeast');

hpos = hds.ax(3).Position;
hcb = colorbar('peer',hds.ax(3));
hds.ax(3).Position = hpos;
hcb.Position(2) = hds.ax(6).Position(2);
hcb.Position(4) = 2*hcb.Position(4);
hcb.YLabel.String = 'log_{10}f(v_x,v_y)';