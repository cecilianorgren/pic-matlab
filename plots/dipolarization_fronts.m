no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

%% plot_line, horizontal
pic = no02m;
comp = 'x';
twpe = 23000;
xlim = pic.xi([1 end])+[0 -0]';
zlim = 1*[-0.5 0.5]+0;
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
varstrs = {{'Bz','ni'};{'Bz','By','Bz'};{'vix';'viy';'viz'};{'Ex','Ey','Ez'}};


h = pic.plot_line(comp,varstrs,'smooth',10); %

%% movie_line
comp = 'x';
twpe = 15000:100:25000;
zlim = 0+1*[-0.5 0.5];
xlim = no02m.xi([1 end])+[50 -50]';

pic = no02m.twpelim(twpe,'exact').xlim(xlim).zlim(zlim);
%pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {{'ni','ne'};{'Ex'};{'Ez'};{'txx([4 6])','tyy([4 6])','tzz([4 6])'};{'txx([3 5])','tyy([3 5])','tzz([3 5])'}};
varstrs = {{'ni','ne','n([1])','n([3 5])','n([4 6])','n([4 6])'};{'Ey'};{'Ez'};{'txx([2 4 6])','tyy([2 4 6])','tzz([2 4 6])'};{'txx([1 3 5])','tyy([1 3 5])','tzz([1 3 5])'}};
varstrs = {{'Bz','By'};{'vix','vex','vExBx'};{'Jx','Jz'};{'Ex','Ez'}};
varstrs = {{'Bz','By'};{'vix','vex','vExBx'};{'Jx','Jz'};{'viy','vey'}};
varstrs = {{'Bx','By','Bz'};{'vix','viy','viz'};{'vex','vey','vez'}};
varstrs = {{'Bz','ni','vix','Ey'}};
%ylim = {[-1 1]*0.99;[-2 2]*0.99;[-5 5]*0.99};

h = pic.movie_line(comp,varstrs,'ylim',{[-1.5 1.5]},'filename',[printpath 'ni_Bz_vix_Ey_zoom']);

%% plot_line, vertical
pic = no02m;
comp = 'x';
twpe = 0;
xlim = [10 11];
zlim = [-10 10];
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
varstrs = {{'Bz','ni'};{'Bz','By','Bz'};{'vix';'viy';'viz'};{'Ex','Ey','Ez'}};


h = pic.plot_line(comp,varstrs,'smooth',10); %


