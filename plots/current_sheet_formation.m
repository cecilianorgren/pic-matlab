% Plot pressure components etcduring early stages.
twci = [20 80];
xc = (no02m.xi(end)-no02m.xi(1))/2;
xlim = xc + 0.5*[-1 1];
zlim = [-3 3];
varstrs = {'pxx([3 5])','pyy([3 5])','pzz([3 5])',...
           'pxx([4 6])','pyy([4 6])','pzz([4 6])'}';

varstrs = {'pxx([1])','pyy([1])','pzz([1])',...
           'pxx([2])','pyy([2])','pzz([2])'}';

varstrs = {'pxy([1])','pyz([1])','pxz([1])',...
           'pxy([3 5])','pyz([3 5])','pxz([3 5])'}';
         

pic = no02m.twcilim(twci).xlim(xlim).zlim(zlim);         
h = pic.plot_timemap('tz',varstrs,'A',0.5);