pic = no02m;
trs = tr100;
filePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories.h5';

for itr = 1:trs.ntr
  tic
  Atmp = no02m.interp(trs(itr).x,trs(itr).z,trs(itr).t,'A');
  tr_id = trs(itr).id;
  h5write_trajs_add(filePath,tr_id,'Ay',Atmp)
  trsA(itr).A = Atmp;
  toc
end

