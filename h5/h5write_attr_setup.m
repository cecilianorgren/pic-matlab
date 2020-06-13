%df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');

pic = bs;
%% Energy partitioning, UB, UK, UT
times = pic.twci;
for it = 1:pic.nt
  pic_tmp = pic.twcilim(times(it));
  Bx = pic_tmp.Bx;
  By = pic_tmp.By;
  Bz = pic_tmp.Bz;
  Babs = sqrt(Bx.^2+By.^2+Bz.^2);
  UB(it) = sum(Babs(:).^2)/2;  
  disp(sprintf('%g',it))
end

%% X line position, Ey at X line, A at X line

times = pic.twci;
for it = 1:pic.nt
  pic_tmp = pic.twcilim(times(it));
  Bx = pic_tmp.Bx;
  Bz = pic_tmp.Bz;
  A = vector_potential(pic_tmp.xi,pic_tmp.zi,Bx,Bz);
  %A = pic_tmp.A;
  [Ainds,Avals] = saddle(A,'sort');
  xXline(it) = pic_tmp.xi(Ainds(1,1));
  zXline(it) = pic_tmp.zi(Ainds(1,2));
  Aval(it) = Avals(1);
  EyXline(it) = mean(mean(pic_tmp.xlim(xXline(it)+[-0.1 0.1]).zlim(zXline(it)+[-0.1 0.1]).Ey));
  disp(sprintf('%g',it))
end

dA = diff(Aval);
dt = diff(times);
dAdt_ = dA./dt;
dAdt = interp1(times(1:end-1)+0.5*dt,dAdt_,times);

% for it = 1:pic.nt
%   pic_tmp = pic.twcilim(times(it));
%   EyXline05(it) = mean(mean(pic_tmp.xlim(xXline(it)+0.5*[-1 1]).zlim(zXline(it)+0.5*[-1 1]).Ey));
% end

%% Write attributes
h5write_attr(pic,times,'RE',EyXline)
h5write_attr(pic,times,'RA',dAdt)
h5write_attr(pic,times,'UB',UB)
h5write_attr(pic,times,'xline_position',[xXline' zXline'])

%% Write ancillary data (not attributes), for example A
pic = bs;
timesteps = pic.twpe;
for time = timesteps(2:end)
  pic_tmp = pic.twpelim(time);
  Bx = pic_tmp.Bx;
  Bz = pic_tmp.Bz;
  A = vector_potential(pic_tmp.xi,pic_tmp.zi,Bx,Bz);
  imagesc(A'); colorbar; pause(0.1)
  h5write_fields_ancillary(pic_tmp,pic_tmp.twpe,'A',A)
  
end