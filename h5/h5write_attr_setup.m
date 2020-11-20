%df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
pic = no02m;
missingAttr = pic.get_missing_attributes('RE');
%pic = pic.twcilim(pic.twci(missingAttr(1:end-2)),'exact');
%pic = pic.twcilim(pic.twci(missingAttr),'exact');
pic = pic.twcilim(pic.twci([missingAttr(1)-2:missingAttr(1)+2]),'exact');
%pic = no02m;
%pic = pic.twpelim(200:200:6000);
%pic = turb.twcilim(1:1:30,'exact');

%% Energy partitioning, UB, UK, UT
times = pic.twci;
for it = 1:pic.nt
  pic_tmp = pic.twcilim(times(it));
  Bx = pic_tmp.Bx;
  By = pic_tmp.By;
  Bz = pic_tmp.Bz;
  Babs = sqrt(Bx.^2+By.^2+Bz.^2);
  UB(it) = sum(Babs(:).^2)/2;  
  disp(sprintf('%g/%g',it,pic.nt))
end

%% X line position, Ey at X line, A at X line

xXlineAll = [];
zXlineAll = [];
EyXlineAll = [];
xXline = [];
zXline = [];
Aval = [];
Aall = nan(pic.nx,pic.nz,pic.nt);
plot(nan,nan); hold on
times = pic.twci;

for it = 1:pic.nt
  pic_tmp = pic.twcilim(times(it));
  if 0
    Bx = pic_tmp.Bx;
    Bz = pic_tmp.Bz;
    if 1 % also magnetic pressure
      By = pic_tmp.By;
      Babs = sqrt(Bx.^2+By.^2+Bz.^2);
      UB(it) = sum(Babs(:).^2)/2;  
    end
    A = vector_potential(pic_tmp.xi,pic_tmp.zi,Bx,Bz);
    try
    h5write_fields_ancillary(pic_tmp,pic_tmp.twpe,'A',A)  
    catch
      disp('Could not write A.')
    end
  else
    A = pic_tmp.A;
  end
  Aall(:,:,it) = A;
  
  [Ainds,Avals] = saddle(A,'sort');
  xXline(it) = pic_tmp.xi(Ainds(1,1));
  zXline(it) = pic_tmp.zi(Ainds(1,2));  
  Aval(it) = Avals(1);
  EyXline(it) = mean(mean(pic_tmp.xlim(xXline(it)+[-0.1 0.1]).zlim(zXline(it)+[-0.1 0.1]).Ey));
  for iX = 1:numel(Avals)
    xXlineAll = [xXlineAll pic_tmp.xi(Ainds(iX,1))];
    zXlineAll = [zXlineAll pic_tmp.zi(Ainds(iX,2))];
    EyXlineAll = [EyXlineAll mean(mean(pic_tmp.xlim(xXlineAll(end)+[-0.1 0.1]).zlim(zXlineAll(end)+[-0.1 0.1]).Ey))];
    scatter(pic_tmp.twci,xXlineAll(end),abs(EyXlineAll(end))*100+1,iX)
    drawnow
  end
  disp(sprintf('%g',it))
end

dA = diff(Aval(1:pic.nt));
dt = diff(times(1:pic.nt));
dAdt_ = dA./dt;
dAdt = interp1(times(1:end-1)+0.5*dt,dAdt_,times);

% for it = 1:pic.nt
%   pic_tmp = pic.twcilim(times(it));
%   EyXline05(it) = mean(mean(pic_tmp.xlim(xXline(it)+0.5*[-1 1]).zlim(zXline(it)+0.5*[-1 1]).Ey));
% end

%% UB
times = pic.twci;
clear UB
for it = 1:pic.nt
  pic_tmp = pic.twcilim(times(it));

  Bx = pic_tmp.Bx;
  Bz = pic_tmp.Bz;
  By = pic_tmp.By;
  Babs = sqrt(Bx.^2+By.^2+Bz.^2);
  UB(it) = sum(Babs(:).^2)/2;  
  disp(sprintf('%g',it))
end


%% Write attributes
indsave = 1:(pic.nt-3);
indsave = (pic.nt-2):pic.nt;
indsave = 1:pic.nt;
indsave = 2:(pic.nt-1);
%indsave = 1:pic.nt;
h5write_attr(pic.subset('t',indsave),times(indsave),'RE',EyXline(indsave))
%h5write_attr(pic.subset('t',indsave),times(indsave),'RA',dAdt(indsave))
h5write_attr(pic.subset('t',indsave),times(indsave),'UB',UB(indsave))
h5write_attr(pic.subset('t',indsave),times(indsave),'xline_position',[xXline(indsave)' zXline(indsave)'])

%% Write ancillary data (not attributes), for example A
pic = no02m;
timesteps = pic.twpelim([15100:100:15900 16100:100:16900],'exact').twpe;

for time = timesteps
  pic_tmp = pic.twpelim(time);
  Bx = pic_tmp.Bx;
  Bz = pic_tmp.Bz;
  A = vector_potential(pic_tmp.xi,pic_tmp.zi,Bx,Bz);
  imagesc(A'); colorbar; pause(0.1)
  h5write_fields_ancillary(pic_tmp,pic_tmp.twpe,'A',A)  
end