%% Harris plasma sheet + cold inflow, without hot uniform background
% Harris plasma sheet + cold inflow, without hot uniform background
% Without the hot background, there is vaccuum in the center with this 
% setup. Should there be a hot background that is limited in z? Let's see
% how it evolves...
% test3: The temperature is correct now, but the pressure seems off.
% test4: Might have been a factor 1/4 missing in the perturbations, both
%        Bx, Bz, and Jy

% nobg = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_test/data_h5/fields.h5');
xlim = [100 102];
zlim = [-10 10];  
pic = nobg.twpelim(440).xlim(xlim).zlim(zlim);
nrows = 4;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;

if 0 % ni
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.zi,squeeze(mean(pic.ni,1)));
  hca.YLabel.String = 'n_i';
  hca.XLabel.String = 'z';
end
if 1 % Bx
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.zi,squeeze(mean(pic.Bx,1)));
  hca.YLabel.String = 'B_x';
  hca.XLabel.String = 'z';
end
if 0 % p
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.zi,squeeze(mean(pic.p([1 3 5]),1)));
  hca.YLabel.String = 'p';
  hca.XLabel.String = 'z';
end
if 1 % T
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.zi,squeeze(mean(pic.t([1]),1)),pic.zi,squeeze(mean(pic.t([2]),1)),pic.zi,squeeze(mean(pic.t([1]),1))+squeeze(mean(pic.t([2]),1)));
  hca.YLabel.String = 't';
  hca.XLabel.String = 'z';
end
if 1 % T
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.zi,squeeze(mean(pic.t([3]),1)),pic.zi,squeeze(mean(pic.t([4]),1)));
  hca.YLabel.String = 't';
  hca.XLabel.String = 'z';
end
if 1 % pT and pB
  hca = h(isub); isub = isub + 1;
  plot(hca,pic.zi,squeeze(mean(pic.p([1]),1))+squeeze(mean(pic.p([2]),1)),'-',pic.zi,squeeze(mean(pic.PB,1)),'--');
  hca.YLabel.String = 'p';
  hca.XLabel.String = 'z';
end

hlinks = linkprop(h,{'XLim'});
