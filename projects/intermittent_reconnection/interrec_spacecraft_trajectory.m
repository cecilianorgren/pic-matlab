%no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');


t = no02m.twci;
x = 125*ones(numel(t),1);
z = 3*ones(numel(t),1);

varstrs = {'Ez','vix','ni','n(1)','n([3 5])'};
for ivar = 1:numel(varstrs)
  varstr = varstrs{ivar};
  vars{ivar} = eval('no02m.interp(x,z,t,varstr);');
end
hca = subplot(1,1,1);
plot(hca,t,cat(2,vars{:}))
hca.XLabel.String = 't\omega_{ci}';
legend(hca,varstrs,'location','best')

% Ez = no02m.interp(x,z,t,'Ez');
% ni = no02m.interp(x,z,t,'ni');
% vix = no02m.interp(x,z,t,'vix');
% nih = no02m.interp(x,z,t,'n(1)');
% nic = no02m.interp(x,z,t,'n([3 5])');
% plot(t,Ez,t,ni,t,vix,t,nih,t,nic)

%%