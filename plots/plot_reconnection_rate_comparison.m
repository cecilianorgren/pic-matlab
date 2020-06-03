%reconnection_rate_comparison
%df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
%df08 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n08/data_h5/fields.h5');


hca = subplot(1,1,1);


dA = diff(recflux);
dt = diff(time);

RA = dA./dt;
time_ = time(2:end)-0.5*dt(1);

plot(hca,df08.twci,df08.RA,df04.twci,df04.RA,time_,RA,'linewidth',1.5)

hca.YLim(1) = 0;
legend('n_c = 0.8 n_0','n_c = 0.4 n_0','n_c = 0.4 n_0 (initialized further out)','location','best')
grid on
hca.XLabel.String = 'time (\omega_{ci}^{-1})';
hca.YLabel.String = 'Reconnection rate (v_AB_0)';
hca.FontSize = 12;