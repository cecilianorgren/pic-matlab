%% Write fields
h5write_fields('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data/','/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields_test.h5',[23000 24000],6)

%% Load pic object
pic = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields_test.h5');

%% Write complementary information to h5 file
h5write_fields_complement(pic)

%% Test what was written
% Load pic object again
pic = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields_test.h5');

%% Plot
npanels = 4;
for ipanel = 1:npanels
  h(ipanel) = subplot(npanels,1,ipanel);
end
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,pic.twci,pic.UB,pic.twci,pic.Uki,pic.twci,pic.Uti,pic.twci,pic.Uke,pic.twci,pic.Ute)
hca.YLabel.String = 'Energies';
legend(hca,{'U_B','U_{ki}','U_{ti}','U_{ke}','U_{te}'})


hca = h(isub); isub = isub + 1;
plot(hca,pic.twci,pic.Axline)
hca.YLabel.String = 'X line A value';
legend(hca,{'A_X'})

hca = h(isub); isub = isub + 1;
plot(hca,pic.twci,pic.xline_position)
hca.YLabel.String = 'X line position';
legend(hca,{'x_X','z_X'})

hca = h(isub); isub = isub + 1;
dA = diff(pic.Axline);
dt = diff(pic.twci);
RA = -dA./dt;
plot(hca,pic.twci,pic.RE,pic.twci(1:end-1)+0.5*dt,RA,'*')
hca.YLabel.String = 'R';
legend(hca,{'R_E','R_A'})