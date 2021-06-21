% compare function_load_data to PIC
twpe = 20000;

no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
pic = no02m.twpelim(twpe);

filepath = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data/fields-%5.0f.dat',twpe);
data = function_load_data(filepath);

%% Compare data from PIC and function_load_data

h = setup_subplots(4,2);
isub = 1;

if 1 % n pic
  hca = h(isub); isub = isub + 1;
  imagesc(hca,no02m.)
end

if 1 % n pic

end
