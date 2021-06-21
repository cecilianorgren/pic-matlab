% compare function_load_data to PIC
twpe = 20000;

no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
pic = no02m.twpelim(twpe);

filepath = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data/fields-%5.0f.dat',twpe);
data = function_load_data(filepath);

%% Compare data from PIC and function_load_data

h = setup_subplots(5,2);
isub = 1;

if 1 % Bz pic
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.Bz');
  hca = colorbar('peer',hca);
end
if 1 % Bz dataverse
  hca = h(isub); isub = isub + 1;
  imagesc(hca,data.x_di,data.z_di,data.Bz')
  hca = colorbar('peer',hca);
end

if 1 % Ey pic
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.Ey');
  hca = colorbar('peer',hca);
end
if 1 % Ey dataverse
  hca = h(isub); isub = isub + 1;
  imagesc(hca,data.x_di,data.z_di,data.Ey')
  hca = colorbar('peer',hca);
end

if 1 % n1 pic
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.n(1)');
  hca = colorbar('peer',hca);
end
if 1 % n1 dataverse
  hca = h(isub); isub = isub + 1;
  imagesc(hca,data.x_di,data.z_di,data.hot_ion.n')
  hca = colorbar('peer',hca);
end

if 0 % vx(1) pic
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.vx(1)');
  hca = colorbar('peer',hca);
end
if 0 % vx(1) dataverse
  hca = h(isub); isub = isub + 1;
  imagesc(hca,data.x_di,data.z_di,data.hot_ion.vx')
  hca = colorbar('peer',hca);
end
if 1 % vx(4) pic
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.vx(4)');
  hca = colorbar('peer',hca);
end
if 1 % vx(4) dataverse
  hca = h(isub); isub = isub + 1;
  imagesc(hca,data.x_di,data.z_di,data.cold_ele_top.vx')
  hca = colorbar('peer',hca);
end

if 0 % Pxx(1) pic
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.pxx(1)');
  hca = colorbar('peer',hca);
end
if 0 % Pxx(1) dataverse
  hca = h(isub); isub = isub + 1;
  imagesc(hca,data.x_di,data.z_di,data.hot_ion.pxx')
  hca = colorbar('peer',hca);
end

if 1 % Pxx(4) pic
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.pxx(4)');
  hca = colorbar('peer',hca);
end
if 1 % Pxx(4) dataverse
  hca = h(isub); isub = isub + 1;
  imagesc(hca,data.x_di,data.z_di,data.cold_ele_top.pxx')
  hca = colorbar('peer',hca);
end

