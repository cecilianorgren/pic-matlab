filepath = '/Users/cecilia/Discs/tesla/software/Smilei4.5/Smilei/GEMchallenge/Fields0.h5';
namelist = '/Users/cecilia/Discs/tesla/software/Smilei4.5/Smilei/GEMchallenge/Harris_new.py';

filepath = '/Users/cecilia/Discs/tesla/software/Smilei4.5/Smilei/3dtest/Fields0.h5';
namelist = '/Users/cecilia/Discs/tesla/software/Smilei4.5/Smilei/GEMchallenge/Harris3d.py';


pic = PIC3D(filepath,namelist);

%%
ix = 550:650;
iy = 350:450;
iz = 5;
Bz = pic.xgrid(ix).ygrid(iy).zgrid(iz).twpelim(pic(1).twpe).Bz;