filepath = '/Users/cecilia/Discs/tesla/software/Smilei4.5/Smilei/GEMchallenge/Fields0.h5';
namelist = '/Users/cecilia/Discs/tesla/software/Smilei4.5/Smilei/GEMchallenge/Harris_new.py';

filepath = '/Users/cecilia/Discs/tesla/software/Smilei4.5/Smilei/3dtest/Fields0.h5';
namelist = '/Users/cecilia/Discs/tesla/software/Smilei4.5/Smilei/3dtest/Harris3d.py';

filepath = '/Users/cno062/tesla/software/Smilei4.5/Smilei/3dtest/Fields0.h5';
namelist = '/Users/cno062/tesla/software/Smilei4.5/Smilei/3dtest/Harris3d.py';


%pic = PIC3D(filepath,namelist); %test
sm = PIC3D(filepath,namelist); %test

%% Examine struture of h5 file
info = h5info(filepath);
%%
ix = 390:410;
iy = 290:310;
iz = 100;
%Bz = pic.xgrid(ix).ygrid(iy).zgrid(iz).twpelim(pic(1).twpe).Bz;
tic
pic_tmp = pic.xgrid(ix).ygrid(iy).zgrid(iz).twpelim(pic(1).twpe);

%% 1D cut
ix = fix(pic.nx/2);
iy = 1:pic.ny;
iz = 100;
%pic_tmp = pic.xgrid(ix).ygrid(iy).zgrid(iz);
pic_tmp = pic.xgrid(ix).zgrid(iz);
tic;
Bx = pic_tmp.Bx;
By = pic_tmp.By;
Bz = pic_tmp.Bz;
Jx = pic_tmp.Jx;
Jy = pic_tmp.Jy;
Jz = pic_tmp.Jz;
Ex = pic_tmp.Ex;
Ey = pic_tmp.Ey;
Ez = pic_tmp.Ez;
toc
%By = pic.zgrid(iz).By;
plot(pic_tmp.yi,squeeze(Bx),pic_tmp.yi,squeeze(By),pic_tmp.yi,squeeze(Bz),...
     pic_tmp.yi,squeeze(Jx),pic_tmp.yi,squeeze(Jy),pic_tmp.yi,squeeze(Jz))
legend({'B_x','B_y','B_z','J_x','J_y','J_z'})

xlabel('y')
%ylabel('B_x')
title(sprintf('x = %g, z = %g',pic_tmp.xi,pic_tmp.zi))
grid on
%imagesc(pic_tmp.xi,pic_tmp.zi,squeeze(Bx)')
%colorbar

%% 2D plot, reconnection plane
ix = 1:sm.nx;
iy = 1:sm.ny;
iz = 100;
%pic_tmp = pic.xgrid(ix).ygrid(iy).zgrid(iz);
pic = sm.xgrid(ix).ygrid(iy).zgrid(iz).tgrid(2);

tic;
Bx = pic.Bx;
By = pic.By;
Bz = pic.Bz;
Jx = pic.Jx;
Jy = pic.Jy;
Jz = pic.Jz;
Ex = pic.Ex;
Ey = pic.Ey;
Ez = pic.Ez;
toc
 
%% 2D plot, equatorial plane
ix = 1:sm.nx;
iy = fix(sm.ny/2);
iz = 1:sm.nz;
it = 2;
%pic_tmp = pic.xgrid(ix).ygrid(iy).zgrid(iz);
pic = sm.xgrid(ix).ygrid(iy).zgrid(iz).tgrid(it);
%%
tic;
Bx = pic.Bx;
By = pic.By;
Bz = pic.Bz;
Jx = pic.Jx;
Jy = pic.Jy;
Jz = pic.Jz;
Ex = pic.Ex;
Ey = pic.Ey;
Ez = pic.Ez;
toc