fid = H5F.open(filePath,'H5F_ACC_RDWR','H5P_DEFAULT');
H5L.delete(fid,'/data/0000019200/dt','H5P_DEFAULT');
H5F.close(fid);
