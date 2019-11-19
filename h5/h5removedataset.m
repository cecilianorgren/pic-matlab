filePath = sim.file;
fid = H5F.open(filePath,'H5F_ACC_RDWR','H5P_DEFAULT');
H5L.delete(fid,'/scalar_timeseries/U/T/1','H5P_DEFAULT');
H5F.close(fid);
