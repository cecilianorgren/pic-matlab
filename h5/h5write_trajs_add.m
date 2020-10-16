function h5write_trajs_add(filePath,tr_id,field,data)
% H5WRITE_TRAJS_ADD Update part of data h5 file.
% H5WRITE_TRAJS_ADD(h5FilePath,tr_id,field,data)
% 
% dirData - directory of data
% h5FilePath - directory and file name
% 

iTr_str = sprintf('%06.0f',tr_id);
group_name = ['/traj/' iTr_str '/'];   
h5write(filePath, [group_name field], data)
    