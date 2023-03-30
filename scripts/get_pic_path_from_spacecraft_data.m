%% Load spacecraft data
units = irf_units;
irf.log('critical')
ic = 1;

localuser = datastore('local','user');
%mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
mms.db_init('local_file_db',['/Volumes/DataRaid/MMS']);
db_info = datastore('mms_db');   

% Time from time interval
tint = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:01:43.00Z');
tint_action = irf.tint('2017-07-25T22:09:30.00Z/2017-07-25T22:11:00.00Z');

c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);

%% Load PIC object
no02m = PIC('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/fields.h5');

%% Cost function
pic_params = {'Ez','By','ne'};
mms_params = {'Ez','By','ne'};