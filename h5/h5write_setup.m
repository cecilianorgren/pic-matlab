%% Write basic data
%h5filepath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5';
%datapath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data/';
h5filepath = '/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/fields.h5';
datapath = '/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data/';
nSpecies = 6;

h5filepath = '/Volumes/DataRaid/cno062/rec_onset_4/data_h5/fields_E01.h5';
datapath = '/Volumes/DataRaid/cno062/rec_onset_4/data_F1_E01/';
nSpecies = 4;

h5filepath = '/Volumes/DataRaid/cno062/susanne_varying_varying/fields.h5';
datapath = '/Volumes/DataRaid/cno062/susanne_varying_varying/data/';
nSpecies = 4;


h5filepath = '/Volumes/DataRaid/cno062/susanne_varying_varying/fields_varying.h5';
datapath = '/Volumes/DataRaid/cno062/susanne_varying_varying/data/';
nSpecies = 4;


h5filepath = '/Volumes/DataRaid/cno062/susanne_varying_varying/fields_varying.h5';
datapath = '/Volumes/DataRaid/cno062/susanne_varying_varying/data/';
nSpecies = 4;

h5filepath = '/Users/cecilia/Discs/spis/data/varying_tite/tite_05/data_h5/fields4.h5';
datapath = '/Users/cecilia/Discs/spis/data/varying_tite/tite_05/data/';
nSpecies = 4;

h5filepath = '/Users/cecilia/Data/PIC/varying_tite/data_F1_E05/fields_newdata.h5';
datapath = '/Users/cecilia/Data/PIC/varying_tite/data_F1_E05/';
nSpecies = 4;


% Susannes rec-onset
% h5filepath = '/Volumes/DataRaid/Susanne-onset/data_h5/fields.h5';
% datapath = '/Volumes/DataRaid/Susanne-onset/data/';
% nSpecies = 2;



%h5filepath = '/Volumes/Fountain/Data/PIC/varying_tite/E05/fields.h5';
%datapath = '/Volumes/DataRaid/cno062/rec_onset_4/data_F1_E05';
%nSpecies = 4;





h5filepath = '/Users/cno062/Data/PIC/rec_onset_4/data_h5/fields_F025_E005_TITE10.h5';
datapath = '/Users/cno062/Data/PIC/rec_onset_4/data_F025_E005_TITE10/';
nSpecies = 4;


h5filepath = '/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5';
datapath = '/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data/';
nSpecies = 6;

%h5filepath = '/Users/cno062/Data/PIC/varying_tite/tite_05/fields.h5';
%datapath = '/Users/cno062/Data/PIC/rec_onset_4/data_F025_E005_TITE10/';
%nSpecies = 4;


%h5write_fields(datapath,h5filepath,3000:1000:20000,nSpecies)
%h5write_fields(datapath,h5filepath,10700:100:10900,nSpecies)
h5write_fields(datapath,h5filepath,16000,nSpecies)

% If you tried at first and it failed (for example I didn't have the 
% external harddrive connected, it might have made a new empty file, which 
% will give a error:
% ...
% Dot indexing is not supported for variables of this type.
% Error in PIC/get_charge (line 2965)
%         iGroup = find(contains({fileInfo.Groups.Name},'/data'));
% ...
% the next time you try, so first remove this file.

% Add complementing data: 
%     UB - Total magnetic energy
%     A_xline - A value at main X line
%     x_xline - x location of main X line
%     z_xline - z location of main X line
%     Ey_xline - Ey at main X line (reconnection electric field)
%%

% First read object
pic = PIC(h5filepath); % If you have many times saved, this can take up to a minute
%pic = pic.twpelim([2200:400:4999],'exact');
h5write_fields_complement(pic)

% To work with the new data, the PIC object needs to be reloaded. Susanne, 
% if you use this hdf5 format, stop to use clear all (otherwise you need to 
% reload the object every time).

%% Add additional times
timesteps = 10000:1000:14000;
h5write_fields(datapath,h5filepath,timesteps,nSpecies)
pic = PIC(h5filepath); % If you have may times saved, this can take up to a minute
h5write_fields_complement(pic.twpelim(timesteps,'exact'))

