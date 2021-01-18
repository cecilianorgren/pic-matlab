%% Write basic data
h5filepath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields_new.h5';
datapath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data/';
nSpecies = 6;
h5write_fields(datapath,h5filepath,15000:1000:16000,nSpecies)
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

% First read object
pic = PIC(h5filepath); % If you have many times saved, this can take up to a minute
h5write_fields_complement(pic)

% To work with the new data, the PIC object needs to be reloaded. Susanne, 
% if you use this hdf5 format, stop to use clear all (otherwise you need to 
% reload the object every time).

%% Add additional times
timesteps = 17000:1000:18000;
h5write_fields(datapath,h5filepath,timesteps,nSpecies)
pic = PIC(h5filepath); % If you have may times saved, this can take up to a minute
h5write_fields_complement(pic.twpelim(timesteps,'exact'))

