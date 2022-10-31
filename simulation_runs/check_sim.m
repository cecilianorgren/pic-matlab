h5filepath = '/Users/Cecilia/Data/PIC/rec_onset_1/data_h5/fields.h5';
datapath = '/Users/Cecilia/Data/PIC/rec_onset_1/data/';
nSpecies = 2;
h5write_fields(datapath,h5filepath,20,nSpecies)

%%
pic = PIC(h5filepath);