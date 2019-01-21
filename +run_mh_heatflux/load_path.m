localuser = datastore('local','user');
%pathFields = ['/Users/' localuser '/PIC/Code/reduced/fields/'];
pathFields = '/Volumes/Fountain/Data/PIC/run_michael_heat_flux/fields/';
pathParts = '/Volumes/Fountain/Data/PIC/run_michael_heat_flux/parts/';
pathSave = ['/Users/' localuser '/Research/Simulations/run_mh_heatflux/'];

fprintf('Fields (%s):',pathFields)
dir([pathFields '/*.dat'])
fprintf('Distributions (%s):',pathParts)
dir([pathParts '/*.dat'])