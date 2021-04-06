timestep = 10000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_test/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_test/data_h5/dists.h5';
distIndRead = 463:667;
distIndRead = 668:918;
nSpecies = 6;
iteration = nobg.twpelim(timestep).iteration;
mass = [25 1 25 1 25 1];
charge = [1 -1 1 -1 1 -1];
tag = 'line4';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tag)
%%
timestep = 5000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/df_cold_protons_04_new_boundary/data_h5/dists.h5';
distIndRead = 131:243;
nSpecies = 6;
iteration = df04n.twpelim(timestep).iteration;
mass = [25  1 25 1 25 1];
charge = [1 -1 1 -1 1 -1];
tag = 'idr vertical';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tag)

%% 
timestep = 24000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 990:910;
tags = arrayfun(@(s)sprintf('A=%g',s),all_tags,'UniformOutput',false);
nSpecies = 6;
iteration = no02m.twpelim(timestep).iteration;
mass = [100 1 100 1 100 1];
charge = [1 -1 1 -1 1 -1];
tag = 'line vertical';
for ic = 1:2000, tags{ic} = tag; end
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tags)

%% 
timestep = 20000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 301:418;
tags = arrayfun(@(s)sprintf('A=%g',s),all_tags,'UniformOutput',false);
nSpecies = 6;
iteration = no02m.twpelim(timestep).iteration;
mass = [100 1 100 1 100 1];
charge = [1 -1 1 -1 1 -1];
%tag = 'idr vertical';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tags)
%% 
timestep = 24000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 1601:1738;1445:1600;1738;%1738
tags = arrayfun(@(s)sprintf('A=7.5',s),1:2000,'UniformOutput',false);
%tags = arrayfun(@(s)sprintf('',s),1:2000,'UniformOutput',false);
% for itag = 1:numel(all_tags)
%   tags{distIndRead(itag)} = sprintf('A=%g',all_tags(itag));
% end

%tags = arrayfun(@(s)sprintf('A=%g',s),all_tags,'UniformOutput',false);
nSpecies = 6;
iteration = no02m.twpelim(timestep).iteration;
mass = [100 1 100 1 100 1];
charge = [1 -1 1 -1 1 -1];
%tag = 'idr vertical';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tags)

%% 
timestep = 23000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 1:280;%1738
tags = arrayfun(@(s)sprintf('A=7.5',s),1:280,'UniformOutput',false);
%tags = arrayfun(@(s)sprintf('',s),1:2000,'UniformOutput',false);
% for itag = 1:numel(all_tags)
%   tags{distIndRead(itag)} = sprintf('A=%g',all_tags(itag));
% end

%tags = arrayfun(@(s)sprintf('A=%g',s),all_tags,'UniformOutput',false);
nSpecies = 6;
iteration = no02m.twpelim(timestep).iteration;
mass = [100 1 100 1 100 1];
charge = [1 -1 1 -1 1 -1];
%tag = 'idr vertical';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tags)


%% From the start, add/overwrite tags later
timestep = 23000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 1:280;
%tags = arrayfun(@(s)sprintf('A=%g',s),all_tags,'UniformOutput',false);
nSpecies = 6;
iteration = no02m.twpelim(timestep).iteration;
mass = [100 1 100 1 100 1];
charge = [1 -1 1 -1 1 -1];
tag = '';
%for ic = 1:2000, tags{ic} = tag; end
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tags)

h5writeatt(h5FilePath, ['/'],'nSpecies', 6);

%% Write tags, setup
%A = no02m.twpelim(timestep).A;
hca = subplot(1,1,1);
imagesc(hca,no02m.xi,no02m.zi,A')
colormap(hca,pic_colors('pasteljet'))
hcb = colorbar('peer',hca);
hca.YDir = 'normal';
hold(hca,'on')
[cc,hh] = contour(hca,no02m.xi,no02m.zi,A',0:0.5:13,'k');
clabel(cc,hh)
hold(hca,'off')

hold(hca,'on')
%hds = ds100.twpelim(twpe).plot_boxes(hca);
hold(hca,'off')

% 1:23 A=10
% 24:47 A=9.5
% 47:74 A=9
% 74:101 A=8.5
% 102:131 A=8
% 132:163 A=7.5
% 164:197 A=7
% 198:234 A=6.5
% 235:275 A=6
% 276:321 A=5.5
% 322:377 A=5
% 378:438 A=4.5
% 439:449 A=4.5 correct - just an island 
% 450:542 A=4
% --
% 543:691 horizontal line z=0
% 692:840 horizontal line z=2 ? check
% 841:989 horizontal line z=4 ? check
% --
% 990:1444 vertical lines
% --
% 1445:1738 A=7.5

%% Write tags, write
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
twpe = 23000;
tag = 'A=7.5';
for idist = 1:280
  dataset = ['/data/' sprintf('%010.0f',2*twpe) '/' sprintf('%05.0f',idist) '/'];  
  h5writeatt(h5FilePath,dataset,'tag',tag)
end

%% Write dervied data, for example f(vpar)
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
twpe = 23000;
iteration = twpe*2;
str_iteration = sprintf('%010.0f',iteration);
dspart = ds100.twpelim(twpe).findtag({'A=7.5'}).dxlim([0 0.3]);
mass_all = no02m.mass;
charge_all = no02m.charge;
iSpecies = 3;

for distnumber = 1%:ds.nd{1}
  ds = dspart.update_inds({distnumber});
  xlo = ds.xi1{1};
  xhi = ds.xi2{1};
  zlo = ds.zi1{1};
  zhi = ds.zi2{1};
  
  xdist = (ds.xi1{1}+ds.xi2{1})/2;
  zdist = (ds.zi1{1}+ds.zi2{1})/2;
  dxdist = ds.xi1{1}-ds.xi2{1};
  dzdist = ds.zi1{1}-ds.zi2{1};
  tdist = repmat(twpe,size(xdist));
  Bx = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'Bx');
  By = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'By');
  Bz = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'Bz');
  Bx = 1;
  By = 0;
  Bz = 0;
  
  % Calculate reduced distribtuions
  for iSpecies = 2;1:6
    fred_tmp = ds.reduce_1d_new('x',[iSpecies],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});
  end
  
  % Write data
  dataset_name_base = ['/data/' str_iteration '/' num2str(distnumber,'%05.0f'),'/fxyz'];
  dataset_name = ['/data/' str_iteration '/' num2str(distnumber,'%05.0f'),'/fvpar/',num2str(iSpecies,'%1.0f')];
  try
    h5create(h5FilePath, dataset_name, size(fxyz));
  catch
    warning('h5 structure %s already exists, overwriting.',dataset_name)      
  end
  h5write(h5FilePath, dataset_name, fxyz);
  h5writeatt(h5FilePath, dataset_name,'x', [xlo xhi]);
  h5writeatt(h5FilePath, dataset_name,'z', [zlo zhi]);
  h5writeatt(h5FilePath, dataset_name,'ic', ic);
  h5writeatt(h5FilePath, dataset_name,'vxa', vxa);
  h5writeatt(h5FilePath, dataset_name,'vya', vya);
  h5writeatt(h5FilePath, dataset_name,'vza', vza);
  h5writeatt(h5FilePath, dataset_name,'axes', axes);
  h5writeatt(h5FilePath, dataset_name,'mass', mass);
  h5writeatt(h5FilePath, dataset_name,'charge', charge);
  
  
end
