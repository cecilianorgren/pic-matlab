%% no02m
timestep = 24000;
dirData = '/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/distributions/24000/';
h5FilePath = '/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 1501:1500;
nSpecies = 6;
iteration = timestep*2;
mass = [100 1 100 1 100 1];
charge = [1 -1 1 -1 1 -1];
tag = '';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tag)


%% Reconnection onset, susannes
timestep = 2400;
dirData = '/Users/cecilia/Data/PIC/rec_onset/distributions/t2400/';
h5FilePath = '/Users/cecilia/Data/PIC/rec_onset/h5data/dists.h5';
distIndRead = 1:117;
nSpecies = 2;
iteration = timestep*2;
mass = [25 1];
charge = [1 -1];
tag = 'map_recsite';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tag)

%%
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
timestep = 20000;
dirData = sprintf('/Users/cno062/Data/PIC/no_hot_bg_n02_m100/distributions/%05.0f/',timestep);
h5FilePath = '/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 1:20;%1738
%tags = arrayfun(@(s)sprintf('A=7.5',s),1:2000,'UniformOutput',false);
tags = arrayfun(@(s)sprintf('',s),1:2000,'UniformOutput',false);
% for itag = 1:numel(all_tags)
%   tags{distIndRead(itag)} = sprintf('A=%g',all_tags(itag));
% end

%tags = arrayfun(@(s)sprintf('A=%g',s),all_tags,'UniformOutput',false);
nSpecies = 6;
%iteration = no02m.twpelim(timestep).iteration;
iteration = 1;
mass = [100 1 100 1 100 1];
charge = [1 -1 1 -1 1 -1];
%tag = 'idr vertical';
h5write_dists(dirData,h5FilePath,distIndRead,nSpecies,mass,charge,timestep,iteration,tags)

%% 
timestep = 24000;
dirData = sprintf('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/distributions/%05.0f/',timestep);
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
distIndRead = 2051:2243;%1738
tags = arrayfun(@(s)sprintf('A=8.0',s),1:3000,'UniformOutput',false);
%tags = arrayfun(@(s)sprintf('',s),1:2500,'UniformOutput',false);
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
twpe = 21000;
tag = '';
%tag = 'A=6.0';
for idist = 1:ds100.twpelim(twpe).nd{1}%1740%:2049%1:280
  dataset = ['/data/' sprintf('%010.0f',2*twpe) '/' sprintf('%05.0f',idist) '/'];    
  h5writeatt(h5FilePath,dataset,'tag',tag)
  dataset = ['/data/' sprintf('%010.0f',2*twpe) '/' sprintf('%05.0f',idist) '/fxyz'];  
  h5writeatt(h5FilePath,dataset,'tag',tag)
end

%% Write derived data, for example f(vpar)
h5FilePath = '/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5';
twpe = 24000;
iteration = twpe*2;
str_iteration = sprintf('%010.0f',iteration);
dspart = ds100.twpelim(twpe).findtag({'A=8.0'}).dxlim([0 0.3]);
%dspart = ds100.twpelim(twpe);
mass_all = no02m.mass;
charge_all = no02m.charge;
%iSpecies = 3;
nInterp = 0;
doInterp = 0;
% twpe24000 is done until 484
% need to redo first 311 indices, because they were overwritten when I had
% the distnumber wrong.
it = 1;
for distcount = 1:dspart.nd{1}
  %disp(sprintf('id = %g/%g',distnumber,dspart.nd{1}))
  ds = dspart.update_inds({distcount});
  distnumber = ds.indices{1};
  dataset_base = ['/data/' str_iteration '/' num2str(distnumber,'%05.0f')];
  disp(sprintf('%s',dataset_base))
  
  
  xlo = ds.xi1{1};
  xhi = ds.xi2{1};
  zlo = ds.zi1{1};
  zhi = ds.zi2{1};  
  xdist = (ds.xi1{1}+ds.xi2{1})/2;
  zdist = (ds.zi1{1}+ds.zi2{1})/2;
  dxdist = ds.xi2{1}-ds.xi1{1};
  dzdist = ds.zi2{1}-ds.zi1{1};
  tdist = repmat(twpe,size(xdist));
  
  if 0 % major error was not here, but in how i did dxdist, but using xlo and xhi direclty is better I think
    % They still differ
    Bx = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'Bx');
    By = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'By');
    Bz = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'Bz');
    Ex = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'Ex');
    Ey = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'Ey');
    Ez = no02m.twpelim(twpe).get_points(xdist,zdist,tdist,dxdist*0.5*[-1 1],'Ez');
  else
    Bx = mean(mean(no02m.twpelim(twpe).xlim([xlo xhi]).zlim([zlo zhi]).Bx));
    By = mean(mean(no02m.twpelim(twpe).xlim([xlo xhi]).zlim([zlo zhi]).By));
    Bz = mean(mean(no02m.twpelim(twpe).xlim([xlo xhi]).zlim([zlo zhi]).Bz));
    Ex = mean(mean(no02m.twpelim(twpe).xlim([xlo xhi]).zlim([zlo zhi]).Ex));
    Ey = mean(mean(no02m.twpelim(twpe).xlim([xlo xhi]).zlim([zlo zhi]).Ey));
    Ez = mean(mean(no02m.twpelim(twpe).xlim([xlo xhi]).zlim([zlo zhi]).Ez));  
  end
  %disp(sprintf('B = [%.2f,%.2f,%.1f], E = [%.2f,%.2f,%.1f]',Bx,By,Bz,Ex,Ey,Ez))
  %Bx = 1;
  %By = 0;
  %Bz = 0;
  
  
  doInterp = 0;
  if 1 % Calculate and write reduced distributions
    fpar = zeros(101,6);
    vpar_center = zeros(101,6);  
    vpar_edges = zeros(102,6); 
    fpar_nointerp = zeros(101,6);
    vpar_center_nointerp = zeros(101,6);
    vpar_edges_nointerp = zeros(102,6);
    tic % check time  
    for iSpecies = 1:6
      fred_tmp = ds.reduce_1d_new('x',[iSpecies],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz},'interp',nInterp);
      fpar(:,iSpecies) = fred_tmp.fvpar;
      vpar_center(:,iSpecies) = fred_tmp.vpar_center;
      vpar_edges(:,iSpecies) = fred_tmp.vpar_edges;      
  %     fred_tmp_nointerp = ds.reduce_1d_new('x',[iSpecies],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz});   
  %     fpar_nointerp(:,iSpecies) = fred_tmp_nointerp.fvpar;
  %     vpar_center_nointerp(:,iSpecies) = fred_tmp_nointerp.vpar_center;
  %     vpar_edges_nointerp(:,iSpecies) = fred_tmp_nointerp.vpar_edges;
    end
    %figure(77); plot(vpar_center,fpar)
    toc    
    if 1 % Write data
      dataset_name_base = [dataset_base,'/fxyz'];
      dataset_name = [dataset_base,'/fpar'];
      try
        h5create(h5FilePath, dataset_name, size(fpar));
      catch
        warning('h5 structure %s already exists, overwriting.',dataset_name)      
      end    
      h5write(h5FilePath, dataset_name, fpar);
      if doInterp == 1, h5writeatt(h5FilePath, dataset_name,'interp', nInterp);
      else doInterp == 1, h5writeatt(h5FilePath, dataset_name,'interp', 0); end    
      h5writeatt(h5FilePath, dataset_name,'x', [xlo xhi]);
      h5writeatt(h5FilePath, dataset_name,'z', [zlo zhi]);
      h5writeatt(h5FilePath, dataset_name,'axes', vpar_center);
      h5writeatt(h5FilePath, dataset_name,'axes_edges', vpar_edges);
      h5writeatt(h5FilePath, dataset_name,'charge', h5readatt(h5FilePath,dataset_name_base,'charge'));
      h5writeatt(h5FilePath, dataset_name,'mass', h5readatt(h5FilePath,dataset_name_base,'mass'));
      %h5writeatt(h5FilePath, dataset_name,'B', [Bx By Bz]);
      %h5writeatt(h5FilePath, dataset_name,'E', [Bx By Bz]);
      h5disp(h5FilePath,dataset_name)
    end
  end
  if 1 % Calculate and write reduced fx,fy,fz
    tic % check time 
    
    fx = zeros(101,6);
    fy = zeros(101,6);
    fz = zeros(101,6);
    v_center = zeros(101,6);  
    v_edges = zeros(101,6);  
    for iSpecies = 1:6
      fred_tmp = ds.f(it,1,iSpecies); % ds is already the picked subindex
      %fred_tmp = ds.f(it,distnumber,iSpecies);
      fx(:,iSpecies) = fred_tmp.fx;
      fy(:,iSpecies) = fred_tmp.fy;
      fz(:,iSpecies) = fred_tmp.fz;
      v_center(:,iSpecies) = fred_tmp.v;      
    end
    if 1 % Write data
      % fx
      dataset_name = [dataset_base,'/fx'];
      try h5create(h5FilePath, dataset_name, size(fx));
      catch warning('h5 structure %s already exists, overwriting.',dataset_name)      
      end  
      h5write(h5FilePath, dataset_name, fx);
      h5writeatt(h5FilePath, dataset_name,'axes', v_center);
      % fy
      dataset_name = [dataset_base,'/fy'];
      try h5create(h5FilePath, dataset_name, size(fy));
      catch warning('h5 structure %s already exists, overwriting.',dataset_name)      
      end  
      h5write(h5FilePath, dataset_name, fy);
      h5writeatt(h5FilePath, dataset_name,'axes', v_center);
      % fz      
      dataset_name = [dataset_base,'/fz'];
      try h5create(h5FilePath, dataset_name, size(fz));
      catch warning('h5 structure %s already exists, overwriting.',dataset_name)      
      end  
      h5write(h5FilePath, dataset_name, fz);
      h5writeatt(h5FilePath, dataset_name,'axes', v_center);
    end
    toc
  end
  
  if 1 % write E, B
    dataset_name = [dataset_base];
    h5writeatt(h5FilePath, dataset_name,'B', [Bx By Bz]);
    h5writeatt(h5FilePath, dataset_name,'E', [Ex Ey Ez]);
    %h5disp(h5FilePath,dataset_name)
  end
  
end
