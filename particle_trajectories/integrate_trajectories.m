%% Get r0,v0 from moments, suitable for inflow at early times
pic = no02m.twpelim(15000);
A = squeeze(pic.A);
z_center = 3:1:7;
x_center = 80:5:125;

[ZC,XC] = meshgrid(z_center,x_center);
XC = reshape(XC,prod(size(XC)),1);
ZC = reshape(ZC,prod(size(ZC)),1);
nboxes = numel(XC);
clear fpeaks
if 0 % A
  Alim = [-24 0];
  ind_keep = zeros(nboxes,1);
  for ibox = 1:nboxes
    xind = find(abs(pic.xi-XC(ibox))==min(abs(pic.xi-XC(ibox))));
    zind = find(abs(pic.zi-ZC(ibox))==min(abs(pic.zi-ZC(ibox))));

    if A(xind,zind)>Alim(1) && A(xind,zind)<Alim(2) % keep
      ind_keep(ibox) = 1;
    else
      ind_keep(ibox) = 0;
    end
  end
  keep_boxes = all_boxes(find(ind_keep==1),:);
  n_boxes = size(keep_boxes,1);
  x0 = XC(find(ind_keep));
  z0 = ZC(find(ind_keep));
else % n
  Alim = [0 8.2];
  vlim = [0 1];
  tlim = [0 0.0005];
  A = pic.A;
  t = pic.t(4);
  n = pic.n(4);
  vx = pic.vx(4);
  vy = pic.vy(4);
  vz = pic.vz(4);
  vabs = sqrt(vx.^2+vy.^2+vz.^2);
  nlim = 0.02;
  ind_keep = zeros(nboxes,1);
  for ibox = 1:nboxes
    xind = find(abs(pic.xi-XC(ibox))==min(abs(pic.xi-XC(ibox))));
    zind = find(abs(pic.zi-ZC(ibox))==min(abs(pic.zi-ZC(ibox))));

    x0(ibox) = XC(ibox);
    y0(ibox) = 0;
    z0(ibox) = ZC(ibox);
    vx0(ibox) = vx(xind,zind);
    vy0(ibox) = vy(xind,zind);
    vz0(ibox) = vz(xind,zind);
    
%     if n(xind,zind)>nlim(1) % keep
%       ind_keep(ibox) = 1;
%     else
%       ind_keep(ibox) = 0;
%     end
    if A(xind,zind)>Alim(1) && A(xind,zind)<Alim(2) % keep
      ind_keep(ibox) = 1;
%       if t(xind,zind)>tlim(2) % too fast discard
%         ind_keep(ibox) = 0;
%       else
%         ind_keep(ibox) = 1;
%       end
    else
      ind_keep(ibox) = 0;
    end
  end
  %x0 = XC(find(ind_keep));
  %z0 = ZC(find(ind_keep));
  x0 = x0(find(ind_keep));
  y0 = y0(find(ind_keep));
  z0 = z0(find(ind_keep));
  vx0 = vx0(find(ind_keep));
  vy0 = vy0(find(ind_keep));
  vz0 = vz0(find(ind_keep));
  for ii = 1:numel(x0)
    fpeaks(ii).x = x0(ii);
    fpeaks(ii).y = y0(ii);
    fpeaks(ii).z = z0(ii);
    fpeaks(ii).vx = vx0(ii);
    fpeaks(ii).vy = vy0(ii);
    fpeaks(ii).vz = vz0(ii);
  end
  hca = subplot(1,1,1);
  imagesc(hca,pic.xi,pic.zi,n');
  colormap(hca,pic_colors('blue_red'))
  colorbar('peer',hca)
  hca.YDir = 'normal';
  hca.Title.String = sprintf('number of position = %g',numel(fpeaks));
  hold(hca,'on')
  plot(hca,[fpeaks.x],[fpeaks.z],'.k')
  hold(hca,'off')
end

%% Integrate trajectories based on fpeaks (which may be based or moments)
%trp = PICTraj('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories_paul.h5');
h5path_traj = '/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories.h5';
%trp = PICTraj();
pic = no02m;
tspan = [75,125];

m = 1/100; 
q = -1;
istart = 28;
ntr_pre = 0;%trp.ntr-(istart-1); % the +(istart-1) is there becuse i already did the 1 and added it to the datafile
%ntr_pre = 0;

for iTr = istart:numel(fpeaks)
  tic
  fprintf('iTr = %g/%g\n',iTr,numel(fpeaks))
  r0 = [fpeaks(iTr).x, fpeaks(iTr).y, fpeaks(iTr).z];
  v0 = [fpeaks(iTr).vx, fpeaks(iTr).vy, fpeaks(iTr).vz];
  tr_tmp = no02m.integrate_trajectory(r0,v0,tspan,m,q);
  
  hca = subplot(2,1,1);
  plot(hca,tr_tmp.x,tr_tmp.z,tr_tmp.x0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca = subplot(2,1,2);
  plot(hca,tr_tmp.y,tr_tmp.z,tr_tmp.y0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  drawnow
%   [Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(tr_tmp.x,tr_tmp.z,tr_tmp.t);  % interpolate
% 
%   tr_tmp.t0 = t0;
%   tr_tmp.Ex = Ex;
%   tr_tmp.Ey = Ey;
%   tr_tmp.Ez = Ez;
%   tr_tmp.Bx = Bx;
%   tr_tmp.By = By;
%   tr_tmp.Bz = Bz;
  
  h5write_trajs(h5path_traj,tr_tmp,'id',ntr_pre+iTr)
  %h5write_trajs('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr_tmp)
  %tr(iPeak,id) = tr_tmp;
  toc
  %catch
  %  continue
  %end
end

%% Make fake trajectories prior to twci = 5, IONS
localuser = datastore('local','user');
loadPath = ['/Users/' localuser '/Data/PIC/no02m_ti_A_preonset/'];
% keep them on the same magnetic field line
twpe_pre = [-10000:400:800];
twci_pre = twpe_pre/200;
nt = numel(twpe_pre);
x0 = 90:5:115; % particles moves down, with small constant vx
% possibility to add some dx or so later
vx0 = -0.025+1*0.05*rand(size(x0));
vx0 = [0.0044, 0.0234, -0.0093, -0.0030, -0.0156, 0.0044];
vy0 = 0.1;

z_part = zeros(nt,numel(x0));
x_part = repmat(x0,nt,1);
z_phase = rand(1,numel(x0));
t_part = 0.9+0.4*rand(1,numel(x0));
vx_part = repmat(vx0,nt,1);

x_part = x_part + cumsum(vx_part,1);


fci = 1;
fce = fci*200;
rci = 0.1;
rce = rci/100;
dtwci = unique(diff(twci_pre));

for it = 1:nt
  twpe = twpe_pre(it);
  
  % Load data  
  data = load([loadPath sprintf('ti_A_twpe=%06.0f.mat',twpe)]);
  A = data.A;
  t = data.ti;
  xi = data.xi;
  zi = data.zi;  
  twpe = data.twpe;
  twci = data.twci;  
  nx = numel(xi);
  nz = numel(zi);
  
  % Find the correct field line and assign the z-position
  S = contourcs(xi,zi,A',Alev*[1 1]);
  S = S(1);
  
  for ix = 1:numel(x0)
    
    ii = find(x>x0(ix),1,'first');
    z_part(it,ix) = S.Y(ii)+rci*t_part(ix)*sin(twpe/200+z_phase(ix));
  end
  plot(x_part,z_part,'.')
  drawnow
end

vx_part = [zeros(size(x0)); diff(x_part,1)./dtwci];
vz_part = [zeros(size(z0)); diff(z_part,1)./dtwci];

%% Pick manual r0 (based on A, with white line), v0
Alev = 7;
twpe = 1000;
A0 = no02m.twpelim(twpe).A;
x = no02m.xi;
z = no02m.zi;
S = contourcs(x,z,A0',Alev*[1 1]);
S = S(1);
x0 = 90:5:115;
x0 = x_part(end,:); % from movie_open_script_rec_onset.m
for ix = 1:numel(x0)
  ii = find(x>x0(ix),1,'first');
  z0(ix) = S.Y(ii);
end

y0 = zeros(size(x0));
%vx0 = zeros(size(x0));
vy0 = zeros(size(x0));
%vz0 = zeros(size(x0));
vx0 = vx_part(end,:);
vz0 = vz_part(end,:);

tspan = [no02m.twpelim(twpe).twci,no02m.twci(end)];

%tspan = [5 10];

m = 1; 
q = 1;
istart = 1;
ntr_pre = trp.ntr-(istart-1); % the +(istart-1) is there becuse i already did the 1 and added it to the datafile
%ntr_pre = 0;

for iTr = istart:numel(x0)
  tic
  fprintf('iTr = %g/%g\n',iTr,numel(x0))
  r0 = [x0(iTr), y0(iTr), z0(iTr)];
  v0 = [vx0(iTr), vy0(iTr), vz0(iTr)];
  tr_tmp = no02m.integrate_trajectory(r0,v0,tspan,m,q);
  
  
  hca = subplot(2,1,1);
  plot(hca,tr_tmp.x,tr_tmp.z,tr_tmp.x0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca = subplot(2,1,2);
  plot(hca,tr_tmp.y,tr_tmp.z,tr_tmp.y0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  drawnow
%   [Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(tr_tmp.x,tr_tmp.z,tr_tmp.t);  % interpolate
% 
%   tr_tmp.t0 = t0;
%   tr_tmp.Ex = Ex;
%   tr_tmp.Ey = Ey;
%   tr_tmp.Ez = Ez;
%   tr_tmp.Bx = Bx;
%   tr_tmp.By = By;
%   tr_tmp.Bz = Bz;
  
  h5write_trajs('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories_paul_t0=5_test.h5',tr_tmp,'id',ntr_pre+iTr)
  %h5write_trajs('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr_tmp)
  %tr(iPeak,id) = tr_tmp;
  toc
  %catch
  %  continue
  %end
end

%% Merge fake and real trajectories
trp = PICTraj('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories_paul_t0=5_test.h5');

ntr = numel(trp);
for itr = 1:ntr
  trall(itr).t = [twci_pre'; trp(itr).t];
  trall(itr).twpe = [twpe_pre'; trp(itr).t*200];
  trall(itr).x = [x_part(:,itr); trp(itr).x];
  trall(itr).z = [z_part(:,itr); trp(itr).z];  
end

%% Make fake trajectories prior to twci = 5, ELECTRONS
localuser = datastore('local','user');
loadPath = ['/Users/' localuser '/Data/PIC/no02m_ti_A_preonset/'];
% keep them on the same magnetic field line
twpe_pre = [-10000:400:800];
twci_pre = twpe_pre/200;
nt = numel(twpe_pre);
x0 = 90:5:115; % particles moves down, with smallconstant vx
% possibility to add some dx or so later
vx0 = -0.025+1*0.05*rand(size(x0));
vx0 = [0.0044, 0.0234, -0.0093, -0.0030, -0.0156, 0.0044];
vy0 = 0.1;

z_part = zeros(nt,numel(x0));
x_part = repmat(x0,nt,1);
z_phase = rand(1,numel(x0));
t_part = 0.9+0.4*rand(1,numel(x0));
vx_part = repmat(vx0,nt,1);

x_part = x_part + cumsum(vx_part,1);


fci = 1;
fce = fci*200;
rci = 0.1;
rce = rci/100;
dtwci = unique(diff(twci_pre));

for it = 1:nt
  twpe = twpe_pre(it);
  
  % Load data  
  data = load([loadPath sprintf('ti_A_twpe=%06.0f.mat',twpe)]);
  A = data.A;
  t = data.ti;
  xi = data.xi;
  zi = data.zi;  
  twpe = data.twpe;
  twci = data.twci;  
  nx = numel(xi);
  nz = numel(zi);
  
  % Find the correct field line and assign the z-position
  S = contourcs(xi,zi,A',Alev*[1 1]);
  S = S(1);
  
  for ix = 1:numel(x0)
    
    ii = find(x>x0(ix),1,'first');
    z_part(it,ix) = S.Y(ii)+rci*t_part(ix)*sin(twpe/200+z_phase(ix));
  end
  plot(x_part,z_part,'.')
  drawnow
end

vx_part = [zeros(size(x0)); diff(x_part,1)./dtwci];
vz_part = [zeros(size(z0)); diff(z_part,1)./dtwci];

%% Pick manual r0 (based on A, with white line), v0
Alev = 7;
twpe = 1000;
A0 = no02m.twpelim(twpe).A;
x = no02m.xi;
z = no02m.zi;
S = contourcs(x,z,A0',Alev*[1 1]);
S = S(1);
x0 = 90:5:115;
x0 = x_part(end,:); % from movie_open_script_rec_onset.m
for ix = 1:numel(x0)
  ii = find(x>x0(ix),1,'first');
  z0(ix) = S.Y(ii);
end

y0 = zeros(size(x0));
vx0 = zeros(size(x0));
vy0 = zeros(size(x0));
vz0 = zeros(size(x0));
vx0 = vx_part(end,:);
vz0 = vz_part(end,:);

tspan = [no02m.twpelim(twpe).twci,no02m.twci(end)];

tspan = [5 50];

m = 1/100; 
q = -1;
istart = 1;
ntr_pre = trp.ntr-(istart-1); % the +(istart-1) is there becuse i already did the 1 and added it to the datafile
%ntr_pre = 0;

for iTr = istart%:numel(x0)
  tic
  fprintf('iTr = %g/%g\n',iTr,numel(x0))
  r0 = [x0(iTr), y0(iTr), z0(iTr)];
  v0 = [vx0(iTr), vy0(iTr), vz0(iTr)];
  tr_tmp = no02m.integrate_trajectory(r0,v0,tspan,m,q);
  
  
  hca = subplot(2,1,1);
  plot(hca,tr_tmp.x,tr_tmp.z,tr_tmp.x0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca = subplot(2,1,2);
  plot(hca,tr_tmp.y,tr_tmp.z,tr_tmp.y0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  drawnow
%   [Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(tr_tmp.x,tr_tmp.z,tr_tmp.t);  % interpolate
% 
%   tr_tmp.t0 = t0;
%   tr_tmp.Ex = Ex;
%   tr_tmp.Ey = Ey;
%   tr_tmp.Ez = Ez;
%   tr_tmp.Bx = Bx;
%   tr_tmp.By = By;
%   tr_tmp.Bz = Bz;
  
  h5write_trajs('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories_paul_t0=5_ele.h5',tr_tmp,'id',ntr_pre+iTr)
  %h5write_trajs('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr_tmp)
  %tr(iPeak,id) = tr_tmp;
  toc
  %catch
  %  continue
  %end
end

%% Merge fake and real trajectories
trp = PICTraj('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories_paul_t0=5_test.h5');

ntr = numel(trp);
for itr = 1:ntr
  trall(itr).t = [twci_pre'; trp(itr).t];
  trall(itr).twpe = [twpe_pre'; trp(itr).t*200];
  trall(itr).x = [x_part(:,itr); trp(itr).x];
  trall(itr).z = [z_part(:,itr); trp(itr).z];  
end

%% Pick manual r0 (based on A, with white line), v0, ELECTRONS
Alev = 7;
twpe = 1000;
A0 = no02m.twpelim(twpe).A;
x = no02m.xi;
z = no02m.zi;
S = contourcs(x,z,A0',Alev*[1 1]);
S = S(1);
x0 = (90+2.5):5:115;
%x0 = x_part(end,:); % from movie_open_script_rec_onset.m
for ix = 1:numel(x0)
  ii = find(x>x0(ix),1,'first');
  z0(ix) = S.Y(ii);
end

y0 = zeros(size(x0));
vx0 = zeros(size(x0));
vy0 = zeros(size(x0));
vz0 = zeros(size(x0));
%vx0 = vx_part(end,:);
%vz0 = vz_part(end,:);
%%
twpe = 1000;
tspan = [no02m.twpelim(twpe).twci,no02m.twci(end)];

(no02m.xi(end)-no02m.xi(1))/2;
x0 = [101 102 103];
x0 = [101 102]+0.5;
z0 = [1 1 1];
y0 = [0 0 0];

vz0 = zeros(size(x0));
vx0 = zeros(size(x0));
vz0 = zeros(size(x0));

%tspan = [5 10];

m = 1/100; 
q = -1;
istart = 1;
ntr_pre = trp.ntr-(istart-1); % the +(istart-1) is there becuse i already did the 1 and added it to the datafile
%ntr_pre = 0;

for iTr = istart:numel(x0)
  tic
  fprintf('iTr = %g/%g\n',iTr,numel(x0))
  r0 = [x0(iTr), y0(iTr), z0(iTr)];
  v0 = [vx0(iTr), vy0(iTr), vz0(iTr)];
  tr_tmp = no02m.integrate_trajectory(r0,v0,tspan,m,q);
  
  
  hca = subplot(2,1,1);
  plot(hca,tr_tmp.x,tr_tmp.z,tr_tmp.x0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hca = subplot(2,1,2);
  plot(hca,tr_tmp.y,tr_tmp.z,tr_tmp.y0,tr_tmp.z0,'o')  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
  drawnow
%   [Ex,Ey,Ez,Bx,By,Bz] = df04.interp_EB3(tr_tmp.x,tr_tmp.z,tr_tmp.t);  % interpolate
% 
%   tr_tmp.t0 = t0;
%   tr_tmp.Ex = Ex;
%   tr_tmp.Ey = Ey;
%   tr_tmp.Ez = Ez;
%   tr_tmp.Bx = Bx;
%   tr_tmp.By = By;
%   tr_tmp.Bz = Bz;
  
  h5write_trajs('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories_paul_t0=5.h5',tr_tmp,'id',ntr_pre+iTr)
  %h5write_trajs('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr_tmp)
  %tr(iPeak,id) = tr_tmp;
  toc
  %catch
  %  continue
  %end
end



