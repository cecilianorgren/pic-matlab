%% See what distributions we have already
timestep = 07000;
root_dir = '/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/';
txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/distributions/%.0f.dat',distnumber);

%% Pattern governed by A, timestep 08000
x_center = 150:1:205;
z_center = 0:1:10;
dx_box = 0.5;
dz_box = 0.5;
% x_center = 150:0.5:205;
% z_center = 0:0.5:10;
% dx_box = 0.25;
% dz_box = 0.25;

xlow = x_center-dx_box;
xhigh = x_center+dx_box;
zlow = z_center-dz_box;
zhigh = z_center+dz_box;

[XC,ZC] = meshgrid(x_center,z_center);

XC = reshape(XC,prod(size(XC)),1);
ZC = reshape(ZC,prod(size(ZC)),1);

nboxes = numel(XC);
Alim = [-24 -17];
Alim = [-24 -18];
Alim = [-20 -16];
ind_keep = zeros(nboxes,1);
for ibox = 1:nboxes
  xind = find(abs(x-XC(ibox))==min(abs(x-XC(ibox))));
  zind = find(abs(z-ZC(ibox))==min(abs(z-ZC(ibox))));
  A(xind,zind);
  if A(xind,zind)>Alim(1) && A(xind,zind)<Alim(2) % keep
    ind_keep(ibox) = 1;
  else
    ind_keep(ibox) = 0;
  end
end

all_boxes = [XC-dx_box XC+dx_box ZC-dz_box ZC+dz_box];
keep_boxes = all_boxes(find(ind_keep==1),:);
n_boxes = size(keep_boxes,1);

figure(401)
hca = subplot(1,1,1);
imagesc(hca,x,z,A')
%imagesc(hca,x,z,pi1.scalar')
hca.XLim = [100 210];
hca.YLim = [-12 12];
hca.Title.String = sprintf('n_boxes = %g',n_boxes);
hca.Title.Interpreter = 'none';
hcb = colorbar('peer',hca);


hold(hca,'on')
for ibox = 1:n_boxes      
  %hstar = plot(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'*k');
  %hpatch = patch(hca,[keep_boxes(ibox,1) keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2)],[keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2) keep_boxes(ibox,1)],'w');
  hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
  hpatch.FaceAlpha = 0;
  hpatch.LineWidth = 1;   
  
end
hold(hca,'off')   

%% Pattern governed by A, timestep 05000
x_center = 184:0.5:205;
z_center = 0:0.5:5;
dx_box = 0.25;
dz_box = 0.25;
% x_center = 150:0.5:205;
% z_center = 0:0.5:10;
% dx_box = 0.25;
% dz_box = 0.25;

xlow = x_center-dx_box;
xhigh = x_center+dx_box;
zlow = z_center-dz_box;
zhigh = z_center+dz_box;

[XC,ZC] = meshgrid(x_center,z_center);

XC = reshape(XC,prod(size(XC)),1);
ZC = reshape(ZC,prod(size(ZC)),1);

nboxes = numel(XC);
Alim = [-24 -17];
Alim = [-23.3 -21.7];
ind_keep = zeros(nboxes,1);
for ibox = 1:nboxes
  xind = find(abs(x-XC(ibox))==min(abs(x-XC(ibox))));
  zind = find(abs(z-ZC(ibox))==min(abs(z-ZC(ibox))));
  A(xind,zind);
  if A(xind,zind)>Alim(1) && A(xind,zind)<Alim(2) % keep
    ind_keep(ibox) = 1;
  else
    ind_keep(ibox) = 0;
  end
end

all_boxes = [XC-dx_box XC+dx_box ZC-dz_box ZC+dz_box];
keep_boxes = all_boxes(find(ind_keep==1),:);
n_boxes = size(keep_boxes,1);

figure(401)
hca = subplot(1,1,1);
imagesc(hca,x,z,jtot.y')
%imagesc(hca,x,z,pi1.scalar')
hca.YDir = 'normal';
hca.XLim = [100 210];
hca.YLim = [-12 12];
hca.Title.String = sprintf('n_boxes = %g',n_boxes);
hca.Title.Interpreter = 'none';
hcb = colorbar('peer',hca);


hold(hca,'on')
for ibox = 1:n_boxes      
  %hstar = plot(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'*k');
  %hpatch = patch(hca,[keep_boxes(ibox,1) keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2)],[keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2) keep_boxes(ibox,1)],'w');
  hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
  hpatch.FaceAlpha = 0;
  hpatch.LineWidth = 1;   
  
end
hold(hca,'off')   

%% Pattern governed by A, timestep 06000
x_center = 180:0.5:205;
z_center = 0:0.5:5;
dx_box = 0.25;
dz_box = 0.25;
% x_center = 150:0.5:205;
% z_center = 0:0.5:10;
% dx_box = 0.25;
% dz_box = 0.25;

xlow = x_center-dx_box;
xhigh = x_center+dx_box;
zlow = z_center-dz_box;
zhigh = z_center+dz_box;

[XC,ZC] = meshgrid(x_center,z_center);

XC = reshape(XC,prod(size(XC)),1);
ZC = reshape(ZC,prod(size(ZC)),1);

nboxes = numel(XC);
Alim = [-24 -17];
Alim = [-23.3 -20.5];
ind_keep = zeros(nboxes,1);
for ibox = 1:nboxes
  xind = find(abs(x-XC(ibox))==min(abs(x-XC(ibox))));
  zind = find(abs(z-ZC(ibox))==min(abs(z-ZC(ibox))));
  A(xind,zind);
  if A(xind,zind)>Alim(1) && A(xind,zind)<Alim(2) % keep
    ind_keep(ibox) = 1;
  else
    ind_keep(ibox) = 0;
  end
end

all_boxes = [XC-dx_box XC+dx_box ZC-dz_box ZC+dz_box];
keep_boxes = all_boxes(find(ind_keep==1),:);
n_boxes = size(keep_boxes,1);

% Plot 
figure(401)
hca = subplot(2,1,1);
imagesc(hca,x,z,B.y')
%imagesc(hca,x,z,pi1.scalar')
hca.YDir = 'normal';
hca.XLim = [140 210];
hca.YLim = [-12 12];
hca.Title.String = sprintf('n_boxes = %g',n_boxes);
hca.Title.Interpreter = 'none';
hcb = colorbar('peer',hca);

hold(hca,'on')
for ibox = 1:n_boxes      
  %hstar = plot(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'*k');
  %hpatch = patch(hca,[keep_boxes(ibox,1) keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2)],[keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2) keep_boxes(ibox,1)],'w');
  hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
  hpatch.FaceAlpha = 0;
  hpatch.LineWidth = 1;   
  
end
hold(hca,'off')   

hca = subplot(2,1,2);
imagesc(hca,x,z,ExB')
%imagesc(hca,x,z,pi1.scalar')
hca.YDir = 'normal';
hca.XLim = [140 210];
hca.YLim = [-12 12];
hca.Title.String = sprintf('n_boxes = %g',n_boxes);
hca.Title.Interpreter = 'none';
hcb = colorbar('peer',hca);


hold(hca,'on')
for ibox = 1:n_boxes      
  %hstar = plot(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'*k');
  %hpatch = patch(hca,[keep_boxes(ibox,1) keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2)],[keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2) keep_boxes(ibox,1)],'w');
  hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
  hpatch.FaceAlpha = 0;
  hpatch.LineWidth = 1;   
  
end
hold(hca,'off')   


%% Tight smaller boxes along separatrix
x_center = 150:0.5:205;
z_center = 0:0.5:10;
dx_box = 0.25;
dz_box = 0.25;

xlow = x_center-dx_box;
xhigh = x_center+dx_box;
zlow = z_center-dz_box;
zhigh = z_center+dz_box;

[XC,ZC] = meshgrid(x_center,z_center);

XC = reshape(XC,prod(size(XC)),1);
ZC = reshape(ZC,prod(size(ZC)),1);

nboxes = numel(XC);
Alim = [-24 -17];
[saddle_locations,saddle_values] = saddle(A);
dA_minus = 1.0;
dA_plus = 0.5;
Alim = max(saddle_values) + [-dA_minus dA_plus];
ind_keep = zeros(nboxes,1);
for ibox = 1:nboxes
  xind = find(abs(x-XC(ibox))==min(abs(x-XC(ibox))));
  zind = find(abs(z-ZC(ibox))==min(abs(z-ZC(ibox))));
  A(xind,zind);
  if A(xind,zind)>Alim(1) && A(xind,zind)<Alim(2) % keep
    ind_keep(ibox) = 1;
  else
    ind_keep(ibox) = 0;
  end
end

all_boxes = [XC-dx_box XC+dx_box ZC-dz_box ZC+dz_box];
keep_boxes = all_boxes(find(ind_keep==1),:);
n_boxes = size(keep_boxes,1);

figure(401)
hca = subplot(1,1,1);
imagesc(hca,x,z,ve1.x')
hca.XLim = [100 210];
hca.YLim = [-12 12];
hca.Title.String = sprintf('n_boxes = %g',n_boxes);
hca.Title.Interpreter = 'none';
hcb = colorbar('peer',hca);


hold(hca,'on')
for ibox = 1:n_boxes      
  %hstar = plot(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'*k');
  %hpatch = patch(hca,[keep_boxes(ibox,1) keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2)],[keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2) keep_boxes(ibox,1)],'w');
  hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
  hpatch.FaceAlpha = 0;
  hpatch.LineWidth = 1;   
  
end
hold(hca,'off')   

%disp(sprintf('%.3f %.3f %.3f %.3f \n'),keep_boxes)

%% Tight smaller boxes around largest electron flow
x_center = 140:0.5:205;
z_center = 0:0.5:10;
dx_box = 0.25;
dz_box = 0.25;

xlow = x_center-dx_box;
xhigh = x_center+dx_box;
zlow = z_center-dz_box;
zhigh = z_center+dz_box;

[XC,ZC] = meshgrid(x_center,z_center);

XC = reshape(XC,prod(size(XC)),1);
ZC = reshape(ZC,prod(size(ZC)),1);

nboxes = numel(XC);
var_lim = [-24 -17];
var = abs(ve1.par);
var_center = max(abs(var(:)));
dvar_minus = 1.2;
dvar_plus = 0.0;
var_lim = var_center + [-dvar_minus dvar_plus];
ind_keep = zeros(nboxes,1);
for ibox = 1:nboxes
  xind = find(abs(x-XC(ibox))==min(abs(x-XC(ibox))));
  zind = find(abs(z-ZC(ibox))==min(abs(z-ZC(ibox))));
  
  if var(xind,zind)>var_lim(1) && var(xind,zind)<var_lim(2) % keep
    ind_keep(ibox) = 1;
  else
    ind_keep(ibox) = 0;
  end
end

all_boxes = [XC-dx_box XC+dx_box ZC-dz_box ZC+dz_box];
keep_boxes = all_boxes(find(ind_keep==1),:);
n_boxes = size(keep_boxes,1);

figure(401)
hca = subplot(1,1,1);
imagesc(hca,x,z,var')
hca.XLim = [100 210];
hca.YLim = [-12 12];
hca.Title.String = sprintf('n_boxes = %g',n_boxes);
hca.Title.Interpreter = 'none';
hcb = colorbar('peer',hca);


hold(hca,'on')
for ibox = 1:n_boxes      
  %hstar = plot(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'*k');
  %hpatch = patch(hca,[keep_boxes(ibox,1) keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2)],[keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2) keep_boxes(ibox,1)],'w');
  hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
  hpatch.FaceAlpha = 0;
  hpatch.LineWidth = 1;   
  
end
hold(hca,'off')   

%% Tight smaller boxes around largest Bz
x_center = 140:0.5:205;
z_center = 0:0.5:10;
dx_box = 0.25;
dz_box = 0.25;

xlow = x_center-dx_box;
xhigh = x_center+dx_box;
zlow = z_center-dz_box;
zhigh = z_center+dz_box;

[XC,ZC] = meshgrid(x_center,z_center);

XC = reshape(XC,prod(size(XC)),1);
ZC = reshape(ZC,prod(size(ZC)),1);

nboxes = numel(XC);

var = abs(B.abs);
var_center = max(abs(var(:)));
dvar_minus = 0.4;
dvar_plus = 0.0;
var_lim = var_center + [-dvar_minus dvar_plus];
ind_keep = zeros(nboxes,1);
for ibox = 1:nboxes
  xind = find(abs(x-XC(ibox))==min(abs(x-XC(ibox))));
  zind = find(abs(z-ZC(ibox))==min(abs(z-ZC(ibox))));
  
  if var(xind,zind)>var_lim(1) && var(xind,zind)<var_lim(2) % keep
    ind_keep(ibox) = 1;
  else
    ind_keep(ibox) = 0;
  end
end

all_boxes = [XC-dx_box XC+dx_box ZC-dz_box ZC+dz_box];
keep_boxes = all_boxes(find(ind_keep==1),:);
n_boxes = size(keep_boxes,1);

figure(401)
hca = subplot(1,1,1);
imagesc(hca,x,z,var')
hca.XLim = [100 210];
hca.YLim = [-12 12];
hca.Title.String = sprintf('n_boxes = %g',n_boxes);
hca.Title.Interpreter = 'none';
hcb = colorbar('peer',hca);


hold(hca,'on')
for ibox = 1:n_boxes      
  %hstar = plot(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'*k');
  %hpatch = patch(hca,[keep_boxes(ibox,1) keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2)],[keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2) keep_boxes(ibox,1)],'w');
  hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
  hpatch.FaceAlpha = 0;
  hpatch.LineWidth = 1;   
  
end
hold(hca,'off')   

%% Pattern governed by A, PIC object
pic = df04;
ind = df04.twpelim(5000).it;
%ind = 26;
A = squeeze(pic(ind).A);
Bz = squeeze(pic(ind).Bz);
ni = squeeze(pic(ind).ni);
vex = squeeze(pic(ind).vex);


[saddle_locations,saddle_values] = saddle(A,'sort');

x_center = (192:0.2:205);
z_center = [0 0.2];
dx_box = 0.1;
dz_box = 0.1;
% x_center = 150:0.5:205;
% z_center = 0:0.5:10;
% dx_box = 0.25;
% dz_box = 0.25;

xlow = x_center-dx_box;
xhigh = x_center+dx_box;
zlow = z_center-dz_box;
zhigh = z_center+dz_box;

%[XC,ZC] = meshgrid(x_center,z_center);
[ZC,XC] = meshgrid(z_center,x_center);

XC = reshape(XC,prod(size(XC)),1);
ZC = reshape(ZC,prod(size(ZC)),1);

nboxes = numel(XC);
Alim = [-24 -17];
Alim = [-24 -18];
%Alim = [-20 -16];
Alim = [-24 saddle_values(1)*0.99];
Alim = [-25 saddle_values(1)*0.99];
Alim = [-23.8 0];
ind_keep = zeros(nboxes,1);
for ibox = 1:nboxes
  xind = find(abs(pic.xi-XC(ibox))==min(abs(pic.xi-XC(ibox))));
  zind = find(abs(pic.zi-ZC(ibox))==min(abs(pic.zi-ZC(ibox))));
  A(xind,zind);
  if A(xind,zind)>Alim(1) && A(xind,zind)<Alim(2) % keep
    ind_keep(ibox) = 1;
  else
    ind_keep(ibox) = 0;
  end
end

all_boxes = [XC-dx_box XC+dx_box ZC-dz_box ZC+dz_box];
keep_boxes = all_boxes(find(ind_keep==1),:);
n_boxes = size(keep_boxes,1);

figure(401)
hca = subplot(1,1,1);
imagesc(hca,pic.xi,pic.zi,ni')
%imagesc(hca,x,z,pi1.scalar')
hca.XLim = [130 240];
hca.YLim = [-10 10];
hca.Title.String = sprintf('twci = %g, twpe = %g, n_boxes = %g',pic.twci(ind),pic.twpe(ind),n_boxes);
hca.Title.Interpreter = 'none';
hcb = colorbar('peer',hca);


hold(hca,'on')
for ibox = 1:n_boxes      
  %hstar = plot(hca,[xlo xlo xhi xhi],[zlo zhi zhi zlo],'*k');
  %hpatch = patch(hca,[keep_boxes(ibox,1) keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2)],[keep_boxes(ibox,1) keep_boxes(ibox,2) keep_boxes(ibox,2) keep_boxes(ibox,1)],'w');
  hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
  hpatch.FaceAlpha = 0;
  hpatch.LineWidth = 1;   
  
end
hold(hca,'off')   