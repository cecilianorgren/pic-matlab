%% See what distributions we have already
timestep = 03000;
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
pic = df04n.twpelim(5000);
%ind = gf05.twpelim(4000).it;
ind = pic.it;
%ind = 26;
ind = 1;
A = squeeze(pic(ind).A);
Bz = squeeze(pic(ind).Bz);
ni = squeeze(pic(ind).ni);
vex = squeeze(pic(ind).vex);
viz = squeeze(pic(ind).viz);
viycold = squeeze(pic(ind).vy([3 5]));
Ez = squeeze(pic(ind).Ez);



[saddle_locations,saddle_values] = saddle(A,'sort');
xx = pic(pic.nt).x_xline;
%%
x_center = (140:0.2:199);
x_center = fix(xx) + [-8 -4 0];
%x_center = (90:0.2:170);
z_center = [4 6]; % 0.4 0.8
z_center = [-4:0.1:4];
dx_box = 0.05;
dz_box = 0.05;
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
Alim = [-24 -18.5];
Alim = [-24 0];
Alim = [-26 0];
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
hca = subplot(2,1,1);
imagesc(hca,pic.xi,pic.zi,viycold')
%imagesc(hca,x,z,pi1.scalar')
hca.XLim = [140 210];
hca.YLim = [-15 15];
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

hca = subplot(2,1,2);
imagesc(hca,pic.xi,pic.zi,Ez')
%imagesc(hca,x,z,pi1.scalar')
hca.XLim = [140 210];
hca.YLim = [-15 15];
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

%% Pattern governed by A, distributions along a field line, PIC object
xlim = [150 200];
zlim = [-1 10];
twpe = 7000;
pic = nobg.twpelim(10000).xlim(xlim).zlim(zlim);
%ind = 26;
A = squeeze(pic.A);
Bz = squeeze(pic.Bz);
ni = squeeze(pic.ni);
vex = squeeze(pic.vex);
viz = squeeze(pic.viz);

Avals = -20;
if numel(Avals) == 1  
  S = contourcs(pic.xi,pic.zi,A',Avals*[1 1]);
else
  S = contourcs(pic.xi,pic.zi,A',Avals);
end

% interpolate/downpolate to fewer point
nq = 200;
d_arc_x = diff(S.X);
d_arc_y = diff(S.Y);
d_arc_distance = sqrt(d_arc_x.^2 + d_arc_y.^2);
arc_distance = [0 cumsum(d_arc_distance)];
d_arcdist = 0.1;
new_arc_distance = arc_distance(1):d_arcdist:arc_distance(end); % try to get atleast one box at start
%rel_arc_distance = arc_distance/arc_distance(end);
%new_rel_arc_distance = linspace(rel_arc_distance(1),rel_arc_distance(end),nq);
% First resample new_arc_distance so that one instance of it appears at z =
% 0 (so that it's top-bot symmetrical).
%

x_new = interp1(rel_arc_distance,S.X,new_rel_arc_distance);
z_new = interp1(rel_arc_distance,S.Y,new_rel_arc_distance);
    

x_center = x_new;
z_center = z_new;

dx_box = 0.1;
dz_box = 0.1;

xlow = x_center-dx_box;
xhigh = x_center+dx_box;
zlow = z_center-dz_box;
zhigh = z_center+dz_box;


n_boxes = numel(x_center);
ind_keep = ones(nboxes,1);
%
if 0
for ibox = 1:nboxes
  xind = find(abs(pic.xi-XC(ibox))==min(abs(pic.xi-XC(ibox))));
  zind = find(abs(pic.zi-ZC(ibox))==min(abs(pic.zi-ZC(ibox))));
  A(xind,zind);
%   if A(xind,zind)>Alim(1) && A(xind,zind)<Alim(2) % keep
%     ind_keep(ibox) = 1;
%   else
%     ind_keep(ibox) = 0;
%   end
end
end

all_boxes = [tocolumn(x_center)-dx_box tocolumn(x_center)+dx_box tocolumn(z_center)-dz_box tocolumn(z_center)+dz_box];
%keep_boxes = all_boxes(find(ind_keep==1),:);
keep_boxes = all_boxes;
n_boxes = size(keep_boxes,1);%size(keep_boxes,1);
%
figure(401)
h = setup_subplots(2,1);
isub = 1;

hca = h(isub); isub = isub + 1;
imagesc(hca,pic.xi,pic.zi,ni')
%imagesc(hca,x,z,pi1.scalar')
%hca.XLim = [60 240];
%hca.YLim = [-10 10];
hca.Title.String = sprintf('twci = %g, twpe = %g, n_boxes = %g',pic.twci,pic.twpe,n_boxes);
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

hca = h(isub); isub = isub + 1;
imagesc(hca,pic.xi,pic.zi,vex')
%imagesc(hca,x,z,pi1.scalar')
%hca.XLim = [60 240];
%hca.YLim = [-10 10];
hca.Title.String = sprintf('twci = %g, twpe = %g, n_boxes = %g',pic.twci,pic.twpe,n_boxes);
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

hlinks = linkprop(h,{'XLim','YLim'});

%% Along aline, typically a field line
pic = nobg.twpelim(10000).xlim([129 210]).zlim([-15 20]);
% Load data
x = pic.xi;
z = pic.zi;
%Bx = pic.Bx;
%By = pic.By;
%Bz = pic.Bz;
A = pic.A; 
%Epar = pic.Epar;
%Ez = pic.Ez;


x0 = 150;
z0 = 14.5;
nsteps = 50.01;

x0 = 150;
z0 = 13.0;
nsteps = 50.01;

%x0 = 140;
%z0 = 12;

%x0 = 150;
%z0 = 10;
%nsteps = 40.99; % if nsteps is not an even number, integration will stop 
                % when the fieldline arclength is above nsteps
                
dx = 0.05;
dy = 0.05;
dz = 0.05;
%nsteps = 40.99; % if nsteps is not an even number, integration will stop 
                % when the fieldline arclength is above nsteps
                
if 1 % field line
  tic;
  [linearc,linex,liney,linez,linebx,lineby,linebz] = fieldline(x0,z0,x,z,Bx,By,Bz,dx,dy,dz,nsteps);
  toc;
  % Interpolate line to desired values
  linearc_new = linearc(1):0.2:linearc(end);
  x_center = interp1(linearc,linex,linearc_new);
  z_center = interp1(linearc,linez,linearc_new);
else % manual
  x_center = 154.6;
  z_center = 4*zeros(numel(x_center));
end
dx_box = 0.1;
dz_box = 0.1;

xlow  = x_center - dx_box;
xhigh = x_center + dx_box;
zlow  = z_center - dz_box;
zhigh = z_center + dz_box;

keep_boxes = [xlow' xhigh' zlow' zhigh'];
keep_boxes(x_center<x0,:) = [];
n_boxes = size(keep_boxes,1);
 
% Plot
figure(402)
h = setup_subplots(3,1);
isub = 1;
doA = 1;

if 1 % Epar
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.Ez')
  hca.Title.String = sprintf('twci = %g, twpe = %g, n_boxes = %g',pic.twci,pic.twpe,n_boxes);
  hca.Title.Interpreter = 'none';
  hcb = colorbar('peer',hca);
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.YDir = 'normal';

  hold(hca,'on')
  for ibox = 1:n_boxes
    hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
    hpatch.FaceAlpha = 0;
    hpatch.LineWidth = 1;   
  end
  if doA    
    hold(hca,'on')
    clim = hca.CLim;
    levA = floor(min(A(:))):1:ceil(max(A(:)));
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;
    contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(hca,'off')
    hca.CLim = clim;
  end
  hold(hca,'off')   
end
if 1 % viy
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.viy')
  hca.Title.String = sprintf('twci = %g, twpe = %g, n_boxes = %g',pic.twci,pic.twpe,n_boxes);
  hca.Title.Interpreter = 'none';
  hcb = colorbar('peer',hca);
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.YDir = 'normal';

  hold(hca,'on')
  for ibox = 1:n_boxes
    hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
    hpatch.FaceAlpha = 0;
    hpatch.LineWidth = 1;   
  end
  if doA    
    hold(hca,'on')
    clim = hca.CLim;
    levA = floor(min(A(:))):1:ceil(max(A(:)));
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;
    contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(hca,'off')
    hca.CLim = clim;
  end
  hold(hca,'off')   
end
if 1 % ni
  hca = h(isub); isub = isub + 1;
  imagesc(hca,pic.xi,pic.zi,pic.ni')
  hca.Title.String = sprintf('twci = %g, twpe = %g, n_boxes = %g',pic.twci,pic.twpe,n_boxes);
  hca.Title.Interpreter = 'none';
  hcb = colorbar('peer',hca);
  colormap(hca,pic_colors('waterfall'))
  hca.CLim = [0 0.2];
  hca.YDir = 'normal';

  hold(hca,'on')
  for ibox = 1:n_boxes
    hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
    hpatch.FaceAlpha = 0;
    hpatch.LineWidth = 1;   
  end
  if doA    
    hold(hca,'on')
    clim = hca.CLim;
    levA = floor(min(A(:))):1:ceil(max(A(:)));
    iAx = 1:4:pic.nx;
    iAz = 1:4:pic.nz;
    contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',levA,'k');
    hold(hca,'off')
    hca.CLim = clim;
  end
  hold(hca,'off')   
end

%% Along a line, interpolated from clicked points
%% plotmap
twpe = 9000;
xlim = [80 210];
zlim = [0 15];
pic = nobg.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {'ni';'Epar'};
clims = {[0 1],[-1 1]}';
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmaps = {cmapwa,cmapbr}';

h = pic.plot_map(varstrs,'A',1,'clim',clims,'cmap',cmaps);

[x,z] = ginput();

arcline = [0; cumsum(sqrt(diff(x).^2 + diff(z).^2))];
darc = 0.2;
arcline_new = arcline(1):darc:arcline(end);

method = 'spline';
x_center = interp1(arcline,x,arcline_new,method);
z_center = interp1(arcline,z,arcline_new,method);

%x_center = interp1(linearc,linex,linearc_new);
%z_center = interp1(linearc,linez,linearc_new);
dx_box = 0.1;
dz_box = 0.1;

xlow  = x_center - dx_box;
xhigh = x_center + dx_box;
zlow  = z_center - dz_box;
zhigh = z_center + dz_box;

keep_boxes = [xlow' xhigh' zlow' zhigh'];
%keep_boxes(x_center<x0,:) = [];
n_boxes = size(keep_boxes,1);

for ip = 1:numel(h)
  hca = h(ip);
  hold(hca,'on')
  for ibox = 1:n_boxes      
    hpatch = patch(hca,keep_boxes(ibox,[1 1 2 2]),keep_boxes(ibox,[3 4 4 3]),'w');
    hpatch.FaceAlpha = 0;
    hpatch.LineWidth = 1;   
  end  
  irf_legend(hca,sprintf('nboxes = %g',n_boxes),[0.98 0.98])
end