df04 = PIC('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/fields.h5');
ds04 = PICDist('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/dists.h5');
tr04 = PICTraj('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5');


%%


%% W (vdotE)
twci0 = 140;
limUstart = 1.5;
%tr = tr04.pass('mass',[0.5 1.5]).pass('x0',[0 220]).pass('t0',twci0+[-1 1]);
%tr = trl(find(trl.Ustart<0.5));
tr = trs(intersect(find([trs.t0] == 160),find(trs.Ustart < limUstart)));
%tr = trs(find([trs.t0] == 160));
%tr = trs;
tr = trs.find([trs.t0] == 160,trs.Ustart < limUstart);
tr = trs.find([trs.t0] == 160);
%tr = trs.find([trs.t0]==160,[trs.x0]>192,[trs.vy0]<0);
tr = trs.find([trs.t0]==160,[trs.vy0]<-0.0,[trs.vz0]<0.0,[trs.x0]>175);

W = zeros(tr.ntr,1); 
Wx = zeros(tr.ntr,1); 
Wy = zeros(tr.ntr,1);
Wz = zeros(tr.ntr,1);
vEx = zeros(tr.ntr,1); 
vEy = zeros(tr.ntr,1);
vEz = zeros(tr.ntr,1);
sumEx = zeros(tr.ntr,1);
sumEy = zeros(tr.ntr,1);
sumEz = zeros(tr.ntr,1);
ncross_pre = zeros(tr.ntr,1);
ncross_post = zeros(tr.ntr,1);

for itr = 1:numel(tr)
  % only include backwards in time from t0
  it_pre = find(tr(itr).t < tr(itr).t0);
  it_post = find(tr(itr).t > tr(itr).t0);
  itime = find(tr(itr).t > 0);
  itime = it_pre;
  nind(itr) = numel(itime);
  
  ncross_pre(itr) = tr(itr).select_inds(it_pre).ncross;
  ncross_post(itr) = tr(itr).select_inds(it_post).ncross;
  
  
  vEx(itr) = nansum(tr(itr).vx(itime).*tr(itr).Ex(itime));
  vEy(itr) = nansum(tr(itr).vy(itime).*tr(itr).Ey(itime));
  vEz(itr) = nansum(tr(itr).vz(itime).*tr(itr).Ez(itime));
  W_ = tr(itr).W; W(itr) = nansum(W_(itime));
  Wx_ = tr(itr).Wx; Wx(itr) = nansum(Wx_(itime));
  Wy_ = tr(itr).Wy; Wy(itr) = nansum(Wy_(itime));
  Wz_ = tr(itr).Wz; Wz(itr) = nansum(Wz_(itime));
  sumEx(itr) = nansum(tr(itr).Ex(itime));
  sumEy(itr) = nansum(tr(itr).Ey(itime));
  sumEz(itr) = nansum(tr(itr).Ez(itime));
end

ncrosses = tr.ncross;
ikeep = find(not(vEx==0));
tr = tr(ikeep);
vEx = vEx(ikeep);
vEy = vEy(ikeep);
vEz = vEz(ikeep);
vE = vEx + vEy + vEz;
% scale to total
vEx = vEx./vE;
vEy = vEy./vE;
vEz = vEz./vE;
%vE(vE==0) = [];
sc = 1e1;
sW = 3e1;
MarkerEdgeColor = 'none';
MarkerFaceColor = 'flat';
MarkerFaceAlpha = 0.7;

nrows = 5;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

isW = [];
if 1 % id of particle, number of bounces
  hca = h(isub); isub = isub + 1;  
  x0 = [tr.x0];
  vz0 = [tr.vz0];
  scatter(hca,x0,vz0,sW,ncrosses,'o','MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha)    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}^{z=0} (t>0)';
  hca.CLim = [0 10];
  hca.YLabel.String = 'v_{z} (t=t_0)';
  hold(hca,'on')
  for itr = 1:tr.ntr    
    text(hca,x0(itr),vz0(itr),sprintf('%g',tr(itr).id),'color',[0 0 0])
  end
  hold(hca,'off')
end
if 1 % id of particle, number of bounces, Nc(x,vy)
  hca = h(isub); isub = isub + 1;  
  x0 = [tr.x0];
  vy0 = [tr.vy0];
  scatter(hca,x0,vy0,sW,ncrosses,'o','MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha)  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}^{z=0} (t>0)';
  hca.CLim = [0 10];
  hca.YLabel.String = 'v_{y} (t=t_0)';
  hold(hca,'on')
  for itr = 1:tr.ntr    
    text(hca,x0(itr),vy0(itr),sprintf('%g',tr(itr).id),'color',[0 0 0])
  end
  hold(hca,'off')
end
if 1 % number of bounces
  hca = h(isub); isub = isub + 1;  
  x0 = [tr.x0];
  vz0 = [tr.vz0];
  scatter(hca,x0,vz0,sW,ncrosses,'o','MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha)  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}^{z=0} (t>0)';
  hca.CLim = [0 10];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 1 % number of bounces, Nc(x,vy)
  hca = h(isub); isub = isub + 1;  
  x0 = [tr.x0];
  vz0 = [tr.vy0];
  scatter(hca,x0,vz0,sW,ncrosses,'o','MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha)  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}^{z=0} (t>0)';
  hca.CLim = [0 10];
  hca.YLabel.String = 'v_{y} (t=t_0)';
end
if 0 % number of bounces t<t0
  hca = h(isub); isub = isub + 1;  
  x0 = [tr.x0];
  vz0 = [tr.vz0];
  scatter(hca,x0,vz0,sW,ncross_pre,'o')  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}^{z=0} (t<t_0)';
  hca.CLim = [0 10];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 0 % number of bounces t>t0
  hca = h(isub); isub = isub + 1;  
  x0 = [tr.x0];
  vz0 = [tr.vz0];
  scatter(hca,x0,vz0,sW,ncross_post,'o')  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}^{z=0} (t>t_0)';
  hca.CLim = [0 10];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 0 % Ustart, mostly to see if initial velocity is too high  
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],sW,tr.Ustart,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'U (t=t_{start})';
  hca.CLim = limUstart*[-1 1];
  hca.YLabel.String = 'v_{z}^{t=t_0}';
end
if 0 % number of bounces
  hca = h(isub); isub = isub + 1;
  ic1 = find(ncrosses == 1);
  x0 = [tr.x0];
  vz0 = [tr.z0];
  scatter(hca,x0(ic1),vz0(ic1),ncrosses(ic1)+10,ncrosses(ic1),'o')
  hold(hca,'on')
  ic1 = find(ncrosses > 1);
  x0 = [tr.x0];
  vz0 = [tr.z0];
  scatter(hca,x0(ic1),vz0(ic1),ncrosses(ic1)+10,ncrosses(ic1),'^')
  hold(hca,'off')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N z=0 crossings';
  hca.CLim = [0 10];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 0 % vxEx + vyEy + vzEz
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],abs(vE)*sc,vE)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vE (t<t_0)';
  hca.CLim = prctile(vE,95)*[0 1];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 0 % vxEx + vyEy + vzEz
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vy0],abs(vE)*sc,vE)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'vE (t<t_0)';
  hca.CLim = prctile(vE,95)*[0 1];
  hca.YLabel.String = 'v_{y} (t=t_0)';
end
if 0 % W(x,vz)
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],sW,W,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha);
  %hs.MarkerFaceColor = 'flat';
  %  hs.MarkerEdgeColor = [0 0 0];
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W (t<t_0)';
  hca.CLim = prctile(abs(W),95)*[-1 1];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 0 % W(x,vy)
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vy0],sW,W,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha);
  %hs.MarkerFaceColor = 'flat';
  %  hs.MarkerEdgeColor = [0 0 0];
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W (t<t_0)';
  hca.CLim = prctile(abs(W),95)*[-1 1];
  hca.YLabel.String = 'v_{y} (t=t_0)';
end
if 1 % Wx(x,vz)
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],sW,Wx,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W_x (t<t_0)';
  hca.CLim = prctile(abs(Wx),95)*[-1 1];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 1 % Wx(x,vy)
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vy0],sW,Wx,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W_x (t<t_0)';
  hca.CLim = prctile(abs(Wx),95)*[-1 1];
  hca.YLabel.String = 'v_{y} (t=t_0)';
end
if 1 % Wy(x,vz)
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],sW,Wy,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W_y (t<t_0)';
  hca.CLim = prctile(abs(Wy),95)*[-1 1];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 1 % Wy(x,vy)
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vy0],sW,Wy,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W_y (t<t_0)';
  hca.CLim = prctile(abs(Wy),95)*[-1 1];
  hca.YLabel.String = 'v_{y} (t=t_0)';
end
if 1 % Wz(x,vz)
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],sW,Wz,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W_z (t<t_0)';
  hca.CLim = prctile(abs(Wz),95)*[-1 1];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 1 % Wz(x,vy)
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vy0],sW,Wz,'MarkerFaceColor',MarkerFaceColor,'MarkerEdgeColor',MarkerEdgeColor,'MarkerFaceAlpha',MarkerFaceAlpha)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W_z (t<t_0)';
  hca.CLim = prctile(abs(Wz),95)*[-1 1];
  hca.YLabel.String = 'v_{y} (t=t_0)';
end
if 0 % Wx/Wy
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],sW,log10(Wx./Wy))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(W_x/W_y) (t<t_0)';
  hca.CLim = [-1 1];
  hca.YLabel.String = 'v_{z0} (t=t_0)';
end
if 0 % Wz/Wy
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],sW,log10(Wz./Wy))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(W_z/W_y) (t<t_0)';
  hca.CLim = [-1 1];
  hca.YLabel.String = 'v_{z0} (t=t_0)';
end
if 0 % Wz/Wx
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],sW,log10(Wz./Wx))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(W_z/W_x) (t<t_0)';
  hca.CLim = [-1 1];
  hca.YLabel.String = 'v_{z0} (t=t_0)';
end
if 0
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],abs(vEx)*sc,vEx)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_xE_x (t<t_0)';
  hca.CLim = prctile(vEx,95)*[-1 1];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 0
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],abs(vEy)*sc,vEy)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_yE_y (t<t_0)';
  hca.CLim = prctile(vEy,95)*[-1 1];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 0
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],abs(vEz)*sc,vEz)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_zE_z (t<t_0)';
  hca.CLim = prctile(vEz,95)*[-1 1];
  hca.YLabel.String = 'v_{z} (t=t_0)';
end
if 0 % vEz/vEy
  hca = h(isub); isub = isub + 1;
  scatter(hca,[tr.x0],[tr.vz0],abs(vEz)*sc,log10(vEz./vEy))
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(v_zE_z/v_yE_y) (t<t_0)';
  hca.CLim = [-1 1];
  hca.YLabel.String = 'v_{z0} (t=t_0)';
end

for ip = 1:npanels
  h(ip).Box = 'on';
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).YLim = [-1 1];
  if 0 % black background
    h(ip).Color = [0 0 0]+0.1;
    h(ip).GridColor = [0 0 0]+0.8;
  end
end
hlinks = linkprop(h,{'XLim','YLim'});
hlinks.Targets(1).YLim = 0.99*[-1 1];
hlinksW = linkprop(h(isW),{'CLim'});
%hlinks_c = linkprop(h(2:4),{'CLim'});
drawnow
colormap(pic_colors('waterfall'))
%colormap(pic_colors('candy'))
compact_panels(0.01)

%%
limUstart = 0.5;
tr = trs.find([trs.t0]==160).pass('mass',[0.5 1.5]).pass('x0',[0 220]);
tr = trl(find(trl.Ustart>0.5));
tr = trs.find([trs.t0] == 160,trs.Ustart > limUstart,trs.ncross == 1);
tr = trs.find([trs.t0] == 140,trs.Ustart < limUstart,trs.ncross < 3,trs.Ustop<0.7);
tr = trs.find([trs.t0] == 160,trs.Ustart < limUstart,trs.ncross < 3,trs.Ustart<0.7);
tr = trs.find([trs.t0] == 160,trs.Ustart < limUstart);
%tr = tr(1:5);

%tr = trs(find([trs.t0] == 160));
%tr = trs;

% for itr = 1:numel(tr)
%   % only include backwards in time from t0
%   it_pre = find(tr(itr).t < tr(itr).t0);
%   it_post = find(tr(itr).t > tr(itr).t0);
%   itime = find(tr(itr).t > 0);
%   itime = it_pre;
%   nind(itr) = numel(itime);
%   
%   ncross_pre(itr) = tr(itr).select_inds(it_pre).ncross;
%   ncross_post(itr) = tr(itr).select_inds(it_post).ncross;
%   
%   
%   vEx(itr) = nansum(tr(itr).vx(itime).*tr(itr).Ex(itime));
%   vEy(itr) = nansum(tr(itr).vy(itime).*tr(itr).Ey(itime));
%   vEz(itr) = nansum(tr(itr).vz(itime).*tr(itr).Ez(itime));
%   W_ = tr(itr).W; W(itr) = nansum(W_(itime));
%   Wx_ = tr(itr).Wx; Wx(itr) = nansum(Wx_(itime));
%   Wy_ = tr(itr).Wy; Wy(itr) = nansum(Wy_(itime));
%   Wz_ = tr(itr).Wz; Wz(itr) = nansum(Wz_(itime));
%   sumEx(itr) = nansum(tr(itr).Ex(itime));
%   sumEy(itr) = nansum(tr(itr).Ey(itime));
%   sumEz(itr) = nansum(tr(itr).Ez(itime));
% end
%
sp = 5;
nrows = 4;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

xlim = [120 220];
ylim = [-10 10];
clim = 3.3e-3*[-1 1];
isW = [];
if 0 % Wx vs x
  hca = h(isub); isub = isub + 1;  
  for itr  = 1:tr.ntr
    Wx = tr(itr).Wx;
    x = tr(itr).x;
    t = tr(itr).t;
    scatter(hca,x,Wx,sp,t,'o')
    if itr == 1, hold(hca,'on'); end
  end
  hold(hca,'off');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 't';
  hca.CLim = [60 210];
  %hca.YLabel.String = 'v_{z} (t=t_0)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % Wx(x,z)
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;  
  for itr  = 1:tr.ntr
    W = tr(itr).Wx;
    x = tr(itr).x;
    z = tr(itr).z;
    scatter(hca,x,z,sp,W,'o')
    if itr == 1, hold(hca,'on'); end
  end
  hold(hca,'off');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W_x';
  hca.CLim = prctile(abs(W),95)*[-1 1];
  %hca.YLabel.String = 'v_{z} (t=t_0)';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
end
if 1 % Wy(x,z)
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;  
  for itr  = 1:tr.ntr
    W = tr(itr).Wy;
    x = tr(itr).x;
    z = tr(itr).z;
    scatter(hca,x,z,sp,W,'o')
    if itr == 1, hold(hca,'on'); end
  end
  hold(hca,'off');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W_y';
  hca.CLim = prctile(abs(W),95)*[-1 1];
  %hca.YLabel.String = 'v_{z} (t=t_0)';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
end
if 1 % Wz vs x
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;  
  for itr  = 1:tr.ntr
    W = tr(itr).Wz;
    x = tr(itr).x;
    z = tr(itr).z;
    scatter(hca,x,z,sp,W,'o')
    if itr == 1, hold(hca,'on'); end
  end
  hold(hca,'off');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W_z';
  hca.CLim = prctile(abs(W),95)*[-1 1];
  %hca.YLabel.String = 'v_{z} (t=t_0)';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
end
if 1 % W vs x
  isW(end+1) = isub;
  hca = h(isub); isub = isub + 1;  
  for itr  = 1:tr.ntr
    W = tr(itr).W;
    x = tr(itr).x;
    z = tr(itr).z;
    scatter(hca,x,z,sp,W,'o')
    if itr == 1, hold(hca,'on'); end
  end
  hold(hca,'off');
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W';
  hca.CLim = prctile(abs(W),95)*[-1 1];
  %hca.YLabel.String = 'v_{z} (t=t_0)';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
end

for ip = 1:npanels
  h(ip).Box = 'on';
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).XLim = xlim;
  h(ip).YLim = ylim;
  h(ip).CLim = clim;
end
%hlinksW = linkprop(h(isW),{'XLim','YLim','CLim'});
h(1).Title.String = sprintf('n = %g',tr.ntr);
hlinks = linkprop(h,{'XLim','YLim'});
hlinksW = linkprop(h(isW),{'CLim'});
drawnow
colormap(pic_colors('waterfall'))
compact_panels(0.01)

%% Some scatter properties ()
limUstart = 0.4;
t0 = 160;
tr = trs.find([trs.t0] == t0,trs.Ustart < limUstart,[trs.x0]<175);
tr = trs.find([trs.t0] == t0,trs.Ustart < limUstart);
tinterp = 160;

sp = 5;
nrows = 2;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

ULim = [-2 2];
isUU = [];
doA = 0; stepA = 1; pic = df04; tstart = 60; tstop = 210;

isXZ = [];
markerStop = 's';
markerStart = 'o';
markerInter = '^';
if 0
  hca = h(isub); isub = isub + 1;
  scatter(hca,tr.ncross,abs(tr.coordgen('z',1)),30,tr.Ustop-tr.Ustart);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'N_{cross}';
  hca.YLabel.String = 'x_{t=60}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'U(t=210)-U(t=60)';
  hca.CLim = [0 2];
  colormap(hca,pic_colors('waterfall'))
end
if 0 % deltaU(Ncr,zstart)
  hca = h(isub); isub = isub + 1;
  tr_ = tr;
  %scatter(hca,tr_.ncross,abs(tr_.coordgen('z',1)),30,tr_.Ustop-tr_.Ustart,'o');
  hold(hca,'on')
  limUstartEdges = [0 0.1 0.2 0.3 0.4 0.5];
  %limUstartEdges = [0.25 0.5];
  legs = cell(numel(limUstartEdges)-1,1);
  markers = {'o','s','^','d','x'};
  for ilim = 1:numel(limUstartEdges)-1
    tr_ = tr.find(trs.Ustart > limUstartEdges(ilim),trs.Ustart < limUstartEdges(ilim+1));
    scatter(hca,tr_.ncross,abs(tr_.coordgen('z',1)),50,tr_.Ustop-tr_.Ustart,'Marker',markers{ilim});    
    legs{ilim} = sprintf('%.2f<U_{start}<%.2f',limUstartEdges(ilim),limUstartEdges(ilim+1));
  end
  legend(hca,legs,'location','best')
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'N_{cross}';
  hca.YLabel.String = 'z_{t=60}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'U(t=210) - U(t=60)';
  hca.CLim = [0 2];
  colormap(hca,pic_colors('waterfall'))
  hca.Box = 'on';
end
if 0 % Ncr(x,z), with interpolated
  isXZ(end+1) = isub;
  ms = 40; % marker size
  hca = h(isub); isub = isub + 1;
  scatter(hca,tr.xstart,tr.zstart,ms,tr.ncross,'Marker',markerStart)
  hold(hca,'on')
  %scatter(hca,tr.interp('x',tinterp),tr.interp('z',tinterp),ms,tr.ncross,'Marker',markerStop)
  scatter(hca,tr.interp('x',tinterp),tr.interp('z',tinterp),ms,tr.ncross,'Marker',markerInter)
  %scatter(hca,[tr.x0],[tr.z0],ms,tr.ncross,'Marker',markerInter)
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Box = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  hca.CLim = [0 10];
  colormap(hca,pic_colors('waterfall'))
  if doA % t0, separatrix
    clim = hca.CLim;
    pic_tmp = pic.twcilim(t0).xlim(hca.XLim).zlim(hca.YLim);
    [~,~,~,Axline] = pic_tmp.saddle;
    A = pic_tmp.A;
    %stepA = 1;
    levA = floor(min(A(:))/stepA)*stepA:stepA:ceil(max(A(:))/stepA)*stepA;
    levA = Axline*[1 1];      
    iAx = 1:5:pic_tmp.nx;
    iAz = 1:5:pic_tmp.nz;
    hold(hca,'on')
    contour(hca,pic_tmp.xi(iAx),pic_tmp.zi(iAz),A(iAx,iAz)',levA,'k')
    hold(hca,'off')
    if 1 % tstart
      pic_tmp = pic.twcilim(tstart).xlim(hca.XLim).zlim(hca.YLim);
      A = pic_tmp.A;
      iAx = 1:5:pic_tmp.nx;
      iAz = 1:5:pic_tmp.nz;
      hold(hca,'on')
      contour(hca,pic_tmp.xi(iAx),pic_tmp.zi(iAz),A(iAx,iAz)',levA,'k')
      hold(hca,'off')
    end
    if 1 % tstop
      pic_tmp = pic.twcilim(tstop).xlim(hca.XLim).zlim(hca.YLim);
      A = pic_tmp.A;
      iAx = 1:5:pic_tmp.nx;
      iAz = 1:5:pic_tmp.nz;
      hold(hca,'on')
      contour(hca,pic_tmp.xi(iAx),pic_tmp.zi(iAz),A(iAx,iAz)',levA,'k')
      hold(hca,'off')
    end
    hca.CLim = clim; 
  end  
end
if 0 % Ncr(x,z)
  isXZ(end+1) = isub;
  ms = 40; % marker size
  hca = h(isub); isub = isub + 1;
  scatter(hca,tr.xstart,tr.zstart,ms,tr.ncross,'Marker',markerStart)
  hold(hca,'on')
  scatter(hca,tr.xstop,tr.zstop,ms,tr.ncross,'Marker',markerStop)
  scatter(hca,[tr.x0],[tr.z0],ms,tr.ncross,'Marker',markerInter)
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Box = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  hca.CLim = [0 10];
  colormap(hca,pic_colors('waterfall'))
  if doA % t0, separatrix
    clim = hca.CLim;
    pic_tmp = pic.twcilim(t0).xlim(hca.XLim).zlim(hca.YLim);
    [~,~,~,Axline] = pic_tmp.saddle;
    A = pic_tmp.A;
    %stepA = 1;
    levA = floor(min(A(:))/stepA)*stepA:stepA:ceil(max(A(:))/stepA)*stepA;
    levA = Axline*[1 1];      
    iAx = 1:5:pic_tmp.nx;
    iAz = 1:5:pic_tmp.nz;
    hold(hca,'on')
    contour(hca,pic_tmp.xi(iAx),pic_tmp.zi(iAz),A(iAx,iAz)',levA,'k')
    hold(hca,'off')
    if 1 % tstart
      pic_tmp = pic.twcilim(tstart).xlim(hca.XLim).zlim(hca.YLim);
      A = pic_tmp.A;
      iAx = 1:5:pic_tmp.nx;
      iAz = 1:5:pic_tmp.nz;
      hold(hca,'on')
      contour(hca,pic_tmp.xi(iAx),pic_tmp.zi(iAz),A(iAx,iAz)',levA,'k')
      hold(hca,'off')
    end
    if 1 % tstop
      pic_tmp = pic.twcilim(tstop).xlim(hca.XLim).zlim(hca.YLim);
      A = pic_tmp.A;
      iAx = 1:5:pic_tmp.nx;
      iAz = 1:5:pic_tmp.nz;
      hold(hca,'on')
      contour(hca,pic_tmp.xi(iAx),pic_tmp.zi(iAz),A(iAx,iAz)',levA,'k')
      hold(hca,'off')
    end
    hca.CLim = clim; 
  end  
end
if 0 % delta U(x,z)
  isXZ(end+1) = isub;
  ms = 40; % marker size
  hca = h(isub); isub = isub + 1;
  scatter(hca,tr.xstart,tr.zstart,ms,tr.dU,'Marker',markerStart)
  hold(hca,'on')
  scatter(hca,tr.xstop,tr.zstop,ms,tr.Ustop-tr.Ustart,'Marker',markerStop)
  scatter(hca,[tr.x0],[tr.z0],ms,tr.Ustop-tr.Ustart,'Marker',markerInter)
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Box = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = '\Delta U';
  hca.CLim = [0 1.5];
  colormap(hca,pic_colors('waterfall'))
end
if 0 % Ustart
  isXZ(end+1) = isub;
  ms = 40; % marker size
  hca = h(isub); isub = isub + 1;
  scatter(hca,tr.xstart,tr.zstart,ms,tr.Ustart,'Marker',markerStart)
  hold(hca,'on')
  scatter(hca,tr.xstop,tr.zstop,ms,tr.Ustart,'Marker',markerStop)
  scatter(hca,[tr.x0],[tr.z0],ms,tr.Ustart,'Marker',markerInter)
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Box = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'U_{start}';
  hca.CLim = ceil(prctile(tr.Ustart,80)*10)/10*[0 1];
  colormap(hca,pic_colors('waterfall'))
end
if 0 % Ustop
  isXZ(end+1) = isub;
  ms = 40; % marker size
  hca = h(isub); isub = isub + 1;
  scatter(hca,tr.xstart,tr.zstart,ms,tr.Ustop,'Marker',markerStart)
  hold(hca,'on')
  scatter(hca,tr.xstop,tr.zstop,ms,tr.Ustop,'Marker',markerStop)
  scatter(hca,[tr.x0],[tr.z0],ms,tr.Ustop,'Marker',markerInter)
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Box = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'z';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'U_{stop}';
  hca.CLim = [0 2];
  colormap(hca,pic_colors('waterfall'))
end
if 0 % Ustart, vs Ustop  
  ms = 40; % marker size
  hca = h(isub); isub = isub + 1;
  scatter(hca,tr.Ustart,tr.Ustop,ms,tr.ncross,'Marker','o')
  %old(hca,'on')
  %scatter(hca,tr.xstop,tr.zstop,ms,tr.Ustop,'Marker',markerStop)
  %scatter(hca,[tr.x0],[tr.z0],ms,tr.Ustop,'Marker',markerInter)
  %hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Box = 'on';
  hca.XLabel.String = 'U_{start}';
  hca.YLabel.String = 'U_{start}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  hca.CLim = [0 20];
  colormap(hca,pic_colors('waterfall'))
end
if 1 % dU vs Wsum
  isUU(end+1) = isub;
  ms = 40; % marker size
  hca = h(isub); isub = isub + 1;
  scatter(hca,tr.Ustop-tr.Ustart,tr.Wsum,ms,tr.ncross,'Marker','o')
  h%old(hca,'on')
  %scatter(hca,tr.xstop,tr.zstop,ms,tr.Ustop,'Marker',markerStop)
  %scatter(hca,[tr.x0],[tr.z0],ms,tr.Ustop,'Marker',markerInter)
  %hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Box = 'on';
  hca.XLabel.String = '\Delta U';
  hca.YLabel.String = 'sum W';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  hca.CLim = [0 20];
  colormap(hca,pic_colors('waterfall'))
end
if 1 % dUx vs Wxsum
  isUU(end+1) = isub;
  ms = 40; % marker size
  hca = h(isub); isub = isub + 1;
  scatter(hca,tr.Uxstop-tr.Uxstart,tr.Wxsum,ms,tr.ncross,'Marker','o')
  %old(hca,'on')
  %scatter(hca,tr.xstop,tr.zstop,ms,tr.Ustop,'Marker',markerStop)
  %scatter(hca,[tr.x0],[tr.z0],ms,tr.Ustop,'Marker',markerInter)
  %hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Box = 'on';
  hca.XLabel.String = '\Delta U_x';
  hca.YLabel.String = 'sum W_x';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  hca.CLim = [0 20];
  colormap(hca,pic_colors('waterfall'))
end
if 1 % dUy vs Wysum
  isUU(end+1) = isub;
  ms = 40; % marker size
  hca = h(isub); isub = isub + 1;
  scatter(hca,tr.Uystop-tr.Uystart,tr.Wysum,ms,tr.ncross,'Marker','o')
  %old(hca,'on')
  %scatter(hca,tr.xstop,tr.zstop,ms,tr.Ustop,'Marker',markerStop)
  %scatter(hca,[tr.x0],[tr.z0],ms,tr.Ustop,'Marker',markerInter)
  %hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Box = 'on';
  hca.XLabel.String = '\Delta U_y';
  hca.YLabel.String = 'sum W_y';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  hca.CLim = [0 20];
  colormap(hca,pic_colors('waterfall'))
end
if 1 % dUz vs Wzsum
  isUU(end+1) = isub;
  ms = 40; % marker size
  hca = h(isub); isub = isub + 1;
  scatter(hca,tr.Uzstop-tr.Uzstart,tr.Wzsum,ms,tr.ncross,'Marker','o')
  %old(hca,'on')
  %scatter(hca,tr.xstop,tr.zstop,ms,tr.Ustop,'Marker',markerStop)
  %scatter(hca,[tr.x0],[tr.z0],ms,tr.Ustop,'Marker',markerInter)
  %hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Box = 'on';
  hca.XLabel.String = '\Delta U_z';
  hca.YLabel.String = 'sum W_z';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  hca.CLim = [0 20];
  colormap(hca,pic_colors('waterfall'))
end


hlinksXZ = linkprop(h(isXZ),{'XLim','YLim'});
hlinksUU = linkprop(h(isUU),{'XLim','YLim'}); axis(h(isUU),'equal')

%% Some timelines with multiple particles
limUstart = 0.5;
tr = trs.find([trs.t0] == 60,trs.Ustart < limUstart);

newt = 60:10:210; % for downsampling

sp = 5;
nrows = 4;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % U(t)
  hca = h(isub); isub = isub + 1;
  for itr = 1:tr.ntr
    plot(hca,tr(itr).t,tr(itr).U);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 't';
  hca.YLabel.String = 'U';
end
if 0 % U(t) downsampled
  hca = h(isub); isub = isub + 1;
  for itr = 1:tr.ntr  
    newvar = resample_average(tr(itr).t,tr(itr).U,newt);
    plot(hca,newt,newvar);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 't';
  hca.YLabel.String = 'U';
end
if 0 % dU(t), this makes no sense, because the time step is not uniform
  hca = h(isub); isub = isub + 1;
  for itr = 1:tr.ntr
    plot(hca,tr(itr).t(1:end-1),diff(tr(itr).U));
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 't';
  hca.YLabel.String = '\Delta U';
end
if 0 % dU(t) downsampled, this makes no sense, because the time step is not uniform
  hca = h(isub); isub = isub + 1;
  for itr = 1:tr.ntr  
    newvar = resample_average(tr(itr).t,tr(itr).U,newt);
    plot(hca,newt(1:end-1)+0.5*diff(newt),diff(newvar));
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 't';
  hca.YLabel.String = '\Delta U';
end
if 1 % Wx(t)
  hca = h(isub); isub = isub + 1;
  for itr = 1:tr.ntr
    plot(hca,tr(itr).t,tr(itr).Wx);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 't';
  hca.YLabel.String = 'W_x';
end
if 1 % Wy(t)
  hca = h(isub); isub = isub + 1;
  for itr = 1:tr.ntr
    plot(hca,tr(itr).t,tr(itr).Wy);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 't';
  hca.YLabel.String = 'W_y';
end
if 1 % Wz(t)
  hca = h(isub); isub = isub + 1;
  for itr = 1:tr.ntr
    plot(hca,tr(itr).t,tr(itr).Wz);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 't';
  hca.YLabel.String = 'W_z';
end

hlinks = linkprop(h,{'XLim'});
compact_panels(0.01)

%% Make reduced velocity distributions
ncross_all = trs.ncross;
tr = trs.find([trs.t0] == 160).lim('z',2+[-0.2 0.2],'t',170+10*[-1 1]);
tr = trs.find([trs.t0] == 160,[trs.z0] == 2).lim('z',2+[-0.2 0.2]);
tr = trs.find([trs.t0]==160,[trs.z0]>1,trs.ncross<10).lim('z',2+[-0.2 0.2]);

ncrlim = [0 5];

sp = 5;
nrows = 4;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 0 % x vz
  hca = h(isub); isub = isub + 1;
  for itr = 1:tr.ntr
    plot(hca,tr(itr).x,tr(itr).vz,'.');
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_z';
end
if 0 % x vy, plot
  hca = h(isub); isub = isub + 1;
  for itr = 1:tr.ntr
    plot(hca,tr(itr).x,tr(itr).vy,'.');
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_y';
end
if 1 % Ncr(x,vz), all ncross, not just in selected subset
  hca = h(isub); isub = isub + 1; 
  crange = unique([tr.z0]); crange = [min(crange) max(crange)];
  crange = ncrlim;
  cmap = pic_colors('waterfall');
  for itr = 1:tr.ntr    
    color = cmap2color(cmap,crange,trs(tr(itr).id).ncross);
    scatter(hca,tr(itr).x,tr(itr).vz,sp,color);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  colormap(hca,cmap)
  hca.CLim = crange;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_z';
end
if 0 % x vz, scatter
  hca = h(isub); isub = isub + 1; 
  crange = unique([tr.z0]); crange = [min(crange) max(crange)];
  crange = ncrlim;
  cmap = pic_colors('waterfall');
  for itr = 1:tr.ntr
    color = cmap2color(cmap,crange,tr(itr).ncross);
    scatter(hca,tr(itr).x,tr(itr).vz,sp,color);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  colormap(hca,cmap)
  hca.CLim = crange;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_z';
end
if 1 % x vy, scatter
  hca = h(isub); isub = isub + 1; 
  crange = unique([tr.z0]); crange = [min(crange) max(crange)];
  crange = ncrlim;
  cmap = pic_colors('waterfall');
  for itr = 1:tr.ntr
    color = cmap2color(cmap,crange,tr(itr).ncross);
    scatter(hca,tr(itr).x,tr(itr).vy,sp,color);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  colormap(hca,cmap)
  hca.CLim = crange;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_y';
end
if 0 % z0(x,vy), scatter
  hca = h(isub); isub = isub + 1; 
  crange = unique([tr.z0]); crange = [min(crange) max(crange)];
  crange = [0 10];
  cmap = pic_colors('waterfall');
  for itr = 1:tr.ntr
    color = cmap2color(cmap,crange,tr(itr).z0);
    scatter(hca,tr(itr).x,tr(itr).vy,sp,color);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'z_{0}';
  colormap(hca,cmap)
  hca.CLim = crange;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_y';
end
if 0 % xstop(x,vy), scatter
  hca = h(isub); isub = isub + 1; 
  crange = unique([tr.xstop]); crange = [min(crange) max(crange)];
  %crange = [0 10];
  cmap = pic_colors('waterfall');
  for itr = 1:tr.ntr
    color = cmap2color(cmap,crange,tr(itr).xstop);
    scatter(hca,tr(itr).x,tr(itr).vy,sp,color);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'x_{stop}';
  colormap(hca,cmap)
  hca.CLim = crange;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_y';
end

%% Make reduced velocity distributions
ncross_all = trs.ncross;
tr = trs.find([trs.t0] == 160).lim('z',2+[-0.2 0.2],'t',170+10*[-1 1]);
tr = trs.find([trs.t0] == 160,[trs.z0] == 2).lim('z',2+[-0.2 0.2]);
tr = trs.find([trs.t0]==160,[trs.z0]>1,trs.ncross<10).lim('z',2+[-0.2 0.2]);
tr = trs.find([trs.t0]==160,[trs.z0]>1,trs.ncross<10).lim('z',2+[-0.1 0.1]).lim('x',170+0.3*[-1 1]);

ncrlim = [0 5];

sp = 5;
nrows = 2;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 0 % x vz
  hca = h(isub); isub = isub + 1;
  for itr = 1:tr.ntr
    plot(hca,tr(itr).vx,tr(itr).vy,'.');
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_z';
end
if 0 % x vy, plot
  hca = h(isub); isub = isub + 1;
  for itr = 1:tr.ntr
    plot(hca,tr(itr).x,tr(itr).vy,'.');
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_y';
end
if 1 % Ncr(vx,vz), all ncross, not just in selected subset,scatter
  hca = h(isub); isub = isub + 1; 
  crange = unique([tr.z0]); crange = [min(crange) max(crange)];
  crange = ncrlim;
  cmap = pic_colors('waterfall');
  for itr = 1:tr.ntr    
    color = cmap2color(cmap,crange,trs(tr(itr).id).ncross);
    scatter(hca,tr(itr).vx,tr(itr).vz,sp,color);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  colormap(hca,cmap)
  hca.CLim = crange;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
end
if 1 % Ncr(vx,vz), all ncross, not just in selected subset,scatter
  hca = h(isub); isub = isub + 1; 
  crange = unique([tr.z0]); crange = [min(crange) max(crange)];
  crange = ncrlim;
  cmap = pic_colors('waterfall');
  for itr = 1:tr.ntr    
    color = cmap2color(cmap,crange,trs(tr(itr).id).ncross);
    crange = [0 tr.ntr];
    color = cmap2color(cmap,crange,itr);
    scatter(hca,tr(itr).vy,tr(itr).vz,sp,color);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  colormap(hca,cmap)
  hca.CLim = crange;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'v_y';
  hca.YLabel.String = 'v_z';
end
if 0 % x vz, scatter
  hca = h(isub); isub = isub + 1; 
  crange = unique([tr.z0]); crange = [min(crange) max(crange)];
  crange = ncrlim;
  cmap = pic_colors('waterfall');
  for itr = 1:tr.ntr
    color = cmap2color(cmap,crange,tr(itr).ncross);
    scatter(hca,tr(itr).x,tr(itr).vz,sp,color);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  colormap(hca,cmap)
  hca.CLim = crange;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_z';
end
if 0 % x vy, scatter
  hca = h(isub); isub = isub + 1; 
  crange = unique([tr.z0]); crange = [min(crange) max(crange)];
  crange = ncrlim;
  cmap = pic_colors('waterfall');
  for itr = 1:tr.ntr
    color = cmap2color(cmap,crange,tr(itr).ncross);
    scatter(hca,tr(itr).x,tr(itr).vy,sp,color);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'N_{cross}';
  colormap(hca,cmap)
  hca.CLim = crange;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_y';
end
if 0 % z0(x,vy), scatter
  hca = h(isub); isub = isub + 1; 
  crange = unique([tr.z0]); crange = [min(crange) max(crange)];
  crange = [0 10];
  cmap = pic_colors('waterfall');
  for itr = 1:tr.ntr
    color = cmap2color(cmap,crange,tr(itr).z0);
    scatter(hca,tr(itr).x,tr(itr).vy,sp,color);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'z_{0}';
  colormap(hca,cmap)
  hca.CLim = crange;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_y';
end
if 0 % xstop(x,vy), scatter
  hca = h(isub); isub = isub + 1; 
  crange = unique([tr.xstop]); crange = [min(crange) max(crange)];
  %crange = [0 10];
  cmap = pic_colors('waterfall');
  for itr = 1:tr.ntr
    color = cmap2color(cmap,crange,tr(itr).xstop);
    scatter(hca,tr(itr).x,tr(itr).vy,sp,color);
    if itr == 1; hold(hca,'on'); end
  end
  hold(hca,'off');   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'x_{stop}';
  colormap(hca,cmap)
  hca.CLim = crange;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v_y';
end



