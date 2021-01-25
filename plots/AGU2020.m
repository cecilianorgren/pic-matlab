%% Initial conditions
pic = no02m;
comp = 'z';
twpe = 18000;
xlim = 53+0.5*[-1 1];
zlim = pic.zi([1 end]);
pic = pic.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');

varstrs = {{'n([1])','n([3 5])','n([1 3 5])'}};

colors = pic_colors('matlab');

h = pic.plot_line(comp,varstrs,'smooth',20);

h.Children(1).Color = colors(3,:);
h.Children(2).Color = colors(1,:);
h.Children(3).Color = colors(2,:);
h.XLim = 10*[-1 1];
legend(h,{'Hot','Cold','Total'},'location','best')
h.Position = [0.10 0.2 0.8 0.7];
h.YLabel.String = 'n (n_0)';
h.FontSize = 14;
h.Title.String = 'Initial density configuration';

%% Movie, ncold/ntot
twpe = 15000:100:24000;
xlim = no02m.xi([1 end])+[40 -40]';
zlim = [-8 8];
pic = no02m.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');

varstrs = {'n([3 5])./n([1 3 5])'}';
clims = {[0 1]};
cmaps = {flipdim(cmapbr,1)};
cbarlabels = {'n_{cold}/n_{tot}'};
filename = [printpath 'no02m_ncold_blue'];

pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'cbarlabels',cbarlabels,'filename',filename);

%% Movie, ncold, Tcold
twpe = 17000;%:100:24000;
xlim = no02m.xi([1 end])+[40 -40]';
zlim = [-8 8];
pic = no02m.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');
cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmapth = pic_colors('thermal');

varstrs = {'n([3 5])','t([3 5])'}';
clims = {[0 0.5]*0.99,[0 0.5]*0.99};
cmaps = {cmapth,cmapth};
cbarlabels = {'n_{i,cold}','T_{i,cold}'};
filename = [printpath 'no02m_ncold_tcold_'];

pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'cbarlabels',cbarlabels,'filename',filename);


%% Some overview
xlim = no02m.xi([1 end])+[40 -40]';
zlim = [-10 10]*0.99;
pic = no02m.twcilim(120).xlim(xlim).zlim(zlim);

xlim_plot = no02m.xi([1 end])+[65 -65]';
ylim_n = [0 1];
cmapbr = pic_colors('blue_red');
cmapth = pic_colors('thermal');
gray = [0.5 0.5 0.5];

varstrs = {'Bz','n([3 5])./n([1 3 5])','n([3])','t(3)'}';
cmaps = {cmapbr,flipdim(cmapbr,1),[cmapth],[cmapth],[gray; cmapbr; gray]};
clims = {[-1 1]*0.99,[0 1]*0.99,[0 0.5]*0.99,[0 0.5]*0.99,[-1 1]*0.99};
cbarlabels = {'B_z','n_{cold}/n_{tot}','n_{cold}^{top}','T_{cold}^{top}','v_{y,cold}^{top}'};

varstrs = {'Bz','n([3 5])./n([1 3 5])','p([3 5])./p([1 3 5])','n([3])','t(3)'}';
cmaps = {cmapbr,flipdim(cmapbr,1),flipdim(cmapbr,1),[cmapth],[cmapth],[gray; cmapbr; gray]};
clims = {[-1 1]*0.99,[0 1]*0.99,[0 1]*0.99,[0 0.5]*0.99,[0 0.5]*0.99,[-1 1]*0.99};
cbarlabels = {'B_z','n_{cold}/n_{tot}','p_{cold}/p_{tot}','n_{cold}^{top}','T_{cold}^{top}','v_{y,cold}^{top}'};


h = pic.plot_map(varstrs,'A',0.5,'cmap',cmaps,'clim',clims,'cbarlabels',cbarlabels,'smooth',10);  

fontsize = 12;
nvars = numel(varstrs);
c_eval('h(?).Position(3) = 0.7;',1:nvars)
c_eval('h(?).FontSize = fontsize;',1:nvars)

%% Some overview
xlim = no02m.xi([1 end])+[40 -40]';
zlim = [-10 10];
pic = no02m.twcilim(120).xlim(xlim).zlim(zlim);

xlim_plot = no02m.xi([1 end])+[65 -65]';
ylim_n = [0 1];
cmapbr = pic_colors('blue_red');
cmapth = pic_colors('thermal');

nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % Bx
  hca = h(isub); isub = isub + 1;
  h_ = pic.plot_map(hca,{'Bz'},'A',1,'cmap',{cmapbr},'cbarlabels',{'B_z'});  
end
if 1 % ncold/ntot
  hca = h(isub); isub = isub + 1;
  h_ = pic.plot_map(hca,{'n([3 5])./n([1 3 5])'},'A',1,...
    'cbarlabels',{'n_{cold}/n_{tot}'},...
    'clim',{[0 1]},'cmap',{flipdim(cmapbr,1)});
end
if 1 % ncold top
  hca = h(isub); isub = isub + 1;
  h_ = pic.plot_map(hca,{'n([3])'},'A',1,...
    'cbarlabels',{'n_{cold}^{top}'},...
    'clim',{[0 0.5]},'cmap',{cmapth});
end


%% Reduced distributions.
xlim = no02m.xi([1 end])+[40 -40]';
zlim = [-10 10];
pic = no02m.twcilim(120).xlim(xlim).zlim(zlim);

xlim_plot = no02m.xi([1 end])+[65 -65]';
ylim_n = [0 1];
cmapbr = pic_colors('blue_red');
cmapth = pic_colors('thermal');

nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % Bx
  hca = h(isub); isub = isub + 1;
  h_ = pic.plot_map(hca,{'Bz'},'A',1,'cmap',{cmapbr},'cbarlabels',{'B_z'});  
end
if 1 % ncold/ntot
  hca = h(isub); isub = isub + 1;
  h_ = pic.plot_map(hca,{'n([3 5])./n([1 3 5])'},'A',1,...
    'cbarlabels',{'n_{cold}/n_{tot}'},...
    'clim',{[0 1]},'cmap',{flipdim(cmapbr,1)});
end
if 1 % ncold top
  hca = h(isub); isub = isub + 1;
  h_ = pic.plot_map(hca,{'n([3])'},'A',1,...
    'cbarlabels',{'n_{cold}^{top}'},...
    'clim',{[0 0.5]},'cmap',{cmapth});
end
%%
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 3;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 1; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.0;
doPhi = 0; colorPhi = [0.5 0.5 0];
doSep = 1; %sep = no02m.twpelim(twpe).separatrix_location;

cmap_dist = pic_colors('waterfall');
legs = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)'};

freds = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
%freds = {fred35_z0,fred35_z0,fred35_z0};
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.fv' labstr]);
    pcolor(hca,fred.x,fred.v,log10(fredplot)')
    shading(hca,'flat')
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy4')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('%s  z = %g',legs{isub-1},unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',16)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 16;
    
    if doSep
      hold(hca,'on')      
      zz = unique(fred.z);
      %[~,iz] = min(abs(sep.z-zz));
      [PKS,LOCS] = findpeaks(-abs(sep.z-zz),sep.x,'minpeakprominence',1);
      %xx = sep.x(iz);
      for iloc = 1:numel(LOCS)
        hSep = plot(hca,LOCS(iloc)*[1 1],hca.YLim,'color',0*colorExB,'linewidth',1,'linestyle',':');
      end
      hold(hca,'off')      
    end
    if doE
      hold(hca,'on')
      xx = eval(['x_z' num2str(unique(fred.z))]);
      vv = eval(['E' labstr '_z' num2str(unique(fred.z))]);
      plot(hca,xx,vv,'color',colorE,'linewidth',0.5,'linestyle','-')            
      hold(hca,'off')
    end
    if 0*doV
      hold(hca,'on')
      plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
      hold(hca,'off')
    end
    if doExB % doExB
      hold(hca,'on')
      xx = eval(['x_z' num2str(unique(fred.z))]);
      vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
      plot(hca,xx,vv,'color',colorExB,'linewidth',0.5,'linestyle','-')
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end
  end
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align

ip = 1;
hpos = h(ip).Position;
hb(ip) = colorbar('peer',h(ip),'location','northoutside');
hb(ip).YLabel.String = 'log_{10}f_{ic}(x,v_x)';
h(ip).Position = hpos;
hb(ip).FontSize = 16;

ip = 2;
hpos = h(ip).Position;
hb(ip) = colorbar('peer',h(ip),'location','northoutside');
hb(ip).YLabel.String = 'log_{10}f_{ic}(x,v_y)';
h(ip).Position = hpos;
hb(ip).FontSize = 16;

ip = 3;
hpos = h(ip).Position;
hb(ip) = colorbar('peer',h(ip),'location','northoutside');
hb(ip).YLabel.String = 'log_{10}f_{ic}(x,v_z)';
h(ip).Position = hpos;
hb(ip).FontSize = 16;

for ip = 1:npanels
  h(ip).Position(2) = h(ip).Position(2) - 0.04;
  h(ip).Position(3) = 0.238;
  h(ip).FontSize = 16;
end
h(1).CLim = 0.99*[-4 2];
h(1).YLim = 0.99*4*[-1 1];

%% Reduced distribution plot, with trajectory initial positions

% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

tr = tr100(783:917);
tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0);

nrows = 3;
ncols = 1;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 1; colorE = [0 0.8 0.8];
doV = 1; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];
doTraj = 0;

cmap_dist = pic_colors('waterfall');

%freds = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
freds = {fred3_z0,fred3_z0,fred3_z0};
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.fv' labstr]);
    pcolor(hca,fred.x,fred.v,log10(fredplot)')
    shading(hca,'flat')
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy4')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 12;
    if doTraj
      hold(hca,'on')
      xx = [tr.x0];
      vv = eval(['[tr.v' labstr '0]']);
      scatter(hca,xx,vv,20,0*[1 1 1])
      hold(hca,'off')
    end    
  end
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
h(1).CLim = 0.99*[-4 2];
h(1).YLim = 0.99*4*[-1 1];

%% Location of integrated particles at different times
times = [90 100 110 120 125];
%times = [90 100];
xlim = [65 100];
zlim = [-6 9];
fontsize = 14;

colors_matlab = pic_colors('matlab');
colors_tr = [colors_matlab(5,:); 1 1 1;colors_matlab(3,:); 0 0 0];

tr = tr100(783:917);
tr1 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]>85);
tr = tr100(783:917);
tr2 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]<0);
tr = tr100(783:917);
tr3 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]<=75);
tr = tr100(783:917);
tr4 = tr.find([tr.Ustart]<0.25,[tr.zstart]>0,[tr.vy0]>0,[tr.x0]>75,[tr.x0]<=85);

tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0);

trajcolordot = nan(tr100.ntr,3);
trajcolordot([tr1.id],:) = repmat(colors(1,:),tr1.ntr,1);
trajcolordot([tr2.id],:) = repmat(colors(2,:),tr2.ntr,1);
trajcolordot([tr3.id],:) = repmat(colors(3,:),tr3.ntr,1);
trajcolordot([tr4.id],:) = repmat(colors(4,:),tr4.ntr,1);
trajcolordot(isnan(trajcolordot)) = [];
trajcolordot = reshape(trajcolordot,numel(trajcolordot)/3,3);


climEz = [-1 1];

nrows = 6;
ncols = 2;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;


tr = tr100(783:917);
tr = tr.find([tr.Ustart]<0.25,[tr.zstart]>0);
doTraj = 1;

if 1 % f(x,vy,t=120), with t0 location of particles  
  freds = {fred35_z0};  
  labstrs = {'y','z','x','y','z','x','y','z'};

  for ifred = 1:numel(freds)  
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.fv' labstr]);
    pcolor(hca,fred.x,fred.v,log10(fredplot)')
    shading(hca,'flat')
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy4')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    hcb = colorbar('peer',hca,'location','northoutside');  
    hcb.YLabel.String = sprintf('log_{10}f_{ic}(x,v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = fontsize;
    if doTraj
      hold(hca,'on')
      for itr_ = 1:4
        xx = eval(sprintf('[tr%g.x0]',itr_));
        vv = eval(sprintf('[tr%g.v%s0]',itr_,labstr));
        hs = scatter(hca,xx,vv,20,[0 0 0],'o');
        hs.MarkerFaceColor = colors_tr(itr_,:);
        %hs. = scatter(hca,xx,vv,20,'MarkerFaceColor',trajcolordot,'Marker','o','MarkerEdgeColor',[0 0 0]);
      end
      hold(hca,'off')
    end
  end
  hca.YLim = 0.99*[-3 3];
  hca.CLim = 0.99*[-4 1.7];  
end

for itime = 1:numel(times) % Ez
  hca = h(isub); isub = isub + 1;
  time = times(itime);
  pic = no02m.twcilim(time).xlim(xlim).zlim(zlim);
  imagesc(hca,pic.xi,pic.zi,pic.Ez')
  colormap(hca,pic_colors('blue_red'));  
  hca.YDir = 'normal';
  if itime == 1
    hcb = colorbar('peer',hca,'location','northoutside');  
    hcb.YLabel.String = 'E_z';
  end
  if 1 % A
    A = pic.A;
    hold(hca,'on')
    clim = hca.CLim;
    hc = contour(hca,pic.xi,pic.zi,A',[0:1:25],'k');
    %hc.Layer = 'top';
    hca.CLim = clim;  
    hold(hca,'off')
  end
  if doTraj
    hold(hca,'on')
    for itr_ = 1:4
      tr_tmp = eval(sprintf('tr%g',itr_));
      tr_tmp = tr_tmp.resample(time);
      xx = [tr_tmp.x];
      zz = [tr_tmp.z];
      %vv = eval(sprintf('[tr_tmp.v%s0]',labstr));
      hs = scatter(hca,xx,zz,20,[0 0 0],'o');
      hs.MarkerFaceColor = colors_tr(itr_,:);
      %hs. = scatter(hca,xx,vv,20,'MarkerFaceColor',trajcolordot,'Marker','o','MarkerEdgeColor',[0 0 0]);
    end
    hold(hca,'off')
  end
  hca.CLim = climEz;[0 0.5];
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  irf_legend(hca,{sprintf('t = %g',time)},[0.02 0.98],'color','k','fontsize',fontsize)
end

if 1 % f(x,vy,t=120), with t0 location of particles  
  freds = {fred35_z0};  
  labstrs = {'z','x','y','z','x','y','z'};

  for ifred = 1:numel(freds)  
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.fv' labstr]);
    pcolor(hca,fred.x,fred.v,log10(fredplot)')
    shading(hca,'flat')
    hca.XLabel.String = 'x (d_i)';
    hca.YLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy4')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    hcb = colorbar('peer',hca,'location','northoutside');  
    hcb.YLabel.String = sprintf('log_{10}f_{ic}(x,v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = fontsize;
    if doTraj
      hold(hca,'on')
      for itr_ = 1:4
        xx = eval(sprintf('[tr%g.x0]',itr_));
        vv = eval(sprintf('[tr%g.v%s0]',itr_,labstr));
        hs = scatter(hca,xx,vv,20,[0 0 0],'o');
        hs.MarkerFaceColor = colors_tr(itr_,:);
        %hs. = scatter(hca,xx,vv,20,'MarkerFaceColor',trajcolordot,'Marker','o','MarkerEdgeColor',[0 0 0]);
      end
      hold(hca,'off')
    end
  end
  hca.YLim = 0.99*[-3 3];
  hca.CLim = 0.99*[-4 1.7];  
end

for itime = 1:numel(times)
  hca = h(isub); isub = isub + 1;
  time = times(itime);
  pic = no02m.twcilim(time).xlim(xlim).zlim(zlim);
  imagesc(hca,pic.xi,pic.zi,pic.n(3)')  
  colormap(hca,pic_colors('thermal'));
  hca.YDir = 'normal';
  if itime == 1
    hcb = colorbar('peer',hca,'location','northoutside');  
    hcb.YLabel.String = 'n_{cold}^{top}';
  end
  if 1 % A
    A = pic.A;
    hold(hca,'on')
    clim = hca.CLim;
    hc = contour(hca,pic.xi(1:2:end),pic.zi(1:2:end),A(1:2:end,1:2:end)',[0:1:25],'k');
    %hc.Layer = 'top';
    hca.CLim = clim;  
    hold(hca,'off')
  end
  if doTraj
    hold(hca,'on')
    for itr_ = 1:4
      tr_tmp = eval(sprintf('tr%g',itr_));
      tr_tmp = tr_tmp.resample(time);
      xx = [tr_tmp.x];
      zz = [tr_tmp.z];
      %vv = eval(sprintf('[tr_tmp.v%s0]',labstr));
      hs = scatter(hca,xx,zz,20,[0 0 0],'o');
      hs.MarkerFaceColor = colors_tr(itr_,:);
      %hs. = scatter(hca,xx,vv,20,'MarkerFaceColor',trajcolordot,'Marker','o','MarkerEdgeColor',[0 0 0]);
    end
    hold(hca,'off')
  end
  hca.CLim = [0 0.5];
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  irf_legend(hca,{sprintf('t = %g',time)},[0.02 0.98],'color',[1 1 1],'fontsize',fontsize)
end

hlinks = linkprop(h,{'XLim'});
h(1).XLim = [65 100];
%compact_panels(h,0.01)

h(nrows+1).Position([2 3 4]) = h(1).Position([2 3 4]);

for ip = 1:npanels
  h(ip).FontSize = fontsize;
%   h(ip).Position(3) = 0.4;
%   h(ip).Position(4) = 0.2;
end
for ip = (nrows+1):npanels
  h(ip).YLabel.String = '';  
  h(ip).YTickLabel = [];  
end
compact_panels(h([2:nrows (nrows+2):npanels]),0.01,0.01)
compact_panels(h([1 nrows+1]),0.01,0.01)
drawnow

c_eval('h(?).Position(4) = h(2).Position(4);',[1 nrows+1])
c_eval('h(?).Position(2) = h(?).Position(2)-0.08;',1:npanels)
c_eval('h(?).Position(2) = h(?).Position(2)+0.09;',[1 nrows+1])
%
% hpos_y_fields = h(2).Position(2);
% c_eval('h(?).Position(2) = hpos_y_fields-0.1;',[2 nrows+2])
% compact_panels(h([2:nrows (nrows+2):npanels]),0.01,0.01)
% compact_panels(h([1 nrows+1]),0.01,0.01)
% c_eval('h(?).Position(2) = h(?).Position(2)-0.05;',[2:nrows (nrows+2):npanels])

%% Reduced distributions, prepare horizontal
twpe = 24000; xlim = [50 155]; zlim = [-15 15];

for zpick = [0 2 4]
  ds = ds100.twpelim(twpe).zfind(zpick).xlim(xlim).findtag({'line horizontal'});
  for id = 1:ds.nd{1}
    xlim_ = [ds.xi1{1}(id) ds.xi2{1}(id)];
    zlim_ = [ds.zi1{1}(id) ds.zi2{1}(id)];    
    pic_lim = no02m.xlim(xlim_).zlim(zlim_).twpelim(twpe);
    %pic = no02m.twpelim(twpe);
    %Bx_ = pic.Bx;
    %By_ = pic.By;
    %Bz_ = pic.Bz;
    eval(sprintf('x_z%g(id) = mean(xlim_);',zpick))
    eval(sprintf('z_z%g(id) = mean(zlim_);',zpick))
%     eval(sprintf('vExBx_z%g(id) = mean(mean(pic_lim.vExBx));',zpick))
%     eval(sprintf('vExBy_z%g(id) = mean(mean(pic_lim.vExBy));',zpick))
%     eval(sprintf('vExBz_z%g(id) = mean(mean(pic_lim.vExBz));',zpick))
    eval(sprintf('Ex_z%g(id) = mean(mean(pic_lim.Ex));',zpick))
    eval(sprintf('Ey_z%g(id) = mean(mean(pic_lim.Ey));',zpick))
    eval(sprintf('Ez_z%g(id) = mean(mean(pic_lim.Ez));',zpick))
    
  end
end

%% Reduced distributions, prepare vertical
twpe = 24000; xlim = [50 155]; zlim = [-15 15];

for xpick = [75 85 95]
  ds = ds100.twpelim(twpe).xfind(xpick).zlim(zlim).findtag({'line vertical'});
  for id = 1:ds.nd{1}
    xlim_ = [ds.xi1{1}(id) ds.xi2{1}(id)];
    zlim_ = [ds.zi1{1}(id) ds.zi2{1}(id)];    
    pic_lim = no02m.xlim(xlim_).zlim(zlim_).twpelim(twpe);
    %pic = no02m.twpelim(twpe);
    %Bx_ = pic.Bx;
    %By_ = pic.By;
    %Bz_ = pic.Bz;
    eval(sprintf('x_x%g(id) = mean(xlim_);',xpick))
    eval(sprintf('z_x%g(id) = mean(zlim_);',xpick))
    eval(sprintf('vExBx_x%g(id) = mean(mean(pic_lim.vExBx));',xpick))
    eval(sprintf('vExBy_x%g(id) = mean(mean(pic_lim.vExBy));',xpick))
    eval(sprintf('vExBz_x%g(id) = mean(mean(pic_lim.vExBz));',xpick))
    eval(sprintf('Ex_x%g(id) = mean(mean(pic_lim.Ex));',xpick))
    eval(sprintf('Ey_x%g(id) = mean(mean(pic_lim.Ey));',xpick))
    eval(sprintf('Ez_x%g(id) = mean(mean(pic_lim.Ez));',xpick))
    eval(sprintf('phiz_x%g = -cumtrapz(z_x%g,Ez_x%g);',xpick,xpick,xpick))
  end
end

%% Reduced distribution plot, vertical
% What to include
% - overview of where boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];

nrows = 1;
ncols = 3;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1; 
doE = 0; colorE = [0 0.8 0.8];
doV = 0; colorV = 0*[1 1 1];
doVe = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];
doSep = 1;
% if doSep
%   xline_pos = no02m.twpelim(twpe).xline_position;
%   sep = no02m.twpelim(twpe).separatrix_location;
% end
hleg = gobjects(0);

cmap_dist = pic_colors('waterfall');

%freds = {fred35_z4,fred35_z4,fred35_z4,fred35_z2,fred35_z2,fred35_z2,fred35_z0,fred35_z0,fred35_z0};
%freds = {fred35_x75,fred35_x80,fred35_x85,fred35_x90,fred35_x95};
%freds = {fred35_x75,fred35_x85,fred35_x95,fred35_x75,fred35_x85,fred35_x95,fred35_x75,fred35_x85,fred35_x95};
freds = {fred3_x75,fred3_x85,fred3_x95};
%freds = {fred3_x85};
%freds = {fred5_x75,fred5_x85,fred5_x95};
%freds = {fred46_x75,fred46_x85,fred46_x95};
%freds = {fred46_x75,fred46_x85,fred46_x95,fred46_x75,fred46_x85,fred46_x95,fred46_x75,fred46_x85,fred46_x95};
%freds = {fred46_x75,fred46_x80,fred46_x85,fred46_x90,fred46_x95};
%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','x','x','x','x','x','x','x','x'};
labstrs = {'y','y','y','y','y','y','y','y','y'};
labstrs = {'z','z','z','z','z','z','z','z','z'};
labstrs = {'z','z','z','y','y','y','x','x','x'};
zrefs = [7.3 4 3.5];
legs = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)'};

for ifred = 1:numel(freds)
  if 1 % fi(v_) z_
    hca = h(isub); isub = isub + 1;
    fred = freds{ifred};
    labstr = labstrs{ifred};
    fredplot = eval(['fred.fv' labstr]);
    pcolor(hca,fred.v,fred.z,log10(fredplot))
    shading(hca,'flat')
    hca.YLabel.String = 'z (d_i)';
    hca.XLabel.String = sprintf('v_{%s}',labstr);
    colormap(hca,pic_colors('candy4')) 
    %irf_legend(hca,{sprintf('f_{cold}(x,v_%s)',labstr)},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{legs{isub-1}},[0.02 0.98],'color',[0 0 0],'fontsize',14)
    irf_legend(hca,{sprintf('x = %g',unique(fred.x))},[0.02 0.9],'color',[0 0 0],'fontsize',14)
    %hcb = colorbar('peer',hca);  
    %hcb.YLabel.String = sprintf('f_{i,cold}(l_{||},v_{%s})',labstr);
    %hca.CLim(2) = prctile(fred.fvx(:),99);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Layer = 'top';
    hca.FontSize = 12;
    if 1*doE
      hold(hca,'on')
      %plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['E' labstr '_x' num2str(unique(fred.x))]);
      hE = plot(hca,smooth(vv,50),zz,'color',0*colorE,'linewidth',1,'linestyle',':');      
      hold(hca,'off')
    end
    if doSep
      hold(hca,'on')      
      xx = unique(fred.x);
      [~,ix] = min(abs(sep.x-xx));
      zz = sep.z(ix);
      hSep = plot(hca,hca.XLim,zz*[1 1],'color',0*colorExB,'linewidth',1,'linestyle',':');      
      hold(hca,'off')      
    end
    if doV
      hold(hca,'on')
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['vi' labstr '_x' num2str(unique(fred.x))]);      
      hv = plot(hca,vv,zz,'color',colorV,'linewidth',1.5);
      hold(hca,'off')
      hleg()
    end
    if doVe
      hold(hca,'on')
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['ve' labstr '_x' num2str(unique(fred.x))]);      
      hv = plot(hca,vv,zz,'color',colorV,'linewidth',1.5);
      hold(hca,'off')
      hleg()
    end
    if doExB
      hold(hca,'on')
      unique(fred.x);
      zz = eval(['z_x' num2str(unique(fred.x))]);
      vv = eval(['vExB' labstr '_x' num2str(unique(fred.x))]);
      hExB = plot(hca,smooth(vv,50),zz,'color',0*colorExB,'linewidth',1,'linestyle','-');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      hold(hca,'off')
    end  
    if doPhi
      hold(hca,'on')
      unique(fred.x);
      %if 1
      try
      zz = eval(['z_x' num2str(unique(fred.x))]);
      pp = eval(['phi' labstr '_x' num2str(unique(fred.x))]);
      vv0 = eval(['vExB' labstr '_x' num2str(unique(fred.x))]);
      vv0 = vv0(end);
      [~,iref]=min(abs(zz-zrefs(ifred)));
      %iref = fix(numel(pp)/2);
      vv_phi = sqrt(2*(abs(pp-pp(iref)))).*sign(pp-pp(iref));
      %(smooth(vv,10)+vv0-vv(end));
      hphi = plot(hca,vv_phi,zz,'color',0*[1 1 1],'linewidth',1,'linestyle','--');
      %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
      end
      hold(hca,'off')
    end   
  end
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim','YLim','CLim'});
position_hcb_peer = h(end).Position;
hcb = colorbar('peer',h(end));
h(end).Position = position_hcb_peer;
hcb.YLabel.String = sprintf('f_{i,cold}(z,v_{%s})',labstr);
hcb.YLabel.FontSize = 14;
% if all([doV doExB])
%   legend([hv,hExB],{sprintf('v_{i,%s}',labstr),sprintf('v_{ExB,%s}',labstr)})
% elseif doV
%   legend([hv],{sprintf('v_{i,%s}',labstr)})
% elseif doExB
%   legend([hExB],{sprintf('v_{ExB,%s}',labstr)})
% end
hh = findobj(h(3),'type','line');
legend(hh([2 1 3]),{'v_{ExB}','v_{\phi}','separarix'},'location','best','edgecolor',[1 1 1],'box','on')
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
compact_panels(0.005,0.01)
c_eval('h(?).Position(1) = h(?).Position(1)-0.05;',1:npanels)
h(1).CLim = 0.99*[-4 2];
h(1).CLim = [-3.5 1.8];
h(1).CLim = [-4 2];
h(1).XLim = 0.99*[-2.5 2.5];
%h(1).CLim = 0.99*[-6 0]; % electrons
%h(1).XLim = 0.99*3*[-1 1];
for ip = 2:npanels
  h(ip).YTickLabels = '';
  h(ip).YLabel.String = '';  
end
for ip = 1:npanels
  h(ip).FontSize = 14;
  h(ip).Position(2) = 0.17;
end

%% Reduced distributions f(x,vy) with pyr