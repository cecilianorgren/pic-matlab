% py = m*vy + q*Ay
% py0 = m*vy0 + q*Ay0
% vy = py0/m - q*Ay/m = py0/m - Bref*x^2/2

f_vy_A = @(py0,A) py0 - A;
f_vy_B = @(py0,Bref,x) py0 - Bref*x^2/2;

f_A0_vy_py0 = @(A,vy) vy + A;

%% Horizontal cut of Bz, A, and -Bref*(x-xref)^2/2
twpe = 24000;
pic = no02m.twpelim(twpe);
comp = 'x';
xlim = pic.xi([1 end])+[70 -70]';
zlim = 1*[-0.5 0.5];
pic = pic.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');

varstrs = {{'Bz'};{'A','-4-0.01.*(xmesh-98).^2/2'}};
h = pic.plot_line(comp,varstrs);

%% Reduced distributions, make

twpe = 24000; xlim = [50 155]; zlim = [-15 15];

for zpick = [0]
  ds = ds100.twpelim(twpe).zfind(zpick).xlim(xlim).findtag({'line horizontal'});
  xdist = (ds.xi1{1}+ds.xi2{1})/2;
  zdist = (ds.zi1{1}+ds.zi2{1})/2;
  pic_lim = no02m.xlim(xlim).zlim(zlim).twpelim(twpe);
  pic = no02m.twpelim(twpe);
  Bx_ = pic.Bx;
  By_ = pic.By;
  Bz_ = pic.Bz;
  Bx = interpfield(pic.xi,pic.zi,Bx_,xdist,zdist); 
  By = interpfield(pic.xi,pic.zi,By_,xdist,zdist); 
  Bz = interpfield(pic.xi,pic.zi,Bz_,xdist,zdist); 
  fred5_tmp = ds.reduce_1d_new('x',[5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred5_z%g = fred5_tmp;',zpick))
  fred3_tmp = ds.reduce_1d_new('x',[3],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred3_z%g = fred3_tmp;',zpick))
  fred35_tmp = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred35_z%g = fred35_tmp;',zpick))
  %fred46_tmp = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx,By,Bz},'pitch',{Bx,By,Bz}); eval(sprintf('fred46_z%g = fred46_tmp;',zpick))    
end

%% Reduced distributions, plot
twpe = 24000;
pic = no02m.twpelim(twpe);
fred = fred35_z0;
twpe = 24000;
fontsize = 12;
% What to include
% - overview of whxx  ere boxes are
fi_clim = [0 0.13];
fe_clim = [0 0.02];
fe_clim = [-4 log10(0.02)];
zlim_line = [-0.5 0.5];
xlim_line = [min(fred3_z0.x) max(fred3_z0.x)];

py0 = [8 7 6.5 5.8 5.2 4.6];
py0 = [8 6.9 6.3 5.6 5.1 4.6];
py0 = [8 6.9 6.3 5.6 5.1 4.6];
py0 = [8 7 6.5 5.8 5.4 4.6];
xmin = [70 70 71 73 80 82];
xmax = [80 85 87 93 95 100];


py0 = [8.5 8.1 7.7 7.4 6.9 6.3 5.8 5.2 4.6];
xmin = [71 71 71 71 71 71.5 73 80 83];
xmax = [74 77 79 81 84 87 92 96 97]+00;
npy = numel(py0);

nrows = 1;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
doE = 1; colorE = [0 0.8 0.8];
doV = 1; colorV = 0*[1 1 1];
doN = 1; colorN = [0 0 0];
doExB = 1; colorExB = 0*[1 1 1]+0.5;
doPhi = 1; colorPhi = [0.5 0.5 0];
doTraj = 1; colorTraj = [0 0 0];
trs = tr100.find([tr100.z0] == 0,[tr100.vy0] > 0.5,[tr100.x0] > 73,[tr100.xstart] < 100);
trs = tr100.find([tr100.z0] == 0,[tr100.vy0] > 0.5).lim('t',[23000 25000]/200);

cmap_dist = pic_colors('waterfall');

%freds = {fred3_z4,fred3_z4,fred3_z4,fred3_z2,fred3_z2,fred3_z2,fred3_z0,fred3_z0,fred3_z0};
labstrs = {'x','y','z','x','y','z','x','y','z'};

if 0 % A(x),
  hca = h(isub); isub = isub + 1; 
  pic_tmp = pic.xlim(xlim_line).zlim(zlim_line);
  plot(hca,pic_tmp.xi,mean(pic_tmp.A,2))  
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'A_{y}';    
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % fi(v_y)
  hca = h(isub); isub = isub + 1;
  fred = fred35_z0;    
  pcolor(hca,fred.x,fred.v,log10(fred.fvy)')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('candy4'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'f_{i,cold}(x,v_{y})';  
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  hca.CLim = 0.99*[-4 2];
  hca.CLim = 0.99*[-4 2];
  hca.YLim = 0.99*4*[-1 1];
  
  for ipy = 1:npy
    hold(hca,'on')
    pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
    A = squeeze(mean(pic_tmp.A,2));    
    xx = pic_tmp.xi;
    yy = f_vy_A(py0(ipy),A);
    plot(hca,xx,yy,':k','linewidth',0.5)
    hold(hca,'off')
    ht = text(hca,xx(end),yy(end),sprintf(' %.1f',py0(ipy)),'fontsize',fontsize);
    ht.HorizontalAlignment = 'left';
    dy = (yy(end)-yy(end-1))./diff(hca.YLim);
    dx = (xx(end)-xx(end-1))./diff(hca.XLim);
    ht.Rotation = atand(dy/dx);
    
  end
  if 0*doE
    hold(hca,'on')
    plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
    hold(hca,'off')
  end
  if 0*doV
    hold(hca,'on')
    plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if 0 % doExB
    hold(hca,'on')
    xx = eval(['x_z' num2str(unique(fred.z))]);
    vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
    plot(hca,xx,vv,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % fi(v_y)
  hca = h(isub); isub = isub + 1;
  fred = fred3_z0;    
  pcolor(hca,fred.x,fred.v,log10(fred.fvy)')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('candy4'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'f_{i,cold}(x,v_{y})';  
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  hca.CLim = 0.99*[-4 2];
  hca.CLim = 0.99*[-2 1];
  hca.YLim = 0.99*4*[-1 1];
  
%   for ipy = 1:npy
%     hold(hca,'on')
%     pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
%     A = squeeze(mean(pic_tmp.A,2));
%     plot(hca,pic_tmp.xi,f_vy_A(py0(ipy),A),'--k')
%     hold(hca,'off')
%   end

   for itr = 1:trs.ntr
      hold(hca,'on')
      plot(hca,trs(itr).x,trs(itr).vy,'k')
      hold(hca,'off')
   end
end
if 0 % fi_top(v_y)/fi_tot(v_y)
  hca = h(isub); isub = isub + 1;
  fred = fred3_z0;    
  pcolor(hca,fred.x,fred.v,fred3_z0.fvy'./fred35_z0.fvy')
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('waterfall'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'f_{i,top}/f_{i,tot}(x,v_{y})';  
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  hca.CLim = 1*[0 1];  
  hca.YLim = 0.99*4*[-1 1];
  
%   for ipy = 1:npy
%     hold(hca,'on')
%     pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
%     A = squeeze(mean(pic_tmp.A,2));
%     plot(hca,pic_tmp.xi,f_vy_A(py0(ipy),A),'--k')
%     hold(hca,'off')
%   end

   for itr = 1:trs.ntr
      hold(hca,'on')
      plot(hca,trs(itr).x,trs(itr).vy,'k')
      hold(hca,'off')
   end
end
if 0 % A0(v_y,A)
  hca = h(isub); isub = isub + 1; 
  A_ = no02m.interp(fred.x,fred.z,no02m.twpelim(twpe).twci,'A');
  [VY,A] = meshgrid(fred.v,A_);  
  A0map = f_A0_vy_py0(A,VY);
  A0map(fred.fvy==0) = NaN;
  %[Ccont,hcont] = contourf(hca,fred.x,fred.v,A0map',-30:1:30);
  [Ccont,hcont] = contourf(hca,fred.x,fred.v,A0map',[0 sort(py0)]);
  clabel(Ccont,hcont,'LabelSpacing',72,'Color','k','FontWeight','bold');
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{y}';
  colormap(hca,pic_colors('waterfall'))   
  irf_legend(hca,{sprintf('z = %g',unique(fred.z))},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  hcb = colorbar('peer',hca,'fontsize',14);  
  hcb.YLabel.String = 'A_{y0}(A_{y,loc},v_y,v_{y0}=0)';  
  irf_legend(hca,{'A_{y0} = A_{y,loc} - v_{y,loc}';'assuming v_{y0}=0'},[0.98 0.1],'color',[0 0 0])
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;
  %hca.CLim = 0.99*[-4 2];
  %hca.CLim = 0.99*[-2 1];
  hca.YLim = 0.99*4*[-1 1];
  
  for ipy = 1:npy
    hold(hca,'on')
    pic_tmp = pic.xlim([xmin(ipy) xmax(ipy)]).zlim(zlim_line);
    A = squeeze(mean(pic_tmp.A,2));
    xx = pic_tmp.xi;
    yy = f_vy_A(py0(ipy),A);
    plot(hca,xx,yy,'--k')
    hold(hca,'off')
    ht = text(hca,xx(end),yy(end),sprintf('p_y = %.1f',py0(ipy)),'fontsize',fontsize);
  end
  if 0*doE
    hold(hca,'on')
    plot(hca,arclength,Ey*max(abs(hca.YLim))/max(abs(Ey)),'color',colorE)
    hold(hca,'off')
  end
  if 0*doV
    hold(hca,'on')
    plot(hca,arclength,viy,'color',colorV,'linewidth',1.5)
    hold(hca,'off')
  end
  if 0 % doExB
    hold(hca,'on')
    xx = eval(['x_z' num2str(unique(fred.z))]);
    vv = eval(['vExB' labstr '_z' num2str(unique(fred.z))]);
    plot(hca,xx,vv,'color',colorExB,'linewidth',1.5)
    %plot(hca,pic_sm.zi,mean(pic_sm.vExBy,1),'color',colorExB,'linewidth',1.5)
    hold(hca,'off')
  end
end
if 0 % ExB_x
  hca = h(isub); isub = isub + 1; 
  pic_tmp = pic.xlim(xlim_line).zlim(zlim_line);
  plot(hca,pic_tmp.xi,mean(pic_tmp.vExBx,2),pic_tmp.xi,mean(pic_tmp.vix,2))
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{x}';    
  legend(hca,{'v_{ExB}','v_i'},'location','best')
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 0 % A(x), for different times  
  hca = h(isub); isub = isub + 1; 
  pic_tmp = no02m.twpelim(19000:1000:24000,'exact').xlim(xlim_line).zlim(zlim_line);
  plot(hca,pic_tmp.xi,squeeze(mean(pic_tmp.A,2)))  
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'A_{y}';    
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
drawnow
compact_panels(h(1:end),0.01,0.05)
%h(1).Title.String = ['nobg, t\omega_{pe} = ' num2str(twpe,'%05.0f')];
fig = gcf; h_ = findobj(fig.Children,'type','axes');
hlinks = linkprop(h(1:end),{'XLim'});
%hlinks.Targets(1).XLim = arclength([1 end]);
%irf_plot_axis_align
for ip = 1:nrows*ncols
  axwidth(ip) = h(ip).Position(3);
  h(ip).FontSize = 14;
end
for ip = 1:nrows*ncols
  h(ip).Position(3) = 0.7;%min(axwidth);
end
% for ip = 1:nrows*ncols
%   h(ip).Position(2) = h(ip).Position(2)-0.05;
% end

%c_eval('h(?).YTickLabel = []; h(?).YLabel = [];',[2 3 5 6 8 9])
%c_eval('h(?).YTick = -10:1:10;',1)

%% Calculate py along trajectory

%trs;
%trs = tr100; 
pic = no02m;
for itr = 1:trs.ntr
  tic
  trsA(itr).A = no02m.interp(trs(itr).x,trs(itr).z,trs(itr).t,'A');
  toc
end
%trsAall=trsA;
%trsA = trsAall([trs.id]);

%% Plot A along trajectory
nrows = 3; ncols = 3;
h = setup_subplots(nrows,ncols);
isub = 1;


if 1 % x,y
  hca = h(isub); isub = isub + 1; 
  
  for itr = 1:trs.ntr
    hold(hca,'on')
    plot(hca,trs(itr).x,trs(itr).y);
    hold(hca,'off')
  end
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'y (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % x,vy
  hca = h(isub); isub = isub + 1; 
  
  for itr = 1:trs.ntr
    hold(hca,'on')
    plot(hca,trs(itr).x,trs(itr).vy);
    hold(hca,'off')
  end
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_y (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % y,vy
  hca = h(isub); isub = isub + 1; 
  
  for itr = 1:trs.ntr
    hold(hca,'on')
    plot(hca,trs(itr).y,trs(itr).vy);
    hold(hca,'off')
  end
  hca.XLabel.String = 'y (d_i)';
  hca.YLabel.String = 'v_y (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % x,y, A as color
  hca = h(isub); isub = isub + 1; 
  
  for itr = 1:trs.ntr
    hold(hca,'on')
    % remember A has wrong sign, change if it is changed at some point
    scatter(hca,trs(itr).x,trs(itr).y,5,-squeeze(trsA(itr).A));
    hold(hca,'off')
  end
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'A';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'y (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % x,y, vy as color
  hca = h(isub); isub = isub + 1; 
  
  for itr = 1:trs.ntr
    hold(hca,'on')
    scatter(hca,trs(itr).x,trs(itr).y,5,trs(itr).vy);
    hold(hca,'off')
  end
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_y';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'y (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % x,y, py = vy+A as color
  hca = h(isub); isub = isub + 1; 
  
  for itr = 1:trs.ntr
    hold(hca,'on')
    py = trs(itr).vy+-squeeze(trsA(itr).A);
    scatter(hca,trs(itr).x,trs(itr).y,5,py);
    hold(hca,'off')
  end
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'py';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'y (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % x,y, py = vy+A, py/py(1) as color
  hca = h(isub); isub = isub + 1; 
  
  for itr = 1:trs.ntr
    hold(hca,'on')
    py = trs(itr).vy+-squeeze(trsA(itr).A);
    scatter(hca,trs(itr).x,trs(itr).y,5,py/py(1));
    hold(hca,'off')
  end
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'py/py(1)';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'y (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % x,y,vx as color
  hca = h(isub); isub = isub + 1; 
  
  for itr = 1:trs.ntr
    hold(hca,'on')
    
    scatter(hca,trs(itr).x,trs(itr).y,5,trs(itr).vx);
    hold(hca,'off')
  end
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_x';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'y (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % x,y, energy as color
  hca = h(isub); isub = isub + 1; 
  
  for itr = 1:trs.ntr
    hold(hca,'on')
    energy = (trs(itr).vx.^2 + trs(itr).vy.^2 + trs(itr).vz.^2)/2;
    scatter(hca,trs(itr).x,trs(itr).y,5,energy/energy(1));
    hold(hca,'off')
  end
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'energy';
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'y (d_i)';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end

%% Plot against time
nrows = 4; ncols = 3;
h = setup_subplots(nrows,ncols);
isub = 1;

itrs = 1:25;%:trs.ntr;

if 1 % Ukin(t)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    ek = (trs(itr).vx.^2 + trs(itr).vy.^2 + trs(itr).vz.^2)/2;
    plot(hca,trs(itr).t,ek);
    hold(hca,'off')
  end  
  hca.YLabel.String = 'U_{kin}';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % Ukin(t)/Ukin(0)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    ek = (trs(itr).vx.^2 + trs(itr).vy.^2 + trs(itr).vz.^2)/2;
    plot(hca,trs(itr).t,ek-ek(1));
    hold(hca,'off')
  end  
  hca.YLabel.String = 'U_{kin}/U_{kin}(1)';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % Upot(t)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    pot = -cumtrapz(trs(itr).x,trs(itr).Ex) + ...
          -cumtrapz(trs(itr).y,trs(itr).Ey) + ...
          -cumtrapz(trs(itr).z,trs(itr).Ez);
    plot(hca,trs(itr).t,pot);
    hold(hca,'off')
  end  
  hca.YLabel.String = 'U_{pot}';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % Upot(t) + Ukin
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    pot = -cumtrapz(trs(itr).x,trs(itr).Ex) + ...
          -cumtrapz(trs(itr).y,trs(itr).Ey) + ...
          -cumtrapz(trs(itr).z,trs(itr).Ez);
    ek = sqrt(trs(itr).vx.^2 + trs(itr).vy.^2 + trs(itr).vz.^2);
    plot(hca,trs(itr).t,pot+ek);
    hold(hca,'off')
  end  
  hca.YLabel.String = 'U_{kin} + U_{pot}';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % Upot/Upot(1) + Ukin/Ukin(1)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    pot = -cumtrapz(trs(itr).x,trs(itr).Ex) + ...
          -cumtrapz(trs(itr).y,trs(itr).Ey) + ...
          -cumtrapz(trs(itr).z,trs(itr).Ez);
    ek = sqrt(trs(itr).vx.^2 + trs(itr).vy.^2 + trs(itr).vz.^2);
    plot(hca,trs(itr).t,pot+ek/ek(1));
    hold(hca,'off')
  end  
  hca.YLabel.String = 'U_{kin}/U_{kin}(1) + U_{pot}/U_{pot}(1)';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % A(t) as color
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    % remember A has wrong sign, change if it is changed at some point
    plot(hca,trs(itr).t,squeeze(-trsA(itr).A));
    hold(hca,'off')
  end
  
  hca.YLabel.String = 'A';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 0 % A(t)/A(t0) as color
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    % remember A has wrong sign, change if it is changed at some point
    plot(hca,trs(itr).t,squeeze(-trsA(itr).A)/squeeze(-trsA(itr).A(1)));
    hold(hca,'off')
  end
  
  hca.YLabel.String = 'A/A_0';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % vy(t)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    plot(hca,trs(itr).t,trs(itr).vy);
    hold(hca,'off')
  end  
  hca.YLabel.String = 'v_y';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 0 % vy(t)/vy(t0)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    plot(hca,trs(itr).t,trs(itr).vy/trs(itr).vy(1));
    hold(hca,'off')
  end  
  hca.YLabel.String = 'v_y/v_{y0}';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % py(t) = vy+A
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    py = trs(itr).vy+-squeeze(trsA(itr).A);
    plot(hca,trs(itr).t,py,trs(itr).t(1),py(1),'k.');
    hold(hca,'off')
  end
  hca.YLabel.String = 'p_y =v_y+A';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 0 % py(t)/py(t_0)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    py = trs(itr).vy+-squeeze(trsA(itr).A);
    plot(hca,trs(itr).t,py/py(1));
    hold(hca,'off')
  end  
  hca.YLabel.String = 'p_y/p_y(1)';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % Ey(t)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    plot(hca,trs(itr).t,trs(itr).Ey);
    hold(hca,'off')
  end  
  hca.YLabel.String = 'E_y';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 0 % dpy/dt
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    py = trs(itr).vy+-squeeze(trsA(itr).A);
    plot(hca,trs(itr).t,gradient(py,trs(itr).t));
    hold(hca,'off')
  end
  hca.YLabel.String = 'dp_y/dt';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % dpy/dt smoothed
  hca = h(isub); isub = isub + 1; 
  np_smooth = 100;
  for itr = itrs
    hold(hca,'on')
    py = trs(itr).vy+-squeeze(trsA(itr).A);
    plot(hca,trs(itr).t,smooth(gradient(py,trs(itr).t),np_smooth));
    hold(hca,'off')
  end
  hca.YLabel.String = {'dp_y/dt',sprintf('(smoothed, np = %g)',np_smooth)};
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % vx(t)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    plot(hca,trs(itr).t,trs(itr).vx);
    hold(hca,'off')
  end  
  hca.YLabel.String = 'v_x';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 1 % y(t)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    plot(hca,trs(itr).t,trs(itr).y);
    hold(hca,'off')
  end  
  hca.YLabel.String = 'y';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 0 % z(t)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    plot(hca,trs(itr).t,trs(itr).z);
    hold(hca,'off')
  end  
  hca.YLabel.String = 'z';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
if 0 % int Ey(t)
  hca = h(isub); isub = isub + 1; 
  
  for itr = itrs
    hold(hca,'on')
    phi = cumtrapz(trs(itr).y,trs(itr).Ey);
    plot(hca,trs(itr).t,phi);
    hold(hca,'off')
  end  
  hca.YLabel.String = '\int E_ydy';
  hca.XLabel.String = 't\omega_{ci}';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end
hlinks = linkprop(h,{'XLim'});
compact_panels(0.01)

%% Scatter plots
nrows = 1; ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % A, vy
  hca = h(isub); isub = isub + 1; 
  
  for itr = 1:trs.ntr
    hold(hca,'on')
    % remember A has wrong sign, change if it is changed at some point
    plot(hca,trs(itr).vy,squeeze(-trsA(itr).A));
    %plot(hca,trs(itr).vy/trs(itr).vy(1),squeeze(trsA(itr).A)/trsA(itr).A(1));
    hold(hca,'off')
  end
  axis(hca,'equal')
  hca.YLabel.String = 'A';
  hca.XLabel.String = 'v_y';  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.FontSize = 12;  
end

%% Horizontal cut of Bz, Ey
twpe = 24000;
pic = no02m.twpelim(twpe);
comp = 'x';
xlim = pic.xi([1 end]);%+[100 -100]';
zlim = 0+0.2*[-1 1];
pic = pic.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');

varstrs = {{'Bz','sqrt(Ey)','-sqrt(Ey)'},{'abs(Bz./Ey)'}}';
h = pic.plot_line(comp,varstrs,'smooth',1);
h(2).YLim = [0 2]; 
h(1).XLim = [60 140];

%% Original position of particle, based on py
py0 = [8.5 8.1 7.7 7.4 6.9 6.3 5.8 5.2 4.6];
np = numel(py0);
pic = no02m(1).xlim([10 11]).zlim([0 15]);
A = mean(pic.A,1);
for ip = 1:np
  [a,ind] = min(abs(A-py0(ip)));
  z0(ip) = pic.zi(ind);
end