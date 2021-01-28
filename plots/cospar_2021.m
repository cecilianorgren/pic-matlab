%% Movie of vex, Epar, log10(ne)

twpe = 15000:100:25000;
xlim = no02m.xi([1 end])+[60 -60]';
xlim = [60 110];
zlim = [-1 8];
pic = no02m.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');

cmapbr = pic_colors('blue_red');
cmapwa = pic_colors('waterfall');
cmapth = pic_colors('thermal');
%cmapjet = colormap('jet');


varstrs = {'ni','log10(ni)'}';
varstrs = {'vex','Epar','log10(ne)'}';
clims = {[-6 6],[-0.5 0.5],[-2 0.5]};

% varstrs = {'n([3 5])','t([3 5])'}';
% clims = {[0 0.5],[0 0.5]};
cmaps = {cmapbr,cmapbr,cmapth};
cbarlabels = {'v_{ex}','E_{||}','log_{10}(n_e)'};
filename = [printpath 'cospar_intro_video_test5'];
pic.movie(varstrs,'A',1,'cmap',cmaps,'clim',clims,'cbarlabels',cbarlabels,'filename',filename,'smooth',3);

%% Movie density hot cold line in z
comp = 'z';
twpe = [1000:1000:24000];
xlim = diff(no02m.xi([1 end]))/2+[-1 1];
zlim = 10*[-1 1];

pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
%pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);
varstrs = {{'ni','ne'};{'Ex'};{'Ez'};{'txx([4 6])','tyy([4 6])','tzz([4 6])'};{'txx([3 5])','tyy([3 5])','tzz([3 5])'}};
varstrs = {{'ni','ne','n([1])','n([3 5])','n([4 6])','n([4 6])'};{'Ey'};{'Ez'};{'txx([2 4 6])','tyy([2 4 6])','tzz([2 4 6])'};{'txx([1 3 5])','tyy([1 3 5])','tzz([1 3 5])'}};
varstrs = {{'Bz','By'};{'vix','vex','vExBx'};{'Jx','Jz'};{'Ex','Ez'}};
varstrs = {{'Bz','By'};{'vix','vex','vExBx'};{'Jx','Jz'};{'viy','vey'}};
varstrs = {{'Bx','By','Bz'};{'vix','viy','viz'};{'vex','vey','vez'}};
varstrs = {{'n(1)','n([3 5])'}};

pic.movie_line(comp,varstrs,'ylim',{[0 1]},'filename',[printpath 'n_cold_hot']);

%% Initial conditions
pic = no02m;
comp = 'z';
twpe = 12000;
pic = no02m.twpelim(twpe);
xlim = pic.x_xline+0.5*[-1 1]; % % at xline
%xlim = 53 + [-1 1]; % unperturbed part of current sheet
zlim = pic.zi([1 end]);
pic = pic.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');

varstrs = {{'n([1])','n([3 5])','n([1 3 5])'}};

colors = pic_colors('matlab');

h = pic.plot_line(comp,varstrs,'smooth',20);

h.YLim = [0 1];
h.Children(1).Color = colors(3,:);
h.Children(2).Color = colors(1,:);
h.Children(3).Color = colors(2,:);
h.XLim = 10*[-1 1];
legend(h,{'Hot','Cold','Total'},'location','best')
h.Position = [0.10 0.2 0.8 0.7];
h.YLabel.String = 'n (n_0)';
h.FontSize = 14;
h.Title.String = 'Initial density configuration';

%% Transition from hot dense to cold tenuous inflow
pic = no02m;
if 0 % Calculate density at and above xline
  pic_Bxline_z05 = pic.interp(pic.x_xline,pic.z_xline+0.5,pic.twci,'Bx');
  pic_nxline_z05 = pic.interp(pic.x_xline,pic.z_xline+0.5,pic.twci,'ne');
  pic_Bxline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'Bx');
  pic_nxline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'ne');
  pic_nxline_z0 = pic.interp(pic.x_xline,pic.z_xline+0,pic.twci,'ne');
  pic_tixline_z1 = pic.interp(pic.x_xline,pic.z_xline+1,pic.twci,'ti');
  pic_vA_z1 = squeeze(pic_Bxline_z1./sqrt(pic_nxline_z1));
  pic_vA_z05 = squeeze(pic_Bxline_z05./sqrt(pic_nxline_z05));
end
% Figure
zlim = [-0.5 0.5];
xlim = [100 300];
pic = pic.zlim(zlim).xlim(xlim);

%pic1 = df04.zlim(zlim).xlim(xlim);
%pic2 = df04n.zlim(zlim).xlim(xlim);
%A1 = squeeze(mean(pic1.A,2));
%A2 = squeeze(mean(pic2.A,2));
%R1 = reconnection_rate(timesteps/200,'A',A_ts,'E',E_ts(:,:,:,2));
Alev = -25:1:25;
doA = 0;
istepA = 3;

nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
hb = gobjects(0);


if 1 % ni, Bx, vA
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic_nxline_z1,pic.twci,pic_Bxline_z1);
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'n, B_x';
  %hca.XLim = [0 0.4];
  legend(hca,{'n(x_x,z_x+1)','B_x(x_x,z_x+1)'},'location','southwest','box','off')
end
if 1 % vA
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic_vA_z1);
  hca.YLim(1) = 0;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'v_A';
  legend(hca,{'v_A(x_x,z_x+1)'},'location','northwest','box','off')
end
if 1 % ER
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic.RE);
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'E_R';
  %hca.XLim = [0 0.4];  
  legend(hca,{'E_y(x_x,z_x)'},'location','northwest','box','off')
end
compact_panels(h,0.01)
for ip = 1:numel(h)
  h(ip).Position(2) = h(ip).Position(2) + 0.05;
%   h(ip).Position(4) = 0.5;
%   h(ip).Position(3) = 0.15;
  %h(ip).XTick = 0:50:500;
  h(ip).XLim = [0 pic.twci(end)];
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).FontSize = 14;
  %hb(ip).FontSize = 14;
end

%% Transition from hot dense to cold tenuous inflow

pic = no02m;
% Figure
zlim = [-0.1 0.1];
xlim = diff(pic.xi([1 end]))/2+[-45 45];
pic = pic.zlim(zlim).xlim(xlim);

%pic1 = df04.zlim(zlim).xlim(xlim);
%pic2 = df04n.zlim(zlim).xlim(xlim);
%A1 = squeeze(mean(pic1.A,2));
%A2 = squeeze(mean(pic2.A,2));
%R1 = reconnection_rate(timesteps/200,'A',A_ts,'E',E_ts(:,:,:,2));
Alev = -25:0.5:25;
doA = 1;
istepA = 3;
if doA
  plA = squeeze(mean(pic.A,2));
end

doXline = 1;

nrows = 3;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
hb = gobjects(0);


if 1 % ER
  hca = h(isub); isub = isub + 1;  
  h_ = plot(hca,pic.twci,pic.RE);
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'E_R';
  %hca.XLim = [0 0.4];  
  legend(hca,{'E_y(x_x,z_x)'},'location','northwest','box','off')
end
if 0 % Ey1
  hca = h(isub); isub = isub + 1;
  varplot = squeeze(mean(pic.Ey,2));
  pcolor(hca,pic.twci,pic.xi,varplot)
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'E_y';
  hca.YLabel.String = 't\omega_{ci}';
  hca.XLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')
    plA = pic.A;
    %contour(hca,pic.xi(iAx),pic.zi(iAz),A(iAx,iAz)',Alev,'k');        
    contour(hca,pic.xi(1:istepA:end),pic.twci,plA(1:istepA:end,:)',Alev,'color',[0 0 0])
    hold(hca,'off')
  end
end
if 1 % Bz1
  hca = h(isub); isub = isub + 1;
  varplot = squeeze(mean(pic.Bz,2));
  pcolor(hca,pic.twci,pic.xi,varplot)
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  %hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'B_z';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  if doA
    hold(hca,'on')        
    [~,hA] = contour(hca,pic.twci,pic.xi(1:istepA:end),plA(1:istepA:end,:),Alev,'color',[0 0 0]);
    hold(hca,'off')
    legend(hA,{'A_y'},'box','off','location','northwest');
  end
  if doXline
    hold(hca,'on')        
    plot(hca,pic.twci,pic.x_xline,'linewidth',1,'color',[0 0 0])
    hold(hca,'off')    
  end
  hca.CLim = [-0.7 0.7];
  colormap(hca,pic_colors('blue_red'))
end
if 1 % n ratio
  hca = h(isub); isub = isub + 1;
  ncold = squeeze(mean(pic.n([3 5]),2));
  nall = squeeze(mean(pic.ni,2));
  varplot = ncold./nall;
  pcolor(hca,pic.twci,pic.xi,varplot)
  shading(hca,'flat')
  colormap(hca,pic_colors('blue_red'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hcb = colorbar('peer',hca);
  %hcb.Location = 'northoutside';
  hb(isub-1) = hcb;
  hcb.YLabel.String = 'n_{cold}/n_{tot}';
  hca.XLabel.String = 't\omega_{ci}';
  hca.YLabel.String = 'x/d_i';
  if doA   
    hold(hca,'on') 
    [~,hA] = contour(hca,pic.twci,pic.xi(1:istepA:end),plA(1:istepA:end,:),Alev,'color',[0 0 0]);
    hold(hca,'off')
    %legend(hA,{'A_y'},'box','off','location','northwest');
  end  
  if doXline
    hold(hca,'on')        
    plot(hca,pic.twci,pic.x_xline,'linewidth',1,'color',[0 0 0])
    hold(hca,'off')    
  end
  hca.CLim = [0 1];
  %colormap(hca,pic_colors('candy'))
  colormap(hca,flipdim(pic_colors('blue_red'),1))
end

compact_panels(0.01);
for ip = 1:numel(h)
  h(ip).Position(2) = h(ip).Position(2) + 0.05;
%   h(ip).Position(4) = 0.5;
  h(ip).Position(3) = 0.7;
  %h(ip).XTick = 0:50:500;
  h(ip).XLim = [10 pic.twci(end)];
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).FontSize = 14;
  %hb(ip).FontSize = 14;
end

%% Illustrating thermalization vs. superposition
colors = pic_colors('matlab');
vt = 1;
n = 1;
v = linspace(-20,30,1000);

f = @(v,n,vd,vt) n./sqrt(vt).*exp(-(v-vd).^2./vt.^2);


if 1
  hca = subplot(2,1,1);
  %set(hca,'ColorOrder',colors); hold(hca,'on')
  plot(hca,v,f(v,1,0,1),'-','color',[0 0 0])
  hold(hca,'on')  
  plot(hca,v,f(v,1,0,5),'--','color',[0 0 0])
  hold(hca,'off')
  hca.Visible = 'off';
  hca.YLim = [0 1];
end
if 1
  hca = subplot(2,1,2);
  vt = 5;
  vdstep = 5;
  plot(hca,v,f(v,1,0*vdstep,vt),'--',v,f(v,1,1*vdstep,vt*1.2),'--',v,f(v,1,2*vdstep,vt*1.4),'--',v,f(v,1,3*vdstep,vt*1.6),'--','color',[0 0 0])
  hold(hca,'on')
  plot(hca,v,f(v,1,0*vdstep,vt)+f(v,1,1*vdstep,vt*1.2)+f(v,1,2*vdstep,vt*1.4)+f(v,1,3*vdstep,vt*1.6),':','color',[0 0 0],'linewidth',1)
  hold(hca,'off')
  %hca.Box = 'off';
  %hca.XTick = [];
  %hca.YTick = [];
  hca.Visible = 'off';
  hca.YLim = [0 1];
end


