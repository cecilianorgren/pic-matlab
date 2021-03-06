% plot_timeseries_adaptive
% Can plot combinations of single valued time series and stacked time
% plots.


%% generic plot with varstrs
varstrs = {'E_ts(:,:,:,2)','B_ts(:,:,:,3)','ve12_ts(:,:,:,1)','vi12_ts(:,:,:,1)','ve12_ts(:,:,:,2)','vi12_ts(:,:,:,2)','ni12_ts'};
clim = {[-1 1],[-1 1],2*[-1 1],2*[-1 1],2*[-1 1],2*[-1 1],[0 1],[]};
varstrs = {'E_ts(:,:,:,2)','B_ts(:,:,:,1)','ve12_ts(:,:,:,2)','vi12_ts(:,:,:,2)','ve12_ts(:,:,:,3)','vi12_ts(:,:,:,3)','ni12_ts','A_ts'};
varstrs = {'E_ts(:,:,:,2)','E_ts(:,:,:,3)','vi12_ts(:,:,:,2)','vi12_ts(:,:,:,3)','vi12_ts(:,:,:,3)-ve12_ts(:,:,:,3)','ni12_ts'}; % 'ni12_ts.*vi12_ts(:,:,:,3)-ne12_ts.*ve12_ts(:,:,:,3)'
clim = {[-0.8 0.8],[-0.8 0.8],[-1 1],[-1 1],[-0.5 0.5],[0 0.11],[]};


varstrs = {'E_ts(:,:,:,2)','E_ts(:,:,:,3)','vi12_ts(:,:,:,3)','ni12_ts'}; % 'ni12_ts.*vi12_ts(:,:,:,3)-ne12_ts.*ve12_ts(:,:,:,3)'
varstrs = {'E_ts(:,:,:,2)','E_ts(:,:,:,3)','viz_ts','ni_ts'}; % 'ni12_ts.*vi12_ts(:,:,:,3)-ne12_ts.*ve12_ts(:,:,:,3)'
clim = {[-0.5 0.5],[-0.2 0.2],[-0.1 0.1],[0 0.6],[]};
cmaps = {'blue_red','blue_red','blue_red','candy','blue_red','blue_red','candy'};

varstrs = {'vAz','vi12_ts(:,:,:,3)','ve12_ts(:,:,:,3)',...
           'vAz-vi12_ts(:,:,:,3)','vAz-ve12_ts(:,:,:,3)',...
           'ni12_ts/ni12_ts(1,100,100)',...
           'B_ts(:,:,:,1).^2'};
         
varstrs = {'vAz','vi12_ts(:,:,:,3)',...
           'vAz-vi12_ts(:,:,:,3)',...
           'ni12_ts/ni12_ts(1,100,100)',...
           'B_ts(:,:,:,1).^2'};
clim = {[-1 1],[-1 1],[-1 1],1+0.1*[-1 1],1+0.3*[-1 1]};
cmaps = {'blue_red','blue_red','blue_red','blue_red','blue_red','blue_red','blue_red','candy2'};


varstrs = {'vi12_ts(:,:,:,3)',...
           've12_ts(:,:,:,3)',...     
           'vi12_ts(:,:,:,1)',...
           've12_ts(:,:,:,1)',...           
           'ni12_ts',...
           'B_ts(:,:,:,1).^2'};
 clim = {[-1 1],[-1 1],0.1*[-1 1],0.1*[-1 1],0.11+[-0.02 0.02],1+0.3*[-1 1]};

varstrs = {'E_ts(:,:,:,2)',...
           'E_ts(:,:,:,3)',...
           };
%varstrs = {'A_ts'};
%clim = {[]};

%cmaps = {'candy','candy','candy','candy','candy','candy','candy','candy'};
npanels = numel(varstrs);
nvars = cellfun(@numel, varstrs);

timeaxis = 'x'; % which plot axis should have the time, x or y
plotaxis = 'z'; % 'x' for horizontal cut, 'z' for vertical cut
if strcmp(plotaxis,'z'), meandirection = 1+1;
elseif strcmp(plotaxis,'x'), meandirection = 2+1; end
zpick = 10+[-0.1 0.1];
%zpick = [0];
xpick = -18.5+0.2*[-1 1];
xpick = 24+0.2*[-1 1];
xpick = 15.5+0.1*[-1 1];
%xpick = xLineXY(1)+0.2*[-1 1];
xpick = 6+0.1*[-1 1];
%xpick = 20-5.6+0.2*[-1 1];
%xpick = 32+0.2*[-1 1];
zind = find_closest_ind(z,zpick);
xind = find_closest_ind(x,xpick);

%xlim = x([1 end/2])'+[100 -100];
xlim = [150 x(fix(end/2))];
xlim = [x(fix(end/2)) x(end)-150];
xlim = [-10 100];%x([1 end])'+[150 -150];
xlim = [x(1) x(end)];
zlim = [-5 5];
zlim = [z(1) z(end)] + 0*[1 -1];
%zlim = [-20 20];

ipx1 = find(x>=xlim(1),1,'first');
ipx2 = find(x<=xlim(2),1,'last');
ipz1 = find(z>=zlim(1),1,'first');
ipz2 = find(z<=zlim(2),1,'last');


% Flux function
doAx = 1; % plot separatrix
doA = 1;
cA = 0*[0.8 0.8 0.8];
nA = 20;
nA = [0:-2:-30];
ipxA = ipx1:20:ipx2;
ipzA = ipz1:20:ipz2;

switch plotaxis % the spatial axis
  case 'x'
    ipx = ipx1:1:ipx2;
    ipz = zind(1):zind(end);
    plot_dep_space = x(ipx);
    plot_dep_space_A = x(ipxA);
    ipzA = find(abs(z-mean(zpick))==min(abs(z-mean(zpick))));
    plotlim = xlim;
    pickind = zpick;
    pickval = z(ipz);
    pickstr = 'z';
    plot_dep_space_str = 'x (d_i)';    
  case 'z'
    ipz = ipz1:1:ipz2;
    ipx = xind(1):xind(end);
    plot_dep_space = z(ipz);    
    plot_dep_space_A = z(ipzA);
    ipxA = find(abs(x-mean(xpick))==min(abs(x-mean(xpick))));
    plotlim = zlim;
    pickind = xpick;
    pickval = x(ipx);
    pickstr = 'x';
    plot_dep_space_str = 'z (d_i)';    
end



linewidth = 1.5;
fontsize = 12;

times = timesteps/200;

% Initialize figure
nvars = numel(varstrs);
npanels = nvars;
maxrows = 6;
nrows = min([npanels,maxrows]);
ncols = ceil(npanels/nrows);
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');
linkaxes(h);

% Panels
isub = 1;
tic;
doColor = 0;
for ipanel = 1:npanels
  hca = h(isub); isub = isub + 1;  
  ivar = ipanel;
  varstr = varstrs{ivar};
  variable = eval(varstr);  
  if strcmp(timeaxis,'x')
    dep_x = timesteps;
    dep_y = plot_dep_space;
    dep_xA = timesteps;
    dep_yA = plot_dep_space_A;
    permute_order = [1 2];
    plot_dep_x_str = 'time (\omega_{ce}^{-1})';
    plot_dep_y_str = plot_dep_space_str;
  elseif strcmp(timeaxis,'y')
    dep_y = timesteps;
    dep_x = plot_dep_space;
    dep_yA = timesteps;
    dep_xA = plot_dep_space_A;
    plot_dep_x_str = plot_dep_space_str;
    plot_dep_y_str = 'time (\omega_{ce}^{-1})';
    permute_order = [2 1];
  end
  plot_data = permute(squeeze(mean(variable(:,ipx,ipz),meandirection)),permute_order);
  himag = imagesc(hca,dep_x,dep_y,plot_data'); 
  hca.Box = 'on';
  hca.XLabel.String = plot_dep_x_str;
  hca.YLabel.String = plot_dep_y_str;
  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = varstr;
  hcb.YLabel.Interpreter = 'none';
  hb(ipanel) = hcb;
  
      
  if doA
    hold(hca,'on')
    plot_data = permute(squeeze(mean(A_ts(:,ipxA,ipzA),meandirection)),permute_order);
    hcont = contour(hca,dep_xA,dep_yA,plot_data',nA,'color',cA,'linewidth',0.5,'displayname','A'); 
    hold(hca,'off')  
  end
  if 0*doAx
    hold(hca,'on')
    [saddle_locations,saddle_values] = saddle(A,'sort');
    hcont = contour(hca,x(ipx),z(ipz),squeeze(A(ipx,ipz))',saddle_values(1)*[1 1],'color',cA,'linewidth',2,'displayname','A_X','linestyle','-'); 
    hold(hca,'off')  
  end
  
  if numel(clim)>=ipanel && not(isempty(clim{ipanel}))
    if iscell(clim)
      h(ipanel).CLim = clim{ipanel}; 
    elseif isnumeric(clim) 
      h(ipanel).CLim = clim; 
    end
  end
  if numel(cmaps)>=ipanel && not(isempty(cmaps{ipanel}))
    colormap(hca,pic_colors(cmaps{ipanel}));
  end  

  hca.FontSize = fontsize;
  hold(hca,'off')
  hca.XGrid = 'off';
  hca.YGrid = 'off';
  if strcmp(plotaxis,'x'), hca.XDir = 'reverse'; end
  hca.YDir = 'normal';
  
end
h(1).Title.String = sprintf('%s = [%.2f,%.2f] (d_i)',pickstr,pickval(1),pickval(end));

drawnow
compact_panels(0.010)
for iPanel = 1:(npanels-1)
  h(iPanel).XLabel.String = [];
end
for iPanel = 1:npanels
  h(iPanel).YDir = 'normal';
end



%%
hca = subplot(1,1,1);
ntimes = numel(timesteps);
for itime = 1:ntimes
  imagesc(x,z,squeeze(A_ts(itime,:,:))')
  hold(hca,'on')
  plot(hca,x(R.Xx(itime)),z(R.Xy(itime)),'ko')
  hold(hca,'off')
  pause(0.1)
end

%%

% use load_resaved_data.m
t = timesteps/50;
x = sim_info.x-mean(sim_info.x);
z = sim_info.z;
[Tx,X] = ndgrid(t,x);
[Tz,Z] = ndgrid(t,z);
% make_time_series

lA = -25:0;

%%
xpick = 0;
ix = find_closest_ind(x,xpick);
rx = 0;
ix = ix + [-rx:1:rx];
iz = 1:numel(z);

h = setup_subplots(3,1);
isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,t,z,squeeze(mean(A(:,ix,iz),2))');
shading(hca,'flat')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'z (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'A';

hold(hca,'on')
C = squeeze(mean(A(:,ix,iz),2));
contour(hca,Tz,Z,C,lA,'k')
hold(hca,'off')
hca.CLim = [-25 0];

colormap(hca,pic_colors('candy'))

hca = h(isub); isub = isub + 1;
icomp = 1;
imagesc(hca,t,z,squeeze(mean(B(:,ix,iz,icomp),2))')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'z (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Bx';

hold(hca,'on')
C = squeeze(mean(A(:,ix,iz),2));
contour(hca,Tz,Z,C,lA,'k')
hold(hca,'off')

hca.CLim = 1*[-1 1];
colormap(hca,pic_colors('blue_red'))

hca = h(isub); isub = isub + 1;
icomp = 2;
imagesc(hca,t,z,squeeze(mean(E(:,ix,iz,icomp),2))')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'z (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Ey';

hold(hca,'on')
C = squeeze(mean(A(:,ix,iz),2));
contour(hca,Tz,Z,C,lA,'k')
hold(hca,'off')

hca.CLim = 0.2*[-1 1];
colormap(hca,pic_colors('blue_red'))

%%
%zlim = [-10 10];
%ipx1 = find(x>xlim(1),1,'first');
%ipx2 = find(x<xlim(2),1,'last');
%ipz1 = find(z>zlim(1),1,'first');
%ipz2 = find(z<zlim(2),1,'last');

zpick = 0;
iz = find_closest_ind(z,zpick);
ix = 1:numel(x);

h = setup_subplots(3,1);
isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,t,x(ix),squeeze(A(:,ix,iz))')
shading(hca,'flat')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'x (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'A';

hold(hca,'on')
C = squeeze(A(:,ix,iz));
contour(hca,Tx,X,C,lA,'k')
hold(hca,'off')

colormap(hca,pic_colors('candy'))
hca.CLim = [-25 0];

hca = h(isub); isub = isub + 1;
icomp = 3;
imagesc(hca,t,x,squeeze(B(:,ix,iz,icomp))')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'x (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Bz';

hold(hca,'on')
C = squeeze(A(:,ix,iz));
contour(hca,Tx,X,C,lA,'k')
hold(hca,'off')

hca.CLim = 0.5*[-1 1];
colormap(hca,pic_colors('blue_red'))

hca = h(isub); isub = isub + 1;
icomp = 2;
imagesc(hca,t,x,squeeze(E(:,ix,iz,icomp))')
hca.XLabel.String = 't (wci^{-1})';
hca.YLabel.String = 'x (di)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Ey';

hold(hca,'on')
C = squeeze(A(:,ix,iz));
contour(hca,Tx,X,C,lA,'k')
hold(hca,'off')

hca.CLim = 0.2*[-1 1];
colormap(hca,pic_colors('blue_red'))

hlink = linkprop(h,{'XLim','YLim'});
hlink.Targets(1).YLim = 60*[-1 1];
