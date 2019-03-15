%% Load data
for timestep = [4200:200:4800 5200:200:5800]
%timestep = 8000;
txtfile = sprintf('/Volumes/Fountain/Data/PIC/df_cold_protons_1/data/fields-%05.0f.dat',timestep); % michael's perturbation

tic; [x,z,E,B,...
  ni1,ne1,ni2,ne2,...
  vi1,ve1,vi2,ve2,...
  ji1,je1,ji2,je2,...
  pi1,pe1,pi2,pe2,...
  ti1,te1,ti2,te2,...
  dfac,teti,nnx,nnz,wpewce,mass,it,time,dt,xmax,zmax,q] = read_fields(txtfile); toc
x0 = mean(x); x = x - x0;

% Calculate auxillary quantities
A = vector_potential(x,z,B.x,B.z); % vector potential
[saddle_locations,saddle_values] = saddle(A);
pic_calc_script

%% Line plot at given x or z, define variable in cell array
if 0
ylims = {[],[],[],[],[],[],[],[]};
% Define what variables to plot
if 1 % Components of combined equation of motion of all species
  varstrs = {{'B.x','B.y','B.z'},...
             {'ne','ni'},...
             {'J.x','J.y','J.z'},...
             {'ji1.x','je1.x','ji2.x','je2.x'},...
             {'pi1.scalar','pe1.scalar','pi2.scalar','pe2.scalar','pe1.scalar+pe2.scalar+pi1.scalar+pi2.scalar'},...
             {'nmvvi1.xx','nmvve1.xx','nmvvi2.xx','nmvve2.xx','nmvve1.xx+nmvve2.xx+nmvvi1.xx+nmvvi2.xx'}...
             ...%{'0.5*B.abs.^2','BB.xx'},...
             ...%{'pe1.scalar+pe2.scalar+pi1.scalar+pi2.scalar','nmvve1.xx+nmvve2.xx+nmvvi1.xx+nmvvi2.xx','0.5*B.abs.^2','BB.xx'}...
             }; vallim = [];
  colororders = {'xyz','g1','xyz','abcd','abcdg','abcdg','','','',''};
  ylims = {[-0.2 0.8],[0.5 3.5],[-2 2],[0 1.5],[0 1.2],[0 1],[],[]};
end
if 0 % Components of equation of motion for the individual species
  varstrs = {{'B.z'},...
             {'ne1','ne2','ni1','ni2'},...
             {'ve1.x','ve2.x','vi1.x','vi2.x'},...
             {'pe1.xx','pe2.xx','pi1.xx','pi2.xx'},...
             {'nmvve1.xx','nmvve2.xx','nmvvi1.xx','nmvvi2.xx'},...             
             }; vallim = [];
  colororders = {'1','bdacm','bdac','bdac','bdac','cd','','',''};
end
if 0
  varstrs = {{'B.x','E.z'},...
             {'ne1','ne2','ni1','ni2'},...
             {'ve1.x','ve2.x','vi1.x','vi2.x'},...
             {'ve1.y','ve2.y','vi1.y','vi2.y'},...
             {'pe1.scalar','pe2.scalar','pi1.scalar','pi2.scalar'},...
             ...%{'pe1.xx+pe1.xy+pe1.xz','pe2.xx+pe2.xy+pe2.xz','pi1.xx+pi1.xy+pi1.xz','pi2.xx+pi2.xy+pi2.xz'},...
             ...%{'-pB','-0.5*(BB.xx+BB.xy+BB.xz)'},...
             {'0.5*B.abs.^2','pe1.scalar+pe2.scalar+pi1.scalar+pi2.scalar','0.5*B.abs.^2+pe1.scalar+pe2.scalar+pi1.scalar+pi2.scalar'}
             }; vallim = [];
  colororders = {'1y','bdacm','bdac','bdac','bdac','1ao','bdac','','',''};
end
if 0 % equilibrium
  varstrs = {{'B.x'},...
             {'ne1','ne2','ni1','ni2'},...
             {'pe1.scalar','pe2.scalar','pi1.scalar','pi2.scalar'},...                          
             {'0.5*B.abs.^2','pe1.scalar+pe2.scalar+pi1.scalar+pi2.scalar','0.5*B.abs.^2+pe1.scalar+pe2.scalar+pi1.scalar+pi2.scalar'}
             }; vallim = [];
  colororders = {'1y','bdacm','bdac','1ao','bdac','','',''};
end
doAddSpeciesExplanation = 0;
speciesIdentification = {'i1','e1','i2','e2'};
speciesExplanation = {'hot ions','hot electrons','cold ions','cold electrons'};

doDiff = [];%[5 6 7];
var_operation = 'diff';
%doDiff = 1;

npanels = numel(varstrs);
nvars = cellfun(@numel, varstrs);

plotaxis = 'x'; % 'x' for horizontal cut, 'z' for vertical cut
zpick = 0;
xpick = 200;
zind = find_closest_ind(z,zpick);
xind = find_closest_ind(x,xpick);

%xlim = x([1 end/2])'+[100 -100];
xlim = [150 x(fix(end/2))];
xlim = [x(fix(end/2)) x(end)-150];
xlim = x([1 end])'+[150 -150];
xlim = [-0 100];0
zlim = [-10 10];

ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
%ipz1:2:ipz2;
%[X,Z] = ndgrid(x,z);
%Xp = X(ipx,ipz);
%Zp = Z(ipx,ipz);
switch plotaxis
  case 'x'
    ipx = ipx1:2:ipx2;
    ipz = zind;
    plot_dep = x(ipx);
    plotlim = xlim;
    pickind = zpick;
    pickval = z(ipz);
    pickstr = 'z';
  case 'z'
    ipz = ipz1:2:ipz2;
    ipx = xind;
    plot_dep = z(ipz);
    plotlim = zlim;
    pickind = xpick;
    pickval = x(ipx);
    pickstr = 'x';
end

linewidth = 1.5;
fontsize = 12;

% Initialize figure
npanels = npanels;
nrows = npanels;
ncols = ceil(npanels/nrows);
npanels = nrows*ncols;
isub = 1; 
clear h hleg;
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1; 
  hold(h(isub-1),'off')

  %h(isub).Box = 'on';
end

    
% Panels
isub = 1;
tic;
doColor = 0;
for ipanel = 1:npanels
  hca = h(isub); isub = isub + 1;
  colors = pic_colors(colororders{ipanel});
  if size(colors,1)>=(nvars(ipanel)), doColor = 1; else doColor = 0; end
  hplot = [];
  for ivar = 1:nvars(ipanel)  
    if ivar == 1, hold(hca,'off');
    elseif ivar == 2,  hold(hca,'on'); end 
    varstr = varstrs{ipanel}{ivar};
    variable = eval(varstr);
    
    if intersect(doDiff,ipanel)
      diff_variable = [0; diff(variable(ipx,ipz))];
      hplot_tmp = plot(hca,plot_dep,diff_variable,'LineWidth',linewidth);
      varstrs{ipanel}{ivar} = sprintf('diff(%s)',varstrs{ipanel}{ivar});      
    else
      hplot_tmp = plot(hca,plot_dep,variable(ipx,ipz),'LineWidth',linewidth);      
    end
    if doColor, hplot_tmp.Color = colors(ivar,:); end
    hplot{ivar} = hplot_tmp;
    hca.XLabel.String = sprintf('%s (d_i)',plotaxis);
    %hca.YLabel.String = varstr;
    %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
    %hca.Title.String = sprintf('%s',varstr); 
    %hca.YLabel.Interpreter = 'none';        
  end
  leg_str_tmp = varstrs{ipanel};
  if doAddSpeciesExplanation
    for ivar = 1:nvars(ipanel)  
      for ispecies = 1:numel(speciesIdentification)
        if strfind(leg_str_tmp{ivar},speciesIdentification{ispecies})
          leg_str_tmp{ivar} = sprintf('%s (%s)',leg_str_tmp{ivar},speciesExplanation{ispecies});
        end
      end    
    end
  end
  hleg_tmp = legend(hca,leg_str_tmp,'box','off','fontsize',fontsize,'location','eastoutside');
  hleg(ipanel) = hleg_tmp;
  hca.FontSize = fontsize;
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XDir = 'reverse';
  if not(isempty(ylims{ipanel}))
    hca.YLim = ylims{ipanel};
  end
end
%hlink = linkprop(h(4:6),'YLim');
hlink = linkprop(h,'XLim'); 
if strcmp(plotaxis,'x'), h(1).XLim = xlim;
elseif strcmp(plotaxis,'z'), h(1).XLim = zlim;
end
%hlink = linkprop(hleg,'Position');
%addprop(hlink,'PlotBoxAspectRatio')
h(1).Title.String = sprintf('%s = %.2f (d_i), time = %g (1/wci) = %g (1/wpe)',pickstr,pickval,time,timestep);
arrayfun(@(x)eval(sprintf('x.Position(3) = 0.6;'),x),h)
drawnow

if ncols == 1
  compact_panels(0.005);
  arrayfun(@(x)eval(sprintf('x.XLabel.String = [];'),x),h(1:end-1)); 
end
%cn.print(sprintf('momentum_equation_2_x_z0_twpe%05.0f',timestep),'path',savedir_root)
%end
end

%% Force terms of different populations separately
if 0
ylims = {[],[],[],[],[],[],[],[],[],[],[]};
% Define what variables to plot
if 1 % Components of combined equation of motion of all species
  varstrs = {{'B.x','B.y','B.z'},...
             {'ji2.x'},...
             {'nmvvi2.xx'},...
             {'pi2.scalar'},...             
             {'ni2.*E.x','ni2.*E_smooth.x'},...
             {'ni2.*vi2xB.x','ni2.*vi2xB.x_yz','ni2.*vi2xB.x_zy'},...
             {'-gradpi2.x','-gradpi2_smooth.x'},...
             {'gradx_nmvvi2xx.x','div_nmvvi2.x'},...
             {'ni2.*E_smooth.x+ni2.*vi2xB.x-gradpi2_smooth.x'}
             }; vallim = [];
  colororders = {'xyz','cyz','c','c','gc','c','gc','c','',''};
  %ylims = {[-0.2 0.8],[0.5 3.5],[-2 2],[0 1.5],[0 1.2],[0 1],[],[]};
end
if 0 % Components of combined equation of motion of all species
  varstrs = {{'B.x','B.y','B.z'},...
             {'ji1.x'},...
             {'nmvvi1.xx'},...
             {'pi1.scalar'},...             
             {'ni1.*E.x','ni1.*E_smooth.x'},...
             {'ni1.*vi1xB.x','ni1.*vi1xB.x_yz','ni1.*vi1xB.x_zy'},...
             {'-gradpi1.x','-gradpi1_smooth.x'},...
             {'gradx_nmvvi1xx.x','div_nmvvi1.x'},...
             {'ni1.*E_smooth.x+ni1.*vi1xB.x-gradpi1_smooth.x'}
             }; vallim = [];
  colororders = {'xyz','c','cyz','c','c','gc','c','gc','c','',''};
  %ylims = {[-0.2 0.8],[0.5 3.5],[-2 2],[0 1.5],[0 1.2],[0 1],[],[]};
end

doAddSpeciesExplanation = 0;
speciesIdentification = {'i1','e1','i2','e2'};
speciesExplanation = {'hot ions','hot electrons','cold ions','cold electrons'};

doDiff = [];%[5 6 7];
var_operation = 'diff';
%doDiff = 1;

npanels = numel(varstrs);
nvars = cellfun(@numel, varstrs);

plotaxis = 'x'; % 'x' for horizontal cut, 'z' for vertical cut
zpick = 0;
xpick = 200;
zind = find_closest_ind(z,zpick);
xind = find_closest_ind(x,xpick);

%xlim = x([1 end/2])'+[100 -100];
xlim = [150 x(fix(end/2))];
xlim = [x(fix(end/2)) x(end)-150];
xlim = x([1 end])'+[150 -150];
xlim = [-0 80];
zlim = [-10 10];

ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
%ipz1:2:ipz2;
%[X,Z] = ndgrid(x,z);
%Xp = X(ipx,ipz);
%Zp = Z(ipx,ipz);
switch plotaxis
  case 'x'
    ipx = ipx1:2:ipx2;
    ipz = zind;
    plot_dep = x(ipx);
    plotlim = xlim;
    pickind = zpick;
    pickval = z(ipz);
    pickstr = 'z';
  case 'z'
    ipz = ipz1:2:ipz2;
    ipx = xind;
    plot_dep = z(ipz);
    plotlim = zlim;
    pickind = xpick;
    pickval = x(ipx);
    pickstr = 'x';
end

linewidth = 1.5;
fontsize = 12;

% Initialize figure
npanels = npanels;
nrows = npanels;
ncols = ceil(npanels/nrows);
npanels = nrows*ncols;
isub = 1; 
clear h hleg;
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1; 
  hold(h(isub-1),'off')

  %h(isub).Box = 'on';
end

    
% Panels
isub = 1;
tic;
doColor = 0;
for ipanel = 1:npanels
  hca = h(isub); isub = isub + 1;
  colors = pic_colors(colororders{ipanel});
  if size(colors,1)>=(nvars(ipanel)), doColor = 1; else doColor = 0; end
  hplot = [];
  for ivar = 1:nvars(ipanel)  
    if ivar == 1, hold(hca,'off');
    elseif ivar == 2,  hold(hca,'on'); end 
    varstr = varstrs{ipanel}{ivar};
    variable = eval(varstr);
    
    if intersect(doDiff,ipanel)
      diff_variable = [0; diff(variable(ipx,ipz))];
      hplot_tmp = plot(hca,plot_dep,diff_variable,'LineWidth',linewidth);
      varstrs{ipanel}{ivar} = sprintf('diff(%s)',varstrs{ipanel}{ivar});      
    else
      hplot_tmp = plot(hca,plot_dep,variable(ipx,ipz),'LineWidth',linewidth);      
    end
    if doColor, hplot_tmp.Color = colors(ivar,:); end
    hplot{ivar} = hplot_tmp;
    hca.XLabel.String = sprintf('%s (d_i)',plotaxis);
    %hca.YLabel.String = varstr;
    %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
    %hca.Title.String = sprintf('%s',varstr); 
    %hca.YLabel.Interpreter = 'none';        
  end
  leg_str_tmp = varstrs{ipanel};
  if doAddSpeciesExplanation
    for ivar = 1:nvars(ipanel)  
      for ispecies = 1:numel(speciesIdentification)
        if strfind(leg_str_tmp{ivar},speciesIdentification{ispecies})
          leg_str_tmp{ivar} = sprintf('%s (%s)',leg_str_tmp{ivar},speciesExplanation{ispecies});
        end
      end    
    end
  end
  hleg_tmp = legend(hca,leg_str_tmp,'box','off','fontsize',fontsize,'location','eastoutside');
  hleg(ipanel) = hleg_tmp;
  hca.FontSize = fontsize;
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XDir = 'reverse';
  if not(isempty(ylims{ipanel}))
    hca.YLim = ylims{ipanel};
  end
end
%hlink = linkprop(h(4:6),'YLim');
hlink = linkprop(h,'XLim'); 
if strcmp(plotaxis,'x'), h(1).XLim = xlim;
elseif strcmp(plotaxis,'z'), h(1).XLim = zlim;
end
%hlink = linkprop(hleg,'Position');
%addprop(hlink,'PlotBoxAspectRatio')
h(1).Title.String = sprintf('%s = %.2f (d_i), time = %g (1/wci) = %g (1/wpe)',pickstr,pickval,time,timestep);
arrayfun(@(x)eval(sprintf('x.Position(3) = 0.6;'),x),h)
drawnow
arrayfun(@(x)eval(sprintf('x.Interpreter = ''none'';'),x),hleg); 

if ncols == 1
  compact_panels(0.005);
  arrayfun(@(x)eval(sprintf('x.XLabel.String = [];'),x),h(1:end-1)); 
end
%cn.print(sprintf('momentum_equation_2_x_z0_twpe%05.0f',timestep),'path',savedir_root)
%end
end

%% Force terms of different populations separately
ylims = {[],[],[],[],[],[],[],[],[],[],[]};
% Define what variables to plot

if 1 % Components of combined equation of motion of all species
  varstrs = {{'B.x','B.y','B.z'},...
             {'E.x','E.y','E.z','E_smooth.x','E_smooth.y','E_smooth.z'},...
             {'ni1','ni2'},...
             {'pi1.scalar','pi2.scalar'},...
             {'vi1xB.x','vi2xB.x'},...
             {'ni1.*vi1xB.x','ni2.*vi2xB.x'},...
             {'ni2.*E.x','ni2.*E_smooth.x'},...
             {'ni2.*vi2xB.x','ni2.*vi2xB.x_yz','ni2.*vi2xB.x_zy'},...
             {'-gradpi1.x','-gradpi1_smooth.x','-gradpi2.x','-gradpi2_smooth.x'},...
             {'gradx_nmvvi2xx.x','-div_nmvvi2.x'},...
             {'ni2.*E_smooth.x+ni2.*vi2xB.x-gradpi2_smooth.x'}
             }; vallim = [];
  colororders = {'xyz','gggxyz','ac','ac','ac','ac','gc','c','gagc','ac','ac','ac','','',''};
  %ylims = {[-0.2 0.8],[0.5 3.5],[-2 2],[0 1.5],[0 1.2],[0 1],[],[]};
end
if 0 % Components of combined equation of motion of all species
  varstrs = {{'B.x','B.y','B.z'},...
             {'ji1.x'},...
             {'nmvvi1.xx'},...
             {'pi1.scalar'},...             
             {'ni1.*E.x','ni1.*E_smooth.x'},...
             {'ni1.*vi1xB.x','ni1.*vi1xB.x_yz','ni1.*vi1xB.x_zy'},...
             {'-gradpi1.x','-gradpi1_smooth.x'},...
             {'gradx_nmvvi1xx.x','div_nmvvi1.x'},...
             {'ni1.*E_smooth.x+ni1.*vi1xB.x-gradpi1_smooth.x'}
             }; vallim = [];
  colororders = {'xyz','c','cyz','c','c','gc','c','gc','c','',''};
  %ylims = {[-0.2 0.8],[0.5 3.5],[-2 2],[0 1.5],[0 1.2],[0 1],[],[]};
end
if 1 % Compare to Yin's plot
  varstrs = {{'B.z'},...
             {'ni1','ni2','ni1+ni2'},...
             {'ti1.scalar','ti2.scalar'},...
             {'pi1.scalar','pi2.scalar'},...
             {'nmvvi1.xx','nmvvi2.xx','nmvvi1.xx+nmvvi2.xx'},...
             {'nmvvi1.xx','pi1.scalar','pB',...
             'nmvvi1.xx+pi1.scalar+pB'},...
             {'nmvvi1.xx+nmvvi2.xx','pi1.scalar+pi2.scalar','pB',...             
             'nmvvi1.xx+nmvvi2.xx+pi1.scalar+pi2.scalar+pB'}...
             %{'-BB.xx','-BB.xy','-BB.xz','-BB.yy','-BB.yz','-BB.zz'}...
             }; vallim = [];
  colororders = {'1','ac1','ac','ac','ac1','yzx1','yzx1','ac','ac','ac','','',''};
  %ylims = {[-0.2 0.8],[0.5 3.5],[-2 2],[0 1.5],[0 1.2],[0 1],[],[]};
end

doAddSpeciesExplanation = 1;
speciesIdentification = {'i1','e1','i2','e2'};
speciesExplanation = {'hot ions','hot electrons','cold ions','cold electrons'};

doDiff = [];%[5 6 7];
var_operation = 'diff';
%doDiff = 1;

npanels = numel(varstrs);
nvars = cellfun(@numel, varstrs);

plotaxis = 'x'; % 'x' for horizontal cut, 'z' for vertical cut
zpick = 0;
xpick = 200;
zind = find_closest_ind(z,zpick);
xind = find_closest_ind(x,xpick);

%xlim = x([1 end/2])'+[100 -100];
xlim = [100 x(fix(end/2))];
xlim = [x(fix(end/2)) x(end)-150];
xlim = x([1 end])'+[150 -150];
xlim = [-00 50];
zlim = [-10 10];

ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
%ipz1:2:ipz2;
%[X,Z] = ndgrid(x,z);
%Xp = X(ipx,ipz);
%Zp = Z(ipx,ipz);
switch plotaxis
  case 'x'
    ipx = ipx1:2:ipx2;
    ipz = zind;
    plot_dep = x(ipx);
    plotlim = xlim;
    pickind = zpick;
    pickval = z(ipz);
    pickstr = 'z';
  case 'z'
    ipz = ipz1:2:ipz2;
    ipx = xind;
    plot_dep = z(ipz);
    plotlim = zlim;
    pickind = xpick;
    pickval = x(ipx);
    pickstr = 'x';
end

linewidth = 1.5;
fontsize = 12;

% Initialize figure
npanels = npanels;
nrows = npanels;
ncols = ceil(npanels/nrows);
npanels = nrows*ncols;
isub = 1; 
clear h hleg;
for ipanel = 1:npanels  
  h(isub) = subplot(nrows,ncols,ipanel); isub = isub + 1; 
  hold(h(isub-1),'off')

  %h(isub).Box = 'on';
end

    
% Panels
isub = 1;
tic;
doColor = 0;
for ipanel = 1:npanels
  hca = h(isub); isub = isub + 1;
  colors = pic_colors(colororders{ipanel});
  if size(colors,1)>=(nvars(ipanel)), doColor = 1; else doColor = 0; end
  hplot = [];
  for ivar = 1:nvars(ipanel)  
    if ivar == 1, hold(hca,'off');
    elseif ivar == 2,  hold(hca,'on'); end 
    varstr = varstrs{ipanel}{ivar};
    variable = eval(varstr);
    
    if intersect(doDiff,ipanel)
      diff_variable = [0; diff(variable(ipx,ipz))];
      hplot_tmp = plot(hca,plot_dep,diff_variable,'LineWidth',linewidth);
      varstrs{ipanel}{ivar} = sprintf('diff(%s)',varstrs{ipanel}{ivar});      
    else
      hplot_tmp = plot(hca,plot_dep,variable(ipx,ipz),'LineWidth',linewidth);      
    end
    if doColor, hplot_tmp.Color = colors(ivar,:); end
    hplot{ivar} = hplot_tmp;
    hca.XLabel.String = sprintf('%s (d_i)',plotaxis);
    %hca.YLabel.String = varstr;
    %hca.Title.String = sprintf('%s, sum(%s) = %g',varstr,varstr,sum(variable(:))); 
    %hca.Title.String = sprintf('%s',varstr); 
    %hca.YLabel.Interpreter = 'none';        
  end
  leg_str_tmp = varstrs{ipanel};
  if doAddSpeciesExplanation
    for ivar = 1:nvars(ipanel)  
      for ispecies = 1:numel(speciesIdentification)
        if strfind(leg_str_tmp{ivar},speciesIdentification{ispecies})
          leg_str_tmp{ivar} = sprintf('%s (%s)',leg_str_tmp{ivar},speciesExplanation{ispecies});
        end
      end    
    end
  end
  hleg_tmp = legend(hca,leg_str_tmp,'box','off','fontsize',fontsize,'location','eastoutside');
  hleg(ipanel) = hleg_tmp;
  hca.FontSize = fontsize;
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XDir = 'reverse';
  if not(isempty(ylims{ipanel}))
    hca.YLim = ylims{ipanel};
  end
end
%hlink = linkprop(h(4:6),'YLim');
hlink = linkprop(h,'XLim'); 
if strcmp(plotaxis,'x'), h(1).XLim = xlim;
elseif strcmp(plotaxis,'z'), h(1).XLim = zlim;
end
%hlink = linkprop(hleg,'Position');
%addprop(hlink,'PlotBoxAspectRatio')
h(1).Title.String = sprintf('%s = %.2f (d_i), time = %g (1/wci) = %g (1/wpe)',pickstr,pickval,time,timestep);
arrayfun(@(x)eval(sprintf('x.Position(3) = 0.6;'),x),h)
drawnow
arrayfun(@(x)eval(sprintf('x.Interpreter = ''none'';'),x),hleg); 

if ncols == 1
  compact_panels(0.005);
  arrayfun(@(x)eval(sprintf('x.XLabel.String = [];'),x),h(1:end-1)); 
end
h(1).YLim = [0 0.7];
h(2).YLim = [0 3.5];
h(3).YLim = [0 0.7];
h(4).YLim = [0 1];
h(5).YLim = [0 1];
h(6).YLim = [0 1.5];
h(7).YLim = [0 1.5];
cn.print(sprintf('mhd_mom_eq_i_x_z0_twpe%05.0f',timestep),'path',savedir_root)
end
