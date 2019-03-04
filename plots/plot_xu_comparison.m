% need to run plot_energy_spectrum.m first

varstrs = {...%{'B.x','B.y','B.z'},...
            {'B.z'},...
           {'ni1','ni2','ni1+ni2'},...
           ...%{'ve1.x','ve2.x','vi1.x','vi2.x'},...
           ...% {'-pe1.xx','-pe2.xx','-pi1.xx','-pi2.xx'},...
            {'ti1.scalar','ti2.scalar'},...  
            {'pi1.scalar','pi2.scalar'}...  
            }; 
ylabelstrs = {'Magnetic field','Density','Temperature','Pressure','Energy'};
vallim = [];
colororders = {'1xyz','aco','acm','acm','bdac','bdac','cd','bdac','','',''};

doAddSpeciesExplanation = 1;
speciesIdentification = {'i1','e1','i2','e2'};
speciesExplanation = {'hot ions','hot electrons','cold ions','cold electrons'};

doDiff = [];%[5 6 7];
var_operation = 'diff';
%doDiff = 1;

npanels = numel(varstrs)+1;
nvars = cellfun(@numel, varstrs);

plotaxis = 'x'; % 'x' for horizontal cut, 'z' for vertical cut
zpick = 1;
xpick = 200;
zind = find_closest_ind(z,zpick);
xind = find_closest_ind(x,xpick);

%xlim = x([1 end/2])'+[100 -100];
xlim = [150 x(fix(end/2))];
xlim = [x(fix(end/2)) x(end)-150];
xlim = x([1 end])'+[150 -150];
zlim = [-10 10];

ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');

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
for ipanel = 1:(npanels-1)
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
  hleg_tmp = legend(hca,leg_str_tmp,'box','off','fontsize',fontsize,'location','northeast');%,'eastoutside');
  hleg(ipanel) = hleg_tmp;
  hca.FontSize = fontsize;
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XDir = 'reverse';
  hca.YLabel.String = ylabelstrs{ipanel};
end

if 1 % plot energy spectrogram
  hca = h(isub); isub = isub + 1;
  % Plot spectrogram
  zz = zpick;    
  ndists = size(f_dist_all.x,1);
  save_ind = find(f_dist_all.z == zz);
  f_dist_fields = fields(f_dist_all);

  for ifields = 1:numel(f_dist_fields)
    eval(['f_dist_z0.' f_dist_fields{ifields} ' = f_dist_all.' f_dist_fields{ifields} '(save_ind,:,:,:,:);'])
  end

  ispecies = [1 3];
  xplot = f_dist_z0.x;
  yplot = squeeze(f_dist_z0.f_energy_centers(1,ispecies(1),:));
  cplot = squeeze(sum(f_dist_z0.f_dist_mean(:,ispecies,:),2));
  cplot = squeeze(sum(f_dist_z0.f_dist_mean(:,ispecies,:),2).*f_dist_z0.f_energy_centers(:,ispecies(1),:).^2);


  pcolor(hca,-xplot,yplot,log10(cplot)'); shading flat;
  %hca.YScale = 'log';
  hcb = colorbar('peer',hca);
  hca.CLim = [-6 -4.2];
  hca.YLabel.String = 'Energy';
  hca.FontSize = fontsize;
  colormap(hca,'jet')
end

%hlink = linkprop(h(4:6),'YLim');
hlink = linkprop(h,'XLim'); 
if strcmp(plotaxis,'x'), h(1).XLim = xlim;
elseif strcmp(plotaxis,'z'), h(1).XLim = zlim;
end
%hlink = linkprop(hleg,'Position');
%addprop(hlink,'PlotBoxAspectRatio')
h(1).Title.String = sprintf('%s = %.0f (d_i), time = %g (1/wci) = %g (1/wpe)',pickstr,pickval,time,timestep);
arrayfun(@(x)eval(sprintf('x.Position(3) = 0.6;'),x),h)
arrayfun(@(x)eval(sprintf('x.XDir = ''reverse'';'),x),h)
arrayfun(@(x)eval(sprintf('x.XLim(1) = 0;'),x),h)
drawnow
compact_panels(0.012)