
pathStat = '/Users/cecilia/MATLAB/pic-matlab/citation_statistics/magnetic_reconnection/Years/';

dirStat = dir([pathStat '*.csv']);
nStat = numel(dirStat);

data = [];
for iStat = 1:nStat  
  tmpData = csvread([dirStat(iStat).folder '/' dirStat(iStat).name],1,0);
  strspl = strsplit(dirStat(iStat).name,{'_','.'});  
  data(iStat).name = strjoin(strspl(1:end-1),' ');
  data(iStat).year = tmpData(:,1);
  data(iStat).ncite = tmpData(:,2);
end

% p
hca = subplot(1,1,1);

set(hca,'ColorOrder',mms_colors('matlab'))  
set(hca,'LineStyleOrder',{'-','-.'},'nextplot','add')
for iStat = 1:nStat
  plot(hca,data(iStat).year,data(iStat).ncite)
  if iStat == 1
    hold(hca,'on')
  elseif iStat == nStat
    hold(hca,'off')
  end   
end

hleg = legend(hca,data.name,'location','northwest');
hleg.Title.String = 'Search term';
hca.XLim = [1970 2020];
hca.YLabel.String = '# Papers';
hca.XLabel.String = 'Year';
hca.Title.String = 'Peer-reviewed papers mentioning ''magnetic reconnection'' + other terms';
irf_legend(hca,{'All search terms include';'magnetic reconnection'},[0.02 0.3],'color',[0 0 0])