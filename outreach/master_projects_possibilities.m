% new master projects
work = {'Spacecraft observations','Global MHD modelling','Programming',...
  'Theory','Dissemination'};
nWorkCat = numel(work);


ip = 1;
projects(ip).title = 'Modelling of space weather';
projects(ip).workload = [0 1 1 0.5 0.2];
ip = ip + 1;

projects(ip).title = 'Electron heating in guide field reconnection';
projects(ip).workload = [1 0 3 1 0.5];
ip = ip + 1;

projects(ip).title = 'Electron heating in guide field reconnection';
projects(ip).workload = [1 0 2 1 0.2];
ip = ip + 1;


% Collect data;
all_workload = vertcat(projects.workload)./repmat(sum(vertcat(projects.workload),2),1,nWorkCat);
all_titles = {projects.title};

hca = subplot(1,1,1);
barh(hca,all_workload,'stacked')
hca.YLabel.String = '';
hca.YTickLabel = {projects.title};
hca.FontSize = 14;

legend(hca,work,'location','eastoutside')

hca.Color = 'white';
hca.XTick = [];
hca.XColor = [1 1 1];
