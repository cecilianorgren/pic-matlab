% Initialize  simulation grid
xmin = 0; xmax = 50;
dx = 1;
xg = xmin:dx:xmax;

% Initialize particles
np = 100; % number of particles

% Particle positions
% uniform
xp = xmin + (xmax-xmin)*rand(1,np); 
fp = @(xg) xg*0 + 1/(xmax-xmin);

if 1
% normal distribution
% with the nomral distribution, some particles might end up outside the grid
sigma = (3*dx); % the standard deviation specifies the 'density gradient'
xp = (xmax-xmin)/2 + sigma*randn(1,np); 
fp = @(xg) normpdf(xg,(xmax-xmin)/2,sigma);
end

if 1 % generate particles like in the code
% normal distribution
% with the nomral distribution, some particles might end up outside the grid
sigma = (3*dx); % the standard deviation specifies the 'density gradient'
fp = @(xg) normpdf(xg,(xmax-xmin)/2,sigma);
fp_sum = sum(fp(xg));
np_grid_all = zeros(size(xg));
xp = [];
for ig = 1:numel(xg)
  np_grid = round(np*normpdf(xg(ig),(xmax-xmin)/2,sigma)/fp_sum);
  xp_grid = xg(ig)-0.5*dx + dx*rand(1,np_grid); 
  xp = [xp,xp_grid];
  np_grid_all(ig) = np_grid;
end

end

% bspline order of particle shape function.
% Each entry will give one density profile.
orders = [0 1 10 20]; % l

% Loop through orders and get density
ii = 0;
for io = orders
  ii = ii + 1;
  nstruct(ii).n = n(xg,xp,io);
  nstruct(ii).order = io;
  nstruct(ii).xp = xp;
  nstruct(ii).nsum = sum(nstruct(ii).n);
end

% Collect output into arrays so that one can plot them simultaneously
xpall = cat(1,nstruct.xp);
nall = cat(1,nstruct.n);
oall = cat(1,nstruct.order);
nsumall = cat(1,nstruct.nsum);

% Plot
nrows = 2;
hca = subplot(nrows,1,1);
hl = plot(hca,xg,nall/np,'-',xpall,xp*0,'.');
%hl = findobj(hca,'type','line'); 
colors = cat(1,hl(1:numel(orders)).Color);
hca.XGrid = 'on';
hca.XTick = xg;

% plot 'real' distribution on a finer grid so that it is smooth
hold(hca,'on')
grid_refinement = 10;
xg_ref = xg(1):(dx/grid_refinement):xg(end);
n_ref = fp(xg_ref);
hf = plot(hca,xg_ref,n_ref,'k-','linewidth',1);
hold(hca,'off')
nsum_ref = sum(n_ref)/grid_refinement*np;

legs = legend(hca,arrayfun(@(x) sprintf('%g',x),oall,'UniformOutput',false),'box','off');
legs.Title.String = {'bspline order of','particle shape','function:'};
irf_legend(hf,'Reference distribution',[0.02 0.98],'color',hf.Color,'fontweight','bold');
set(hca,"ColorOrder",[colors; hf.Color])
irf_legend(hca,arrayfun(@(x) sprintf('int(n) = %g',x),[nsumall;nsum_ref],'UniformOutput',false),[0.02 0.5],'fontweight','bold');

hca = subplot(nrows,1,2);
histogram(hca,xp,xg)
hca.XGrid = 'on';
hca.XTick = xg;
hold(hca,'on')
plot(xg,np_grid_all)
hold(hca,'off')
hca.YLabel.String = 'Particles in cell';

linkprop(findobj(gcf,'type','axes'),{'XLim'});

function n = n(xg,xp,l)
% n(xg,xp,l) Calculate density from macroparticles of different shape.
% Returns density on grid.
%   xg - simulation grid
%   xp - particle positions
%   l - bspline order of particle shape function

  % Get size of grid, for now only for uniform grid spacing
  dx = xg(2)-xg(1);

  % Define knots and particle shape from particle shape order.
  % How knots is related to the bspline order was just guess work.
  knots = 0:(l+1);
  curve = bspline(knots-mean(knots));
  
  % Initialize 'density'
  n = zeros(size(xg));
   
  % Loop through particles
  for ip = 1:numel(xp)
    ng = fnval(curve,(xg-xp(ip))/dx);
    ng((xg-xp(ip))/dx<-mean(knots)) = 0;
    ng((xg-xp(ip))/dx>mean(knots)) = 0;
    n = n + ng;
  end
end